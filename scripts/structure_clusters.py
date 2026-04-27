"""
structure_clusters.py

Convert FoldTree's all-vs-all structural alignment output into a Cosmograph-
ready graph, then detect communities the same way `community_detection.py`
does for sequence data: top-K trim, connected components, asynchronous LPA
within each component.

Input is FoldTree's `allvall_1.csv` with columns:
    query, target, fident, alnlen, mismatch, gapopen,
    qstart, qend, tstart, tend, evalue, bits, lddt, lddtfull, alntmscore

We use `alntmscore` (TM-score from the alignment) as the structural similarity
metric. It ranges 0-1 with higher = more similar, so unlike the sequence
pipeline we don't take -log10 — the score is already a usable weight.

Output mirrors community_detection.py:
    structure_nodes.csv: id, connected_component_id, community_group,
                        source_set, subset
    structure_links.csv: source, target, weight

Run as a Snakemake script via:
    rule detect_structure_clusters:
        input:  data = "results/foldtree_run/allvall_1.csv"
        output: nodes = "results/clusters/structure_nodes.csv",
                links = "results/clusters/structure_links.csv"
        script: "scripts/structure_clusters.py"
"""

import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms.community import asyn_lpa_communities


# Snakemake injects these when the rule fires.
input_csv = snakemake.input.data
nodes_out = snakemake.output.nodes
links_out = snakemake.output.links

# Optional config knobs with defaults — no required new keys in config.yaml.
MIN_SCORE = snakemake.config.get("structure_min_score", 0.5)
MAX_OUTBOUND_EDGES = snakemake.config.get("structure_max_outbound_edges", 4)
SCORE_COLUMN = snakemake.config.get("structure_score_column", "alntmscore")

# =============================================================================
# 1) LOAD ALIGNMENT
# =============================================================================
# FoldTree's allvall_1.csv has a 15-column schema, but one column (`lddtfull`)
# is itself a comma-separated list of per-residue scores. Foldseek writes it
# unquoted, so naive CSV parsing sees hundreds of extra fields per row and
# pandas chokes on the ragged shape. We parse manually by position: the first
# 13 fields and the last field are stable; everything between is the lddtfull
# mess we don't need.

import csv

WANTED_COLUMNS = ["query", "target", "fident", "alnlen", "evalue", "bits", "lddt", "alntmscore"]

records = []
with open(input_csv) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)  # discard header row
    for row in reader:
        if len(row) < 15:
            continue
        records.append({
            "query":      row[0],
            "target":     row[1],
            "fident":     row[2],
            "alnlen":     row[3],
            "evalue":     row[10],
            "bits":       row[11],
            "lddt":       row[12],
            "alntmscore": row[-1],
        })

df = pd.DataFrame(records)
print(f"Loaded {len(df):,} pairwise alignments from {input_csv}")

if SCORE_COLUMN not in df.columns:
    raise ValueError(
        f"Score column '{SCORE_COLUMN}' not found in parsed data. "
        f"Available columns: {list(df.columns)}. "
        f"If you set `structure_score_column` to something other than "
        f"the standard Foldseek output names, add it to WANTED_COLUMNS "
        f"in the parser above."
    )

df[SCORE_COLUMN] = pd.to_numeric(df[SCORE_COLUMN], errors="coerce")
df = df.dropna(subset=["query", "target", SCORE_COLUMN])
df = df[df["query"] != df["target"]]

print(f"After dropping self-alignments and NaNs: {len(df):,} edges")

# =============================================================================
# 2) FILTER BY MINIMUM SCORE
# =============================================================================
# FoldTree runs Foldseek with `-e inf` to retain every pair, including very
# weak structural matches. Most of those are noise and would just produce a
# hairball graph. Cut at MIN_SCORE before doing any community detection.

df = df[df[SCORE_COLUMN] >= MIN_SCORE]
print(f"After score filter (>= {MIN_SCORE}): {len(df):,} edges")

if df.empty:
    raise ValueError(
        f"No edges passed the score filter at {MIN_SCORE}. "
        f"Check `structure_min_score` config or the FoldTree output."
    )

# Standardize the weight column name for downstream code, regardless of which
# metric was chosen.
df = df.rename(columns={SCORE_COLUMN: "weight"})

# =============================================================================
# 3) TRIM TO TOP-K OUTBOUND EDGES PER NODE  (matches sequence pipeline)
# =============================================================================
# Same sparsification as community_detection.py. Each query keeps its top-K
# strongest structural matches. Inbound edges are uncapped.

df_trimmed = (
    df.sort_values("weight", ascending=False)
      .groupby("query", as_index=False)
      .head(MAX_OUTBOUND_EDGES)
)

print(
    f"Edges before trim: {len(df):,}  "
    f"after trim: {len(df_trimmed):,}  "
    f"({100 * len(df_trimmed) / max(len(df), 1):.1f}% retained)"
)

# =============================================================================
# 4) BUILD GRAPH
# =============================================================================

G = nx.from_pandas_edgelist(
    df_trimmed,
    source="query",
    target="target",
    edge_attr="weight",
)

print(f"Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

# =============================================================================
# 5) CONNECTED COMPONENTS — broad superclusters
# =============================================================================

connected_components = list(nx.connected_components(G))
print(f"Connected components: {len(connected_components):,}")

# =============================================================================
# 6) WITHIN EACH COMPONENT — async label propagation for sub-communities
# =============================================================================

node_records = []
global_community_id = 0

for component_id, member_nodes in enumerate(connected_components):
    if len(member_nodes) <= 2:
        for node in member_nodes:
            node_records.append({
                "id": node,
                "connected_component_id": component_id,
                "community_group": global_community_id,
            })
        global_community_id += 1
        continue

    subgraph = G.subgraph(member_nodes)
    sub_communities = asyn_lpa_communities(
        subgraph,
        weight="weight",
        seed=42,
    )

    for community_members in sub_communities:
        for node in community_members:
            node_records.append({
                "id": node,
                "connected_component_id": component_id,
                "community_group": global_community_id,
            })
        global_community_id += 1

print(
    f"Communities found: {global_community_id:,}  "
    f"(avg {global_community_id / max(len(connected_components), 1):.2f} "
    f"per connected component)"
)

# =============================================================================
# 7) PARSE NODE METADATA
# =============================================================================
# FoldTree uses bare UniProt accessions as identifiers (e.g., "A0AAJ6H1K8"),
# not the pipe-delimited form from the sequence step. So we can't just split
# on '|' the way community_detection.py does.
#
# To recover source_set / subset for downstream visualization, join against
# the sequence pipeline's seq_nodes.csv if it's available. The accession is
# the middle field of the original ID (dark_xantho_tr|A0AAJ6H1K8|...).

nodes_df = pd.DataFrame(node_records)

# Try to recover source_set and subset by joining with seq_nodes if present.
# This is best-effort; if seq_nodes isn't available the columns are filled
# with "unknown" so the schema stays consistent with the sequence output.
seq_nodes_path = snakemake.config.get(
    "seq_nodes_for_metadata",
    "results/clusters/seq_nodes.csv",
)

try:
    seq_nodes = pd.read_csv(seq_nodes_path)
    # Extract the accession (middle pipe field) as the join key.
    seq_nodes["accession"] = (
        seq_nodes["id"].str.split("|", n=2, expand=True).iloc[:, 1]
    )
    metadata = (
        seq_nodes[["accession", "source_set", "subset"]]
        .drop_duplicates(subset="accession")
        .set_index("accession")
    )
    nodes_df = nodes_df.merge(
        metadata,
        left_on="id",
        right_index=True,
        how="left",
    )
    n_matched = nodes_df["source_set"].notna().sum()
    print(
        f"Matched {n_matched:,} of {len(nodes_df):,} nodes to sequence "
        f"metadata via accession"
    )
    nodes_df["source_set"] = nodes_df["source_set"].fillna("unknown")
    nodes_df["subset"] = nodes_df["subset"].fillna("unknown")
except FileNotFoundError:
    print(
        f"Sequence node metadata not found at {seq_nodes_path}; "
        f"filling source_set/subset with 'unknown'."
    )
    nodes_df["source_set"] = "unknown"
    nodes_df["subset"] = "unknown"

# =============================================================================
# 8) EXPORT
# =============================================================================

nodes_df.to_csv(nodes_out, index=False)

links_df = (
    df_trimmed[["query", "target", "weight"]]
    .rename(columns={"query": "source"})
)
links_df.to_csv(links_out, index=False)

print(
    f"\nWrote {len(nodes_df):,} nodes "
    f"({nodes_df['connected_component_id'].nunique():,} components, "
    f"{nodes_df['community_group'].nunique():,} communities) "
    f"and {len(links_df):,} links."
)