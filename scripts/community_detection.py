"""
community_detection.py

Build the protein similarity graph from MMseqs2 all-vs-all output, then
identify protein communities using a two-step procedure that matches
Durairaj et al. (2023):

    1. Trim the graph by keeping at most the top-4 outbound edges per node
       (sparsifies the graph and emphasizes closest relationships).
    2. Find connected components (broad supergroups of homologous sequences).
    3. Within each connected component, run asynchronous label propagation
       to identify finer-grained communities.

Each node ends up labeled with both:
    - connected_component_id : the broad cluster
    - community_group        : the sub-cluster (globally unique across
                               the whole graph)

Run as a Snakemake script via:
    rule community_detection:
        input:  tsv = "results/all_vs_all_alignment.tsv"
        output: nodes = "results/map_nodes.csv",
                links = "results/map_links.csv"
        script: "scripts/community_detection.py"
"""

import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms.community import asyn_lpa_communities


# Snakemake injects these when the rule fires.
input_tsv = snakemake.input.tsv
nodes_out = snakemake.output.nodes
links_out = snakemake.output.links

# Maximum outbound edges retained per node. Matches paper.
MAX_OUTBOUND_EDGES = 4

# =============================================================================
# 1) LOAD ALIGNMENT — MMseqs2 easy-search default columns
# =============================================================================

columns = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bitscore",
]
df = pd.read_csv(input_tsv, sep="\t", names=columns)

# Drop self-alignments — they distort community layout and topology.
df = df[df["query"] != df["target"]]

# =============================================================================
# 2) WEIGHT EDGES — -log10(E-value), with floor for zero-E-value entries
# =============================================================================

# Some E-values come back as 0.0 from MMseqs2 (very strong matches that
# underflow). Replace with the smallest observed nonzero, or 1e-300 if
# no nonzero entries exist (defensive fallback).
nonzero_evalues = df.loc[df["evalue"] > 0, "evalue"]
min_evalue = nonzero_evalues.min() if not nonzero_evalues.empty else 1e-300
df.loc[df["evalue"] == 0, "evalue"] = min_evalue
df["weight"] = -np.log10(df["evalue"])

# =============================================================================
# 3) TRIM TO TOP-4 OUTBOUND EDGES PER NODE  (paper-faithful sparsification)
# =============================================================================
# For each `query` (source node), keep only its 4 strongest edges by weight.
# This dramatically reduces graph density while preserving each node's
# closest relationships - which is exactly what the community structure
# should reflect.
#
# Note: this does NOT cap edges symmetrically. A target node might receive
# many inbound edges. That's intentional - popular hubs are still hubs.
# The sparsification only limits how many neighbors each node "votes" for
# during propagation.

df_trimmed = (
    df.sort_values("weight", ascending=False)
      .groupby("query", as_index=False)
      .head(MAX_OUTBOUND_EDGES)
)

print(
    f"Edges before trim: {len(df):,}  "
    f"after trim: {len(df_trimmed):,}  "
    f"({100 * len(df_trimmed) / len(df):.1f}% retained)"
)

# =============================================================================
# 4) BUILD GRAPH — undirected, edge weights for label propagation
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
# Each connected component is a maximal set of nodes reachable through
# significant-alignment edges. These are guaranteed-homologous groups: every
# pair of sequences in the same component is connected via a chain of
# alignments above the E-value/coverage thresholds.

connected_components = list(nx.connected_components(G))
print(f"Connected components: {len(connected_components):,}")

# =============================================================================
# 6) WITHIN EACH COMPONENT — async label propagation for sub-communities
# =============================================================================
# For each connected component, build the induced subgraph and run async
# LPA on it. Each community gets a globally unique ID by combining
# component_id with a within-component community index.

node_records = []                     # rows for the nodes CSV
global_community_id = 0               # monotonic counter across all components

for component_id, member_nodes in enumerate(connected_components):
    # Singletons and tiny components: skip the LPA step, treat the whole
    # component as a single community. LPA on 1-2 node graphs is trivial
    # and just adds overhead.
    if len(member_nodes) <= 2:
        for node in member_nodes:
            node_records.append({
                "id": node,
                "connected_component_id": component_id,
                "community_group": global_community_id,
            })
        global_community_id += 1
        continue

    # Induce the subgraph on this component. networkx returns a view, which
    # is cheap. asyn_lpa_communities yields one set of node IDs per community.
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
# 7) EXPORT — nodes (with both labels) and links (Cosmograph-compatible)
# =============================================================================

nodes_df = pd.DataFrame(node_records)

# Nodes csv has headers:    'id', 'connected_component_id', 'community_group'
# Example row:              dark_xantho_tr|A0AAJ6H1K8|A0AAJ6H1K8_9XANT,0,0

# parse node 'id' str to build other data columns for further analysis
# use most efficient method, list comprehension, and only once
nodes_df[['source_set', 'subset']] = nodes_df['id'].str.split('_', n=2, expand=True).iloc[:, :2]

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