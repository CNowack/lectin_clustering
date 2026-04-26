"""
select_seq_clusters.py

Identify candidate novel lectins by finding proteins that bridge
known lectin communities at within-family similarity strength.

Strategy:
  1. Define "known lectin clusters": communities whose lectin-fraction
     meets `min_lectin_fraction` and size meets `min_known_cluster_size`.
  2. For each, compute the median within-cluster edge weight as that
     cluster's admission threshold (scaled by `weight_strictness_multiplier`).
  3. Find candidate proteins (non-lectin, in any community) that connect
     to >= `min_bridge_count` known clusters via edges meeting each
     respective cluster's threshold.
  4. Emit candidates with the clusters they bridge and edge weights, plus
     a summary file for sanity-checking thresholds.

Recall: link weights are -log10(E-value), so HIGHER weight = stronger match.
"""

import pandas as pd
import numpy as np
from collections import defaultdict

# Snakemake injection
nodes_in = snakemake.input.nodes
links_in = snakemake.input.links
candidates_out = snakemake.output.target_seq_clusters
summary_out = snakemake.output.summary

# Config-driven parameters
MIN_LECTIN_FRACTION = snakemake.config["min_lectin_fraction"]
MIN_BRIDGE_COUNT = snakemake.config["min_bridge_count"]
WEIGHT_MULT = snakemake.config["weight_strictness_multiplier"]
LECTIN_LABEL = snakemake.config["lectin_category_label"]
MIN_CLUSTER_SIZE = snakemake.config["min_known_cluster_size"]

# =============================================================================
# 1) LOAD
# =============================================================================

nodes = pd.read_csv(nodes_in)
links = pd.read_csv(links_in)

print(f"Loaded {len(nodes):,} nodes, {len(links):,} links")

# Lookup: node id -> community_group. Used to map endpoints of every edge
# to their communities in O(1).
id_to_community = dict(zip(nodes["id"], nodes["community_group"]))

# =============================================================================
# 2) IDENTIFY KNOWN LECTIN CLUSTERS
# =============================================================================
# A "known lectin cluster" is a community where the lectin fraction meets
# the threshold and the community is large enough to have a meaningful
# internal weight distribution.

community_stats = (
    nodes.assign(is_lectin=(nodes["source_set"] == LECTIN_LABEL))
    .groupby("community_group")
    .agg(size=("id", "count"), lectin_fraction=("is_lectin", "mean"))
)

known_clusters = community_stats[
    (community_stats["lectin_fraction"] >= MIN_LECTIN_FRACTION)
    & (community_stats["size"] >= MIN_CLUSTER_SIZE)
].index

print(
    f"Known lectin clusters: {len(known_clusters):,} "
    f"(of {len(community_stats):,} total communities)"
)

if len(known_clusters) == 0:
    raise ValueError(
        "No communities passed the known-lectin-cluster filters. "
        "Check `min_lectin_fraction`, `min_known_cluster_size`, and the "
        "`source_set` column in your nodes file."
    )

# =============================================================================
# 3) COMPUTE PER-CLUSTER ADMISSION THRESHOLD
# =============================================================================
# For each known cluster, the threshold is the median weight of edges
# whose BOTH endpoints are in that cluster, scaled by WEIGHT_MULT.
# This represents "how similar a protein has to be to count as belonging
# to this family, in this family's own terms".

# Annotate every edge with the community of each endpoint.
links = links.assign(
    source_community=links["source"].map(id_to_community),
    target_community=links["target"].map(id_to_community),
)

# Within-cluster edges: source and target in the same known cluster.
within_cluster = links[
    (links["source_community"] == links["target_community"])
    & (links["source_community"].isin(known_clusters))
]

cluster_thresholds = (
    within_cluster.groupby("source_community")["weight"].median() * WEIGHT_MULT
)

# Some known clusters may have no within-cluster edges retained after the
# top-4 trim. Drop them — we have no basis for a threshold.
clusters_with_threshold = cluster_thresholds.index
dropped = set(known_clusters) - set(clusters_with_threshold)
if dropped:
    print(
        f"Warning: {len(dropped)} known clusters have no internal edges and "
        f"will be skipped."
    )

print(
    f"Threshold weights — min: {cluster_thresholds.min():.2f}, "
    f"median: {cluster_thresholds.median():.2f}, "
    f"max: {cluster_thresholds.max():.2f}"
)

# =============================================================================
# 4) FIND BRIDGE EDGES TO KNOWN CLUSTERS
# =============================================================================
# A "bridge edge" goes from a non-lectin candidate to a member of a known
# lectin cluster, with weight >= that cluster's threshold. Edges in the
# trimmed graph are directionally stored (we kept top-4 outbound per query),
# but biologically the relationship is symmetric — count both directions.

lectin_ids = set(nodes.loc[nodes["source_set"] == LECTIN_LABEL, "id"])

# Build candidate->cluster bridges from both edge directions.
# Each entry: candidate_id -> {cluster_id: best_weight_seen}
bridges = defaultdict(dict)

def record_bridges(edges, candidate_col, lectin_col, cluster_col):
    """
    Iterate edges where `candidate_col` holds the candidate id and
    `lectin_col` holds the known-lectin endpoint, with `cluster_col`
    giving the lectin's cluster.
    """
    for cand, cluster, weight in zip(
        edges[candidate_col], edges[cluster_col], edges["weight"]
    ):
        threshold = cluster_thresholds.get(cluster)
        if threshold is None or weight < threshold:
            continue
        # Keep the strongest weight per (candidate, cluster) pair.
        prev = bridges[cand].get(cluster)
        if prev is None or weight > prev:
            bridges[cand][cluster] = weight

# Direction A: candidate is `source`, lectin is `target`.
edges_a = links[
    (~links["source"].isin(lectin_ids))
    & (links["target"].isin(lectin_ids))
    & (links["target_community"].isin(clusters_with_threshold))
]
record_bridges(edges_a, "source", "target", "target_community")

# Direction B: candidate is `target`, lectin is `source`.
edges_b = links[
    (~links["target"].isin(lectin_ids))
    & (links["source"].isin(lectin_ids))
    & (links["source_community"].isin(clusters_with_threshold))
]
record_bridges(edges_b, "target", "source", "source_community")

print(f"Candidates with at least one qualifying bridge: {len(bridges):,}")

# =============================================================================
# 5) APPLY MIN_BRIDGE_COUNT AND EMIT
# =============================================================================

candidate_records = []
for cand_id, cluster_to_weight in bridges.items():
    if len(cluster_to_weight) < MIN_BRIDGE_COUNT:
        continue
    candidate_records.append({
        "id": cand_id,
        "n_bridged_clusters": len(cluster_to_weight),
        "bridged_clusters": ";".join(map(str, sorted(cluster_to_weight))),
        "max_bridge_weight": max(cluster_to_weight.values()),
        "mean_bridge_weight": np.mean(list(cluster_to_weight.values())),
    })

candidates_df = pd.DataFrame(candidate_records)

# Carry forward useful node metadata so downstream steps don't have to re-join.
if not candidates_df.empty:
    candidates_df = candidates_df.merge(
        nodes[["id", "community_group", "source_set", "subset"]],
        on="id",
        how="left",
    )

candidates_df.to_csv(candidates_out, index=False)

# =============================================================================
# 6) SUMMARY — for sanity-checking thresholds before structural clustering
# =============================================================================

summary_lines = [
    f"Total nodes:                     {len(nodes):,}",
    f"Total communities:               {len(community_stats):,}",
    f"Known lectin clusters:           {len(known_clusters):,}",
    f"  ...with usable thresholds:     {len(clusters_with_threshold):,}",
    f"Candidates passing filters:      {len(candidates_df):,}",
    f"",
    f"Bridge-count distribution among candidates:",
]
if not candidates_df.empty:
    bridge_dist = candidates_df["n_bridged_clusters"].value_counts().sort_index()
    for n, count in bridge_dist.items():
        summary_lines.append(f"  bridges {n}: {count:,} candidates")

summary_lines += [
    "",
    "Per-cluster threshold weights (top 20 by size):",
]
top_clusters = community_stats.loc[clusters_with_threshold].nlargest(20, "size")
for cluster_id, row in top_clusters.iterrows():
    threshold = cluster_thresholds[cluster_id]
    n_cand = sum(
        1 for v in bridges.values() if cluster_id in v
    )
    summary_lines.append(
        f"  cluster {cluster_id}: size={int(row['size'])}, "
        f"threshold={threshold:.2f}, candidates_bridged={n_cand}"
    )

summary_text = "\n".join(summary_lines)
print("\n" + summary_text)

with open(summary_out, "w") as f:
    f.write(summary_text + "\n")