# Walkthrough
## Step 1 
— Load. Read both CSVs. Build an id → community_group dict so we can look up community membership for any edge endpoint without repeated joins.
## Step 2 
— Find known lectin clusters. Group nodes by community_group, compute size and lectin-fraction. A cluster is "known" if its lectin-fraction meets min_lectin_fraction (default 1.0 — purely lectin) and it has at least min_known_cluster_size members (default 3 — needed for a meaningful internal median). Hard-error if zero pass; that almost certainly means a misconfigured category column rather than a real biological signal.
## Step 3 
— Compute per-cluster thresholds. For each known cluster, take all edges where both endpoints are inside the cluster — these are the within-family alignments. The median of their weights is that cluster's "this is what within-family similarity looks like" benchmark. Multiply by weight_strictness_multiplier to get the admission threshold. Tight families get high thresholds; divergent families get lower ones — automatically. Clusters with no surviving within-cluster edges (post top-4 trim) get dropped with a warning, since there's no basis to set a threshold.
## Step 4 
— Find bridge edges. A bridge goes from an unknown protein to a known lectin in a known cluster, with weight ≥ that cluster's threshold. Because the trimmed graph stored edges directionally (top-4 outbound per query), I scan both directions and merge. For each candidate I keep a dict of {bridged_cluster: best_weight_to_that_cluster}, which both deduplicates multiple edges to the same cluster and records the strongest evidence.
## Step 5 
— Apply bridge count and emit. Filter candidates whose bridge dict has at least min_bridge_count distinct clusters. With the default of 1, this just means "any unknown that hits any known lectin cluster at threshold strength." Bumping to 2 reproduces the paper. The output preserves community metadata so structural clustering can ask "which sequence community did this candidate come from" later.
## Step 6 
— Summary file. Counts at each filter, the bridge-count distribution, and the top-20 largest clusters with their thresholds and how many candidates bridged each. This is your sanity check before kicking off FoldTree — if one cluster is responsible for 90% of candidates, or if every threshold is suspiciously low, you'll see it here.

**The two main tuning knobs:**

* min_bridge_count: precision lever. 1 → inclusive (current default). 2 → unifiers across families (paper). Higher → only the most striking bridges.
* weight_strictness_multiplier: stringency lever, orthogonal to bridge count. Below 1.0 admits weaker matches; above 1.0 demands stronger-than-median similarity.

You can also tighten min_lectin_fraction to, say, 0.8 if you want to admit communities that are mostly lectin but include a couple of unknowns (which can happen when an unknown has clustered tightly with lectins — slightly circular but defensible). Default 1.0 is the conservative choice.