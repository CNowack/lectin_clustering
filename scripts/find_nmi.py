import pandas as pd
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import AgglomerativeClustering
import numpy as np

# --- Load links and filter nodes with < 4 connections ---
links = pd.read_csv("results/foldtree_input/cosmo_links.csv")
degree = pd.concat([links["source"], links["target"]]).value_counts()
valid_ids = degree[degree >= 4].index

links = links[links["source"].isin(valid_ids) & links["target"].isin(valid_ids)]

# --- Parse PHYLIP distance matrix ---
with open("results/foldtree_input/foldtree_fastmemat.txt") as f:
    n = int(f.readline().strip())
    ids, rows = [], []
    for line in f:
        parts = line.strip().split()
        ids.append(parts[0])
        rows.append([float(x) for x in parts[1:]])

# Filter distance matrix to valid_ids only
dist_df = pd.DataFrame(rows, index=ids, columns=ids)
dist_df = dist_df.loc[dist_df.index.isin(valid_ids), dist_df.columns.isin(valid_ids)]

# --- Cut into 4 clusters ---
clustering = AgglomerativeClustering(n_clusters=4, metric="precomputed", linkage="average")
labels = clustering.fit_predict(dist_df.values)
fold_df = pd.DataFrame({"id": dist_df.index, "fold_cluster": labels})

# --- Ground truth families ---
nodes = pd.read_csv("results/foldtree_input/cosmo_nodes.csv", usecols=["id", "protein_family"])
nodes = nodes[nodes["protein_family"] != "Other"]

# --- Merge and NMI ---
df = nodes.merge(fold_df, on="id", how="inner")
print(f"Matched proteins: {len(df)}")
print(df["protein_family"].value_counts(), "\n")

nmi = normalized_mutual_info_score(df["protein_family"], df["fold_cluster"])
print(f"NMI foldtree clusters vs protein families: {nmi:.4f}")