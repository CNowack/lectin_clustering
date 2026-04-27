# scripts/filter_cosmograph.py
import pandas as pd

# --- Load ---
links = pd.read_csv("results/foldtree_input/cosmo_links.csv")
nodes = pd.read_csv("results/foldtree_input/cosmo_nodes.csv")

# --- Filter nodes with < 4 connections ---
degree = pd.concat([links["source"], links["target"]]).value_counts()
valid_ids = degree[degree >= 4].index

links = links[links["source"].isin(valid_ids) & links["target"].isin(valid_ids)]
nodes = nodes[nodes["id"].isin(valid_ids)]

print(f"Nodes: {len(nodes)}, Links: {len(links)}")

# --- Save ---
nodes.to_csv("results/foldtree_input/cosmo_nodes_filtered.csv", index=False)
links.to_csv("results/foldtree_input/cosmo_links_filtered.csv", index=False)
print("Saved to results/foldtree_input/cosmo_nodes_filtered.csv and cosmo_links_filtered.csv")