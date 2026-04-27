# scripts/prepare_cosmograph.py
import pandas as pd

# --- Parse family membership from fasta headers ---
def parse_ids(fasta_path):
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # UniProt format: >sp|A0A0P0VIP0|LRSK7_ORYSJ or just >A0A0P0VIP0
                parts = line.strip().lstrip(">").split("|")
                acc = parts[1] if len(parts) > 1 else parts[0].split()[0]
                ids.append(acc)
    return set(ids)

ig_ids        = parse_ids("../data/ig.fasta")
lectin_ids    = parse_ids("../data/lectins.fasta")
ribosomal_ids = parse_ids("../data/ribosomal.fasta")
kinase_ids    = parse_ids("../data/kinases.fasta")

def assign_family(acc):
    if acc in lectin_ids:    return "Lectin"
    if acc in kinase_ids:    return "Kinase"
    if acc in ig_ids:        return "Immunoglobulin"
    if acc in ribosomal_ids: return "Ribosomal"
    return "Other"

# --- Nodes ---
nodes = pd.read_csv("../results/foldtree_input/finalset.csv", usecols=["Entry", "Protein names", "Organism"])
nodes = nodes.rename(columns={"Entry": "id"})
nodes["protein_family"] = nodes["id"].apply(assign_family)

# --- Links ---
links = pd.read_csv(
    "../results/foldtree_input/allvall_1.csv",
    sep="\t", header=None,
    usecols=[0, 1, 12],
    names=["source", "target", "tmscore"]
)

# After links filtering, filter nodes with cluster size >= 4
from collections import Counter

# Filter isolated nodes — keep only those with at least 4 connections
degree = pd.concat([links["source"], links["target"]]).value_counts()
valid_ids = degree[degree >= 4].index
nodes = nodes[nodes["id"].isin(valid_ids)]
links = links[links["source"].isin(valid_ids) & links["target"].isin(valid_ids)]

print(nodes["protein_family"].value_counts())
nodes.to_csv("../results/foldtree_input/cosmo_nodes.csv", index=False)

links = links[links["source"] != links["target"]]
links = links[links["tmscore"] >= 0.5]
links.to_csv("../results/foldtree_input/cosmo_links.csv", index=False)

print(f"Nodes: {len(nodes)}, Links: {len(links)}")