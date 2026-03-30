import pandas as pd
import networkx as nx
import numpy as np
from networkx.algorithms.community import asyn_lpa_communities

# Snakemake passes these variables automatically
input_tsv = snakemake.input.tsv
nodes_out = snakemake.output.nodes
links_out = snakemake.output.links

# Load the MMseqs2 easy-search default format
columns = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen", 
    "qstart", "qend", "tstart", "tend", "evalue", "bitscore"
]
df = pd.read_csv(input_tsv, sep='\t', names=columns)

# Remove self-alignments (a node linking to itself distorts the visualization layout)
df = df[df['query'] != df['target']]

# Calculate edge weights proportional to E-value. 
# We use -log10(evalue). We must handle E-values of 0.0 to prevent infinite weights.
min_evalue = df[df['evalue'] > 0]['evalue'].min()
if pd.isna(min_evalue):
    min_evalue = 1e-300
df['evalue'] = df['evalue'].replace(0, min_evalue)
df['weight'] = -np.log10(df['evalue'])

# Construct the undirected graph
G = nx.from_pandas_edgelist(df, source='query', target='target', edge_attr='weight')

# Execute Asynchronous Label Propagation
communities = asyn_lpa_communities(G, weight='weight')

# Map each protein to its detected community
node_mapping = []
for community_id, comm in enumerate(communities):
    for node in comm:
        node_mapping.append({"id": node, "community_group": community_id})

# Export Nodes for Cosmograph
nodes_df = pd.DataFrame(node_mapping)
nodes_df.to_csv(nodes_out, index=False)

# Export Links (Edges) for Cosmograph
links_df = df[['query', 'target', 'weight']].rename(columns={'query': 'source'})
links_df.to_csv(links_out, index=False)