import pandas as pd
import os
from cosmograph import cosmo, set_api_key

# web app api key
cosmo_api_key = os.environ.get("COSMO_API_KEY")
set_api_key(cosmo_api_key)

# snakemake inputs
nodes = snakemake.input.nodes
links = snakemale.input.links
graph = snakemake.output.graph
cosmo_project = snakemake.input.project

points = pd.read_csv(nodes)
links = pd.read_csv(links)

# Nodes csv has headers:    'id', 'connected_component_id', 'community_group', 'category', 'type'
# Example row:              dark_xantho_tr|A0AAJ6H1K8|A0AAJ6H1K8_9XANT,0,0,dark,xantho

# Links csv has headers:    'source', 'target', 'weight'
# Example row:              lectin_nr_tr|B9GAY9|B9GAY9_ORYSJ,lectin_nr_tr|I1R0K4|I1R0K4_ORYGL,323.3062153431158

g = cosmo(
    points = points,
    links = links,
    point_id_by = 'id',
    point_label_by = 'id',
    point_color_by = 'category',
    point_shape_by = 'type',
    link_source_by = 'source',
    link_target_by = 'target',
    link_strength_by = 'weight',
    # enhancements
    # point_include_columns=['glycan_specificity', 'organism', 'pdb_id'] need to add columns
)

project_id = widget.export_project_by_name("sequence_clusters")
print(f"Project exported! View and download PNG here: https://run.cosmograph.app/project/{project_id}")


# create a dummy output file to force snakemake to check if the script runs correctly
seq_cluster.done