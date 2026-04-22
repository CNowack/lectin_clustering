configfile: "config.yaml"

rule all:
    input:
        "results/map_nodes.csv",
        "results/map_links.csv"

# --- Build Database ---
rule download_database:
    output: 
        directory("data/batches/")
    params: 
        query = config["dataset_query"], prefix = "lectins"
    shell: 
        "python scripts/batch_download.py {params.query:q} {output} {params.prefix}"

# --- MMseq2 Clustering ---
rule sequence_clustering:
    input:
        query = "data/query.fasta"
    output:
        clusters = "results/clusters.tsv"
    conda:
        "envs/mmseqs2.yaml"
    threads: 4
    resources:
        # mem_mb=16000, # request 16GB RAM
        # runtime="04:00:00" # max runtime of 4 hours
    shell:
        # Create a temporary directory for MMseqs2 calculations
        # Execute search using parameters from the paper
        # purge temp files
        """
        mkdir -p results/tmp_search
        mmseqs easy-search {input.query} {input.query} {output.clusters} results/tmp_search \
            --threads {threads} \
            -e 1e-4 \
            -c 0.5 \
            --cov-mode 0
        rm -rf results/tmp_search
        """

# Step 3A: Build the graph and run Asynchronous Label Propagation
rule get_connected_components:
    input:
        edges = "results/all_vs_all_alignment.tsv"
    output:
        components = "results/connected_components.tsv",
        communities = "results/communities.tsv"
    conda:
        "envs/afdb_graph_env.yaml"
    threads: 8
    shell:
        """
        python external/afdb90v4/scripts/get_connected_components.py \
            --input {input.edges} \
            --out_comp {output.components} \
            --out_comm {output.communities} \
            --threads {threads}
        """
# fake the metadata file normally generated from the MongoDB
rule generate_local_metadata:
    input: "data/test_lectins.fasta"
    output: "data/protein_metadata.csv"
    shell: "python scripts/generate_metadata.py {input} {output}"

# Step 3B: Calculate the "Darkness" of each community
rule summarize_communities:
    input:
        communities = "results/communities.tsv",
        metadata = "data/protein_metadata.csv" # See the Accompanying Script below
    output:
        summary = "results/communities_summary.tsv"
    conda:
        "envs/afdb_graph_env.yaml"
    shell:
        """
        python external/afdb90v4/scripts/get_communities_summary.py \
            --communities {input.communities} \
            --metadata {input.metadata} \
            --output {output.summary}
        """

# Step 3C: Format the data for the Cosmograph visualizer
rule make_communities_map:
    input:
        edges = "results/all_vs_all_alignment.tsv",
        communities = "results/communities.tsv",
        summary = "results/communities_summary.tsv"
    output:
        nodes = "results/map_nodes.csv",
        links = "results/map_links.csv"
    conda:
        "envs/afdb_graph_env.yaml"
    shell:
        """
        python external/afdb90v4/scripts/make_communities_map.py \
            --edges {input.edges} \
            --communities {input.communities} \
            --summary {input.summary} \
            --out_nodes {output.nodes} \
            --out_links {output.links}
        """