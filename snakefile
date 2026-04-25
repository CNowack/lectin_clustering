configfile: "config.yaml"

rule all:
    input:
        "results/map_nodes.csv",
        "results/map_links.csv"

# --- MMseqs2 Representative Consolidation ---
# Cluster at 50% identity to collapse 12.7M to something more managable

rule representative_clustering:
    input:
        query = config["input_fasta"]
    output:
        rep_query = "data/query_rep50_rep_seq.fasta"
    params:
        prefix = "data/query_rep50"
    benchmark:
        "results/benchmarks/representative_clustering.tsv"
    conda:
        "envs/mmseqs2.yaml"
    threads: config["threads"]
    shell:
        """
        mkdir -p data/temp_cluster50
        mmseqs easy-cluster {input.query} {params.prefix} data/temp_cluster50 \
            --min-seq-id 0.5 \
            -c 0.8 \
            --cov-mode 0 \
            --threads {threads}
        rm -rf data/temp_cluster50
        """

# --- Add Labeled Controls ---
# Concat the labeled lectins onto the representative sequences fasta

rule add_controls:
    input:
        rep_query = "data/query_rep50_rep_seq.fasta",
        lectins = "data/lectins.fasta.gz"
    output:
        combined = "data/query_all.fasta"
    shell:
        """
        cat {input.rep_query} <(gunzip -c {input.lectins}) > {output.combined}
        """

# --- MMseq2 Clustering ---
# All-vs-All alignment
rule sequence_clustering:
    input:
        query = "data/query_all.fasta"
    output:
        clusters = "results/all_vs_all_alignment.tsv"
    benchmark:
        "results/benchmarks/sequence_clustering.tsv"
    conda:
        "envs/mmseqs2.yaml"
    threads: config["threads"]
    shell:
        # Create a temporary directory for MMseqs2 calculations
        # Execute search using parameters from the paper
        # purge temp files
        """
        mkdir -p results/tmp_search
        mmseqs easy-search {input.query} {input.query} {output.clusters} results/tmp_search \
            --threads {threads} \
            -e {config[mmseqs_evalue]} \
            -c {config[mmseqs_coverage]} \
            --cov-mode 0
        rm -rf results/tmp_search
        """

# --- Detect Clusters ---
rule community_detection:
    input:
        tsv = "results/all_vs_all_alignment.tsv"
    output:
        nodes = "results/clusters/seq_nodes.csv",
        links = "results/clusters/seq_links.csv"
    benchmark:
        "results/benchmarks/community_detection.tsv"
    conda:
        "envs/afdb_graph.yaml"
    script:
        "scripts/community_detection.py"

# --- Sequence Cosmograph ---
rule seq_cosmograph:
    input:
        nodes = "results/clusters/seq_nodes.csv"
        links = "results/clusters/seq_links.csv"
        project = "sequence_clusters"
    output:
        graph = "results/cosmograph/seq_clusters.done"
    benchmark:
        "results/benchmarks/seq_cosmopgraph.tsv"
    conda:
        "envs/cosmo.yaml"
    script:
        "scripts/seq_cosmograph.py"

# --- Select Lectin Clusters ---
rule select_clusters:
    input:
        nodes = "results/clusters/seq_nodes.csv"
        links = "results/clusters/seq_links.csv"
    output:
        target_seq_clusters = "results/clusters/target_seq_clusters.csv"
    conda:
        "envs/base.yaml"
    script:
        "scripts/select_seq_clusters.py"

# --- FoldTree Structural Clustering ---

# need to downweight the sequence similarity parameter

rule structural_clustering:
    input:
        target_seq_clusters = "results/clusters/target_seq_clusters.csv"
    output:
        data = "results/foldtree/structure_alignment.tsv"
    conda:
        "envs/foldtree.yaml"
    script:
        "scripts/cluster_by_structure.py"
    shell:
        """
        foldtree
        """

# --- Detect Structural Clusters ---
rule detect_structure_clusters:
    input:
        data = "results/foldtree/structure_alignment.csv"
    output:
        nodes = "results/clusters/structure_nodes.csv"
        links = "results/clusters/structure_links.csv"
    benchmark:
        "results/benchmarks/detect_structure_clusters.tsv"
    conda:
        "envs/afdb_graph.yaml"
    script:
        "scripts/structure_clusters.py"


# --- Structure Cosmograph ---
rule structure_cosmograph:
    input:
        nodes = "results/clusters/structure_nodes.csv"
        links = "results/clusters/structure_links.csv"
        project = "structure_clusters"
    output:
        graph = "results/cosmograph/structure_clusters.done"
    conda:
        "envs/cosmo.yaml"
    benchmark:
        "results/benchmarks/structure_cosmograph.tsv"
    script:
        "scripts/structure_cosmograph.py"