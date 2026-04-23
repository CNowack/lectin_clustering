configfile: "config.yaml"

rule all:
    input:
        "results/map_nodes.csv",
        "results/map_links.csv"

# --- MMseq2 Clustering ---
rule sequence_clustering:
    input:
        query = config["input_fasta"]
    output:
        clusters = "results/clusters.tsv"
    benchmark:
        "results/benchmarks/sequence_clustering.tsv"
    conda:
        "envs/mmseqs2.yaml"
    threads: config["threads"]
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
            -e {config[mmseqs_evalue]} \
            -c {config[mmseqs_coverage]} \
            --cov-mode 0
        rm -rf results/tmp_search
        """

rule community_detection:
    input:
        tsv = "results/clusters.tsv"
    output:
        nodes = "results/map_nodes.csv",
        links = "results/map_links.csv"
    benchmark:
        "results/benchmarks/community_detection.tsv"
    conda:
        "envs/afdb_graph.yaml"
    script:
        "scripts/community_detection.py"

# fake the metadata file normally generated from the MongoDB
rule generate_local_metadata:
    input: "data/query.fasta"
    output: "data/protein_metadata.csv"
    shell: "python scripts/generate_metadata.py {input} {output}"
