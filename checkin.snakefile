configfile: "config.yaml"

rule all:
    input:
        "results/network_nodes.csv",
        "results/network_links.csv",
        "results/network_nodes_labeled.csv",
        "results/scores.txt"

# --- Download Test Set ---
rule download_lectin_fasta:
    output:
        lectins = "data/lectins.fasta"
    params:
        limit = "200"
    shell:
        """
        mkdir -p data
        curl -k -o {output.lectins} "https://rest.uniprot.org/uniprotkb/search?format=fasta&query=(protein_name:lectin)+AND+(reviewed:true)&size={params.limit}"
        """

# Download 3 other families, limit to 200 of each class
# 1. Immunoglobulins, Should be closest, some bind carbohydrates
# 2. Kinases, Should not be similar
# 3. Ribosomal, Should be very disimilar 
rule download_non_lectin_fasta:
    output:
        ig = "data/ig.fasta",
        kinases = "data/kinases.fasta",
        ribosomal = "data/ribosomal.fasta"
    params:
        limit = "200"
    shell:
        """
        mkdir -p data
        curl -k -o {output.ig} "https://rest.uniprot.org/uniprotkb/search?format=fasta&query=(protein_name:immunoglobulin)+AND+(reviewed:true)&size={params.limit}"
        curl -k -o {output.kinases} "https://rest.uniprot.org/uniprotkb/search?format=fasta&query=(protein_name:kinase)+AND+(reviewed:true)&size={params.limit}"
        curl -k -o {output.ribosomal} "https://rest.uniprot.org/uniprotkb/search?format=fasta&query=(protein_name:ribosomal)+AND+(reviewed:true)&size={params.limit}"
        """

rule concat_fastas:
    input:
        lectins = "data/lectins.fasta",
        ig = "data/ig.fasta",
        kinases = "data/kinases.fasta",
        ribosomal = "data/ribosomal.fasta"
    output:
        combined = "data/combined.fasta"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        mkdir -p data
        cat {input.lectins} {input.ig} {input.kinases} {input.ribosomal} \
            | seqkit rmdup -n > {output.combined}
        echo "Unique sequences: $(grep -c '^>' {output.combined})"
        """

rule sequence_clustering:
    input:
        query = "data/combined.fasta"
    output:
        clusters = "results/clusters.tsv"
    conda:
        "envs/mmseqs2.yaml"
    threads: 4
    resources:
        # mem_mb=16000, # request 16GB RAM
        # runtime="04:00:00" # max runtime of 4 hours
    shell:
        """
        # Create a temporary directory for MMseqs2 calculations
        mkdir -p results/tmp_search

        # Execute search using parameters from the paper
        mmseqs easy-search {input.query} {input.query} {output.clusters} results/tmp_search \
            --threads {threads} \
            -e 1e-4 \
            -c 0.5 \
            --cov-mode 0

        # purge temp files
        rm -rf results/tmp_search
        """

# --- Protein Community Detection ---
rule detect_communities:
    input:
        tsv = "results/clusters.tsv"
    output:
        nodes = "results/network_nodes.csv",
        links = "results/network_links.csv"
    conda:
        "envs/network.yaml"
    script:
        "scripts/community_detection.py"

rule color_by_class:
    # To ensure its part of the DAG
    input:
        nodes = "results/network_nodes.csv",
    output:
        nodes = "results/network_nodes_labeled.csv"
    conda:
        "envs/base.yaml"
    shell:
        """
        python /projectnb/ds596/students/dsch28/lectin_clustering/data/to-label-csv.py {input.nodes} {output.nodes}
        """

rule score_clusters:
    # To ensure its part of the DAG
    input:
        nodes = "results/network_nodes_labeled.csv",
        links = "results/network_links.csv",
    output:
        file = "results/scores.txt",
        nodes = "results/filtered_nodes.csv",
        link = "results/filtered_links.csv",
    conda:
        "envs/base.yaml"
    shell:
        """
        python /projectnb/ds596/students/dsch28/lectin_clustering/results/score-cluster.py {input.nodes} {input.links} {output.file} {output.nodes} {output.link}
        """



## Old project structure 
### -- Identify structural homology search method (Foldseek/FoldTree).
### -- Install and run corresponding software package.
### -- Build initial version of the pipeline.
### -- Run pipeline on a small subset of known lectins (positive controls).
### -- Validate.
### -- Define full scope of dataset (targeting specific folds like the OB-fold).
### -- Build structural down-sampling method.
### -- Run full pipeline.
### -- Cluster by binding domain similarity.
### -- Attempt to classify by tissue tropism and evolutionary relationship.

## New project structure
### -- Cluster UniRef50 using MMseqs2
### -- Cluster Lectins using FoldSeek
### -- Bridge neworks to indentify novel sequence clusters of lectins from structural homolog clusters
### -- Attmept to assign tissure tropism
