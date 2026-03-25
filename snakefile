configfile: "config.yaml"

rule all:
    input:
        "results/seq_clustering.tsv"

rule sequence_clustering:
    input:
        query = "data/query.fasta"
    output:
        clusters = "results/clusters.tsv"
    conda:
        "envs/mmseqs2.yaml"
    script:
        """
        # Create a temporary directory for MMseqs2 calculations
        mkdir -p results/tmp_search

        # Execute search using parameters from the paper
        mmseqs easy-search {input.query} {input.query} {output.results} results/tmp_search \
            --threads {threads} \
            -e 1e-4 \
            -c 0.5 \
            --cov-mode 0

        # Purge temporary files to maintain data sanctity
        rm -rf results/tmp_search
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
