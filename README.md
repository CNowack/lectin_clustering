# lectin_clustering
A project to identify novel clusters of lectin binding proteins through structural homology clustering.

### Current Run Command
```bash
snakemake -c 4 --use-conda
```

## Background & Motivation
Evolutionary relationships between proteins are traditionally inferred from sequence similarity. However, protein structure is more conserved than sequence over long evolutionary timescales—two proteins can share virtually no sequence identity yet retain nearly identical folds. Sequence-based methods miss homologous relationships in highly diverged protein families. We address this with FoldTree, a structural phylogenetics pipeline that replaces sequence alignments with structure-based distances (via Foldseek) to build evolutionary trees.

We are going to implement this methodology to locate and understand candidate glycan-binding proteins (lectins) more broadly. While AB5 toxins are a known test case, many other lectin families (such as the OB-fold) share highly conserved overall folds despite low pairwise sequence identity. This makes them an ideal test case for FoldTree because:
* Known models of lectin folds provide a structural framework for identifying novel systems that share the same overall architecture.
* Identifying these proteins is critical for understanding host-pathogen interactions and how proteins "read" the glycome.
* We can use this pipeline to discover entirely new lectins that sequence-based tools simply cannot see.

## Objectives
1. Use a structural homology search to identify similar proteins in the AlphaFold database (starting from known glycan-binding fold families like the OB-fold).
2. Extend the analysis by clustering the newly identified proteins according to their proposed binding specificities.
3. Once we have those clusters, infer biological understanding of these sites using genus and tissue tropism.

## Why This Is Interesting Regardless of Outcome
* If the binding site is structurally conserved: This provides evidence that it is a functionally meaningful, ancestrally conserved feature—potentially relevant to how these proteins recognize specific sugars across different species.
* If it is not conserved: That is equally informative—it would suggest the binding specificity is a lineage-specific adaptation or driven by convergent evolution, demonstrating the practical boundaries of structure-based subdomain analysis.

# Execution:
This pipeline runs using snakemake to connect the work of Durairaj et al. to the novel discovery of new proteins clusters based on a known protein with known function.

## Protein Test Set for Pipeline Validation
Sourced from UniProtKB, filtering by lectin famliy and review status: 
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(protein_name:lectin)+AND+(reviewed:true)"

- Size and composition:
    840 Lectin family proteins

# Acknowledgements and References
All base work was reproduced from Durairaj et al. 2023.

"Uncovering new families and folds in the natural protein universe"

    Nature (2023). DOI: https://doi.org/10.1038/s41586-023-06622-3

```bibtex
@article{Durairaj2023,
  title={Uncovering new families and folds in the natural protein universe},
  author={Durairaj, Janani and Waterhouse, Andrew M. and Mets, Toomas and others},
  journal={Nature},
  year={2023},
  publisher={Nature Publishing Group},
  doi={10.1038/s41586-023-06622-3},
  url={[https://doi.org/10.1038/s41586-023-06622-3](https://doi.org/10.1038/s41586-023-06622-3)}
}
```

### Author's Code:
* https://github.com/ProteinUniverseAtlas/dbuilder
* https://github.com/TurtleTools/geometricus/tree/master/training
* https://github.com/ProteinUniverseAtlas/AFDB90v4

    The external code provided by the authors are managed via Git submodules in `external/`.
