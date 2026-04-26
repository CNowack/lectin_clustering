# lectin_clustering
A project to identify novel clusters of lectin binding proteins through structural outlier detection.

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
1. Cluster large sequence space using MMseqs2 and loose thresholds
2. Identify cluster containing known lectin family proteins and gather surrounding superset
3. Cluster superset by FoldTree
4. Attempt to identify novel lectin or glycan binding proteins co-localized with known lectins
5. Attempt to assign tropism to structural clusters

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

## Full sequence dataset
**All sequences sourced from UniProtKB using** `download_by_order.py` driven by `batch_download_2.py` to split the 12,735,307 protein sequences within the Proteobacteria class Gammaproteobacteria into groups by taxonomic order to facilitate downloading.
### Parameters
* `REVIEWED = False` - include unreviewed TrEMBL entries
* `TAXONOMY = [1224]` - NCBI taxonomy filter for Proteobacteria phylum
* `LENGTH_RANGE = [100, 600]` - Restrict AA length
* `REQUIRE_AFDB = True ` - Make sure each entry has a Alphafold-predicted structure for FoldTree

All *dark* sequences were download and compiled into a single database `fasta.gz` by `tag_concat.py` that labels each order in the fasta header. This step was repeated for both categories of the lectin sequences. Once the *dark* sequences are collapsed into representative seqeunces (`rule representative_clustering`) the lectins are concatenated to the output (`rule add_controls`) to preserve all controls. The resulting sequence set is passed to the MMseqs2 sequence clustering step.

840         - Reviewed Lectins (lectin_rev)

190563      - Not Reviewed Lectins (lectin_nr)

12543904    - Unreveiwed *"dark proteins"* (dark_<class prefix>)

    aeromonadales:         "dark_aero"
    alteromonadales:       "dark_altero"
    chromatiales:          "dark_chroma"
    enterobacterales:      "dark_entero"
    legionellales:         "dark_legion"
    oceanospirillales:     "dark_oceano"
    pasteurellales:        "dark_pasteur"
    pseudomonadales:       "dark_pseudo"
    thiotrichales:         "dark_thio"
    vibrionales:           "dark_vibrio"
    xanthomonadales:       "dark_xantho"

    Wrote 12,735,307 sequences to data/query.fasta.gz
        size on disk: 2528.6 MB

    Per-file counts:
           840  data/lectins.fasta
       190,563  data/lectins_nr.fasta
       217,155  data/by_order/aeromonadales.fasta.gz
       788,265  data/by_order/alteromonadales.fasta.gz
       285,426  data/by_order/chromatiales.fasta.gz
     5,053,590  data/by_order/enterobacterales.fasta.gz
       185,603  data/by_order/legionellales.fasta.gz
       671,808  data/by_order/oceanospirillales.fasta.gz
       295,444  data/by_order/pasteurellales.fasta.gz
     3,091,323  data/by_order/pseudomonadales.fasta.gz
       161,711  data/by_order/thiotrichales.fasta.gz
       966,591  data/by_order/vibrionales.fasta.gz
       826,988  data/by_order/xanthomonadales.fasta.gz

Duplicate Lectins removed: 

**= 12,542,692  - Total Gammaproteobacteria sequences**

### Future Dataset Expansions (all of Proteobacteria)

Expand to whole Proteobacteria phylum (26,181,537 protein sequences), NCBI ID: 1224

    "data/proteobacteria/acidithiobacillia.fasta.gz":           "dark_p-acidi",
    "data/proteobacteria/alphaproteobacteria.fasta.gz":         "dark_p-alpha",
    "data/proteobacteria/betaproteobacteria.fasta.gz":          "dark_p-beta",
    "data/proteobacteria/deltaproteobacteria.fasta.gz":         "dark_p-delta",
    "data/proteobacteria/environmental.fasta.gz":               "dark_p-env",
    "data/proteobacteria/epsilonproteobacteria.fasta.gz":       "dark_p-epsilon",
    "data/proteobacteria/hydrogenophilalia.fasta.gz":           "dark_p-hydro",
    "data/proteobacteria/incertae_sedis.fasta.gz":              "dark_p-incertae",
    "data/proteobacteria/lambdaproteobacteria.fasta.gz":        "dark_p-lambda",
    "data/proteobacteria/magnetococcia.fasta.gz":               "dark_p-magneto",
    "data/proteobacteria/muproteobacteria.fasta.gz":            "dark_p-mu",
    "data/proteobacteria/oligoflexia.fasta.gz":                 "dark_p-oligo",
    "data/proteobacteria/unclassified.fasta.gz":                "dark_p-unclass",
    "data/proteobacteria/zetaproteobacteria.fasta.gz":          "dark_p-zeta"

Wrote 13,116,527 sequences to data/proteo_query.fasta.gz
  size on disk: 2659.0 MB

Per-file counts:
        48,590  data/proteobacteria/acidithiobacillia.fasta.gz
     7,802,211  data/proteobacteria/alphaproteobacteria.fasta.gz
     4,432,690  data/proteobacteria/betaproteobacteria.fasta.gz
       164,370  data/proteobacteria/deltaproteobacteria.fasta.gz
         1,462  data/proteobacteria/environmental.fasta.gz
       621,554  data/proteobacteria/epsilonproteobacteria.fasta.gz
         6,172  data/proteobacteria/hydrogenophilalia.fasta.gz
            68  data/proteobacteria/incertae_sedis.fasta.gz
         4,343  data/proteobacteria/lambdaproteobacteria.fasta.gz
         8,598  data/proteobacteria/magnetococcia.fasta.gz
        11,656  data/proteobacteria/muproteobacteria.fasta.gz
             1  data/proteobacteria/oligoflexia.fasta.gz
         4,755  data/proteobacteria/unclassified.fasta.gz
        10,057  data/proteobacteria/zetaproteobacteria.fasta.gz

Duplicate Lectins removed: 4,703

Total Preoteobacteria: 13,111,824

Concatenated Proteobacteria and Gammaprotebacteria: 25,654,516

added 3 controls

Grand total: input = 25,654,519

### Dataset Workflow:
1. Download sets of sequeces from UniProt
2. Tag each set of sequences based on their origin (Proteobacteria class or Gammaproteobacteria order) and concatenate into one fasta.gz
  - `python tag_concat.py`
3. Deduplicate (remove unreviewed lectins from dark query)
  - `python deduplicate_lectins.py`
4. **Optional:** Concatenate with other fasta.gz files
  - `cat <1.fasta.gz> <2.fasta.gz> > concatenated.fasta.gz`
5. **Highly Suggested:** Add controls to deduplicated fasta.gz
  - `
  - `gzip -c data/controls.fasta >> data/input_query_dedup.fasta.gz`
6. Deduplicated fasta.gz are snakemake entry point

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
