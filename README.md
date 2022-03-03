# Comparative_snATAC_snRNA
This directory contains analysis scripts for comparative snATAC-seq and snRNA-seq experiments performed on human, chimpanzee and rhesus macaque.

The folders are organized as:

**snATAC_seq: Scripts to preprocess snATAC-seq dataset after peak calling. Scripts are numbered based on their order.
**

01: Finds consensus peaks (aka CREs) across samples per species.

02: Finds common CREs across species (uses liftOver on chimpanzee and macaque to find common coordinates with human).

03: Finds overlap between raw reads and final set of CREs in each species own coordinates.

04: Generates CRE-Cell count matrix.

05: Dimensionality reduction and cell clustering per species.

06: Gene activity scores per cell using Cicero.

07: Broad annotation of cell types per species.

08: Subtype annotation of broadly annotated cell types per species. Accompanying snRNA-seq is used as a reference to achieve this.

09: Differentially accessible region (DAR) analysis between species per cell type.

10: Motif enrichments of DAR and other CREs of interest per cell type per species.

11: Human accelerated region (HAR) enrichments of differentially accessible CREs between species.

12: Ancient human variant enrichments of differentially accessible CREs (DA-CREs) between species and identification of species-specifically regulated CREs.

13: Overlap of DA-CREs and DEGs per species per cell type compared to randomized background.

14: Overlap of CREs in this study with Kozlenkov et al. (https://www.pnas.org/doi/10.1073/pnas.2011884117)

15: Motif enrichment of open-chromatin regions associated with CREs in Kozlenkov et al. (https://www.pnas.org/doi/10.1073/pnas.2011884117)

16: Overlap of CREs in this study with activity regulated CREs in Boulting et al. (https://www.nature.com/articles/s41593-020-00786-1)


**snRNA_seq: Scripts to preprocess snRNA-seq dataset after count matrix generation. Scripts are numbered based on their order.
**

01: Broad annotation of cell types per species.

02: Subtype annotation of broadly annotated cell types after integrating species to identify homologous cell types.

03: Differential gene expression (DEG) analysis between species and identification of species-specifically regulated genes.

04: Gene ontology (GO) enrichment of DEGs between species.

05: Enrichment of DEGs with genes dysregulated in disease based on bulk transcriptomic comparisons (disease dataset: https://www.science.org/doi/10.1126/science.aat8127)

06: Enrichment of DEGs with genes dysregulated in disease based on cell type specific transcriptomic comparisons.


**Supplementary_Functions: Contains functions used in other scripts
**

Please do not hesitate contact me with any questions: emre.caglayan@utsouthwestern.edu
