# Comparative_snATAC_snRNA
This directory contains analysis scripts for comparative snATAC-seq and snRNA-seq experiments performed on human, chimpanzee and rhesus macaque.

The folders are organized as:
\
\
**snATAC_seq: Scripts to preprocess snATAC-seq dataset after peak calling. Scripts are numbered based on their order.**

01: Finds consensus peaks (aka CREs) across samples per species.

02: Finds common CREs across species (uses liftOver on chimpanzee and macaque to find common coordinates with human).

03: Finds overlap between raw reads and final set of CREs in each species own coordinates.

04: Generates CRE-Cell count matrix.

05: Dimensionality reduction and cell clustering per species.

06: Gene activity scores per cell using Cicero.

07: Broad annotation of cell types per species.

08: Subtype annotation of broadly annotated cell types per species. Accompanying snRNA-seq is used as a reference to achieve this.

09: Differentially accessible region analysis between species per cell type.

10: Motif enrichments of CREs with human-specific accessibility per cell type.

11_01: Human accelerated region (HAR) identification on cortical CREs.

11_02_03: Human accelerated region (HAR) enrichments of differentially accessible CREs between species.

12_01: Identify modern human-specific substitutions by comparing to 3 neanderthals and 1 denisovan.

12_02: Enrichment between modern human-specific substitutions and human-specific CRE accessibility changes per cell type.

13: Overlap of differentially accessible CREs and DEGs per species per cell type compared to randomized background.

14: Identification of human-specific substitutions within CREs compared to great ape species.

\
\
\
**snRNA_seq: Scripts to preprocess snRNA-seq dataset after count matrix generation. Scripts are numbered based on their order.**

01: Broad annotation of cell types per species.

02: Subtype annotation of broadly annotated cell types after integrating species to identify homologous cell types.

03: Statistical comparison of cell type ratios between species.

04: Differential gene expression (DEG) analysis between species and identification of species-specifically regulated genes.

05: Gene ontology (GO) enrichment of DEGs between species.

\
\
\
**Supplementary_Functions: Contains functions used in other scripts**
\
\
\
Please do not hesitate contact me with any questions: emre.caglayan@utsouthwestern.edu
