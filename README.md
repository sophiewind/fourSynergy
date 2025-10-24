# fourSynergy: Ensemble based 4C-seq analysis

Circular Chromosome Conformation Capture Sequencing (4C-seq) is a sequencing technique enabling the identification of chromatin interactions. fourSynergy allows to identify interactions using ensemble based methods. It leverages syngeries between the existing 4C-seq analysis tools r3Cseq (https://doi.org/10.1093/nar/gkt373), fourSig (https://doi.org/10.1093/nar/gku156), peakC (https://doi.org/10.1093/nar/gky443) and R.4Cker (https://doi.org/10.1371/journal.pcbi.1004780).

# Requirements

To run fourSynergy, you need R Version 4.5.1 or higher.

# Installation

fourSynergy is available at Bioconductor. The latest version can be installed via:

```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("fourSynergy")
```
The latest version available at github can be installed via:

```
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::install_github("sophiewind/fourSynergy")
```

# Key features

fourSyngery offers a variety of analyses:

- Base tool interaction results
- Visualization of base tool interactions
- Ensemble interaction calling
- Visualization of ensemble interactions
- Differential interaction calling
- Visualization of differential interactions
- Highlighting of genes of interest

# Detailed documentation

A detailed documentation can be found in the vignette available within this repository.

For questions, feature requests and errors, so not hestitate to contact Sophie Wind (sophie.wind@uni-muenster.de).
