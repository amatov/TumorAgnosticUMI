# Dependencies

This is analysis code developed on a private institutional HPC cluster
(Aarhus University's GenomeDK). There is no package manifest or
environment lock file -- this document lists what's needed, gathered by
inspecting `library()`/`require()`/`source()`/`import` calls across the
scripts.

## R packages (CRAN / Bioconductor)

Install via `install.packages()` (CRAN) or `BiocManager::install()`
(Bioconductor):

```
abind, caret, devtools, dplyr, FactoMineR, FNN, GenomicRanges (Bioconductor),
ggplot2, ggpubr, glmnet, gplots, imager, matrixStats, mixtools, precrec,
pvclust, qlcMatrix, RColorBrewer, readxl, reshape2, reticulate, ROCit,
ROCR, seewave, spatstat, stringr, tidyr, transport, vioplot, xlsx
```

The README also describes using Bioconductor's **deepSNV** and
**shearwater** packages for variant calling; they are not directly
`library()`-called in any script in this repository, but are part of the
described analysis pipeline upstream of these scripts.

## Python packages

`workflow.py` and `UMI_seq_fragment_length.py` require:

```
gwf        # workflow engine specific to Aarhus University's GenomeDK cluster
pysam
numpy
```

`gwf` in particular assumes you are running on GenomeDK (or a similarly
configured Slurm-based cluster with `gwf` installed) -- `workflow.py`
will not run as a plain Python script outside that environment.

## Sourced files NOT included in this repository

The following files are referenced via `source()` calls in various
scripts but are not present here. Scripts that call them will fail with
a "cannot open file" error unless you supply your own copies:

- `tools.R`
- `sw_piles.R`
- `duplex_tools.R`
- `cmapply.R`
- `image_plot.R`
- `auc.R`, `recoder.R`, `scaler.R`, `confusion_plot.R` (from a separate,
  unpublished personal `utility_functions` collection)

## Sourced files that ARE included, but referenced by a cluster path

`read_bed.R` and `mutationScore.R` are both present as top-level files in
this repository, but several scripts `source()` them via an absolute
GenomeDK cluster path rather than a relative path within this repo. If
you've cloned this repository and want to use the local copies, change
those `source(...)` lines to point at your local clone instead (each is
flagged with an `# EDIT:` comment).
