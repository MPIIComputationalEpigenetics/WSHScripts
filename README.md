## Introduction

This repository contains all scripts relevant to reproduce the results in Scherer M., et al., *Quantitative comparison of within-sample heterogeneity scores for DNA methylation data* , Nucleic Acids Research, 2020, [10.1093/nar/gkaa120](https://doi.org/10.1093/nar/gkaa120). These scores are qFDRP, FDRP, PDR, MHL, Epipolymorphism and Entropy. The repository contains the following folders:

### [misc](misc/)

Contains scripts that are used in various occassions by the pipelines, for instance to add a methylation call string (similar to bismark) to a bam file.

### [plotting_scripts](plotting_scripts/)

Contains R scripts needed to reproduce the plots/analysis from the publication, both for synthetic and biological data.

### [real_world_applications](real_world_applications/)

Contains the scripts for the calculation on the two RRBS data sets, i.e. the [healthy blood data set](real_world_applications/KielCohort) and the [Ewing sarcoma cancer example](real_world_applications/Ewing).

### [scores](scores/)

Contains the implementation of the scores.

### [simulation](simulation/)

Contains the pipelines used to compare the WSH scores on synthetic data.

## License

Please note that the GNU GPL-3 license applies to all of the code in this repository with the exception of the MHL scripts, which are located in *scores/MHL* and whose license/permission to use is stated in the source files.
