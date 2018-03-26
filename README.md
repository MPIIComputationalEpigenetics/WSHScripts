## Introduction

This repository contains all scripts relevant to reproduce the results in **Quantitative comparison of Intra-Sample Heterogeneity Scores for DNA Methylation Data**. These scores are qFDRP, FDRP, PDR, MHL, Epipolymorphism and Entropy. The repository contains the following folders:

### misc

Contains scripts that are used in various occassions by the pipelines, for instance to add a methylation call string (similar to bismark) to a bam file.

## plotting_scripts

Contains R scripts needed to reproduce the plots/analysis from the publication, both for synthetic and biological data.

### real_world_application

Contains the scripts for the calculation on the two RRBS data sets, which are a healthy blood data set (Kiel Cohort) and a Ewing sarcoma cancer example (Ewing).

### scores

Contains the implementation of the scores.

### simulation

Contains the pipelines used to compare the ISH scores on synthetic data.
