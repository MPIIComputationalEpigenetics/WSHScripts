# Simulation experiments to compare ISH scores

This file should guide you how to use the scripts located in the directory to simulate a specific scenario and then to test the ISH scores FDRP, qFDRP, PDR, Epipolymorphism, Entropy and MHL on the simulated data. First, you should select one of the scenarios you want to simulate. Those are located in *annotation_generator* and briefly described here. 

* create_CT_heterogeneity.R

	Simulates a randomly selected number of cell types (between 2 and 10) and introduces a methylation switch in each of the cell types therefore creating a truly heterogeneous region (THR). With probability 50%, the methylation state does not change, then we call this a negative example. Number of cell types and whether it is a negative control is stored in the file *numCT.txt*, the THR location is stored in *dmr_location.txt* in the cell types folders.

* create_sample_purity.R

	Simulates a contaminating cell type in a population of cells that appears to be pure. The contamination has a different methylation state in a randomly selected region (defined as the THR). This location is stored in *contamination/dmr/dmr_location.txt*; the simulated sample purity (between 50 and 100%) and whether it is a negative control in *purity.txt*.

* create_ASM.R

	Simulates a region exhibiting allele-specific methylation by creating two cell types (here rather alleles), one of them switching the methylation state, while the other remains constant over the region. The location of the ASM region and if its a negative control or not are stored in *asm/asm_location.txt*.

* create_erosion.R

	Simulates DNA methylation erosion by stochastically introducing methylated or non methylated CpGs (by a random parameter **alpha**) in a defined region. This region is then replicated **gamma** times to simulate stochasticity of selecting reads for sequencing. All those parameters, whether the region is a negative control or not, and the location of the THR is stored in the folder.

* create_coverage.R

	Simluates dependecies of the ISH scores on different coverage data sets.

* create_read_length.R

	Investigates possible dependencies of the ISH scores on different read length ranging from 40-150bp.

* create_sequencing_errors.R

	Randomly introduces artificial sequencing errors into the reads to model dependency of the ISH score on this phenomenon. Here, the sequencing error level is stepwise increased from 1 to 10 percent.


After selecting one of the scenarios, you need to create the sample annotation sheet, which is one the inputs to the scripts and also is responsible for (optional) distribution of the individual jobs across a scientific compute cluster. This can be done by employing the corresponding script.

```bash
Rscript your_scenario.R
```

**NOTE:** The RnBeads package (https://www.bioconductor.org/packages/release/bioc/html/RnBeads.html), as well as the annotation package RnBeads.hg38 (https://www.bioconductor.org/packages/release/data/experiment/html/RnBeads.hg38.html), need to be installed in the R installation you are using.

Next, you'll need to add the path to the produced sample annotation sheet in the first line of *simulation_project_config.yaml* as your sample annotation sheet. Please also specify your output folder in this file.

```yaml
metadata:
 sample_annotation: annotation_generator/YOUR_SCENARIO.csv
 output_dir: YOUR_OUTPUT_DIRECTORY
 pipelines_dir: .
```

Since you selected one of the scenarios, you then need to modify the corresponding configuration file (.yaml) in *pipelines*. Here you need to add paths to certain files, for instance to bismark, samtools and methclone (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0472-5). Please install those tools, if you haven't done so yet. 

If you employ your pipeline on a Sun Grid Engine (SGE), you are ready to go. First install looper (http://looper.readthedocs.io) and then just specify the environment variable `LOOPERENV=path_to_compute_config.yaml`. Add the path to looper (normally in `~/.local/bin`) to your PATH variable and then start the pipeline with:

```bash
looper run --compute sge simulation_project_config.yaml
```

Leaving `--compute sge` out, you'll employ the pipelines sequentially on your machine. If you have another scientific compute cluster configuration, please read in the looper documentation (http://looper.readthedocs.io) on how to adjust *compute_config.yaml* to your setting.
