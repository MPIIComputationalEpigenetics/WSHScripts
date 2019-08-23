#####################################################################################################
#' create_rnbSet.R
#----------------------------------------------------------------------------------------------------
#' This scripts creates an RnBSet object from the given cov files for the simualtion experiment, using
#' RnBeads' default parameters and import functions. Furthermore, the region of interest (ROI) for 
#' MHL calculation is inferred and exported from the RnBSet.
suppressPackageStartupMessages(library(RnBeads))
#' the only input is the folder for each of the simulation runs, which contains a (minimal) sample 
#' annotation sheet and the cov-files produced by bismark
cmds <- commandArgs(trailingOnly=TRUE)
folder <- cmds[1]

data.source <- c(file.path(folder,'covs'),file.path(folder,'sample_annotation.csv'),'filename')
rnb.options(identifiers.column='sample_id',
	import.bed.style='bismarkCov',
	assembly='hg38',
	disk.dump.bigff=F)
rnb.set <- rnb.execute.import(data.source=data.source,data.type='bs.bed.dir')
save.rnb.set(rnb.set,file.path(folder,'rnbSet'))
anno <- annotation(rnb.set)
write.table(anno[,1:3],file.path(folder,"roi.bed"),sep="\t",row.names=FALSE,quote=FALSE)
