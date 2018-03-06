#######################################################################################################
#' compute_entropy.R
#' This scripts converts the output of methclone to genomic bins and calculates methylation entropy
#' for those regions
#' It takes the following input parameters:
#'	[1] output file of methclone in the form of a tab-delimited text file, full path required
#'	[2] name of the final output file
#'	[3] output folder, where the result should be written out

#' Load required libraries
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RnBeads))

# Load the input arguments from the command line
cmds <- commandArgs(trailingOnly=TRUE)
methclone.data <- cmds[1]
output.name <- cmds[2]
output.folder <- cmds[3]

#' This function calculates methylation entropy from a single genomic region (line in the methclone output)
calculate.entropy <- function(line){
	percentages <- as.numeric(line[15:30])
	percentages <- percentages[percentages!=0]
	entropy <- -0.25*(sum((percentages/100)*log2(percentages/100)))
	entropy <- round(entropy,4)
	entropy
}

logger.start('Reading in data')
methclone.data <- read.csv(methclone.data,sep='\t')
logger.completed()
output.frame <- data.frame(chromosome=methclone.data$chr,start=methclone.data$start,end=methclone.data$end,strand=methclone.data$strand)
logger.start('Calculating Methylation Entropy')
values <- apply(methclone.data,1,calculate.entropy)
logger.completed()
output.frame <- data.frame(output.frame,Entropy=values)
logger.start('Writing Result')
write.csv(output.frame,file.path(output.folder,paste0(output.name,'.csv')),quote=FALSE,row.names=FALSE)
logger.completed()
