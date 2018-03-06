#######################################################################################################
#' compute_epipolymorphism.R
#' This scripts converts the output of methclone to genomic bins and calculates Epipolymorphism
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

#' This function calculates epipolymorphism from a single genomic region (line in the methclone output)
calculate.epipoly <- function(line){
	percentages <- as.numeric(line[15:30])
	epipoly <- 1-(sum((percentages/100)^2))
	epipoly <- round(epipoly,4)
	epipoly
}

logger.start('Reading in data')
methclone.data <- read.csv(methclone.data,sep='\t')
logger.completed()
output.frame <- data.frame(chromosome=methclone.data$chr,start=methclone.data$start,end=methclone.data$end,strand=methclone.data$strand)
logger.start('Calculating Epipolymorphism')
values <- apply(methclone.data,1,calculate.epipoly)
logger.completed()
output.frame <- data.frame(output.frame,Epipolymorphism=values)
logger.start('Writing Result')
write.csv(output.frame,file.path(output.folder,paste0(output.name,'.csv')),quote=FALSE,row.names=FALSE)
logger.completed()
