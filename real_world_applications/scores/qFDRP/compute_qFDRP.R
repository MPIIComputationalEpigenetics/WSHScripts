######################################################################################################
#' compute_qFDRP.R
#' This script is used to preprocess the data as input to calculate_qFDRP.R
#' The script needs the following command line parameters:
#'	[1]	name of the output file
#'	[2]	folder to write out the calculated qFDRP as a CSV file
#'	[3]	input bam file
#'	[4]	path to a RnBSet object to get the annotation in the human genome

# Load the input commands from the command line
cmds <- commandArgs(trailingOnly=TRUE)
name <- cmds[1]
folder <- cmds[2]
bam <- cmds[3]
rnb.path <- cmds[4]
cores <- as.numeric(cmds[5])

# Load the functions needed for the qFDRP calculation
source('calculate_qFDRP.R')

options(warn=-1)
logger.start("Run")
if(!file.exists(file.path(folder,'annotation.RData'))){	
	logger.start("Loading RnBSet")
	#' loading the RnBSet object for which the analysis should be conducted
	rnbSet <- load.rnb.set(rnb.path)
	logger.completed()
	logger.start("Setting coverage threshold")
	annotation <- annotation(rnbSet)	
	coverage <- covg(rnbSet)
	rm(rnbSet)
	#' removing CpGs located on the gonosomes to avoid biases towards one of the genders
	if(all(c("chrX","chrY")%in%annotation$Chromosome)){
		coverage <- as.matrix(coverage[annotation$Chromosome!="chrX",])
		annotation <- annotation[annotation$Chromosome!="chrX",]
		coverage <- as.matrix(coverage[annotation$Chromosome!="chrY",])
		annotation <- annotation[annotation$Chromosome!="chrY",]
	}
	#' calculating mean coverage for each CpG site in the RnBSet object
	mean_coverage <- rowMeans(coverage,na.rm=TRUE)
	rm(coverage)
	#' setting the coverage threshold and only consider CpGs with a sufficiently large coverage in the dataset
	high_enough <- mean_coverage >= 10
	rm(mean_coverage)
	annotation <- annotation[high_enough,]
	rm(high_enough)
	logger.completed()
	logger.start("Selecting subset of CpGs")
	annotation <- GRanges(Rle(annotation$Chromosome),IRanges(start=annotation$Start,end=annotation$End),annotation$Strand)
	if(!file.exists(folder)){
		dir.create(folder)
	}
	save(annotation,file=file.path(folder,'annotation.RData'))
	#' we only consider unique positions
	names <- paste(as.character(seqnames(annotation)),start(ranges(annotation)),end(ranges(annotation)),sep="_")
	unique <- unique(names)
	match <- match(unique,names)
	rm(unique)
	annotation <- annotation[match]
	logger.completed()
}else{
	load(file.path(folder,'annotation.RData'))
}

logger.start(paste("qFDRP calculation for",name))
calculate.qfdrps(bam,annotation,folder,name,cores=cores,window.size=51)
logger.completed()

options(warn=1)
logger.completed()
