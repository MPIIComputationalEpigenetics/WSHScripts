######################################################################################################
#' compute_PDR.R
#' This script is used to preprocess the data as input to calculate_PDR.R
#' The script needs the following command line parameters:
#'	[1]	name of the output file
#'	[2]	folder to write out the calculated FDRP as a CSV file
#'	[3]	input bam file
#'	[4]	path to a RnBSet object to get the annotation in the human genome

#' Taken from https://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

# Load the input commands from the command line
cmds <- commandArgs(trailingOnly=TRUE)
name <- cmds[1]
folder <- cmds[2]
bam <- cmds[3]
rnb.path <- cmds[4]
cores <- as.numeric(cmds[5])

#' Load the functions needed for calculation
source(file.path(dirname(thisFile()),"calculate_PDR.R"))

options(warn=-1)
if(!file.exists(file.path(folder,'annotation.RData'))){	
	logger.start("Run")
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
	#' setting the coverage threshold and only consider CpGs with a sufficiently large coverage in the 	dataset
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

logger.start(paste("PDR calculation for",name))
calculate.pdrs(bam,annotation,folder,name,cores=cores)
logger.completed()
