######################################################
## This script calculates the PDR values for all   ##
## CpG positions in the given rnbSet with a         ##
######################################################

#' loading the required packages
suppressPackageStartupMessages(library(RnBeads))
#.libPaths(c(.libPaths(),"/TL/deep/projects/work/mage/lib"))
#' suppressPackageStartupMessages(library(RnBeads.hg38))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(doParallel))

convert <- function(y,string){
	if(nchar(string) >= y+2){
		s <- substr(string,y+1,y+2)
		ret <- ifelse(s == 'CG',TRUE,FALSE)
		ret
	}
}

#This function return TRUE, if the read is discordant, FALSE otherwise
classify.read <- function(index,match_read_cpg,starts_cpgs,starts_reads,seqs_reads){
	covered_cpgs <- match_read_cpg[[index]]
	#' This is one of the criteria implemented in the PDR, each read has to contain at least 4 CpGs to
	#' be used for the classification into discordant/concordant
	if(length(covered_cpgs)<4){
		return(NA)
	}
	start_cpgs <- starts_cpgs[covered_cpgs]
	names(start_cpgs) <- start_cpgs
	start_of_read <- starts_reads[index]
	start_cpgs <- start_cpgs-start_of_read
	sequence <- seqs_reads[index]
	representation <- lapply(start_cpgs,convert,sequence)
	representation <- unlist(representation)
	#' A read is only concordant if all of the CpGs show the same methylation status
	concordant <- (all(representation) || all(!representation))
	return(!concordant)
}
		

#' This function calculated the pdr for a given site
#' cpg is the vector of read indeces which contain the given CpG site for which the calculation 
#' should be performed
compute.pdr <- function(cpg,reads){
	#' we only consider the calculation if the cpg site is covered by more than 10 reads
	#' Another requirement stated in the PDR paper
	if(length(cpg)<10){
		return(NA)
	}

	selected <- reads[cpg]
	rm(reads)
	
	#' We select the reads that contain the given CpG site
	values <- values(selected)[,"isDiscordant"]
	rm(selected)
	#' we calculate the PDR as \frac{#discordant reads}{#all reads} that contain the given CpG
	values <- unlist(values)	
	pdr <- mean(values,na.rm=TRUE)
	pdr
}

calculate.pdr.chromosome <- function(bam, anno){
	chromosome <- as.character(seqnames(anno))[1]
	logger.start(paste('Calculation of',chromosome))
	start <- start(ranges(anno)[1])
	end <- end(ranges(anno)[length(anno)])
	if(!(chromosome %in% names(scanBamHeader(bam)[[1]]))){
		chromosome <- unique(na.omit(as.numeric(unlist(strsplit(chromosome,"[^0-9]+")))))
	}
	which <- paste0(chromosome,":",start,"-",end)
	which <- GRanges(which)
	param <- ScanBamParam(which=which,what="seq",mapqFilter=35,flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
	reads <- readGAlignments(bam,param=param)

	range_reads <- GRanges(reads)
	rm(reads)
	newStyle <- mapSeqlevels(seqlevels(range_reads),'UCSC')
	newStyle <- newStyle[!is.na(newStyle)]
	range_reads <- renameSeqlevels(range_reads,newStyle)

	#' we only analyze those CpGs that are covered (on average) by enough reads in the complete dataset
	overlap <- findOverlaps(range_reads,anno,ignore.strand=TRUE)
	query <- queryHits(overlap)
	query <- unique(query)
	range_reads <- range_reads[query]

	####### REPRESENTATION ############
	#'This part clasifies all reads into either discordant or concordant
	logger.start(paste('Representation',chromosome))
	range_cpgs <- ranges(anno)
	starts_cpgs <- start(range_cpgs)
	rm(range_cpgs)
	seqs_reads <- as.character(values(range_reads)$seq)
	starts_reads <- start(ranges(range_reads))
	overlap <- findOverlaps(range_reads,anno,ignore.strand=TRUE)
	match_read_cpg <- as(overlap,"list")
	overlap <- findOverlaps(anno,range_reads,ignore.strand=TRUE)
	pdrs <- as.list(rep(NA,length(anno)))
	rm(anno)
	#' we classify each read into either discordant or concordant
	classified_reads <- as.list(1:length(range_reads))
	classified_reads <- lapply(classified_reads,classify.read,match_read_cpg,starts_cpgs,starts_reads,seqs_reads)
	rm(match_read_cpg)
	rm(starts_reads)
	rm(seqs_reads)
	rm(starts_cpgs)
	values(range_reads) <- DataFrame(cbind('isDiscordant'=classified_reads))
	rm(classified_reads)
	logger.completed()
	
	######### PDR CALCULATION ###########
	#' Starting from the classification of the reads, we now calculate the actual PDR
	#' values for the CpG sites of interest
	logger.start(paste('PDR',chromosome))
	match_cpg_reads <- as(overlap,"list")
	rm(overlap)
	null <- lapply(match_cpg_reads,function(x){length(x)>0})
	null <- unlist(null)
	match_cpg_reads <- match_cpg_reads[null]
	pdrs_actual <- lapply(match_cpg_reads,compute.pdr,range_reads)
	#' when we do not have a read that covers this site, we set the FDRP for this site to NA
	pdrs[null] <- pdrs_actual
	pdrs <- unlist(pdrs)
	logger.completed()
	logger.completed()
	pdrs
}

split.anno <- function(anno){
	start <- start(ranges(anno)[1])
	end <- end(ranges(anno)[length(anno)])
	end <- round(end/2,0)
	part1 <- anno[end(ranges(anno))<=end]
	part2 <- anno[start(ranges(anno))>end]
	return(GRangesList(part1,part2))
}

###### calculate.pdr ############################
#' This function calculates the PDRs for all analyzed CpG sites in 
#' the bam file of the corresponding sample
#'
#' @param bam_file: 	bam-file of the sample to be analyzed
#' @param anno: 	positional annotation that consists of the CpGs to be analyzed
#' @param path:		path to a folder where the PDR CSV file should be written out
#' @param output_name:	intended name of the output file
#' @param cores:	number of cores available for the analysis
#'
#' This function writes the results into the given working directory as a csv-file

calculate.pdrs <- function(bam_file,anno,path,output_name,cores=1){
	bam <- BamFile(bam_file)
	cl <- makeCluster(cores,outfile=file.path(path,paste0('log_',output_name,'.log')))
	registerDoParallel(cl)
	anno <- split(anno,seqnames(anno))
	anno <- anno[lengths(anno)>0]
	if(length(anno)>=22){
		anno <- anno[1:22]
		first <- split.anno(anno[[1]])
		anno <- c(first,anno[2:22])
		second <- split.anno(anno[[3]])
		anno <- c(anno[1:2],second,anno[4:23])
	}
	pdrs <- foreach(chromosome=anno,.combine='c',.packages=c('RnBeads','GenomicAlignments','Rsamtools','rtracklayer'),.export=c('calculate.pdr.chromosome','compute.pdr','classify.read','bam','convert')) %dopar%{
		calculate.pdr.chromosome(bam,chromosome)
	}
	stopCluster(cl)
	pdrs <- unlist(pdrs)
	write.csv(pdrs,paste0(path,output_name,".csv"),row.names=FALSE)
}	
