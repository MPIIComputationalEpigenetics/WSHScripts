########################################################################################################################
## calculate_FDRP.R
## last update: 2018-02-27
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## This script calculates the FDRP values for all CpG positions in the given annotation with a 
## sufficently large coverage		
########################################################################################################################

#' loading the required packages
suppressPackageStartupMessages(library(RnBeads))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(doParallel))

## G L O B A L S #######################################################################################################

WINDOW.SIZE <- 50

## F U N C T I O N S ###################################################################################################

#' convert
#' This function converts a sequence representation as a string into a binary (TRUE/FALSE)
#' representation according to the methylation state of the CpG dinucleotide at the given position.
#' 
#' @param y:		position of the CpG nucleotide in the string
#' @param string: 	sequence of the read as a string
#' @return:	 TRUE, if the CpG dinucleotide was methylated at the given position and FALSE otherwise
#'
#' @author Michael Scherer

convert <- function(y,string){
	if(nchar(string) >= y+2){
		s <- substr(string,y+1,y+2)
		ret <- ifelse(s == 'CG',TRUE,FALSE)
		ret
	}
}

#' toCpGs
#' This function converts the sequence representation as the result of a sequencing technology to a 
#' binary representation with TRUE values for the methylated CpGs and FALSE values for the remaining 
#' CpGs
#' 
#' @param index:			index of the read that should be converted to the custom representation
#' @param match_read_cpg:	list containing the mapping from the reads to the CpG sites covered by 
#'							the read
#' @param starts_cpgs:		vector of positions for the CpG sites
#' @param starts_reads:		starting positions of the reads
#' @param seqs_reads:		raw sequence of the read
#' @return:				representation of the read as a logical vector for each CpG site in the read
#'						where TRUE represents a methylated CpG and FALSE a non-methylated (or a sequence
#'						variation)
#'@details:				This function calls the function convert for each CpG site covered by the read.
#'						Sequence variations and unmethylated CpG sites are treated in the same way.
#'
#' @author Michael Scherer 

toCpGs <- function(index,match_read_cpg,starts_cpgs,starts_reads,seqs_reads){
	covered_cpgs <- match_read_cpg[[index]]
	start_cpgs <- starts_cpgs[covered_cpgs]
	names(start_cpgs) <- start_cpgs
	start_of_read <- starts_reads[index]
	start_cpgs <- start_cpgs-start_of_read
	sequence <- seqs_reads[index]
	representation <- lapply(start_cpgs,convert,sequence)
	representation <- unlist(representation)
	representation
}

#' restrict
#' This function restricts a read to those CpG sites that are at most WINDOW.SIZE away from one
#' another.
#'
#' @param positions:	positions of the CpG sites in the read of interest
#' @param cpg:			position of the CpG site of interest
#' @return:				positions of the CpG sites that are at most WINDOW.SIZE away from one
#'						another
#'
#' @author Michael Scherer

restrict <- function(positions,cpg){
  if(!any(positions == cpg)) return(NA)
  distances <- abs(as.numeric(positions)-as.numeric(cpg))
  positions <- positions[distances<=WINDOW.SIZE]
  all.valid.intervals <- sapply(min(as.numeric(positions)):max(as.numeric(positions)),function(x){
    interval <- x:(x+(WINDOW.SIZE-1))
    if(max(interval)>max(positions)) return(NULL)
    if(!(cpg%in%interval)) return(NULL)
    c(x,x+(WINDOW.SIZE-1),sum(interval%in%positions))
  })
  all.valid.intervals <- all.valid.intervals[sapply(all.valid.intervals,function(x)!is.null(x))]
  if(length(all.valid.intervals)==0) return(positions)
  interval.info <- t(as.data.frame(all.valid.intervals))
  sel.interval <- interval.info[which.max(interval.info[,3]),]
  positions <- positions[as.numeric(positions)>=sel.interval[1]&as.numeric(positions)<=sel.interval[2]]
  return(positions)
}

#' compute.discordant
#' This function decides for a given read pair if it is discordant or not.
#'
#' @param index:	index of the read pair
#' @param read1:	vector of indeces for the query reads for the read pairs
#' @param read2:	vector of indeces for the subject reads for the read pairs
#' @param values:	GRanges object containing the reads from which the read pairs
#'					were calculated
#' @param site:		position of the CpG site of interest
#' @return:			classification of the read pair as discordant or concordant
#'
#' @details:	The given read pair is only concordant, if all CpG sites that overlap between the
#' 				reads reflect the same methylation status. Otherwise the read pair is classified
#' 				as discordant.
#'
#' @author Michael Scherer

compute.discordant <- function(index,read1,read2,values,site){
	#' this is the first read
	q <- read1[index]
	#' this is the second
	s <- read2[index]
	#' the values are stored in the corresponding part of the GRanges object
	v1 <- values[q]
	#' we check if we have any information about the mehtylation state available
	if(length(v1)>0){
		v1 <- v1[[1]]
	}else{
		return(NA)
	}
	v2 <- values[s]
	if(length(v2)>0){
		v2 <- v2[[1]]
	}else{
		return(NA)
	}
	#' now we match the corresponding positions to each other
	names1 <- names(v1)
	names2 <- names(v2)
	both <- intersect(names1,names2)
	both <- restrict(both,site)
    if(length(both)<2){
        return(NA)
    }
	if(any(is.na(both))){
		return(NA)
	}
	intersection1 <- v2[both]
	intersection1 <- intersection1[!is.na(intersection1)]
	intersection2 <- v1[both]
	intersection2 <- intersection2[!is.na(intersection2)]
	#' two reads are only concordant to each other if they reflect the same methylation
	#' state at each position covered by the reads
	discordant <- any(intersection1 != intersection2)
	discordant
}
		
#' compute.fdrp
#' This function compute the FDRP score for a given cpg site.
#'
#' @param pos:		index of the cpg site
#' @param cpg:		list containing the mapping from each cpg site to the reads that
#'					contain this site
#' @param reads:	GRanges object containing the reads needed for the calculation
#'					of the FDRP, already converted into the custom representation	
#' @param site:		position of the CpG site of interest as an integer number
#' @return:			FDRP score for the given CpG site
#'
#' @details:	This function is called by calculate.fdrps.chromosome each CpG present on
#' 				the chormosome and calls compute.discordant for each pair of reads.
#'
#' @author Michael Scherer

compute.fdrp <- function(pos,cpg,reads,site){
	cpg <- cpg[[pos]]
	site <- site[[pos]]
	#' we only consider the calculation if there are more than 2 reads for this CpG site
	if(length(cpg)<3){
		return(NA)
	}
	#' the maximum number of reads to be calculated is 40, which already are
	#' (1/2)*39*40 = 780 read pairs for a single site. Otherwise 40 reads are 
	#' are sampled from the all reads.
	if(length(cpg)>40){
		cpg <- sample(cpg,40)
	}
	selected <- reads[cpg]
	rm(reads)
	
	#' here we calculate all pairs of the reads that the corresponding CpG site covers
	#' we require the overlap to be at least 35 bp long
	overlap <- findOverlaps(selected,selected,minoverlap=35,ignore.strand=TRUE)
	query <- queryHits(overlap)
	if(length(query)==0){
		return(NA)
	}
	subject <- subjectHits(overlap)
	#' we only consider each of the pairs once and the findOverlaps function calculates
	#' all pairs in both directions
	smaller <- query<subject
	query <- query[smaller]
	subject <- subject[smaller]

	#' now we start with the calculation of the FDRP for the sample
	values <- values(selected)[,"CpG"]
	rm(selected)
	#' query is as long as all read pairs we want to consider at this specific position
	ret <- as.list(1:length(query))
	ret <- lapply(ret,compute.discordant,query,subject,values,site)
	ret <- unlist(ret)
	#' we acutally calculate the FDRP as \frac{#discordant read pairs}{#all read pairs}
	fdrp <- sum(ret,na.rm=TRUE)/length(ret)
	fdrp
}

#' calculate.fdrps.chromosome 
#' This function calculates the FDRP scores for the reads given in the bam files
#' for the CpGs present in the annotation.
#'
#' @param bam:	bam file with the reads from a bisulfite sequencing technique
#'				already aligned to a reference genome
#' @param anno:	annotation as a GRanges object with the CpG sites to be analyzed
#' @return:		vector of fdrp scores for the given CpG sites in anno
#'
#' @details:	This function is called by calculate.fdrps for each chromosome
#'				separately and calls toCpG as well as compute.fdrp.
#'
#' @author Michael Scherer

calculate.fdrps.chromosome <- function(bam, anno){
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
	#' This part converts the raw sequencing reads from the alignment into
	#' an representation, where only CpG positions are considered and from whom we 
	#' can infer discordance or concordance of reads
	logger.start(paste('Representation',chromosome))
	range_cpgs <- ranges(anno)
	starts_cpgs <- start(range_cpgs)
	rm(range_cpgs)
	seqs_reads <- as.character(values(range_reads)$seq)
	starts_reads <- start(ranges(range_reads))
	overlap <- findOverlaps(range_reads,anno,ignore.strand=TRUE)
	match_read_cpg <- as(overlap,"list")
	overlap <- findOverlaps(anno,range_reads,ignore.strand=TRUE)
	fdrps <- as.list(rep(NA,length(anno)))
	rm(anno)
	#' for each read we convert the covered CpG sites into a custom representation
	read_representation <- as.list(1:length(range_reads))
	read_representation <- lapply(read_representation,toCpGs,match_read_cpg,starts_cpgs,starts_reads,seqs_reads)
	rm(match_read_cpg)
	rm(starts_reads)
	rm(seqs_reads)
	values(range_reads) <- DataFrame(cbind(CpG=read_representation))
	rm(read_representation)
	logger.completed()
	
	######### FDRP CALCULATION ###########
	#' we only calculate the FDRP for the CpGs that are acutally covered by any read in the
	#' corresponding sample
	logger.start(paste('FDRP',chromosome))
	match_cpg_reads <- as(overlap,"list")
	rm(overlap)
	null <- lapply(match_cpg_reads,function(x){length(x)>0})
	null <- unlist(null)
	match_cpg_reads <- match_cpg_reads[null]
	starts_cpgs <- starts_cpgs[null]
	toApply <- 1:length(starts_cpgs)
	fdrps_actual <- lapply(toApply,compute.fdrp,match_cpg_reads,range_reads,starts_cpgs)
	#' when we do not have a read that covers this site, we set the FDRP for this site to NA
	fdrps[null] <- fdrps_actual
	fdrps <- unlist(fdrps)
	#' here we set the name for the corresponding position in the genome
	logger.completed()
	logger.completed()
	fdrps
}

#' split.anno
#' This function splits an annotation as a GRanges object into a GRangesList 
#' with two entries.
#'
#' @param anno:	annotation as a GRanges object to be split into two annotations
#' @return:		GRangesList object with two more or less equally sized parts
#'
#' @details:	the first with lower starting points as the half of the
#' maximum starting number and the second with the remaining. In the best case,
#' the two annotations are equally sized in the end.
#'
#' @author Michael Scherer

split.anno <- function(anno){
	start <- start(ranges(anno)[1])
	end <- end(ranges(anno)[length(anno)])
	end <- round(end/2,0)
	part1 <- anno[end(ranges(anno))<=end]
	part2 <- anno[start(ranges(anno))>end]
	return(GRangesList(part1,part2))
}

#' calculate.fdrps
#' This function calculates the FDRPs for all CpG sites in 
#' the annotation from the reads provided by the bam file.
#'
#' @param bam_file: 	bam-file of the sample to be analyzed
#' @param anno: 		CpG sites to be analyzed as a GRanges object
#' @param path:			path to a folder where the FDRP CSV file should be written out
#' @param output_name:	name of the output file
#' @param cores:		number of cores available for the analysis
#' @param window.size:	window size used to restrict the concordance/discordance classification of each read pair
#' 						DEFAULT: 50 as the maximum distance
#'
#' This function writes the FDRP scores into the given directory as a csv-file.
#'
#' @author Michael Scherer

calculate.fdrps<- function(bam_file,anno,path,output_name,cores=1,window.size){
	if(!missing(window.size)){
		WINDOW.SIZE <<- window.size
	}
	bam <- BamFile(bam_file)
	if(!file.exists(file.path(path,'log'))){
		dir.create(file.path(path,'log'))
	}
	cl <- makeCluster(cores,outfile=file.path(path,'log',paste0('log_',output_name,'.log')))
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
	fdrps <- foreach(chromosome=anno,.combine='c',.packages=c('RnBeads','GenomicAlignments','Rsamtools','rtracklayer'),.export=c('calculate.fdrps.chromosome','compute.fdrp','compute.discordant','toCpGs','bam','convert','restrict','WINDOW.SIZE')) %dopar%{
		calculate.fdrps.chromosome(bam,chromosome)
	}
	stopCluster(cl)
	fdrps <- unlist(fdrps)
	write.csv(fdrps,paste0(path,output_name,".csv"),row.names=FALSE)
}	
