#####################################################################################################
#' add_bismark_meth_string.R
#' This script adds a methylation call string, as introduced with bismark
#' (http://www.bioinformatics.babraham.ac.uk/projects/bismark/) to any bam file. Currently, this 
#' scripts operates on the human reference genome hg38, but can be easily extended to other 
#' reference genomes.

#' Load required libraries
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(RnBeads))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

cmds <- commandArgs(trailingOnly=TRUE)
bam.file <- cmds[1]
store.path <- cmds[2]
store.name <- cmds[3]
cores <- as.numeric(cmds[4])

### FUNCTIONS ###
#' create.string
#' This function converts a given sequence into the bisulfite conversion representation. Only if a
#' CG is found at a CpG position, was the corresponding cytosine methylated.
#' @param index	Index of the read to be converted
#' @param lens	Lenghts of the reads
#' @param seqs	Sequences of the reads (as a string)
#' @param posis Positions of the CpGs in the corresponding read 
#' @return	The methylation call string according to the one added by bismark
create.string <- function(index,lens,seqs,posis){
  len <- lens[index]
  ret.string <- rep(".",len)
  seq.string <- seqs[index]
  pos <- posis[[index]]
  for(i in pos){
    if(substr(seq.string,i,i+1)=="CG"){
      ret.string[i] <- "Z"
    }else{
      ret.string[i] <- "z"
    }
  }
  paste0(ret.string,collapse = "")
}

### MAIN ###
bam.file <- BamFile(bam.file)
which.all <- c("1:1-248956422","2:1-242193529","3:1-198295559","4:1-190214555","5:1-181538259","6:1-170805979","7:1-159345973","8:1-145138636","9:1-138394717","10:1-133797422","11:1-135086622","12:1-133275309","13:1-114364328","14:1-107043718","15:1-101991189","16:1-90338345","17:1-83257441","18:1-80373285","19:1-58617616","20:1-64444167","21:1-46709983","22:1-50818468")
return.reads <- c()
if(!file.exists(file.path(store.path,"log"))){
	dir.create(file.path(store.path,"log"))
}
cl <- makeCluster(cores,outfile=file.path(store.path,"log","loginfo.log"))
registerDoParallel(cl)
return.reads <- foreach(which=which.all, .combine="c", .packages=c("RnBeads","GenomicAlignments","rtracklayer","Rsamtools"), .export=c("create.string","bam.file")) %dopar%{
  logger.start(substr(which,1,5))
  which <- GRanges(which)
  param <- ScanBamParam(which=which,what=scanBamWhat())
  reads <- readGAlignments(bam.file,param=param,use.names = T)
  cpgs <- unlist(rnb.get.annotation("CpG",assembly="hg38"))
  op <- findOverlaps(reads,cpgs)
  op <- as.list(op)
  starts.cpgs <- start(cpgs)
  starts.cpgs <- lapply(op,function(x,starts){starts[x]},starts.cpgs)
  starts.reads <- start(reads)
  pos.cpgs <- lapply(1:length(starts.reads),function(x,starts.cpgs,starts.reads){starts.cpgs[[x]]-starts.reads[[x]]+1},starts.cpgs,starts.reads)
  lengths <- values(reads)$qwidth
  seqs <- as.character(values(reads)$seq)
  xm.tag.string <- lapply(1:length(starts.reads),create.string,lengths,seqs,pos.cpgs)
  values(reads) <- data.frame(values(reads),"XM"=unlist(xm.tag.string),stringsAsFactors = F)
  logger.completed()
  reads
}
reads <- GAlignmentsList(return.reads)
reads <- unlist(reads)
export(reads,file.path(store.path,store.name))
