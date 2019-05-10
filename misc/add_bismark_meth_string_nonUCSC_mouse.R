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
which.all <- c("1:1-195471971","2:1-182113224","3:1-160039680","4:1-156508116","5:1-151834684","6:1-149736546","7:1-145441459","8:1-129401213","9:1-124595110","10:1-130694993","11:1-122082543","12:1-120129022","13:1-120421639","14:1-124902244","15:1-104043685","16:1-98207768","17:1-94987271","18:1-90702639","19:1-61431566")
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
  newStyle <- mapSeqlevels(seqlevels(reads),'UCSC')
  newStyle <- newStyle[!is.na(newStyle)]
  reads <- renameSeqlevels(reads,newStyle)
  cpgs <- unlist(rnb.get.annotation("CpG",assembly="mm10"))
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
rm(return.reads)
gc()
reads <- unlist(reads)
export(reads,file.path(store.path,store.name))
