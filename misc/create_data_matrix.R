#####################################################################################################
#' create_data_matrix.R
#----------------------------------------------------------------------------------------------------
#' This scripts combines the scores for the individual samples of a ISH run to a combined data matrix.
#' This matrix then can be used to further explore the scores, similar to a methylation matrix.
#'
#' parameters to set
#' pipeline_folder: the results_pipeline folder, which was created in the output folder of your pipeline
pipeline.folder <- "."
#' save.path: path to which the results should be saved. Those include the data matrices and 
#' annotations as .RData files that contain the annotations for further exploration
save.path <- "."

#' We need the RnBeads package for these steps, if not done so, please install RnBeads (https://rnbeads.org)
library(RnBeads)
folders <- list.dirs(pipeline.folder,full.names = T,recursive = F)
sample.names <- list.files(pipeline.folder,full.names = F)
sample.names <- sample.names[!grepl("[.]",sample.names)]

#' This script contains three parts:
#' (i):		First we summarize the values for FDRP, qFDRP and PDR, since those scores produce scores
#'     		for the same annotation and therefore also produce matrices of the same dimension. Each
#' 			of the matrices is stored separately as a CSV file
#' (ii): 	Then we summarize the values for MHL, which also produces scores for single CpGs but 
#' 		 	produces a data matrix of a different dimension than the first. The annotation is not the
#'			same for all samples, thus we select the largest file and use this as the base annotation
#' (iii):	Last, we summarize the values for Epipolymorphism and Entropy, which have the following
#'			drawbacks: no single CpG, but for 4 CpGs, thus overlapping annotation is produced. 
#'          Also the annotation is not the same for all samples, thus we select the largest file and
#'          use this as the base annotation

#' (i) Summarize values for FDRP, qFDRP and PDR
measures <- c("FDRP","qFDRP","PDR")
path.extension <- c("FDRP.csv","qFDRP.csv","PDR.csv")
res <- list()
for(i in 1:length(measures)){
  dat <- c()
  for(j in 1:length(folders)){
    file <- file.path(folders[j],measures[i],path.extension[i])
    temp.data <- read.csv(file)
    temp.data <- unlist(temp.data)
    names(temp.data) <- sample.names[j]
    dat <- cbind(dat,temp.data)
    colnames(dat)[j] <- sample.names[j]
  }
  file <- gzfile(file.path(save.path,paste0(measures[i],".csv.gz")),"w")
  write.csv(dat,file)
  close(file)
}
load(file.path(folders[1],"FDRP","annotation.RData"))
save(annotation,file=file.path(save.path,"FDRP_qFDRP_PDR_annotation.RData"))

#' (ii) Summarize MHL scores
max <- 0
which.max <- ""
for(i in 1:length(sample.names)){
  name <- sample.names[i]
  file <- folders[i]
  size <- file.info(file.path(file,"mhl.txt"))$size
  if(size>max){
    max <- size
    which.max <- name
  }
}

complete.data <- c()
i <- 1
name <- sample.names[i]
file <- folders[i]
mhl <- read.table(file.path(pipeline.folder,which.max,"mhl.txt"),sep="\t",skip=1)
base.anno <- as.character(mhl[,1])
base.anno <- strsplit(base.anno,'[[:punct:]]')
base.anno <- as.data.frame(base.anno)
base.anno <- t(base.anno)
base.anno <- GRanges(Rle(base.anno[,1]),IRanges(as.numeric(base.anno[,2]),as.numeric(base.anno[,3])))
for(i in 1:length(sample.names)){
  name <- sample.names[i]
  file <- folders[i]
  mhl <- read.table(file.path(file,"mhl.txt"),sep="\t",skip=1)
  anno <- as.character(mhl[,1])
  anno <- strsplit(anno,'[[:punct:]]')
  anno <- as.data.frame(anno)
  anno <- t(anno)
  anno <- GRanges(Rle(anno[,1]),IRanges(as.numeric(anno[,2]),as.numeric(anno[,3])))
  op <- findOverlaps(base.anno,anno,type = "equal")
  new.data <- rep(NA,length(base.anno))
  new.data[queryHits(op)] <- mhl[subjectHits(op),2]
  complete.data <- cbind(complete.data,new.data)
  colnames(complete.data)[i] <- name
}
annotation <- base.anno
save(annotation,file=file.path(save.path,"mhl_annotation.RData"))
gzfi <- gzfile(file.path(save.path,"mhl.csv.gz"),"w")
write.csv(complete.data,gzfi)
close(gzfi) 

#' (iii) Summarize Epipolymorphism and Entropy
#' Epipolymorphism
max <- 0
which.max <- ""
for(i in 1:length(sample.names)){
  name <- sample.names[i]
  file <- folders[i]
  size <- file.info(file.path(file,"epipoly.csv"))$size
  if(size>max){
    max <- size
    which.max <- name
  }
}

complete.data <- c()
i <- 1
name <- sample.names[i]
file <- folders[i]
epipoly <- read.csv(file.path(pipeline.folder,which.max,"epipoly.csv"))
base.anno <- GRanges(Rle(epipoly$chromosome),IRanges(start=epipoly$start,end=epipoly$end))
for(i in 1:length(sample.names)){
  name <- sample.names[i]
  file <- folders[i]
  epipoly <- read.csv(file.path(file,"epipoly.csv"))
  anno <- GRanges(Rle(epipoly$chromosome),IRanges(start=epipoly$start,end=epipoly$end))
  op <- findOverlaps(base.anno,anno,type = "equal")
  new.data <- rep(NA,length(base.anno))
  new.data[queryHits(op)] <- epipoly[subjectHits(op),5]
  complete.data <- cbind(complete.data,new.data)
  colnames(complete.data)[i] <- name
}
annotation <- base.anno
save(annotation,file=file.path(save.path,"epipoly_annotation.RData"))
gzfi <- gzfile(file.path(save.path,"epipoly.csv.gz"),"w")
write.csv(complete.data,gzfi)
close(gzfi)

#' Entropy
complete.data <- c()
for(i in 1:length(sample.names)){
  name <- sample.names[i]
  file <- folders[i]
  entropy <- read.csv(file.path(file,"entropy.csv"))
  anno <- GRanges(Rle(entropy$chromosome),IRanges(start=entropy$start,end=entropy$end))
  op <- findOverlaps(base.anno,anno,type = "equal")
  new.data <- rep(NA,length(base.anno))
  new.data[queryHits(op)] <- entropy[subjectHits(op),5]
  complete.data <- cbind(complete.data,new.data)
  colnames(complete.data)[i] <- name
}
annotation <- base.anno
save(annotation,file=file.path(save.path,"entropy_annotation.RData"))
gzfi <- gzfile(file.path(save.path,"entropy.csv.gz"),"w")
write.csv(complete.data,gzfi)
close(gzfi)
