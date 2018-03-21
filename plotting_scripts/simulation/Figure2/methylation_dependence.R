#######################################################################################################
#' methylation_dependence.R
#-----------------------------------------------------------------------------------------------------
#' This scripts produces a plot similar to Figure 2 to visualize the connection of the scores
#' to methylation level. You should specify the folder, where the results of
#' performing the simulation experiment on a swithing methylation state are stored as an argument
#' and the folder where to store the results.
#' folder.path: path to the results_pipeline folder for the simulation experiment
#' plot.path: folder to store the results
#' 
#' A scatterplot for each of the scores, Spearman correlations for dependence on site-wise coverage
#' and CpG density as well as a heatmap visualizing the correlations.
folder.path <- "path_to_results_pipeline"
plot.path <- getwd()
library(reshape2)
library(ggplot2)
library(RnBeads)
result <- as.data.frame(matrix(nrow=6,ncol=1))
colnames(result) <- c('Methylation')
row.names(result) <- c('FDRP','qFDRP','PDR','MHL','Epipoly','Entropy')
folders <- list.files(folder.path,full.names = TRUE)

# FDRP, qFDRP, PDR
df <- c()
for(file in folders){
  fdrp <- read.csv(file.path(file,'FDRP','FDRP.csv'))
  qfdrp <- read.csv(file.path(file,'qFDRP','qFDRP.csv'))
  pdr <- read.csv(file.path(file,'PDR','PDR.csv'))
  data <- cbind(as.numeric(unlist(fdrp)),as.numeric(unlist(qfdrp)),as.numeric(unlist(pdr)))
  load(file.path(file,'FDRP','annotation.RData'))
  anno <- annotation
  rnb.set <- load.rnb.set(file.path(file,'rnbSet.zip'))
  anno_set <- annotation(rnb.set)
  anno_set <- GRanges(Rle(anno_set$Chromosome),IRanges(start = anno_set$Start,end=anno_set$End))
  op <- findOverlaps(anno,anno_set)
  methData <- meth(rnb.set)
  methData <- methData[subjectHits(op)]
  data <- data[queryHits(op),]
  data <- data.frame(data,Methylation=methData)
  colnames(data) <- c('FDRP','qFDRP','PDR','Methylation')
  data <- as.data.frame(data)
  df <- rbind(df,data)
}
df <- as.data.frame(df)
colnames(df) <- c('FDRP','qFDRP','PDR','Methylation')
meth_factor <- rep('0-5 %',dim(df)[1])
for(i in seq(0,95,by=5)){
  meth_factor[df$Methylation>i/100&df$Methylation<=(i+5)/100] <- paste0(i,'-',i+5,' %')
}
meth_factor <- factor(meth_factor)
meth_factor <- factor(meth_factor,levels=c("0-5 %","5-10 %","10-15 %","15-20 %","20-25 %","25-30 %","30-35 %","35-40 %","40-45 %","45-50 %","50-55 %","55-60 %","60-65 %","65-70 %","70-75 %","75-80 %","80-85 %","85-90 %","90-95 %", "95-100 %"))
df$Methylation <- meth_factor
plot <- ggplot(df,aes(x=Methylation,y=FDRP))+geom_boxplot()+theme(text=element_text(size=30,face='bold'),
                                                                  panel.background=element_rect(color='black',fill='white'),
                                                                  panel.grid.major=element_line(color='grey80'),
                                                                  panel.grid.minor=element_line(color='grey95'),
                                                                  axis.text.x=element_text(size=20,angle=45,hjust=1))+
  xlab('Methylation')+ylab('FDRP')+ylim(-0.01,1.01)
ggsave(file.path(plot.path,'MethylationvsFDRP_segmentation_unit.pdf'),plot,device="pdf",height=11,width=8.5,units="in")

plot <- ggplot(df,aes(x=Methylation,y=qFDRP))+geom_boxplot()+theme(text=element_text(size=30,face='bold'),
                                                                  panel.background=element_rect(color='black',fill='white'),
                                                                  panel.grid.major=element_line(color='grey80'),
                                                                  panel.grid.minor=element_line(color='grey95'),
                                                                  axis.text.x=element_text(size=20,angle=45,hjust=1))+
  xlab('Methylation')+ylab('qFDRP')+ylim(-0.01,1.01)
ggsave(file.path(plot.path,'MethylationvsqFDRP_segmentation_unit.pdf'),plot,device="pdf",height=11,width=8.5,units="in")

plot <- ggplot(df,aes(x=Methylation,y=PDR))+geom_boxplot()+theme(text=element_text(size=30,face='bold'),
                                                                  panel.background=element_rect(color='black',fill='white'),
                                                                  panel.grid.major=element_line(color='grey80'),
                                                                  panel.grid.minor=element_line(color='grey95'),
                                                                  axis.text.x=element_text(size=20,angle=45,hjust=1))+
  xlab('Methylation')+ylab('PDR')+ylim(-0.01,1.01)
ggsave(file.path(plot.path,'MethylationvsPDR_segmentation_unit.pdf'),plot,device="pdf",height=11,width=8.5,units="in")

# MHL
df <- c()
for(file in folders){
  mhl <- tryCatch(read.csv(file.path(file,'mhl.txt'),sep='\t',skip=1),error=function(e){e})
  if(inherits(mhl,'error')) next
  anno_mhl <- as.character(mhl[,1])
  anno_mhl <- strsplit(anno_mhl,'[[:punct:]]')
  anno_mhl <- as.data.frame(anno_mhl)
  anno_mhl <- t(anno_mhl)
  anno_mhl <- GRanges(Rle(anno_mhl[,1]),IRanges(as.numeric(anno_mhl[,2]),as.numeric(anno_mhl[,3])))
  rnb.set <- load.rnb.set(file.path(file,'rnbSet.zip'))
  anno_set <- annotation(rnb.set)
  anno_set <- GRanges(Rle(anno_set$Chromosome),IRanges(anno_set$Start,anno_set$End))
  op <- findOverlaps(anno_mhl,anno_set)
  methData <- meth(rnb.set)
  methData <- methData[subjectHits(op)]
  mhl <- as.numeric(mhl[,2])
  mhl <- mhl[queryHits(op)]
  mhl <- aggregate(mhl,by=list(subjectHits(op)),mean,na.rm=TRUE)
  mhl <- mhl[,2]
  data <- data.frame(mhl,Methylation=methData)
  df <- rbind(df,data)
}
colnames(df) <- c('MHL','Methylation')
df <- as.data.frame(df)
meth_factor <- rep('0-5 %',dim(df)[1])
for(i in seq(0,95,by=5)){
  meth_factor[df$Methylation>i/100&df$Methylation<=(i+5)/100] <- paste0(i,'-',i+5,' %')
}
meth_factor <- factor(meth_factor)
meth_factor <- factor(meth_factor,levels=c("0-5 %","5-10 %","10-15 %","15-20 %","20-25 %","25-30 %","30-35 %","35-40 %","40-45 %","45-50 %","50-55 %","55-60 %","60-65 %","65-70 %","70-75 %","75-80 %","80-85 %","85-90 %","90-95 %", "95-100 %"))
df$Methylation <- meth_factor
plot <- ggplot(df,aes(x=Methylation,y=MHL))+geom_boxplot()+theme(text=element_text(size=30,face='bold'),
                                                                  panel.background=element_rect(color='black',fill='white'),
                                                                  panel.grid.major=element_line(color='grey80'),
                                                                  panel.grid.minor=element_line(color='grey95'),
                                                                  axis.text.x=element_text(size=20,angle=45,hjust=1))+
  xlab('Methylation')+ylab('MHL')+ylim(-0.01,1.01)
ggsave(file.path(plot.path,'MethylationvsMHL_segmentation_unit.pdf'),plot,device="pdf",height=11,width=8.5,units="in")

# more CpG scores
df <- c()
for(file in folders){
  epipoly <- read.csv(file.path(file,'epipoly.csv'))
  anno_epi <- GRanges(Rle(epipoly$chromosome),IRanges(epipoly$start,epipoly$end))
  rnb.set <- load.rnb.set(file.path(file,'rnbSet.zip'))
  anno_set <- annotation(rnb.set)
  anno_set <- GRanges(Rle(anno_set$Chromosome),IRanges(anno_set$Start,anno_set$End))
  op <- findOverlaps(anno_set,anno_epi)
  epipoly <- epipoly$Epipolymorphism
  epipoly <- epipoly[subjectHits(op)]
  epipoly <- tryCatch(aggregate(epipoly,by=list(queryHits(op)),mean,na.rm=TRUE),error=function(e){e})
  if(inherits(epipoly,'error')) next
  epipoly <- epipoly[,2]
  methData <- meth(rnb.set)
  methData <- methData[unique(queryHits(op))]
  entropy <- read.csv(file.path(file,'entropy.csv'))
  entropy <- entropy$Entropy
  entropy <- entropy[subjectHits(op)]
  entropy <- tryCatch(aggregate(entropy,by=list(queryHits(op)),mean,na.rm=TRUE),error=function(e){e})
  if(inherits(entropy,'error')) next
  entropy <- entropy[,2]
  data <- data.frame(epipoly,entropy,methData)
  df <- rbind(df,data)
}
colnames(df) <- c('Epipolymorphism','Entropy','Methylation')
data <- as.data.frame(df)
meth_factor <- rep('0-5 %',dim(df)[1])
for(i in seq(0,95,by=5)){
  meth_factor[df$Methylation>i/100&df$Methylation<=(i+5)/100] <- paste0(i,'-',i+5,' %')
}
meth_factor <- factor(meth_factor)
meth_factor <- factor(meth_factor,levels=c("0-5 %","5-10 %","10-15 %","15-20 %","20-25 %","25-30 %","30-35 %","35-40 %","40-45 %","45-50 %","50-55 %","55-60 %","60-65 %","65-70 %","70-75 %","75-80 %","80-85 %","85-90 %","90-95 %", "95-100 %"))
df$Methylation <- meth_factor
plot <- ggplot(df,aes(x=Methylation,y=Epipolymorphism))+geom_boxplot()+theme(text=element_text(size=30,face='bold'),
                                                                  panel.background=element_rect(color='black',fill='white'),
                                                                  panel.grid.major=element_line(color='grey80'),
                                                                  panel.grid.minor=element_line(color='grey95'),
                                                                  axis.text.x=element_text(size=20,angle=45,hjust=1))+
  xlab('Methylation')+ylab('Epipolymorphism')+ylim(-0.01,1.01)
ggsave(file.path(plot.path,'MethylationvsEpipoly_segmentation_unit.pdf'),plot,device="pdf",height=11,width=8.5,units="in")

plot <- ggplot(df,aes(x=Methylation,y=Entropy))+geom_boxplot()+theme(text=element_text(size=30,face='bold'),
                                                                             panel.background=element_rect(color='black',fill='white'),
                                                                             panel.grid.major=element_line(color='grey80'),
                                                                             panel.grid.minor=element_line(color='grey95'),
                                                                             axis.text.x=element_text(size=20,angle=45,hjust=1))+
  xlab('Methylation')+ylab('Entropy')+ylim(-0.01,1.01)
ggsave(file.path(plot.path,'MethylationvsEntropy_segmentation_unit.pdf'),plot,device="pdf",height=11,width=8.5,units="in")
