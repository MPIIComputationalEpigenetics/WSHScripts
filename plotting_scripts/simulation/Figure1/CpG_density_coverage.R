#######################################################################################################
#' coverage_dependence.R
#-----------------------------------------------------------------------------------------------------
#' This scripts produces a plot similar to Figure 1D to visualize the overall dependence of the scores
#' on site-wise coverage and CpG density. You should specify the folder, where the results of
#' performing the simulation experiment on coverage dependence are stored as an argument and the
#' folder where to store the results. This scripts selects from the pipeline folder those cases,
#' where 25,000 reads were simulated.
#' as another.
#' folder.path: path to the results_pipeline folder for the simulation experiment
#' plot.path: folder to store the results
#' 
#' A scatterplot for each of the scores, Spearman correlations for dependence on site-wise coverage
#' and CpG density as well as a heatmap visualizing the correlations.
folder.path <- "path_to_results_pipeline"
plot.path <- getwd()
result <- as.data.frame(matrix(nrow=6,ncol=1))
row.names(result) <- c('FDRP','qFDRP','PDR','MHL','Epipoly','Entropy')
library(RnBeads)

folders <- list.files(folder.path,pattern='25000',full.names=TRUE)
####################################################################################################
#' CpG density
#'--------------------------------------------------------------------------------------------------
data <- c()
#' FDRP
for(case in folders){
  fdrp <- read.csv(file.path(case,'FDRP','FDRP.csv'))
  load(file.path(case,'FDRP','annotation.RData'))
  mini <- min(start(annotation))
  maxi <- max(end(annotation))
  chrom <- as.character(unique(seqnames(annotation)))
  tiles <- seq(mini,maxi,by=50)
  tiles <- paste0(chrom,":",tiles,"-",tiles+50)
  tiles <- GRanges(tiles)
  op <- findOverlaps(tiles,annotation)
  op <- op[match(unique(subjectHits(op)),subjectHits(op))]
  op <- as.list(op)
  cpg_count <- lengths(op)
  op <- findOverlaps(annotation,tiles,select='first')
  fdrp <- aggregate(fdrp,by=list(op),mean,na.rm=TRUE)
  fdrp <- fdrp[,2]
  cpg_count <- cpg_count[cpg_count!=0]
  nas <- is.na(fdrp)
  cpg_count <- cpg_count[!nas]
  fdrp <- fdrp[!nas]
  temp <- cbind(FDRP=fdrp,Count=cpg_count)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('FDRP','Count')
plot <- ggplot(data,aes(x=Count,y=FDRP))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  ylim(-0.01,1.01)+ylab('FDRP')+geom_smooth(method='lm',color='black',se=FALSE)
ggsave(file.path(plot.path,'FDRP.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$FDRP,data$Count,method='spearman')
if(is.na(cori)) cori <- cor(na.omit(data)$FDRP,na.omit(data)$Count,method='spearman')
result['FDRP',1] <- cori

#' qFDRP
data <- c()
for(case in folders){
  fdrp <- read.csv(file.path(case,'qFDRP','qFDRP.csv'))
  load(file.path(case,'qFDRP','annotation.RData'))
  mini <- min(start(annotation))
  maxi <- max(end(annotation))
  chrom <- as.character(unique(seqnames(annotation)))
  tiles <- seq(mini,maxi,by=50)
  tiles <- paste0(chrom,":",tiles,"-",tiles+50)
  tiles <- GRanges(tiles)
  op <- findOverlaps(tiles,annotation)
  op <- op[match(unique(subjectHits(op)),subjectHits(op))]
  op <- as.list(op)
  cpg_count <- lengths(op)
  op <- findOverlaps(annotation,tiles,select='first')
  fdrp <- aggregate(fdrp,by=list(op),mean,na.rm=TRUE)
  fdrp <- fdrp[,2]
  cpg_count <- cpg_count[cpg_count!=0]
  nas <- is.na(fdrp)
  cpg_count <- cpg_count[!nas]
  fdrp <- fdrp[!nas]
  temp <- cbind(FDRP=fdrp,Count=cpg_count)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('qFDRP','Count')
plot <- ggplot(data,aes(x=Count,y=qFDRP))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  ylim(-0.01,1.01)+ylab('qFDRP')+geom_smooth(method='lm',color='black',se=FALSE)
ggsave(file.path(plot.path,'qFDRP.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$qFDRP,data$Count,method='spearman')
if(is.na(cori)) cori <- cor(na.omit(data)$qFDRP,na.omit(data)$Count,method='spearman')
result['qFDRP',1] <- cori

#' PDR
data <- c()
for(case in folders){
  fdrp <- read.csv(file.path(case,'PDR','PDR.csv'))
  load(file.path(case,'PDR','annotation.RData'))
  mini <- min(start(annotation))
  maxi <- max(end(annotation))
  chrom <- as.character(unique(seqnames(annotation)))
  tiles <- seq(mini,maxi,by=50)
  tiles <- paste0(chrom,":",tiles,"-",tiles+50)
  tiles <- GRanges(tiles)
  op <- findOverlaps(tiles,annotation)
  op <- op[match(unique(subjectHits(op)),subjectHits(op))]
  op <- as.list(op)
  cpg_count <- lengths(op)
  op <- findOverlaps(annotation,tiles,select='first')
  fdrp <- aggregate(fdrp,by=list(op),mean,na.rm=TRUE)
  fdrp <- fdrp[,2]
  cpg_count <- cpg_count[cpg_count!=0]
  nas <- is.na(fdrp)
  cpg_count <- cpg_count[!nas]
  fdrp <- fdrp[!nas]
  temp <- cbind(FDRP=fdrp,Count=cpg_count)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('PDR','Count')
plot <- ggplot(data,aes(x=Count,y=PDR))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  ylim(-0.01,1.01)+ylab('PDR')+geom_smooth(method='lm',color='black',se=FALSE)
ggsave(file.path(plot.path,'PDR.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$PDR,data$Count,method='spearman')
if(is.na(cori)) cori <- cor(na.omit(data)$PDR,na.omit(data)$Count,method='spearman')
result['PDR',1] <- cori

#' Epipolymorphism
data <- c()
for(case in folders){
  epipoly <- read.csv(file.path(case,'epipoly.csv'))
  anno_epi <- GRanges(Rle(epipoly$chromosome),IRanges(epipoly$start,epipoly$end))
  if(length(anno_epi)==0) next
  mini <- min(start(anno_epi))
  maxi <- max(end(anno_epi))
  chrom <- as.character(unique(seqnames(anno_epi)))
  tiles <- seq(mini,maxi,by=50)
  tiles <- paste0(chrom,":",tiles,"-",tiles+50)
  tiles <- GRanges(tiles)
  op <- findOverlaps(tiles,anno_epi)
  op <- op[match(unique(subjectHits(op)),subjectHits(op))]
  op <- as.list(op)
  cpg_count <- lengths(op)
  op <- findOverlaps(anno_epi,tiles,select='first')
  epipoly <- epipoly$Epipolymorphism
  epipoly <- aggregate(epipoly,by=list(op),mean,na.rm=TRUE)
  epipoly <- epipoly[,2]
  cpg_count <- cpg_count[cpg_count!=0]
  nas <- is.na(epipoly)
  cpg_count <- cpg_count[!nas]
  epipoly <- epipoly[!nas]
  temp <- cbind(Epipoly=epipoly,Count=cpg_count)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('Epipoly','Count')
plot <- ggplot(data,aes(x=Count,y=Epipoly))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  ylim(-0.01,1.01)+ylab('Epipolymorphism')+geom_smooth(method='lm',color='black',se=FALSE)
ggsave(file.path(plot.path,'Epipoly.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Epipoly,data$Count,method='spearman')
if(is.na(cori)) cori <- cor(na.omit(data)$Epipoly,na.omit(data)$Count,method='spearman')
result['Epipoly',1] <- cori

#' Entropy
data <- c()
for(case in folders){
  epipoly <- read.csv(file.path(case,'entropy.csv'))
  anno_epi <- GRanges(Rle(epipoly$chromosome),IRanges(epipoly$start,epipoly$end))
  if(length(anno_epi)==0) next
  mini <- min(start(anno_epi))
  maxi <- max(end(anno_epi))
  chrom <- as.character(unique(seqnames(anno_epi)))
  tiles <- seq(mini,maxi,by=50)
  tiles <- paste0(chrom,":",tiles,"-",tiles+50)
  tiles <- GRanges(tiles)
  op <- findOverlaps(tiles,anno_epi)
  op <- op[match(unique(subjectHits(op)),subjectHits(op))]
  op <- as.list(op)
  cpg_count <- lengths(op)
  op <- findOverlaps(anno_epi,tiles,select='first')
  epipoly <- epipoly$Entropy
  epipoly <- aggregate(epipoly,by=list(op),mean,na.rm=TRUE)
  epipoly <- epipoly[,2]
  cpg_count <- cpg_count[cpg_count!=0]
  nas <- is.na(epipoly)
  cpg_count <- cpg_count[!nas]
  epipoly <- epipoly[!nas]
  temp <- cbind(Epipoly=epipoly,Count=cpg_count)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('Entropy','Count')
plot <- ggplot(data,aes(x=Count,y=Entropy))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  ylim(-0.01,1.01)+ylab('Entropy')+geom_smooth(method='lm',color='black',se=FALSE)
ggsave(file.path(plot.path,'Entropy.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Entropy,data$Count,method='spearman')
if(is.na(cori)) cori <- cor(na.omit(data)$Entropy,na.omit(data)$Count,method='spearman')
result['Entropy',1] <- cori

#' MHL
data <- c()
for(case in folders){
  mhl <- tryCatch(read.table(file.path(case,'mhl.txt'),sep='\t',skip=1),error=function(e){e})
  if(inherits(mhl,'error')) next
  anno_mhl <- as.character(mhl[,1])
  anno_mhl <- strsplit(anno_mhl,'[[:punct:]]')
  anno_mhl <- as.data.frame(anno_mhl)
  anno_mhl <- t(anno_mhl)
  anno_mhl <- GRanges(Rle(anno_mhl[,1]),IRanges(as.numeric(anno_mhl[,2]),as.numeric(anno_mhl[,3])))
  mini <- min(start(anno_mhl))
  maxi <- max(end(anno_mhl))
  chrom <- as.character(unique(seqnames(anno_mhl)))
  tiles <- seq(mini,maxi,by=50)
  tiles <- paste0(chrom,":",tiles,"-",tiles+50)
  tiles <- GRanges(tiles)
  op <- findOverlaps(tiles,anno_mhl)
  op <- op[match(unique(subjectHits(op)),subjectHits(op))]
  op <- as.list(op)
  cpg_count <- lengths(op)
  op <- findOverlaps(anno_mhl,tiles,select='first')
  mhl <- as.numeric(mhl[,2])
  mhl <- aggregate(mhl,by=list(op),mean,na.rm=TRUE)
  mhl <- mhl[,2]
  cpg_count <- cpg_count[cpg_count!=0]
  nas <- is.na(mhl)
  cpg_count <- cpg_count[!nas]
  mhl <- mhl[!nas]
  temp <- cbind(Epipoly=mhl,Count=cpg_count)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('MHL','Count')
plot <- ggplot(data,aes(x=Count,y=MHL))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  ylim(-0.01,1.01)+ylab('MHL')+geom_smooth(method='lm',color='black',se=FALSE)
ggsave(file.path(plot.path,'MHL.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$MHL,data$Count,method='spearman')
if(is.na(cori)) cori <- cor(na.omit(data)$MHL,na.omit(data)$Count,method='spearman')
result['MHL',1] <- cori

####################################################################################################
#' site-wise coverage
#'--------------------------------------------------------------------------------------------------

#' FDRP, qFDRP, PDR
data <- c()
cors <- matrix(nrow=6,ncol=1)
row.names(cors) <- c("FDRP",'qFDRP','PDR','MHL','Epipolymorphism','Entropy')
for(file in folders){
  rnb.set <- load.rnb.set(file.path(file,'rnbSet.zip'))
  anno_set <- annotation(rnb.set)
  anno_set <- GRanges(Rle(anno_set$Chromosome),IRanges(start=anno_set$Start,end=anno_set$End))
  fdrp <- read.csv(file.path(file,'FDRP','FDRP.csv'))
  load(file.path(file,'FDRP','annotation.RData'))
  op <- findOverlaps(annotation,anno_set)
  fdrp <- fdrp[queryHits(op),]
  covg <- covg(rnb.set)
  covg <- covg[subjectHits(op),]
  qfdrp <- read.csv(file.path(file,'qFDRP','qFDRP.csv'))
  qfdrp <- qfdrp[queryHits(op),]
  pdr <- read.csv(file.path(file,'PDR','PDR.csv'))
  pdr <- pdr[queryHits(op),]
  temp <- cbind(Coverage=covg,FDRP=fdrp,qFDRP=qfdrp,PDR=pdr)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('Coverage','FDRP','qFDRP','PDR')
colors <- rnb.getOption('colors.category')
temp <- colors[2]
colors[2] <- colors[5]
colors[5] <- temp
cors['FDRP',] <- cor(data$Coverage,data$FDRP,method='spearman')
cors['qFDRP',] <- cor(data$Coverage,data$qFDRP,method='spearman')
cors['PDR',] <- cor(na.omit(data)$Coverage,na.omit(data)$FDRP,method='spearman')

#' MHL
data <- c()
for(file in folders){
  rnb.set <- load.rnb.set(file.path(file,'rnbSet.zip'))
  anno_set <- annotation(rnb.set)
  anno_set <- GRanges(Rle(anno_set$Chromosome),IRanges(start=anno_set$Start,end=anno_set$End))
  mhl <- tryCatch(read.table(file.path(file,'mhl.txt'),sep='\t',skip=1),error=function(e){e})
  if(inherits(mhl,'error')) next
  anno_mhl <- as.character(mhl[,1])
  anno_mhl <- strsplit(anno_mhl,'[[:punct:]]')
  anno_mhl <- as.data.frame(anno_mhl)
  anno_mhl <- t(anno_mhl)
  anno_mhl <- GRanges(Rle(anno_mhl[,1]),IRanges(as.numeric(anno_mhl[,2]),as.numeric(anno_mhl[,3])))
  op <- findOverlaps(anno_set,anno_mhl)
  load(file.path(file,'FDRP','annotation.RData'))
  covg <- covg(rnb.set)
  covg <- covg[queryHits(op),]
  mhl <- as.numeric(mhl[subjectHits(op),2])
  temp <- cbind(Coverage=covg,MHL=mhl)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('Coverage','MHL')
cors['MHL',] <- cor(data$Coverage,data$MHL,method='spearman')

#' Epipolymorphism and Entropy
data <- c()
for(file in folders){
  rnb.set <- load.rnb.set(file.path(file,'rnbSet.zip'))
  anno_set <- annotation(rnb.set)
  anno_set <- GRanges(Rle(anno_set$Chromosome),IRanges(start=anno_set$Start,end=anno_set$End))
  epipoly <- read.csv(file.path(file,'epipoly.csv'))
  anno_4 <- GRanges(Rle(epipoly$chromosome),IRanges(start=epipoly$start,end=epipoly$end))
  if(length(anno_4)<1){
    next
  }
  op <- findOverlaps(anno_4,anno_set)
  epipoly <-epipoly[unique(queryHits(op)),5]
  covg <- covg(rnb.set)
  covg <- covg[subjectHits(op),]
  covg <- aggregate(covg,by=list(queryHits(op)),mean,na.rm=TRUE)
  covg <- covg[,2]
  entropy <- read.csv(file.path(file,'entropy.csv'))
  entropy <- entropy[unique(queryHits(op)),5]
  temp <- cbind(Coverage=covg,Epipolymorphism=epipoly,Entropy=entropy)
  data <- rbind(data,temp)
}
colnames(data) <- c('Coverage','Epipolymorphism','Entropy')
data <- as.data.frame(data)
cors["Epipolymorphism",] <- cor(data$Epipolymorphism,data$Coverage,method='spearman')
cors["Entropy",] <- cor(data$Entropy,data$Coverage,method='spearman')

result <- cbind(result,cors)
colnames(result) <- c("Measure","CpG_density","Coverage")

write.csv(result,file.path(plot.path,"spearman_correlations.csv"))

data <- melt(result)
colnames(data) <- c("Measure","Feature","Value")
data$Feature <- as.character(data$Feature)
data$Feature[data$Feature%in%"CpG_density"] <- "CpG Density"
data$Feature <- factor(data$Feature,levels=rev(c('CpG Density','Coverage')))
data$Measure <- as.character(data$Measure)
data$Measure[data$Measure%in%"Epipoly"] <- "Epipolymorphism"
data$Measure <- factor(data$Measure,levels=rev(c("FDRP","qFDRP","PDR","MHL","Epipolymorphism","Entropy")))
plot <- ggplot(data,aes(x=Feature,y=Measure,fill=Value))+geom_tile()+xlab("")+
  ylab("Score")+theme(text=element_text(size=30,face="bold"),
                                                      panel.background = element_rect(fill='white'),
                                                      legend.key.size = unit(10,"mm"),
                                                      legend.text = element_text(size=20),
                                                      legend.title = element_text(size=25,face="bold"),
                                                      legend.key = element_rect(fill='white'),legend.position = 'top')+
  scale_fill_continuous(name="Correlation")+geom_text(aes(y=Measure,label=round(Value,2)),color="white",size=15)
ggsave(file.path(plot.path,"heatmap.pdf"),device="pdf",height=11,width=8.5,units="in")
