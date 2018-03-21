#######################################################################################################
#' score_comparison.R
#-----------------------------------------------------------------------------------------------------
#' This scripts produces a plot similar to Figure 2 to visualize the interconnections of the scores.
#' You should specify the folder, where the results of performing the simulation experiment on
#' coverage dependence are stored as an argument and the folder where to store the results. This
#' scripts selects from the pipeline folder those cases, where 25,000 reads were simulated.
#' folder.path: path to the results_pipeline folder for the simulation experiment
#' plot.path: folder to store the results
#' 
#' A scatterplot for each pair of scores and Spearman correlations for between the values produced
#' by the plots.
folder.path <- "path_to_results_pipeline"
plot.path <- getwd()
result <- as.data.frame(matrix(nrow=6,ncol=5))
colnames(result) <- c('qFDRP','PDR','MHL','Epipoly','Entropy')
row.names(result) <- c('FDRP','qFDRP','PDR','MHL','Epipoly','Entropy')
library(RnBeads)
folders <- list.files(folder.path,pattern = '25000',full.names = TRUE)

# FDRP, qFDRP, PDR
data <- c()
for(case in folders){
  fdrp <- unlist(read.csv(file.path(case,'FDRP','FDRP.csv')))
  qfdrp <- unlist(read.csv(file.path(case,'qFDRP','qFDRP.csv')))
  pdr <- unlist(read.csv(file.path(case,'PDR','PDR.csv')))
  temp <- data.frame(FDRP=fdrp,qFDRP=qfdrp,PDR=pdr)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('FDRP','qFDRP','PDR')
plot <- ggplot(data,aes(x=qFDRP,y=PDR))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'qFDRPvsPDR.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$PDR,data$qFDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$qFDRP,na.omit(data)$PDR,method = 'spearman')
result['qFDRP','PDR'] <- cori

plot <- ggplot(data,aes(x=FDRP,y=PDR))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'FDRPvsPDR.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$PDR,data$FDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$FDRP,na.omit(data)$PDR,method = 'spearman')
result['FDRP','PDR'] <- cori

plot <- ggplot(data,aes(x=FDRP,y=qFDRP))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'FDRPvsqFDRP.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$FDRP,data$qFDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$FDRP,na.omit(data)$qFDRP,method = 'spearman')
result['FDRP','qFDRP'] <- cori


# FDRP, qFDRP, PDR and MHL
data <- c()
for(case in folders){
  fdrp <- unlist(read.csv(file.path(case,'FDRP','FDRP.csv')))
  qfdrp <- unlist(read.csv(file.path(case,'qFDRP','qFDRP.csv')))
  pdr <- unlist(read.csv(file.path(case,'PDR','PDR.csv')))
  temp <- data.frame(FDRP=fdrp,qFDRP=qfdrp,PDR=pdr)
  load(file.path(case,'FDRP','annotation.RData'))
  mhl <- tryCatch(read.table(file.path(case,'mhl.txt'),sep='\t',skip=1),error=function(e){e})
  if(inherits(mhl,"error"))next
  anno_mhl <- as.character(mhl[,1])
  anno_mhl <- strsplit(anno_mhl,'[[:punct:]]')
  anno_mhl <- as.data.frame(anno_mhl)
  anno_mhl <- t(anno_mhl)
  anno_mhl <- GRanges(Rle(anno_mhl[,1]),IRanges(as.numeric(anno_mhl[,2]),as.numeric(anno_mhl[,3])))
  op <- findOverlaps(annotation,anno_mhl)
  temp <- temp[queryHits(op),]
  mhl <- unlist(mhl[,2])[subjectHits(op)]
  temp <- data.frame(temp,MHL=mhl)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('FDRP','qFDRP','PDR','MHL')
plot <- ggplot(data,aes(x=FDRP,y=MHL))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'FDRPvsMHL.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$MHL,data$FDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$FDRP,na.omit(data)$MHL,method = 'spearman')
result['FDRP','MHL'] <- cori

plot <- ggplot(data,aes(x=qFDRP,y=MHL))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'qFDRPvsMHL.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$MHL,data$qFDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$qFDRP,na.omit(data)$MHL,method = 'spearman')
result['qFDRP','MHL'] <- cori

plot <- ggplot(data,aes(x=PDR,y=MHL))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'PDRvsMHL.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$MHL,data$PDR,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$PDR,na.omit(data)$MHL,method = 'spearman')
result['PDR','MHL'] <- cori


# FDRP, qFDRP, PDR and Epipolymorphism, Entropy
data <- c()
for(case in folders){
  fdrp <- unlist(read.csv(file.path(case,'FDRP','FDRP.csv')))
  qfdrp <- unlist(read.csv(file.path(case,'qFDRP','qFDRP.csv')))
  pdr <- unlist(read.csv(file.path(case,'PDR','PDR.csv')))
  temp <- data.frame(FDRP=fdrp,qFDRP=qfdrp,PDR=pdr)
  load(file.path(case,'FDRP','annotation.RData'))
  epipoly <- read.csv(file.path(case,'epipoly.csv'))
  anno_epi <- GRanges(Rle(epipoly$chromosome),IRanges(epipoly$start,epipoly$end))
  entropy <- read.csv(file.path(case,'entropy.csv'))
  epipoly <- epipoly$Epipolymorphism
  entropy <- entropy$Entropy
  op <- findOverlaps(annotation,anno_epi)
  op <- op[match(unique(queryHits(op)),queryHits(op))]
  e_data <- data.frame(Epipolymorphism=epipoly,Entropy=entropy)
  temp <- temp[queryHits(op),]
  temp <- tryCatch(aggregate(temp,by=list(subjectHits(op)),mean,na.rm=TRUE),error=function(e)e)
  if(inherits(temp,'error'))next
  temp <- temp[,c(2,3,4)]
  e_data <- e_data[unique(subjectHits(op)),]
  temp <- data.frame(temp,e_data)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('FDRP','qFDRP','PDR','Epipolymorphism','Entropy')

plot <- ggplot(data,aes(x=FDRP,y=Epipolymorphism))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'FDRPvsEpipoly.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Epipolymorphism,data$FDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$FDRP,na.omit(data)$Epipolymorphism,method = 'spearman')
result['FDRP','Epipoly'] <- cori

plot <- ggplot(data,aes(x=qFDRP,y=Epipolymorphism))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'qFDRPvsEpipoly.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Epipolymorphism,data$qFDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$qFDRP,na.omit(data)$Epipolymorphism,method = 'spearman')
result['qFDRP','Epipoly'] <- cori

plot <- ggplot(data,aes(x=PDR,y=Epipolymorphism))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'PDRvsEpipoly.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Epipolymorphism,data$PDR,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$PDR,na.omit(data)$Epipolymorphism,method = 'spearman')
result['PDR','Epipoly'] <- cori

plot <- ggplot(data,aes(x=FDRP,y=Entropy))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'FDRPvsEntropy.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Entropy,data$FDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$FDRP,na.omit(data)$Entropy,method = 'spearman')
result['FDRP','Entropy'] <- cori

plot <- ggplot(data,aes(x=qFDRP,y=Entropy))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'qFDRPvsEntropy.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Entropy,data$qFDRP,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$qFDRP,na.omit(data)$Entropy,method = 'spearman')
result['qFDRP','Entropy'] <- cori

plot <- ggplot(data,aes(x=PDR,y=Entropy))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'PDRvsEntropy.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$PDR,data$Entropy,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$PDR,na.omit(data)$Entropy,method = 'spearman')
result['PDR','Entropy'] <- cori


# MHL and Entropy, Epipoly
data <- c()
for(case in folders){
  mhl <- tryCatch(read.table(file.path(case,'mhl.txt'),sep='\t',skip=1),error=function(e){e})
  if(inherits(mhl,"error"))next
  anno_mhl <- as.character(mhl[,1])
  anno_mhl <- strsplit(anno_mhl,'[[:punct:]]')
  anno_mhl <- as.data.frame(anno_mhl)
  anno_mhl <- t(anno_mhl)
  anno_mhl <- GRanges(Rle(anno_mhl[,1]),IRanges(as.numeric(anno_mhl[,2]),as.numeric(anno_mhl[,3])))
  epipoly <- read.csv(file.path(case,'epipoly.csv'))
  anno_epi <- GRanges(Rle(epipoly$chromosome),IRanges(epipoly$start,epipoly$end))
  entropy <- read.csv(file.path(case,'entropy.csv'))
  epipoly <- epipoly$Epipolymorphism
  entropy <- entropy$Entropy
  op <- findOverlaps(anno_mhl,anno_epi)
  op <- op[match(unique(queryHits(op)),queryHits(op))]
  e_data <- data.frame(Epipolymorphism=epipoly,Entropy=entropy)
  mhl <- mhl[queryHits(op),]
  mhl <- tryCatch(aggregate(mhl[,2],by=list(subjectHits(op)),mean,na.rm=TRUE),error=function(e)e)
  if(inherits(mhl,'error'))next
  mhl <- mhl[,c(2)]
  e_data <- e_data[unique(subjectHits(op)),]
  temp <- data.frame(MHL=mhl,e_data)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('MHL','Epipolymorphism','Entropy')

plot <- ggplot(data,aes(x=MHL,y=Entropy))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'MHLvsEntropy.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$MHL,data$Entropy,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$MHL,na.omit(data)$Entropy,method = 'spearman')
result['MHL','Entropy'] <- cori

plot <- ggplot(data,aes(x=MHL,y=Epipolymorphism))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'MHLvsEpipoly.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Epipolymorphism,data$MHL,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$MHL,na.omit(data)$Epipolymorphism,method = 'spearman')
result['MHL','Epipoly'] <- cori


# Entropy and Epipoly
data <- c()
for(case in folders){
  epipoly <- read.csv(file.path(case,'epipoly.csv'))
  entropy <- read.csv(file.path(case,'entropy.csv'))
  epipoly <- epipoly$Epipolymorphism
  entropy <- entropy$Entropy
  temp <- data.frame(Epipolymorphism=epipoly,Entropy=entropy)
  data <- rbind(data,temp)
}
data <- as.data.frame(data)
colnames(data) <- c('Epipolymorphism','Entropy')

plot <- ggplot(data,aes(y=Entropy,x=Epipolymorphism))+geom_point(color='firebrick4')+stat_density_2d(aes(fill=..level..),geom='polygon',alpha=0.75)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_gradient(low='dodgerblue4',high='white')+
  xlim(-0.01,1.01)+ylim(-0.01,1.01)+geom_abline(slope=1,intercept=0)+xlab("")+ylab("")
ggsave(file.path(plot.path,'EntropyvsEpipoly.png'),plot,device="png",height=11,width=8.5,units="in",dpi=300)
cori <- cor(data$Epipolymorphism,data$Entropy,method = 'spearman')
if (is.na(cori)) cori <- cor(na.omit(data)$Entropy,na.omit(data)$Epipolymorphism,method = 'spearman')
result['Epipoly','Entropy'] <- cori

write.csv(result,file.path(plot.path,"spearman_correlations.csv"))
