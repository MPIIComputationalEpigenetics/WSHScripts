#######################################################################################################
#' snp_dependence.R
#-----------------------------------------------------------------------------------------------------
#' This scripts produces a plot similar to Figure 1F to visualize the overall dependence of the scores
#' on single nucleotide polymorphisms (SNPs) introduced into the reads. You should specify the folder,
#' where the results of performing the simulation experiment are stored as an argument and the folder
#' where to store the results as another.
#' folder.path: path to the results_pipeline folder for the simulation experiment
#' plot.path: folder to store the results
#' 
#' A single plot is the result.
folder.path <- "path_to_results_pipeline"
plot.path <- getwd()
library(reshape2)
library(ggplot2)
library(RnBeads)
folders <- list.files(folder.path,full.names = T)
data <- data.frame()
for(file in folders){
  snps <- as.numeric(unlist(strsplit(file,"_"))[length(unlist(strsplit(file,"_")))])
  fdrp <- tryCatch(read.csv(file.path(file,'FDRP','FDRP.csv')),error=function(e){c()})
  fdrp <- mean(as.numeric(unlist(fdrp)),na.rm=TRUE)
  qfdrp <- tryCatch(read.csv(file.path(file,'qFDRP','qFDRP.csv')),error=function(e){c()})
  qfdrp <- mean(as.numeric(unlist(qfdrp)),na.rm=TRUE)
  pdr <- tryCatch(read.csv(file.path(file,'PDR','PDR.csv')),error=function(e){c()})
  pdr <- mean(as.numeric(unlist(pdr)),na.rm=TRUE)
  mhl <- tryCatch(read.table(file.path(file,'mhl.txt'),sep='\t',skip = 1),error=function(e)c())
  mhl <- mean(as.numeric(unlist(mhl[,2])),na.rm=TRUE)
  epipoly <- tryCatch(read.csv(file.path(file,'epipoly.csv')),error=function(e){c()})
  epipoly <- mean(as.numeric(unlist(epipoly$Epipolymorphism)),na.rm=TRUE)
  entropy <- tryCatch(read.csv(file.path(file,'entropy.csv')),error=function(e){c()})
  entropy <- mean(as.numeric(unlist(entropy$Entropy)),na.rm=TRUE)
  data <- rbind(data,c(snps,fdrp,qfdrp,pdr,mhl,epipoly,entropy))
}
colnames(data) <- c('SNPs','FDRP','qFDRP','PDR','MHL','Epipolymorphism','Entropy')
agg <- aggregate(data,by=list(data$SNPs),function(x){c(mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm = TRUE))})
agg <- agg[,-c(2)]
melted <- c()
for(i in 2:7){
  frame <- agg[,i]
  add <- data.frame(Measure=colnames(data)[i],SNPs=agg[,1],Mean=frame[,1],SD=frame[,2])
  melted <- rbind(melted,add)
}
colors <- rnb.getOption('colors.category')
temp <- colors[2]
colors[2] <- "#76bf23ff"
colors[1] <- "#00806fff"
colors[5] <- temp
plot <- ggplot(melted,aes(x=SNPs,y=Mean,ymin=Mean-2*SD,ymax=Mean+2*SD,color=Measure,shape=Measure))+geom_point(size=3)+
  geom_line(size=1.2)+theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
                            panel.grid.major = element_line(color='grey80'),panel.grid.minor = element_line(color='grey95'),
                            legend.key = element_rect(fill='white'),legend.position = 'top',legend.text = element_text(size=20))+
  scale_color_manual(values=colors)+xlab('Number of SNPs')+ylim(0,1)+xlim(0,9.9)
ggsave(paste0(plot.path,"snp_dependence.pdf"),plot,device="pdf",height=11,width=8.5,units="in")
