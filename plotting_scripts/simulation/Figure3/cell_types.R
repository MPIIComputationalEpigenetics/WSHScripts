#####################################################################################################
#' cell_types.R
#----------------------------------------------------------------------------------------------------
#' This script produces figures like Figure 3C in the paper: A positional view of the ISH scores
#' together with the methylation level. You need to specify two arguments
#' folder.path: the results_pipeline folder with the simulation results obtained by applying the 
#'          simulation scripts for cell type simulation
#' plot.path: a path, where the corresponding plots should be stored
#' for each simulation run conducted, a plot will be created, where the file name is the starting
#' chrN_abced_abced+length
folder.path <- "path_to_results_pipeline"
plot.path <- "path_to_store_plots" 
library(reshape2)
library(ggplot2)
library(RnBeads)
folders <- list.files(folder.path,full.names = TRUE)
count <- 1
for(file in folders){
  data <- data.frame()
  fdrp <- read.csv(file.path(file,'FDRP','FDRP.csv'))
  qfdrp <- read.csv(file.path(file,'qFDRP','qFDRP.csv'))
  pdr <- read.csv(file.path(file,'PDR','PDR.csv'))
  data <- cbind(as.numeric(unlist(fdrp)),as.numeric(unlist(qfdrp)),as.numeric(unlist(pdr)))
  load(file.path(file,'FDRP','annotation.RData'))
  anno <- annotation
  data <- cbind(data,start(anno))
  mhl <- read.table(file.path(file,'mhl.txt'),sep='\t',skip=1)
  anno_mhl <- as.character(mhl[,1])
  anno_mhl <- strsplit(anno_mhl,'[[:punct:]]')
  anno_mhl <- as.data.frame(anno_mhl)
  anno_mhl <- t(anno_mhl)
  anno_mhl <- GRanges(Rle(anno_mhl[,1]),IRanges(as.numeric(anno_mhl[,2]),as.numeric(anno_mhl[,3])))
  op <- findOverlaps(anno,anno_mhl)
  mhl <- as.numeric(mhl[,2])
  mhl <- mhl[subjectHits(op)]
  mhl <- aggregate(mhl,by=list(queryHits(op)),mean,na.rm=TRUE)
  mhl <- mhl[,2]
  repla <- rep(NA,dim(data)[1])
  repla[unique(queryHits(op))] <- mhl
  data <- data.frame(data,repla)
  epipoly <- read.csv(file.path(file,'epipoly.csv'))
  anno_epi <- GRanges(Rle(epipoly$chromosome),IRanges(epipoly$start,epipoly$end))
  op <- findOverlaps(anno,anno_epi)
  epipoly <- epipoly$Epipolymorphism
  epipoly <- epipoly[subjectHits(op)]
  if(length(op)!=0){
   epipoly <- aggregate(epipoly,by=list(queryHits(op)),mean,na.rm=TRUE)
   epipoly <- epipoly[,2]
  }
  repla <- rep(NA,dim(data)[1])
  repla[unique(queryHits(op))] <- epipoly
  data <- data.frame(data,repla)
  entropy <- read.csv(file.path(file,'entropy.csv'))
  anno_en <- GRanges(Rle(entropy$chromosome),IRanges(entropy$start,entropy$end))
  op <- findOverlaps(anno,anno_en)
  entropy <- entropy$Entropy
  entropy <- entropy[subjectHits(op)]
  if(length(op)!=0){
    entropy <- aggregate(entropy,by=list(queryHits(op)),mean,na.rm=TRUE)
    entropy <- entropy[,2]
  }
  repla <- rep(NA,dim(data)[1])
  repla[unique(queryHits(op))] <- entropy
  data <- data.frame(data,repla)
  rnb.set <- load.rnb.set(file.path(file,'rnbSet.zip'))
  methData <- meth(rnb.set)
  anno_set <- annotation(rnb.set)
  anno_set <- GRanges(Rle(anno_set$Chromosome),IRanges(start=anno_set$Start,end=anno_set$End))
  chr <- seqnames(anno_set)[1]
  min <- min(start(anno_set))
  max <- max(end(anno_set))
  island <- annotation(rnb.set,type='cpgislands')
  island <- island[which.max(island$CpG),]
  start_island <- island$Start
  end_island <- island$End
  op <- findOverlaps(anno,anno_set)
  methData <- methData[subjectHits(op),]
  data <- data.frame(methData,data)
  colnames(data) <- c('Methylation','FDRP','qFDRP','PDR','Position','MHL','Epipolymorphism','Entropy')
  data <- as.data.frame(data)
  melted <- melt(data,id='Position')
  colnames(melted)[2:3] <- c('Measure','Value')
  colors <- rnb.getOption('colors.category')
  temp <- colors[2]
  colors[2] <- colors[5]
  colors[5] <- temp
  plot <- ggplot(melted,aes(x=Position,y=Value,color=Measure,fill=Measure))+facet_grid(Measure~.)+geom_point(data=melted[melted$Measure=='Methylation',],color='black',size=1)+geom_point(data=melted[melted$Measure!='Methylation',],size=0.1)+geom_bar(data=melted[melted$Measure!='Methylation',],stat='identity')+theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=20,face='bold'),panel.grid.major = element_line(color='grey80'),panel.grid.minor = element_line(color='grey95'),legend.key = element_rect(fill='white'),legend.position = 'none',strip.text = element_text(size=12))+scale_color_manual(values=colors,guide=FALSE)+scale_fill_manual(values=colors)+
    xlim(min,max)+scale_y_continuous(breaks=c(0,0.5,1))+ggtitle(chr)
  ggsave(paste0(plot.path,paste(chr,min,max,sep="_"),".pdf"),plot,device="pdf",height=11,width=8.5,units="in")
  count <- count +1
}
