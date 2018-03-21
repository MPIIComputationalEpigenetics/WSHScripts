#######################################################################################################
#' coverage_dependence.R
#-----------------------------------------------------------------------------------------------------
#' This scripts produces a plot similar to Figure 1B to visualize the overall dependence of the scores
#' on sequencing coverage. You should specify the folder, where the results of performing the 
#' simulation experiment are stored as an argument and the folder where to store the results
#' as another.
#' folder.path: path to the results_pipeline folder for the simulation experiment
#' plot.path: folder to store the results
#' 
#' A single plot is the result.
folder.path <- "path_to_results_pipeline"
plot.path <- getwd()
mround <- function(x,base){
        base*round(x/base)
}

library(RnBeads)
folders <- list.dirs(folder.path,full.names = TRUE,recursive=FALSE)
data <- data.frame()
for(file in folders){
  fdrp <- tryCatch(read.csv(file.path(file,'FDRP','FDRP.csv')),error=function(e){c()})
  fdrp <- mean(unlist(fdrp),na.rm=TRUE)
  qfdrp <- tryCatch(read.csv(file.path(file,'qFDRP','qFDRP.csv')),error=function(e){c()})
  qfdrp <- mean(as.numeric(unlist(qfdrp)),na.rm=TRUE)
  pdr <- tryCatch(read.csv(file.path(file,'PDR','PDR.csv')),error=function(e){c()})
  pdr <- mean(as.numeric(unlist(pdr)),na.rm=TRUE)
  mhl <- tryCatch(read.table(file.path(file,'mhl.txt'),sep='\t',skip = 1),error=function(e){c()})
  mhl <- mean(as.numeric(unlist(mhl[,2])),na.rm=TRUE)
  epipoly <- tryCatch(read.csv(file.path(file,'epipoly.csv')),error=function(e){c()})
  epipoly <- mean(as.numeric(unlist(epipoly$Epipolymorphism)),na.rm=TRUE)
  entropy <- tryCatch(read.csv(file.path(file,'entropy.csv')),error=function(e){c()})
  entropy <- mean(as.numeric(unlist(entropy$Entropy)),na.rm=TRUE)
  rnb.set <- tryCatch(load.rnb.set(file.path(file,'rnbSet.zip')),error=function(e){e})
  if(inherits(rnb.set,"error")) next
  cov <- round(mean(covg(rnb.set)),2)
  data <- rbind(data,c(fdrp,qfdrp,pdr,mhl,epipoly,entropy,cov))
}
colnames(data) <- c('FDRP','qFDRP','PDR','MHL','Epipolymorphism','Entropy','Coverage')
data$Coverage <- mround(data$Coverage,5)
data <- aggregate(data,by=list(data$Coverage),mean,na.rm=TRUE)
data <- data[,-1]
melted <- melt(data,id='Coverage')
colnames(melted)[2:3] <- c('Measure','Value')
colors <- rnb.getOption('colors.category')
temp <- colors[2]
colors[2] <- "#76bf23ff"
colors[1] <- "#00806fff"
colors[5] <- temp
plot <- ggplot(melted,aes(x=Coverage,y=Value,color=Measure,shape=Measure))+geom_point(size=3)+geom_line(size=1.2)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=30,face='bold'),
        panel.grid.major = element_line(color='grey80'),panel.grid.minor = element_line(color='grey95'),
        legend.key = element_rect(fill='white'),legend.position = 'top',legend.text = element_text(size=20))+
  scale_color_manual(values=colors)
ggsave(paste0(plot.path,"coverage_dependence.pdf"),plot,device="pdf",height=11,width=8.5,units="in")
