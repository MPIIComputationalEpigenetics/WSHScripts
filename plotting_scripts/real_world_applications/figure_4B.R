#####################################################################################################
#' figure_4B.R
#'---------------------------------------------------------------------------------------------------
#' This scripts produces Figure 4B from the manuscript. You need to specify the folders, where you 
#' stored the data matrices for the Ewing data set. You should have produced the score with the
#' scripts in ISH_scripts/real_world_applications/Ewing and then summarized the indiviudal samples to
#' a data matrix with misc/create_data_matrix.R. For this script three inputs are required.
#' 
#' ewing.path: path to the Ewing data matrix and annotation
#' score: the score you'd like to plot (qFDRP per default)
#' plot.path: folder to store the results in
ewing.path <- "path_to_ewing"
score <- "qfdrp"
plot.path <- getwd()
library(RnBeads)
ish.matrix <- read.csv(file.path(ewing.path,paste0(score,".csv.gz")))
ish.matrix <- ish.matrix[,-1]
mean.ish <- rowMeans(ish.matrix,na.rm=T)
load(file.path(ewing.path,paste0(score,"_annotation.RData")))
rnb.load.annotation.from.db("ensembleRegBuildBPall","hg38")
reg.build <- unlist(rnb.get.annotation("ensembleRegBuildBPall","hg38"))
op <- findOverlaps(annotation,reg.build)
reg.anno <- rep(NA,length(mean.ish))
reg.anno[queryHits(op)] <- values(reg.build)$elementType[subjectHits(op)]
toPlot <- data.frame(ISH=mean.ish,Annotation=reg.anno)
plot <- ggplot(toPlot,aes(x=Annotation,y=ISH))+geom_violin(fill="grey50") + geom_boxplot(fill="grey80",color="black",alpha=0.25) + ylim(-0.01,1.01)+
  theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=20,face='bold'),
        panel.grid.major = element_line(color="grey80"),
        legend.key = element_rect(fill='white'),legend.position = 'none')+scale_fill_manual(values=rnb.getOption("colors.category"),na.value="grey80")+ylab(score)
ggsave(file.path(plot.path,"figure_4B.pdf"),plot,device="pdf",height=11,width=8.5,units="in")
