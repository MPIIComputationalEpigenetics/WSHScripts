#####################################################################################################
#' figure_4A.R
#'---------------------------------------------------------------------------------------------------
#' This scripts produces Figure 4A from the manuscript. You need to specify the folders, where you 
#' stored the data matrices for both real world applications (Kiel Cohort and Ewing). You should
#' have produced the score with the scripts in ISH_scripts/real_world_applications and then
#' summarized the indiviudal samples to a data matrix with misc/create_data_matrix.R. For this
#' script four inputs are required.
#' 
#' ewing.path: path to the Ewing data matrix and annotation
#' kiel.path: path to the Kiel Cohort data matrix and annotation
#' score: the score you'd like to plot (qFDRP per default)
#' plot.path: folder to store the results in
ewing.path <- "path_to_ewing"
kiel.path <- "path_to_kiel"
score <- "qfdrp"
plot.path <- getwd()
library(RnBeads)
qfdrp.ewing <- read.csv(file.path(ewing.path,paste0(score,".csv.gz")))
qfdrp.kiel <- read.csv(file.path(kiel.path,paste0(score,".csv.gz")))
qfdrp.ewing <- qfdrp.ewing[,-1]
qfdrp.kiel <- qfdrp.kiel[,-1]
load(file.path(ewing.path,paste0(score,"_annotation.RData")))
anno.ewing <- annotation
load(file.path(kiel.path,paste0(score,"_annotation.RData")))
anno.kiel <- annotation
op <- findOverlaps(anno.ewing,anno.kiel)
qfdrp.ewing <- qfdrp.ewing[queryHits(op),]
qfdrp.kiel <- qfdrp.kiel[subjectHits(op),]

qfdrp <- data.frame(qfdrp.ewing,qfdrp.kiel)

trait <- read.csv("sample_annotation.csv")
trait <- trait[match(colnames(qfdrp.ewing),trait$sample_id),]
trait <- trait$disease_type_detailed

means.healthy <- rowMeans(qfdrp[,(ncol(qfdrp.ewing)+1):ncol(qfdrp)],na.rm=T)
means.tissue <- rowMeans(qfdrp[,which(trait%in%"Ewing_Tumor")],na.rm=T)
means.cl <- rowMeans(qfdrp[,which(trait%in%"Ewing_Cell_Line")],na.rm=T)
means.eMSC <- rowMeans(qfdrp[,which(trait%in%"MSC_Ewing")],na.rm=T)
means.MSC <- rowMeans(qfdrp[,which(trait%in%"MSC_normal")],na.rm=T)

to.plot <- data.frame(healthy=means.healthy,tissue=means.tissue,CL=means.cl,eMSC=means.eMSC,MSC=means.MSC)
to.plot <- melt(to.plot)
colnames(to.plot) <- c("Sample_Type","Value")

plot <- ggplot(to.plot,aes(x=Sample_Type,y=Value,fill=Sample_Type))+geom_violin()+geom_boxplot(fill="grey80",color="black",alpha=0.25)+theme_bw()+theme(text=element_text(face="bold",size=15),axis.text.x=element_text(angle=45,hjust=1))+scale_fill_manual(values=rnb.getOption("colors.category"))+ylim(-0.01,1.01)
ggsave(file.path(plot.path,"figure_4A.pdf"),plot,device="pdf",height=11,width=8.5,units="in")
