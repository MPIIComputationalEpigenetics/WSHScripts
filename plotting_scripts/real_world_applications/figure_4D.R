#####################################################################################################
#' figure_4D.R
#'---------------------------------------------------------------------------------------------------
#' This scripts produces Figure 4D from the manuscript. You need to specify the folders, where you 
#' stored the data matrix for the Ewing data set. You should have produced the score with the scripts
#' in ISH_scripts/real_world_applications/Ewing and then summarized the indiviudal samples to a data
#' matrix with misc/create_data_matrix.R. For this script three inputs are required.
#' 
#' ewing.path: path to the Ewing data matrix and annotation
#' score: the score you'd like to plot (qFDRP per default)
#' plot.path: folder to store the results in
ewing.path <- "path_to_ewing"
score <- "qfdrp"
plot.path <- getwd()
library(RnBeads)
library(ggfortify)
qfdrp <- read.csv(file.path(ewing.path,paste0(score,".csv.gz")))
qfdrp <- qfdrp[,-1]
qfdrp <- qfdrp[,unlist(lapply(colnames(qfdrp),function(x)grepl("EwS_T",x)))]
load(file.path(ewing.path,paste0(score,"_annotation.RData")))
reg.build <- unlist(rnb.get.annotation("genes","hg38"))
op <- findOverlaps(annotation,reg.build)
to.agg <- qfdrp[queryHits(op),]
agg.by <- paste0(seqnames(reg.build[subjectHits(op),]),"_",start(reg.build[subjectHits(op),]))
to.agg <- aggregate(to.agg,by=list(agg.by),mean,na.rm=T)
name <- to.agg[,1]
to.agg <- to.agg[,-1]
row.names(to.agg) <- name

trait <- read.csv("sample_annotation.csv")
trait <- trait[match(colnames(to.agg),trait$sample_id),]

class <- trait$sample_type
cols <- rep("black",length(class))
new.cols <- c("#7c7b04","#49f9c1")
for(i in 1:length(unique(class))){
	type <- unique(class)[i]
	cols[class%in%type] <- new.cols[i]
}
pca.model <- prcomp(t(na.omit(to.agg)))
pdf(file.path(plot.path,"figure_4D.pdf"))
autoplot(pca.model,label=F,colour=cols)+theme_bw()+theme(text=element_text(face="bold",size=15))
dev.off()

pdf(file.path(plot.path,"figure_4D_legend.pdf"))
legend.text <- unique(class)
legend.cols <- unique(cols)
plot.new()
legend("topright",legend=legend.text,fill=legend.cols)
dev.off()
