#####################################################################################################
#' figure_4B.R
#'---------------------------------------------------------------------------------------------------
#' This scripts produces Figure 4C from the manuscript. You need to specify the folders, where you 
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
ish.matrix <- ish.matrix[,unlist(lapply(colnames(ish.matrix),function(x)grepl("EwS_T",x)))]
load(file.path(ewing.path,paste0(score,"_annotation.RData")))
reg.build <- unlist(rnb.get.annotation("genes","hg38"))
op <- findOverlaps(annotation,reg.build)
to.agg <- ish.matrix[queryHits(op),]
agg.by <- paste0(seqnames(reg.build[subjectHits(op),]),"_",start(reg.build[subjectHits(op),]))
to.agg <- aggregate(to.agg,by=list(agg.by),median,na.rm=T)
name <- to.agg[,1]
to.agg <- to.agg[,-1]
row.names(to.agg) <- name
to.agg <- na.omit((to.agg))
vars <- apply(to.agg, 1, var, na.rm=T)
ordered <- order(vars,decreasing = T)
most.var <- to.agg[ordered,]
most.var <- most.var[1:1000,]

trait <- read.csv("sample_annotation.csv")
trait <- trait[match(colnames(most.var),trait$sample_id),]

class <- trait$sample_type
cols <- rep("black",length(class))
new.cols <- c("#7c7b04","#49f9c1")
for(i in 1:length(unique(class))){
	type <- unique(class)[i]
	cols[class%in%type] <- new.cols[i]
}

png(file.path(plot.path,"figure_4C.png"),height=11,width=8.5,units="in",res=300)
heatmap.2(as.matrix(most.var),trace="none",ColSideColors = cols,breaks=seq(0,0.6,by=0.01))
dev.off()
