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
ewing.path <- "path_to_Ewing_results"
score <- "qfdrp"
plot.path <- getwd()
library(RnBeads)
library(pheatmap)
library(RColorBrewer)
ish.matrix <- read.csv(file.path(ewing.path,paste0(score,".csv.gz")))
ish.matrix <- ish.matrix[,-1]
ish.matrix <- ish.matrix[,unlist(lapply(colnames(ish.matrix),function(x)grepl("EwS_T",x)))]
ffpe <- c("EwS_T122","EwS_T123","EwS_T124",
"EwS_T125","EwS_T126","EwS_T127","EwS_T129","EwS_T130","EwS_T131","EwS_T132","EwS_T133")
ish.matrix <- ish.matrix[,!(colnames(ish.matrix) %in% ffpe)]
load(file.path(ewing.path,paste0(score,"_annotation.RData")))
reg.build <- unlist(rnb.get.annotation("genes","hg38"))
op <- findOverlaps(annotation,reg.build)
to.agg <- ish.matrix[queryHits(op),]
agg.by <- values(reg.build)$symbol[subjectHits(op)]
to.agg <- aggregate(to.agg,by=list(agg.by),median,na.rm=T)
name <- to.agg[,1]
to.agg <- to.agg[,-1]
row.names(to.agg) <- name
to.agg <- na.omit((to.agg))
vars <- apply(to.agg, 1, var, na.rm=T)
ordered <- order(vars,decreasing = T)
most.var <- to.agg[ordered,]
most.var <- most.var[1:1000,]

trait <- read.csv("/DEEP_fhgfs/projects/mscherer/data/RRBS/Ewing_complete/src/first_analysis/sample_annotation.csv")
trait <- trait[match(colnames(most.var),trait$sample_id),]

class <- data.frame(Sex=trait$patient_sex)
ann.cols <- list(Sex=c(f="#f3449e",m="#2da9d9"))
row.names(class) <- trait$sample_id
breaksList <- seq(0,1,by=0.01)

png(file.path(plot.path,"figure_4C.png"),height=160,width=160,unit="mm",res=300)
pheatmap(as.matrix(most.var),breaks=breaksList,annotation_col=class,show_rownames=F,show_colnames=F,annotation_colors=ann.cols,
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))
dev.off()

clusti <- hclust(dist(t(as.matrix(most.var))))
count(cutree(clusti,2))

