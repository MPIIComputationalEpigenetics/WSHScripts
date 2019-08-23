############# figure_5B.R ######################################################
#' This script reproduces Figure 5B from the manuscript. First, one should execute the
#' tumor_purity_prediction.R script to select the sites consistently linked to tumor purity
#' and then specify them here.
#' @param score.path: path to the results in matrix form (e.g. qfdrp.csv.gz)
#' @param score.annotation: the path to the .RData file containing the genomic annotations of the rows
#'           of score.path
#' @param selected.sites: path to a CSV file where the selected sites used for prediction are stored
#' @param save.heatmap: path to a folder where the final heatmap is to be stored
score.path <- "path_to_score.csv.gz"
score.annotation <- "path_to_score_annotation.RData"
selected.sites <- "selected_sites.csv"
save.heatmap <- "heatmap_folder"

sample.annotation <- "../sample_annotation.csv"

#### SCRIPT ########
library(RnBeads)
library(glmnet)
library(pheatmap)
library(RColorBrewer)
qfdrp <- read.csv(score.path)
qfdrp <- qfdrp[,-1]
qfdrp <- qfdrp[,unlist(lapply(colnames(qfdrp),function(x)grepl("EwS_T",x)))]
ffpe <- c("EwS_T122","EwS_T123","EwS_T124",
"EwS_T125","EwS_T126","EwS_T127","EwS_T129","EwS_T130","EwS_T131","EwS_T132","EwS_T133")
qfdrp <- qfdrp[,!(colnames(qfdrp)%in%ffpe)]

ages <- read.csv(sample.annotation)
ages <- ages[unlist(lapply(as.character(ages$sample_id),function(x)grepl("EwS_T",x))),]
ages <- ages[!(ages$sample_id%in%ffpe),]
ages <- ages[ages$sample_id!="EwS_T123",]
ages <- ages[match(colnames(qfdrp),ages$sample_id),]
row.names(ages) <- ages$sample_id
tp <- ages$tumor_purity
tp[tp %in% "not quantified"] <- NA
ages$tumor_purity <- as.numeric(as.character(tp))

load(score.annotation)
sel.sites <- read.csv(selected.sites)
sel.anno <- makeGRangesFromDataFrame(sel.sites)
op <- findOverlaps(sel.anno,annotation,type="equal")

dat <- qfdrp[subjectHits(op),]
res <- ages$tumor_purity
dat <- dat[,!is.na(res)]
res <- res[!is.na(res)]
dat <- as.data.frame(t(dat))
lin.mod <- glm(res~.,data=dat)

all.dat <- qfdrp[subjectHits(op),]
pred.purity <- predict(lin.mod,as.data.frame(t(all.dat)))
ages$predicted_purity <- pred.purity

breaksList <- seq(0,1,by=0.01)
anno.cols <- list("LUMP"=c("#223a88ff","#753a88ff","#b82a8aff","#cc2b5eff"),"tumor_purity"=c("#0284b0ff","#02aab0ff","#00cdacff","#00cd67ff"),"predicted_purity"=c("#0284b0ff","#02aab0ff","#00cdacff","#00cd67ff"))
png(file.path(save.heatmap,"_selected_sites_heatmap.png"),height=160,width=160,unit="mm",res=300)
pheatmap(t(all.dat),breaks=breaksList,annotation_row=ages[,c("tumor_purity","predicted_purity")],annotation_colors=anno.cols,
	show_rownames=F,show_colnames=F,
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))
dev.off()

