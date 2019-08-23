#####################################################################################################
#' figure_4D.R
#'---------------------------------------------------------------------------------------------------
#' This scripts produces Figure 4G from the manuscript. You need to specify the folders, where you 
#' stored the data matrix for the Ewing data set. You should have produced the score with the scripts
#' in ISH_scripts/real_world_applications/Ewing and then summarized the indiviudal samples to a data
#' matrix with misc/create_data_matrix.R. Furthermore, GO and LOLA enrichment results are exported
#' as csv files. For this script four inputs are required.
#' 
#' ewing.path: path to the Ewing data matrix and annotation
#' score: the score you'd like to plot (qFDRP per default)
#' lola.db.path: path to store the LOLA database object
#' plot.path: folder to store the results in
ewing.path <- "path_to_ewing"
score <- "qfdrp"
lola.db.path <- getwd()
plot.path <- getwd()
library(RnBeads)
library(GOstats)
library(LOLA)
qfdrp.ewing <- read.csv(file.path(ewing.path,paste0(score,".csv.gz")))
qfdrp.ewing <- qfdrp.ewing[,-1]
ffpe.samples <- c('EwS_T122','EwS_T123','EwS_T124','EwS_T125','EwS_T126','EwS_T127','EwS_T129',
'EwS_T130','EwS_T131','EwS_T132','EwS_T133')
qfdrp.ewing <- qfdrp.ewing[,-which(colnames(qfdrp.ewing)%in%ffpe.samples)]
load(file.path(ewing.path,paste0(score,"_annotation.RData")))
qfdrp <- data.frame(qfdrp.ewing)

reg.build <- unlist(rnb.get.annotation("genes","hg38"))
op <- findOverlaps(annotation,reg.build)
to.agg <- qfdrp[queryHits(op),]
agg.by <- values(reg.build)$entrezID[subjectHits(op)]
longer <- lengths(strsplit(agg.by,";"))>0
agg.by[longer] <- unlist(lapply(strsplit(agg.by,";"),function(x)x[1]))[longer]
to.agg <- aggregate(to.agg,by=list(agg.by),mean,na.rm=T)
name <- to.agg[,1]
to.agg <- to.agg[,-1]
row.names(to.agg) <- name

trait <- read.csv("sample_annotation.csv")
trait <- trait[match(colnames(to.agg),trait$sample_id),]

trait <- trait$disease_type_detailed
trait <- trait[!is.na(trait)]
ind.g1 <- which(grepl("MSC",trait))
ind.g2 <- which(trait %in% c("Ewing_Tumor"))

p.vals <- limmaP(as.matrix(to.agg),inds.g1=ind.g1,inds.g2=ind.g2)
p.vals.adj.fdr <- p.adjust(p.vals,method="fdr")

means.g1 <- rowMeans(to.agg[,ind.g1],na.rm=T)
means.g2 <- rowMeans(to.agg[,ind.g2],na.rm=T)

universe.ids <- row.names(to.agg)
universe <- reg.build[unique(subjectHits(op))]
db.dir <- lola.db.path
dir <- downloadLolaDbs(db.dir,"LOLACore")
lola.db <- loadRegionDB(dir$hg38)

gene.ids.hyper <- row.names(to.agg)[which(p.vals.adj.fdr<0.01&means.g1>means.g2)]
params <- new("GOHyperGParams",annotation="org.Hs.eg.db",geneIds = gene.ids.hyper, universeGeneIds = universe.ids, ontology = "BP",conditional = TRUE, testDirection = "over")
test.result.hyper <- hyperGTest(params)
res.table.hyper <- summary(test.result.hyper)
res.table.hyper$p.val.adj.fdr <- p.adjust(res.table.hyper$Pvalue,method="fdr",n=length(test.result.hyper@pvalue.order))
write.csv(res.table.hyper,file.path(plot.path,"limma_GO_hyper.csv"),row.names=F)

gene.ids.hypo <- row.names(to.agg)[which(p.vals.adj.fdr<0.01&means.g1<means.g2)]
params <- new("GOHyperGParams",annotation="org.Hs.eg.db",geneIds = gene.ids.hypo, universeGeneIds = universe.ids, ontology = "BP",conditional = TRUE, testDirection = "over")
test.result.hypo <- hyperGTest(params)
res.table.hypo <- summary(test.result.hypo)
res.table.hypo$p.val.adj.fdr <- p.adjust(res.table.hypo$Pvalue,method="fdr",n=length(test.result.hypo@pvalue.order))
write.csv(res.table.hypo,file.path(plot.path,"limma_GO_hypo.csv"),row.names=F)

user.set.hyper <- universe[which(p.vals.adj.fdr<0.01&means.g1>means.g2)]
lola.res.hyper <- runLOLA(userSets=user.set.hyper,universe,lola.db)
write.csv(lola.res.hyper,file.path(plot.path,"limma_LOLA_hyper.csv"))

user.set.hypo <- universe[which(p.vals.adj.fdr<0.01&means.g1<means.g2)]
lola.res.hypo <- runLOLA(userSets=user.set.hypo,universe,lola.db)
write.csv(lola.res.hypo,file.path(plot.path,"limma_LOLA_hyper.csv"))

plot <- lolaBarPlot(lolaDb=lola.db,lolaRes=lola.res.hyper,maxTerms=15)+theme_bw()
ggsave(file.path(plot.path,"lolaBarPlot_hyper_limma_001.pdf"),device="pdf",height=11,width=8.5,unit="in")
