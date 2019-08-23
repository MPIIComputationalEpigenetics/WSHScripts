#####################################################################################################
#' figure_4C.R
#'---------------------------------------------------------------------------------------------------
#' This scripts produces Figure 4C from the manuscript. You need to specify the folders, where you 
#' stored the data matrix for the Ewing data set. You should have produced the scores with the scripts
#' in ISH_scripts/real_world_applications/Ewing and then summarized the individual samples to a data
#' matrix with misc/create_data_matrix.R. For this script three inputs are required.
#' 
#' ewing.path: path to the Ewing data matrix and annotation
#' score: the score you'd like to plot (qFDRP per default)
#' plot.path: folder to store the results in
ewing.path <- "path_to_ewing"
score <- "qfdrp"
plot.path <- getwd()

library(RnBeads)
library(GOstats)

volcano.plot <- function(df2p,rank.cutoff,rerank=F){
	cn.d <- "mean.diff"
	cn.pa <- "p.vals.adj.fdr"
	cn.pl <- "p.vals.adj.fdr"

	if(rerank){
		df2p$rank <- rank(df2p[,cn.pa])
	}

	dont.plot.p.val <- all(is.na(df2p[,cn.pa]))
		
	pp <- ggplot(df2p) + aes_string(cn.d, paste0("-log10(",cn.pl,")")) +
		geom_point(color=ifelse(df2p[,cn.pa]<rank.cutoff,ifelse(df2p[,cn.d]>0,"#a65500ff","#4d2989ff"),"black"))
	
	pp <- pp+theme_bw()+theme(text=element_text(size=15))+xlab("Difference")+ylab("-log10(FDR-adjusted p-value)")
	return(pp)
}

combinedRanking.tab <- function(tt,rerank=FALSE){
	rank.mat <- c()
	for (i in 1:ncol(tt)){
		rrs <- rank(tt[,i],na.last="keep",ties.method="min")
		if (!all(is.na(rrs))) {
			rank.mat <- cbind(rank.mat,rrs)
		}
	}
	if (is.null(rank.mat)){
		logger.warning("Could not compute combined ranking: To few non-NA columns specified")
		return(rep(NA,nrow(tt)))
	}
	res <- rowMaxs(rank.mat,na.rm=FALSE)
	res[res==-Inf] <- NA
	if (rerank) res <- rank(res,na.last="keep",ties.method="min")
	return(res)
}

auto.select.rank.cut <- function(p,r,alpha=0.1){
	res <- 0
	lp <- length(p)
	j <- 1L:lp
	oa <- order(r)
	od <- rev(oa)
	p.oa <- p[oa]
	cmin.d <- cummin(p[od])
	cquant <- cummax(p.oa)
	inds.better.than.tail <- which(cquant<cmin.d)
	if (length(inds.better.than.tail) > 0){
		L <- max(inds.better.than.tail)
		res <- r[oa][L]
	}
	return(res)
}


eps <- 0.01

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
agg.by <- paste0(seqnames(reg.build[subjectHits(op),]),"_",start(reg.build[subjectHits(op),]))
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

to.agg <- na.omit(to.agg)

means.g1 <- rowMeans(to.agg[,ind.g1],na.rm=T)
means.g2 <- rowMeans(to.agg[,ind.g2],na.rm=T)


p.vals <- limmaP(as.matrix(to.agg),inds.g1=ind.g1,inds.g2=ind.g2)
p.vals.adj.fdr <- p.adjust(p.vals,method="fdr")

comp.frame <- data.frame(p.vals=p.vals,p.vals.adj.fdr=p.vals.adj.fdr,mean.diff=means.g1-means.g2,log.quot=log2((means.g1+eps)/(means.g2+eps)))
na.frame <- unlist(apply(comp.frame,1,function(x)any(is.na(x))))
to.agg <- to.agg[!na.frame,]
means.g1 <- means.g1[!na.frame]
means.g2 <- means.g2[!na.frame]
comp.frame <- na.omit(comp.frame)

dm4ranking <- cbind(-abs(comp.frame$mean.diff),-abs(comp.frame$log.quot),comp.frame$p.vals)

ranks <- combinedRanking.tab(dm4ranking)
comp.frame$rank <- ranks

select.cut <- auto.select.rank.cut(comp.frame$p.vals.adj.fdr,comp.frame$rank)

ranks <- combinedRanking.tab(dm4ranking)
comp.frame$rank <- ranks

select.cut <- auto.select.rank.cut(comp.frame$p.vals.adj.fdr,comp.frame$rank)

plot <- volcano.plot(comp.frame,rank.cutoff=0.01,rerank=F)
ggsave(file.path(plot.path,paste0("figure_4C.png")),plot,height=11,width=8.5,unit="in",dpi=300)

