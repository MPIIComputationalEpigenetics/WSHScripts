#######################################################################################################
#' ROC_erosion.R
#-----------------------------------------------------------------------------------------------------
#' This script produces plots similar to Figure 2D: genome browser-like views of the WSH scores,
#' the methylation levels in individual cell types and aggregated over all cell types. Furthermore,
#' ROC curves are produced representing the classification performance of a particular score or the
#' average methylation level. Three paths need to be given to the script.
#' folder.path: path to the results_pipeline folder for the simulation experiment
#' plot.path: folder to store the results
#' roc.path: folder to which the ROC curves are to be stored
#' plot.type: character, representing if pdf or png images are to be created
#' 
library(reshape2)
library(ggplot2)
library(RnBeads)
library(pROC)
library(PRROC)
folder.path <- "path_to_erosion_results_pipeline"
plot.path <- "path_for_pdfs"
roc.path <- "path_for_rocs"
plot.type <- "pdf"

folders <- list.files(folder.path,full.names = TRUE)
count <- 1
roc.data <- c()
for(file in folders){
  data <- data.frame()
  fdrp <- tryCatch(read.csv(file.path(file,'FDRP','FDRP.csv')),error=function(e)e)
  if(inherits(fdrp,"error")){
	next
  }
  qfdrp <- tryCatch(read.csv(file.path(file,'qFDRP','qFDRP.csv')),error=function(e)e)
  if(inherits(qfdrp,"error")){
	next
  }
  pdr <- tryCatch(read.csv(file.path(file,'PDR','PDR.csv')),error=function(e)e)
  if(inherits(pdr,"error")){
	next
  }
  data <- cbind(as.numeric(unlist(fdrp)),as.numeric(unlist(qfdrp)),as.numeric(unlist(pdr)))
  load(file.path(file,'FDRP','annotation.RData'))
  chr <- paste0("chr",substr(unlist(strsplit(file,"chr"))[2],1,2))
  chr <- gsub("_","",chr)
  anno <- annotation
  sel.rows <- seqnames(anno) %in% chr
  anno <- anno[sel.rows]
  data <- data[as.logical(sel.rows),]
  data <- cbind(data,start(anno))
  mhl <- tryCatch(read.table(file.path(file,'mhl.txt'),sep='\t',skip=1),error=function(e){
	data.frame()
  })
  if(nrow(mhl)>0){
  	anno_mhl <- as.character(mhl[,1])
  	anno_mhl <- strsplit(anno_mhl,'[[:punct:]]')
  	anno_mhl <- as.data.frame(anno_mhl)
  	anno_mhl <- t(anno_mhl)
  	anno_mhl <- GRanges(Rle(anno_mhl[,1]),IRanges(as.numeric(anno_mhl[,2]),as.numeric(anno_mhl[,3])))
  	op <- findOverlaps(anno,anno_mhl)
	repla <- rep(NA,dim(data)[1])
	if(length(op)!=0){
	  	mhl <- as.numeric(mhl[,2])
	  	mhl <- mhl[subjectHits(op)]
	  	mhl <- aggregate(mhl,by=list(queryHits(op)),mean,na.rm=TRUE)
	  	mhl <- mhl[,2]
  		repla[unique(queryHits(op))] <- mhl
	}
  	data <- data.frame(data,repla)
  }else{
	data <- data.frame(data,rep(NA,nrow(data)))
  }
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
  colnames(data) <- c('FDRP','qFDRP','PDR','Position','MHL','Epipolymorphism','Entropy')
  cov.files <- list.files(file.path(file,"covs"),full.names=T)
  ct.count <- 1
  for(ct in cov.files){
	 meth.info <- read.table(ct,sep="\t")
	 anno_set <- GRanges(Rle(meth.info[,1]),IRanges(start=as.numeric(meth.info[,2]),end=as.numeric(meth.info[,3])))
	 methData <- c()
	 for(i in 1:length(anno_set)){
		if(i > length(anno_set)){
			break
		}
		dists <- distance(anno_set[i],anno_set)
		dists[i] <- NA
		is.same <- which(dists==0)
		new.meth <- mean(c(as.numeric(meth.info[i,4]),as.numeric(meth.info[is.same,4])))
		methData <- c(methData,new.meth)
		if(length(is.same)>0){
			meth.info <- meth.info[-is.same,]
			anno_set <- anno_set[-is.same]
		}
	 }
	 min <- min(start(anno_set[seqnames(anno_set) %in% chr]))
	 max <- max(end(anno_set[seqnames(anno_set) %in% chr]))
	 op <- findOverlaps(anno,anno_set)
	 add.data <- rep(NA,nrow(data))
	 add.data[queryHits(op)] <- methData[subjectHits(op)]
	 add.data <- add.data/100
	 data <- data.frame(add.data,data)
     if(!grepl("merged",ct)){
     	colnames(data)[1] <- paste("Methylation_CT_",ct.count)
	 	ct.count <- ct.count+1
	 }else{
		colnames(data)[1] <- "Methylation_Mixture"
	 }
  }
  data <- as.data.frame(data)
  melted <- melt(data,id='Position')
  colnames(melted)[2:3] <- c('Measure','Value')
  colors <- rnb.getOption('colors.category')
  temp <- colors[2]
  colors[2] <- colors[5]
  colors[5] <- temp
  colors <- c(colors,rep("black",ct.count))
  plot <- ggplot(melted,aes(x=Position,y=Value,color=Measure,fill=Measure))+facet_grid(Measure~.)+geom_point(data=melted[grepl('Methylation',melted$Measure),],color='black',size=1)+geom_point(data=melted[!grepl('Methylation',melted$Measure),],size=0.1)+geom_bar(data=melted[!grepl('Methylation',melted$Measure),],stat='identity')+theme(panel.background = element_rect(fill='white',color='black'),text=element_text(size=20,face='bold'),panel.grid=element_blank(),legend.key = element_rect(fill='white'),legend.position = 'none',strip.text = element_text(size=12))+scale_color_manual(values=colors,guide=FALSE)+scale_fill_manual(values=colors)+
    scale_y_continuous(breaks=c(0,0.5,1))+ggtitle(chr)
	dmr.info <- readLines(file.path(file,"erosion_score.txt"))
	is.negative <- dmr.info[1] == "Negative control"
	if(is.negative){
		plot.name <- paste0(plot.path,paste(chr,min,max,"negative",sep="_"),".",plot.type)
	}else{
		plot.name <- paste0(plot.path,paste(chr,min,max,sep="_"),".",plot.type)
	}
 	ggsave(plot.name,plot,device=plot.type,height=11,width=8.5,units="in")
	dmr.start <- as.numeric(dmr.info[4])
	dmr.end <- as.numeric(dmr.info[5])
  	dmr.start <- min(start(anno))+dmr.start
	dmr.end <- min(start(anno))+dmr.end
	in.dmr <- start(anno) >= dmr.start  & start(anno) <= dmr.end
	data <- data[,!(grepl("Methylation_CT",colnames(data)))]
	data <- data[,!(grepl("Position",colnames(data)))]
	p.vals <- apply(data,2,function(x){
		outside <- x[!in.dmr]
		inside <- x[in.dmr]
		p.val <- tryCatch(t.test(inside,outside)$p.value,error=function(e)NULL)
		if(is.null(p.val)){
			return(NA)
		}else if(is.na(p.val)){
			p.val <- 1
		}
		p.val
	})
	roc.data <- rbind(roc.data,c(is.negative,p.vals))
  count <- count +1
}
roc.data <- as.data.frame(roc.data)
colnames(roc.data) <- c("IsNegative","Methylation",'FDRP','qFDRP','PDR','MHL','Epipolymorphism','Entropy')
line.types <- c("solid","solid", "22", "42", "44", "13", "1343")
names(line.types) <- c("Methylation",'FDRP','qFDRP','PDR','MHL','Epipolymorphism','Entropy')
for(score in c("Methylation",'FDRP','qFDRP','PDR','MHL','Epipolymorphism','Entropy')){
	roc.obj <- roc(response=!roc.data$IsNegative,predictor=1-roc.data[,score])
	pdf(paste0(roc.path,"/",score,".pdf"),height=11,width=8.5)
	plot(roc.obj,print.auc=T,lty=line.types[score])
	dev.off()
	nas.score <- is.na(roc.data[,score])
	pr.obj <- pr.curve(scores.class0=1-roc.data[!nas.score,score],weights.class0=!roc.data$IsNegative[!nas.score],curve=T)
	pdf(paste0(roc.path,"/",score,"_PR.pdf"),height=11,width=8.5)
	plot(pr.obj,print.auc=T,lty=line.types[score],legend=F,color="black")
	dev.off()
}
write.csv(roc.data,file.path(roc.path,"roc_data.csv"))
sens.cutoff <- 0.01
for(score in c("Methylation",'FDRP','qFDRP','PDR','MHL','Epipolymorphism','Entropy')){
	fps <- sum(roc.data[,score]<sens.cutoff & roc.data$IsNegative,na.rm=T)
	fns <- sum(roc.data[,score]>sens.cutoff & !(roc.data$IsNegative),na.rm=T)
	tps <- sum(roc.data[,score]<sens.cutoff & !(roc.data$IsNegative),na.rm=T)
	tns <- sum(roc.data[,score]>sens.cutoff & roc.data$IsNegative,na.rm=T)
	conf.matrix <- cbind(c("True State/Predicted State","True","False"),c("True",tps,fps),c("False",fns,tns))
	write.csv(conf.matrix,file.path(roc.path,paste0(score,".csv")))
}

