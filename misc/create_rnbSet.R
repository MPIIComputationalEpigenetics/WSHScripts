suppressPackageStartupMessages(library(RnBeads))

cmds <- commandArgs(trailingOnly=TRUE)
folder <- cmds[1]

data.source <- c(file.path(folder,'covs'),file.path(folder,'sample_annotation.csv'),'filename')
rnb.options(identifiers.column='sample_id',
	import.bed.style='bismarkCov',
	assembly='hg38')
rnb.set <- rnb.execute.import(data.source=data.source,data.type='bs.bed.dir')
save.rnb.set(rnb.set,file.path(folder,'rnbSet'))
anno <- annotation(rnb.set)
chr <- as.character(anno$Chromosome[1])
min <- min(anno$Start)
max <- min+50000
range <- GRanges(paste0(chr,':',min,'-',max))
cpgs <- rnb.get.annotation('CpG','hg38')[[chr]]
op <- findOverlaps(range,cpgs)
roi <- cpgs[subjectHits(op)]
roi <- data.frame(Chromosome=seqnames(roi),Start=start(roi),End=end(roi))
roi <- roi[seq(1,dim(roi)[1],by=2),]
write.table(roi,file.path(folder,"roi.bed"),sep="\t",row.names=FALSE,quote=FALSE)
