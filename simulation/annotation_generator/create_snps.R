library(RnBeads)
anno <- rnb.get.annotation('CpG','hg38')
lengths <- unlist(lapply(anno,function(x){max(end(x))}))
starts <- unlist(lapply(anno,function(x){min(start(x))}))
names(starts) <- names(lengths)
data <- matrix(nrow=10000,ncol=5)
lengths <- lengths[-c(22,23,24)]
complete <- sum(as.numeric(lengths))
probs <- lengths/complete
i <- 0
while(i<10000){
  print(i)
  chrom <- sample(names(lengths),1,prob=probs)
  start <- sample((starts[chrom]+100000):lengths[chrom]-100000,1)
  if((start+50000 > lengths[chrom]) || (is.na(start))){
    if(i>=10){
      i <- i-10
    }else{
      i <- 0
    }
    print(i)
    next
  }
  end <- start+50000
  id <- paste(chrom,start,sep='_')
  for(j in seq(1,10,by=1)){
    rid <- paste0(id,"_",j)
    line <- c(ID=rid,Chr=chrom,Start=start,End=end,SNPs=j)
    data[i+j,] <- line
  }
  i <- i+10
}
data <- data.frame(data,rep('SNPs',10000))
colnames(data) <- c('sample_name','chr','start','end','snps','library')
write.csv(data,'sample_annotation_snps.csv',quote=FALSE,row.names=FALSE)
