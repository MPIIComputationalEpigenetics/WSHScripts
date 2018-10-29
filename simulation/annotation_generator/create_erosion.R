library(RnBeads)
anno <- rnb.get.annotation('CpG','hg38')
lengths <- unlist(lapply(anno,function(x){max(end(x))}))
starts <- unlist(lapply(anno,function(x){min(start(x))}))
names(starts) <- names(lengths)
data <- matrix(nrow=1000,ncol=4)
lengths <- lengths[-c(22,23,24)]
complete <- sum(as.numeric(lengths))
probs <- lengths/complete
i <- 0
while(i<1000){
  print(i)
  chrom <- sample(names(lengths),1,prob=probs)
  start <- sample((starts[chrom]+100000):lengths[chrom]-100000,1)
  if((start+50000 > lengths[chrom]) || (is.na(start))){
    if(i>=1){
      i <- i-1
    }else{
      i <- 0
    }
    print(i)
    next
  }
  end <- start+50000
  id <- paste(chrom,start,sep='_')
  line <- c(ID=id,Chr=chrom,Start=start,End=end)
  data[i+1,] <- line
  i <- i+1
}
data <- data.frame(data,rep('EROSION',1000))
colnames(data) <- c('sample_name','chr','start','end','library')
write.csv(data,'sample_annotation_EROSION.csv',quote=FALSE,row.names=FALSE)
