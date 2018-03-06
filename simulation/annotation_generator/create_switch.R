library(RnBeads)
anno <- rnb.get.annotation('CpG','hg38')
lengths <- lengths(anno)
starts <- unlist(lapply(anno,function(x){min(start(x))}))
names(starts) <- names(lengths)
data <- matrix(nrow=100,ncol=4)
lengths <- lengths[-c(22,23,24)]
complete <- sum(lengths)
probs <- lengths/complete
i <- 0
while(i<100){
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
  #j <- sample(20000:50000,1)
  id <- paste(chrom,start,sep='_')
  line <- c(ID=id,Chr=chrom,Start=start,End=end)
  data[i+1,] <- line
  i <- i+1
}
data <- data.frame(data,rep('SWITCH',100))
colnames(data) <- c('sample_name','chr','start','end','library')
write.csv(data,'sample_annotation_switch.csv',quote=FALSE,row.names=FALSE)
