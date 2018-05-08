library(RnBeads)
anno <- rnb.get.annotation('CpG','hg38')
lengths <- unlist(lapply(anno,function(x){max(end(x))}))
starts <- unlist(lapply(anno,function(x){min(start(x))}))
names(starts) <- names(lengths)
data <- matrix(nrow=12000,ncol=6)
lengths <- lengths[-c(22,23,24)]
complete <- sum(as.numeric(lengths))
probs <- lengths/complete
i <- 0
while(i<12000){
  print(i)
  chrom <- sample(names(lengths),1,prob=probs)
  start <- sample((starts[chrom]+100000):lengths[chrom]-100000,1)
  if((start+50000 > lengths[chrom]) || (is.na(start))){
    if(i>=12){
      i <- i-12
    }else{
      i <- 0
    }
    print(i)
    next
  }
  end <- start+50000
  id <- paste(chrom,start,sep='_')
  for(j in seq(40,150,by=10)){
    rid <- paste0(id,"_",j)
    num_reads <- round(250000*(5/j),0)
    line <- c(ID=rid,Chr=chrom,Start=start,End=end,Reads=num_reads,Length=j)
    data[i+(j/10)-3,] <- line
  }
  i <- i+12
}
data <- data.frame(data,rep('Read_Length',12000))
colnames(data) <- c('sample_name','chr','start','end','reads','length','library')
write.csv(data,'sample_annotation_read_length.csv',quote=FALSE,row.names=FALSE)
