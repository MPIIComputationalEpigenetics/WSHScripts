library(RnBeads)
anno <- rnb.get.annotation('CpG','hg38')
lengths <- lengths(anno)
starts <- unlist(lapply(anno,function(x){min(start(x))}))
names(starts) <- names(lengths)
data <- matrix(nrow=1000,ncol=5)
lengths <- lengths[-c(22,23,24)]
complete <- sum(lengths)
probs <- lengths/complete
i <- 0
while(i<1000){
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
  #j <- sample(20000:50000,1)
  id <- paste(chrom,start,sep='_')
  for(j in seq(1,10,by=1)){
    rid <- paste0(id,"_",j)
    line <- c(ID=rid,Chr=chrom,Start=start,End=end,Error=j)
    data[i+j,] <- line
  }
  i <- i+10
}
data <- data.frame(data,rep('Errors',1000))
colnames(data) <- c('sample_name','chr','start','end','error','library')
write.csv(data,'sample_annotation_errors.csv',quote=FALSE,row.names=FALSE)
