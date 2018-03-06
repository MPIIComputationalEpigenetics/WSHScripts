library(RnBeads)
islands <- unlist(rnb.get.annotation('cpgislands','hg38')[-c(22,23,24)])
i <- 1
data <- matrix(nrow = 100,ncol=6)
while(i<=100){
  selected <- sample(islands,1)
  chr <- as.character(seqnames(selected))
  start_island <- start(selected)
  end_island <- end(selected)
  size <- end_island-start_island
  remaining.size <- 50000-size
  add <- round(remaining.size/2,0)
  start <- start_island-add
  if(start<0){
    i <- i-1
    next
  }
  end <- end_island+add
  if(end>seqlengths(islands)[chr]){
    i <- i-1
    next
  }
  reads <- sample(20000:50000,1)
  id <- paste(chr,start_island,end_island,reads,sep='_')
  line <- c(sample_name=id,chr=chr,start=start,end=end,reads=reads,library='CGI')
  data[i,] <- line
  i <- i+1
}
colnames(data) <- c('sample_name','chr','start','end','reads','library')
write.csv(data,'sample_annotation_CGI.csv',quote=FALSE,row.names=FALSE)
