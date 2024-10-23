library(SpiecEasi)
set.seed(0413)

getSEnetwork<-function(string){
  prefix="Generated Tables/SpiecEasi/"
  suffix="_OTU.csv"
  fileName=paste0(prefix,string,suffix,sep="")
  outfile=paste0("Generated Tables/SpiecEasi/RDSobj/",string,"-spieceasi-network-mb.rds", sep="")
  tissue.otu<-as.matrix(read.csv(fileName, row.names = 1))
  
  tissue.net <- spiec.easi(tissue.otu, method = 'mb', lambda.min.ratio=1e-5, nlambda=30, pulsar.params=list(rep.num=999))
  
  saveRDS(tissue.net, file=outfile)
  rm(tissue.otu)
  rm(tissue.net)
}


TissueList<-c("Bud","Flower","Green","Coloring","Blue", "Stem", "Leaf", "PLP20", "PLP21", "SWMREC")

for(i in TissueList){
  getSEnetwork(i)
}