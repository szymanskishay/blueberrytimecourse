if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)
library(SpiecEasi)
library(igraph)

cytoscapePing() #connect to cytoscape
#below script loads SpiecEasi generated networks then prepares them to be used in Cytoscape
doCytoscapes<-function(var){
  varRDS<-paste0("Generated Tables/SpiecEasi/RDSobj/", var, "-spieceasi-network-mb.rds", sep="")
  var.net<-readRDS(varRDS)
  varOTU<-paste0("Generated Tables/SpiecEasi/", var, "_OTUs.csv", sep="")
  var.net<-readRDS(varRDS)
  var.OTU<-as.matrix(read.csv(varOTU, row.names = 1))
  var.se.mat<-symBeta(getOptBeta(var.net), mode='maxabs')
  var.mat<-as.matrix(var.se.mat)
  rownames(var.mat)<-colnames(var.OTU)
  colnames(var.mat)<-colnames(var.OTU)
  net<-graph.adjacency(var.mat, mode='undirected', weighted = TRUE, diag=FALSE)
  createNetworkFromIgraph(net, var)
}

varlist<-c("Bud","Flower","Green","Coloring","Blue", "Stem", "Leaf", "PLP20", "PLP21", "SWMREC")
for(i in varlist){
  doCytoscapes(i)
}

