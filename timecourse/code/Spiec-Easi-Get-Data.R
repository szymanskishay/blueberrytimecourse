library(SpiecEasi)
library(igraph)

#This uses the outputs of the script SpiecEasi
#directory <- 'path/to/timecourse'
#setwd(directory)

getNetworkInfo <- function(var){
  set.seed(0413)
  varRDS<-paste0("Generated Tables/SpiecEasi/RDSobj/", var, "-spieceasi-network-mb-rf.rds", sep="")
  varOTU<-paste0("Generated Tables/SpiecEasi/", var, "_refilter_OTUs.csv", sep="")
  varTax<-paste0("Generated Tables/SpiecEasi/", var, "_refilter_tax.csv", sep="")
  var.net<-readRDS(varRDS)
  var.OTU<-as.matrix(read.csv(varOTU, row.names = 1))
  var.Tax<-as.matrix(read.csv(varTax, row.names = 1))
  
  var.se.mat<-symBeta(getOptBeta(var.net), mode='maxabs')
  var.mat<-as.matrix(var.se.mat)
  rownames(var.mat)<-colnames(var.OTU)
  colnames(var.mat)<-colnames(var.OTU)
  
  net<-graph.adjacency(var.mat, mode='undirected', weighted = TRUE, diag=FALSE)
  V(net)$name <- colnames(var.OTU)
  #betweenness (weights converted to distances)
  if(max(abs(E(net)$weight)) < 1){
    net.dist<-net
    weights.dist <- 1-abs(E(net.dist)$weight)
    E(net.dist)$weight <- weights.dist
  }
  #Weigthed degree 
  net.abs <- net
  E(net.abs)$weight <- abs(E(net.abs)$weight)
  #Calculate centrality metrics
  net.alpha <- alpha.centrality(net)
  net.strength <- strength(net.abs)
  net.bet <- betweenness(net.dist, v=V(net.dist))
  net.cn <- closeness(net.dist)
  net.pr <- page_rank(net)$vector
  net.hs <- hub_score(net)$vector
  #summary of centrality metrics
  central_summary <- as.data.frame(net.alpha)
  colnames(central_summary) <- ("Alpha_centrality")
  rownames(central_summary) <- colnames(var.OTU)
  central_summary$Weighted_Vertex_Degree <- net.strength
  central_summary$Betweenness_centrality <- net.bet
  central_summary$Closeness_centrality <- net.cn
  central_summary$Page_rank <- net.pr
  central_summary$Hub_score <- net.hs
  
  cs_outfile <- paste0("tables/central_summary_",var,".csv", sep="")
  write.csv(central_summary, file=cs_outfile)
  metrics <- central_summary
  
  #cluster nodes into modules
  #use distance conversion to address negative weights 
  net.walktrap <- walktrap.community(net.dist)
  net.multilevel <- multilevel.community(net.dist)
  wt<-net.walktrap$membership
  ml<-net.multilevel$membership
  modules = data.frame(as_ids(V(net)), wt, ml)
  metrics$walktrap<-net.walktrap$membership
  metrics$multilevel<-net.multilevel$membership
  write.csv(metrics, file=paste0("Generated Tables/SpiecEasi/", var, "_metrics_with_clusters.csv"))
  write.csv(modules, file=paste0("Generated Tables/SpiecEasi/",var,"_cluster_nodes.csv"))
  
  node.similarity.jac <- similarity(net, vids=V(net), mode="all", method="jaccard")
  colnames(node.similarity.jac)<-rownames(node.similarity.jac)<-colnames(var.OTU)
  
  write.csv(node.similarity.jac, file=paste0("Generated Tables/SpiecEasi/", var,"node_similarity_jaccard.csv"))
  
  net.AP <- as_ids(articulation.points(net))
  write.csv(net.AP, file=paste0("Generated Tables/SpiecEasi/", var,"_articulation_points.csv"))
}

varlist<-c("PLP20","PLP21","SWMREC","Blue","Coloring","Green","Flower","Bud","Leaf","Stem")

for(i in varlist){
  getNetworkInfo(i)
}