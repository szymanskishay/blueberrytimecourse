
ps.course<-readRDS("pscourse.rds")
Tissue_for_SE <- function(phyloseq, expression, outputheader){
  require(phyloseq)
  require(seqtime)
  
  ps.tissue<-do.call("subset_samples", list(quote(phyloseq), substitute(expression)))
  ps.tissue.otu<-otu_table(ps.tissue)
  ps.tissue.taxa<-tax_table(ps.tissue)
  ps.tissue.f = filterTaxonMatrix(ps.tissue.otu, minocc=length(colnames(otu_table(ps.tissue)))*0.1, keepSum=TRUE, return.filtered.indices = TRUE)
  ps.tissue.otu.f = ps.tissue.f$mat
  ps.tissue.taxa.f = ps.tissue.taxa[setdiff(1:nrow(ps.tissue.taxa), ps.tissue.f$filtered.indices),]
  ps.tissue.taxa.f = ps.tissue.taxa.f[,1:7]
  filteredsumtax=c("Filtered_K", "Filtered_P","Filtered_C","Filtered_O", "Filtered_F","Filtered_G","Filtered_S")
  taxa.f = rbind(ps.tissue.taxa.f, filteredsumtax)
  rownames(taxa.f)[nrow(taxa.f)]="OTU_0"
  rownames(ps.tissue.otu.f)[nrow(ps.tissue.otu.f)]="OTU_0"
  updatedotus=otu_table(ps.tissue.otu.f, taxa_are_rows = TRUE)
  updatedtaxa=tax_table(taxa.f)
  otu_out<-t(otu_table(updatedotus)[which(rowSums(otu_table(updatedotus))>=10),])
  outfile_otu<-paste("Generated Tables/SpiecEasi/", outputheader, "OTUs.csv", sep="_")
  outfile_tax<-paste("Generated Tables/SpiecEasi/", outputheader, "tax.csv", sep="_")
  write.csv(otu_out, file=outfile_otu)
  write.csv(as.matrix(updatedtaxa), file=outfile_tax)
}
for(tissue in c("Bud", "Flower", "Green", "Coloring", "Blue", "Leaf", "Stem", "PLP20", "PLP21", "SWMREC")){
  Tissue_for_SE(ps.course, Tissue == tissue, paste0(tissue, "_", sep=""))  
}
