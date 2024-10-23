library(metagenomeSeq)
ps.course<-readRDS('pscourse.RDS')
set.seed(0413) #seed for reproducability because NMDS being used...
otu_table(ps.course) <- prune_taxa(taxa_sums(ps.course)>10, ps.course) #removes spurious OTUs
#ps.course<- subset_samples(course_true, Tissue%in%c("Bud", "Blue", "Green", "Flower", "Coloring"))

###### Supplemental Figure 1
ps.course.n = phyloseq_to_metagenomeSeq(ps.course) #n for normalized
p_biom_course<-cumNormStat(ps.course.n)
biom_quant_course<-cumNorm(ps.course.n, p=p_biom_course)
normFactors(biom_quant_course)
ps.course.nf <-MRcounts(biom_quant_course, norm=T)
#create physeq object with normalized otu table
otu_table(ps.course) <- otu_table(ps.course.nf, taxa_are_rows = TRUE)
tissue_colors<-c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Black", "Leaf"="darkgreen", "Stem"="Brown")
course.ord = ordinate(ps.course, method ="NMDS", distance="bray") #NMDS to avoid the need for assumptions of homogenity of variance, bray= bray-curtis dissimilarity. Jaccard produces similar results. 
# this function will visualize the ordination along the two axes that account for the most variation
course.pca.SY = plot_ordination(ps.course, course.ord, color="SiteYear") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  scale_color_manual(values=c("PLP20"="red", "SWMREC"="darkgreen", "PLP21" = "blue"))+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1)+
  theme(legend.position = "bottom")+
  labs(tag = "A", title="Site Year")
course.pca.T = plot_ordination(ps.course, course.ord, color="Tissue") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  scale_color_manual(values=tissue_colors)+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1)+
  theme(legend.position = "bottom")+
  labs(tag = "B", title="Tissue")
tiff(filename="Figures/Supplemental_Figure_1.tiff", width = 1000, height = 600, units = "px")
ggarrange(course.pca.SY, course.pca.T)
dev.off()


#####Supplemental Figure 2
ps.course<-readRDS('pscourse.RDS') #reset the ps.course object because it was modified above
ps.course<-subset_samples(ps.course, Tissue%in%c("Bud", "Flower", "Green", "Coloring", "Blue"))
otu_table(ps.course) <- prune_taxa(taxa_sums(ps.course)>10, ps.course) #redo removal of potentially spurrious taxa

ps.course.2020<-subset_samples(ps.course, Year%in%c("2020") & Location%in%c("PLP")) #PLP20
ps.course.2021<-subset_samples(ps.course, Year%in%c("2021")) #PLP21
ps.course.SWMREC<-subset_samples(ps.course, Location%in%c("SWMREC")) #SWMREC
#S2A
ps.course.n = phyloseq_to_metagenomeSeq(ps.course) #n for normalized
p_biom_course<-cumNormStat(ps.course.n)
biom_quant_course<-cumNorm(ps.course.n, p=p_biom_course)
normFactors(biom_quant_course)
ps.course.nf <-MRcounts(biom_quant_course, norm=T)
#create physeq object with normalized otu table
otu_table(ps.course) <- otu_table(ps.course.nf, taxa_are_rows = TRUE)
tissue_colors<-c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Black")
course.ord = ordinate(ps.course, method ="NMDS", distance="bray") #NMDS to avoid the need for assumptions of homogenity of variance, bray= bray-curtis dissimilarity. Jaccard produces similar results. 
# this function will visualize the ordination along the two axes that account for the most variation
course.pca.T2 = plot_ordination(ps.course, course.ord, color="Tissue") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  scale_color_manual(values=tissue_colors)+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1)+
  theme(legend.position = "bottom")+
  labs(tag = "A", title="All Samples")


ps.course.2020.n = phyloseq_to_metagenomeSeq(ps.course.2020) #n for normalized
p_biom_course.2020<-cumNormStat(ps.course.2020.n)
biom_quant_course.2020<-cumNorm(ps.course.2020.n, p=p_biom_course.2020)
normFactors(biom_quant_course.2020)
ps.course.2020.nf <-MRcounts(biom_quant_course.2020, norm=T)
#create physeq object with normalized otu table
otu_table(ps.course.2020) <- otu_table(ps.course.2020.nf, taxa_are_rows = TRUE)

course.2020.ord = ordinate(ps.course.2020, method ="NMDS", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
course.2020.pca = plot_ordination(ps.course.2020, course.2020.ord, color="Tissue") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  scale_color_manual(values=c("Blue"="blue", "Green"="green", "Coloring"="red", "Bud"="orange", "Flower"="Black"))+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1, position="jitter")+
  labs(tag = "B", title="PLP20")
plot(course.2020.pca)


ps.course.2021.n = phyloseq_to_metagenomeSeq(ps.course.2021) #n for normalized
p_biom_course.2021<-cumNormStat(ps.course.2021.n)
biom_quant_course.2021<-cumNorm(ps.course.2021.n, p=p_biom_course.2021)
normFactors(biom_quant_course.2021)
ps.course.2021.nf <-MRcounts(biom_quant_course.2021, norm=T)
#create physeq object with normalized otu table
otu_table(ps.course.2021) <- otu_table(ps.course.2021.nf, taxa_are_rows = TRUE)

course.2021.ord = ordinate(ps.course.2021, method ="NMDS", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
course.2021.pca = plot_ordination(ps.course.2021, course.2021.ord, color="Tissue") + 
  theme_classic()+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  scale_color_manual(values=c("Blue"="blue", "Green"="green", "Coloring"="red", "Bud"="orange", "Flower"="Black"))+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1, position="jitter")+
  labs(tag = "C", title="PLP21")
plot(course.2021.pca)


ps.course.SWMREC.n = phyloseq_to_metagenomeSeq(ps.course.SWMREC) #n for normalized
p_biom_course.SWMREC<-cumNormStat(ps.course.SWMREC.n)
biom_quant_course.SWMREC<-cumNorm(ps.course.SWMREC.n, p=p_biom_course.SWMREC)
normFactors(biom_quant_course.SWMREC)
ps.course.SWMREC.nf <-MRcounts(biom_quant_course.SWMREC, norm=T)
#create physeq object with normalized otu table
otu_table(ps.course.SWMREC) <- otu_table(ps.course.SWMREC.nf, taxa_are_rows = TRUE)

course.SWMREC.ord = ordinate(ps.course.SWMREC, method ="NMDS", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
course.SWMREC.pca = plot_ordination(ps.course.SWMREC, course.SWMREC.ord, color="Tissue") + 
  theme_classic()+
  labs(tag = "D", title="SWMREC")+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18, face="bold"), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  scale_color_manual(values=c("Blue"="blue", "Green"="green", "Coloring"="red", "Bud"="orange", "Flower"="Black"))+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1, position="jitter")
plot(course.SWMREC.pca)

tiff(filename="Figures/Supplemental_Figure_2.tiff", units="px", height = 800, width = 800)
ggarrange(course.pca.T2, course.2020.pca, course.2021.pca, course.SWMREC.pca, common.legend=TRUE, legend="right")

dev.off()
