ps.course<-readRDS(file="pscourse.rds") #read in default unaltered files

#ps.course<-subset_samples(course, Tissue%in%c("Bud", "Flower", "Green", "Coloring", "Blue")) this subsets everything other than stem/leaf 
course<-ps.course
course.2020<-subset_samples(course, Year%in%c("2020")) #subsets 2020 to look at location effects restricted to that year
course.PLP20<-subset_samples(course, Year%in%c("2020") & Location%in%c("PLP"))
course.PLP21<-subset_samples(course, Year%in%c("2021"))
course.SWMREC<-subset_samples(course, Location%in%c("SWMREC"))
library(metagenomeSeq)
course.n = phyloseq_to_metagenomeSeq(course) #n for normalized
p_biom_course<-cumNormStat(course.n) #looking to normalize
biom_quant_course<-cumNorm(course.n, p=p_biom_course) #normalizing
normFactors(biom_quant_course)
course.nf <-MRcounts(biom_quant_course, norm=T) #normalizing
otu_course <- as.data.frame(otu_table(course))
taxa_course <- as.data.frame(as.matrix(tax_table(course)))
metadata_course <- as.data.frame(as.matrix(sample_data(course)))

#adonis
model.matrix(~Tissue, strata="Replicate", data=metadata_course)
library(vegan)
adonis.whole.T<-adonis2(t(otu_course) ~ Tissue, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)
adonis.whole.W<-adonis2(t(otu_course) ~ Week, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)
adonis.whole.F<-adonis2(t(otu_course) ~ Fortnight, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)
adonis.whole.TF<-adonis2(t(otu_course) ~ Tissue*Fortnight, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)
adonis.whole.TW<-adonis2(t(otu_course) ~ Tissue*Week, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)
adonis.whole.Y<-adonis2(t(otu_course) ~ Year, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)
adonis.whole.L<-adonis2(t(otu_course) ~ Location, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)
adonis.whole.SY<-adonis2(t(otu_course) ~ SiteYear, strata=metadata_course$Replicate,  data=metadata_course, permutations=9999)

blank<-c("--","--","--","--","--")
holder<-rbind(adonis.whole.T,blank,adonis.whole.W,blank,adonis.whole.F,blank,adonis.whole.TF,blank,adonis.whole.TW,blank,adonis.whole.Y,blank,adonis.whole.L,blank,adonis.whole.SY)
write.table(holder, file="Generated Tables/adoniswhole-alt.txt") #used to make Table 1

# betadispersion
vegan::vegdist(t(otu_course), method="bray") -> dist_otu_course #make base info
permdisp_otu_course.tissue <- betadisper(dist_otu_course, metadata_course$Tissue)#Bdisp for tissue
betadisp.spray.tissue<-anova(permdisp_otu_course.tissue, permutations = 9999)#set to write

permdisp_otu_course.week<- betadisper(dist_otu_course, metadata_course$Week)#Bdisp Week
betadisp.spray.week<-anova(permdisp_otu_course.week, permutations = 9999)#set to write

permdisp_otu_course.Fortnight<- betadisper(dist_otu_course, metadata_course$Fortnight)#Bdisp Fortnight
betadisp.spray.Fortnight<-anova(permdisp_otu_course.Fortnight, permutations = 9999)#set to write


betaholder<-rbind(betadisp.spray.tissue,betadisp.spray.week,betadisp.spray.Fortnight)
row.names(betaholder)<-c("Tissue", "Residuals_T", "Week", "Residuals_W", "Fortnight", "Residuals_F")
write.table(betaholder, file="Generated Tables/betadisp_whole.txt") #this was used to make Table 1
########
#repeat above for each site year. 
course.PLP20.n = phyloseq_to_metagenomeSeq(course.PLP20) #n for normalized
p_biom_course.PLP20<-cumNormStat(course.PLP20.n)
biom_quant_course.PLP20<-cumNorm(course.PLP20.n, p=p_biom_course.PLP20)
normFactors(biom_quant_course.PLP20)
course.PLP20.nf <-MRcounts(biom_quant_course.PLP20, norm=T) 
otu_course.PLP20 <- as.data.frame(otu_table(course.PLP20))
taxa_course.PLP20 <- as.data.frame(as.matrix(tax_table(course.PLP20)))
metadata_course.PLP20 <- as.data.frame(as.matrix(sample_data(course.PLP20)))

#adonis
adonis.PLP20.T<-adonis2(t(otu_course.PLP20) ~ Tissue, strata=metadata_course.PLP20$Replicate,  data=metadata_course.PLP20, permutations=9999)
adonis.PLP20.W<-adonis2(t(otu_course.PLP20) ~ Week, strata=metadata_course.PLP20$Replicate,  data=metadata_course.PLP20, permutations=9999)
adonis.PLP20.F<-adonis2(t(otu_course.PLP20) ~ Fortnight, strata=metadata_course.PLP20$Replicate,  data=metadata_course.PLP20, permutations=9999)

adonis.PLP20.TW<-adonis2(t(otu_course.PLP20) ~ Tissue*Week, strata=metadata_course.PLP20$Replicate,  data=metadata_course.PLP20, permutations=9999)
adonis.PLP20.TF<-adonis2(t(otu_course.PLP20) ~ Tissue*Fortnight, strata=metadata_course.PLP20$Replicate,  data=metadata_course.PLP20, permutations=9999)

# betadispersion
vegan::vegdist(t(otu_course.PLP20), method="bray") -> dist_otu_course.PLP20 #make base info
permdisp_otu_course.PLP20.tissue <- betadisper(dist_otu_course.PLP20, metadata_course.PLP20$Tissue)#Bdisp for tissue
betadisp.course.tissue.PLP20<-anova(permdisp_otu_course.PLP20.tissue, permutations = 9999)#set to write
print(betadisp.course.tissue.PLP20) #used for Table 1

permdisp_otu_course.PLP20.Week <- betadisper(dist_otu_course.PLP20, metadata_course.PLP20$Week)#Bdisp for Week
betadisp.course.Week.PLP20<-anova(permdisp_otu_course.PLP20.Week, permutations = 9999)#set to write
print(betadisp.course.Week.PLP20) #used for Table 1

permdisp_otu_course.PLP20.Fortnight <- betadisper(dist_otu_course.PLP20, metadata_course.PLP20$Fortnight)#Bdisp for Fortnight
betadisp.course.Fortnight.PLP20<-anova(permdisp_otu_course.PLP20.Fortnight, permutations = 9999)#set to write
print(betadisp.course.Fortnight.PLP20) #used for Table 1


#####

course.PLP21.n = phyloseq_to_metagenomeSeq(course.PLP21) #n for normalized
p_biom_course.PLP21<-cumNormStat(course.PLP21.n)
biom_quant_course.PLP21<-cumNorm(course.PLP21.n, p=p_biom_course.PLP21)
normFactors(biom_quant_course.PLP21)
course.PLP21.nf <-MRcounts(biom_quant_course.PLP21, norm=T) 
otu_course.PLP21 <- as.data.frame(otu_table(course.PLP21))
taxa_course.PLP21 <- as.data.frame(as.matrix(tax_table(course.PLP21)))
metadata_course.PLP21 <- as.data.frame(as.matrix(sample_data(course.PLP21)))

#adonis
adonis.PLP21.T<-adonis2(t(otu_course.PLP21) ~ Tissue, strata=metadata_course.PLP21$Replicate,  data=metadata_course.PLP21, permutations=9999)
adonis.PLP21.W<-adonis2(t(otu_course.PLP21) ~ Week, strata=metadata_course.PLP21$Replicate,  data=metadata_course.PLP21, permutations=9999)
adonis.PLP21.F<-adonis2(t(otu_course.PLP21) ~ Fortnight, strata=metadata_course.PLP21$Replicate,  data=metadata_course.PLP21, permutations=9999)

adonis.PLP21.TW<-adonis2(t(otu_course.PLP21) ~ Tissue*Week, strata=metadata_course.PLP21$Replicate,  data=metadata_course.PLP21, permutations=9999)
adonis.PLP21.TF<-adonis2(t(otu_course.PLP21) ~ Tissue*Fortnight, strata=metadata_course.PLP21$Replicate,  data=metadata_course.PLP21, permutations=9999)

# betadispersion
vegan::vegdist(t(otu_course.PLP21), method="bray") -> dist_otu_course.PLP21 #make base info
permdisp_otu_course.PLP21.tissue <- betadisper(dist_otu_course.PLP21, metadata_course.PLP21$Tissue)#Bdisp for tissue
betadisp.course.tissue.PLP21<-anova(permdisp_otu_course.PLP21.tissue, permutations = 9999)#set to write
print(betadisp.course.tissue.PLP21) #used for Table 1

permdisp_otu_course.PLP21.Week <- betadisper(dist_otu_course.PLP21, metadata_course.PLP21$Week)#Bdisp for Week
betadisp.course.Week.PLP21<-anova(permdisp_otu_course.PLP21.Week, permutations = 9999)#set to write
print(betadisp.course.Week.PLP21) #used for Table 1

permdisp_otu_course.PLP21.Fortnight <- betadisper(dist_otu_course.PLP21, metadata_course.PLP21$Fortnight)#Bdisp for Fortnight
betadisp.course.Fortnight.PLP21<-anova(permdisp_otu_course.PLP21.Fortnight, permutations = 9999)#set to write
print(betadisp.course.Fortnight.PLP21) #used for Table 1





course.SWMREC.n = phyloseq_to_metagenomeSeq(course.SWMREC) #n for normalized
p_biom_course.SWMREC<-cumNormStat(course.SWMREC.n)
biom_quant_course.SWMREC<-cumNorm(course.SWMREC.n, p=p_biom_course.SWMREC)
normFactors(biom_quant_course.SWMREC)
course.SWMREC.nf <-MRcounts(biom_quant_course.SWMREC, norm=T) 
otu_course.SWMREC <- as.data.frame(otu_table(course.SWMREC))
taxa_course.SWMREC <- as.data.frame(as.matrix(tax_table(course.SWMREC)))
metadata_course.SWMREC <- as.data.frame(as.matrix(sample_data(course.SWMREC)))

#adonis
adonis.SWMREC.T<-adonis2(t(otu_course.SWMREC) ~ Tissue, strata=metadata_course.SWMREC$Replicate,  data=metadata_course.SWMREC, permutations=9999)
adonis.SWMREC.W<-adonis2(t(otu_course.SWMREC) ~ Week, strata=metadata_course.SWMREC$Replicate,  data=metadata_course.SWMREC, permutations=9999)
adonis.SWMREC.F<-adonis2(t(otu_course.SWMREC) ~ Fortnight, strata=metadata_course.SWMREC$Replicate,  data=metadata_course.SWMREC, permutations=9999)

adonis.SWMREC.TW<-adonis2(t(otu_course.SWMREC) ~ Tissue*Week, strata=metadata_course.SWMREC$Replicate,  data=metadata_course.SWMREC, permutations=9999)
adonis.SWMREC.TF<-adonis2(t(otu_course.SWMREC) ~ Tissue*Fortnight, strata=metadata_course.SWMREC$Replicate,  data=metadata_course.SWMREC, permutations=9999)

# betadispersion
vegan::vegdist(t(otu_course.SWMREC), method="bray") -> dist_otu_course.SWMREC #make base info
permdisp_otu_course.SWMREC.tissue <- betadisper(dist_otu_course.SWMREC, metadata_course.SWMREC$Tissue)#Bdisp for tissue
betadisp.course.tissue.SWMREC<-anova(permdisp_otu_course.SWMREC.tissue, permutations = 9999)#set to write
print(betadisp.course.tissue.SWMREC) #used for Table 1

permdisp_otu_course.SWMREC.Week <- betadisper(dist_otu_course.SWMREC, metadata_course.SWMREC$Week)#Bdisp for Week
betadisp.course.Week.SWMREC<-anova(permdisp_otu_course.SWMREC.Week, permutations = 9999)#set to write
print(betadisp.course.Week.SWMREC) #used for Table 1

permdisp_otu_course.SWMREC.Fortnight <- betadisper(dist_otu_course.SWMREC, metadata_course.SWMREC$Fortnight)#Bdisp for Fortnight
betadisp.course.Fortnight.SWMREC<-anova(permdisp_otu_course.SWMREC.Fortnight, permutations = 9999)#set to write
print(betadisp.course.Fortnight.SWMREC) #used for Table 1

