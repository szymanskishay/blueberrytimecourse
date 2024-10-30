library(devtools)
BiocManager::install("DESeq2")
BiocManager::install("genefilter")
install.packages("pheatmap")
BiocManager::install("edgeR")
install.packages("RJSONIO")
install.packages('C:/Users/sakuy/Downloads/Tax4Fun_0.3.1.tar.gz', repos=NULL, type='source')
devtools::install_github("joey711/biom")
devtools::install_github("xia-lab/MicrobiomeAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
devtools::install_github("kylebittinger/qiimer")
library(MicrobiomeAnalystR)
library(ggplot2)
library(tidyverse)
directory<-"C:/Users/sakuy/Desktop/timecourse/"
setwd(directory)
ps.course<-readRDS("pscourse.rds")
write.table(otu_table(ps.course), file = "Generated Tables/lefse/otu_table_extra.txt")
write.table(sample_data(ps.course), file = "Generated Tables/lefse/sample_data_extra.txt")
write.table(tax_table(ps.course), file = "Generated Tables/lefse/tax_table_extra.txt")
#The above lines write out and make copies of the processed and pruned OTU table and taxonomy table and sample data to import into MicrobiomeAnalyst
#After writing these out, I added to the first line and space of each file #SAMPLE followed by a space, except for the taxonomy file which has had added #TAXONOMY followed by a space. This allows for import into MicrobiomeAnalystR
mbCourse<-Init.mbSetObj()
mbCourse<-SetModuleType(mbCourse, "mdp")
mbCourse<-ReadSampleTable(mbCourse, "Generated Tables/lefse/sample_data_extra.txt");
mbCourse<-Read16STaxaTable(mbCourse, "Generated Tables/lefse/tax_table_extra.txt");
mbCourse<-Read16SAbundData(mbCourse, "Generated Tables/lefse/otu_table_extra.txt","text","Others/Not_specific","T","false");
mbCourse<-SanityCheckData(mbCourse, "text","sample","true");
mbCourse<-SanityCheckSampleData(mbCourse);
mbCourse<-SetMetaAttributes(mbCourse, "1")
mbCourse<-PlotLibSizeView(mbCourse, "Generated Tables/lefse/norm_libsizes_0","png"); #will create a png of library sizes
mbCourse<-SanityCheckData(mbCourse, "text","count","true");
mbCourse<-SanityCheckSampleData(mbCourse);
mbCourse<-SetMetaAttributes(mbCourse, "1")
mbCourse<-PlotLibSizeView(mbCourse, "Generated Tables/lefse/norm_libsizes_1","png"); #labeled by week
mbCourse<-CreatePhyloseqObj(mbCourse, "text","Others/Not_specific","F" , "false")
mbCourse<-ApplyAbundanceFilter(mbCourse, "prevalence", 5, 0.05); #filtering OTUs by sample prevalence (must be in at least 5% of samples with at least 5 reads)
mbCourse<-ApplyVarianceFilter(mbCourse, "iqr", 0.1); #filtering by variance (filter by low variance)
mbCourse<-PerformNormalization(mbCourse, "none", "none", "none") #The object needs to have a normalization component, but this arrangement removes pre-normalization as the utilized LefSe package already deals with normalization its own ways
mbCourse<-PerformLefseAnal(mbCourse, 0.05, "fdr", 2.0, "Tissue", "F", "NA", "Species")
mbCourse<-PlotLEfSeSummary(mbCourse, 15, "dot", "Generated Tables/lefse/Tissue2","png"); # outputs native LefSe plot, but part of the scripting breaks and doesn't allow you to set width properly so names get cut off
#So now, we need to make our own plots from this same data

mbCourse[["analSet"]][["lefse"]][["resTable"]]->tissuetable

tissue_plot<-tissuetable[order(tissuetable$LDAscore),]
tissue_plot$Name=factor(row.names(tissue_plot), levels=row.names(tissue_plot))
#examining the current state of things, creating a vector species of interest to this study to help streamline presentation of data
ofinterest<-c("Cladosporium_herbarum", "Epicoccum_nigrum","Aureobasidium_pullulans", "Peltaster_fructicola", "Sporobolomyces_roseus", "Taphrina_letifera", "Dioszegia_hungarica", "Peltaster_gemmifer","Filobasidium_floriforme","Colletotrichum_fioriniae", "Filobasidium_oeirense", "Papiliotrema_fusca","Dioszegia_zsoltii","Sporobolomyces_ruberrimus","Botrytis_cinerea")
tissue_plot %>% filter(LDAscore >=2.25) %>% filter(FDR < 0.05) %>% filter(row.names(.) %in% ofinterest) -> tissue_plot
#filters by LDA score and False Discovery Rate for significant ones, also selects for taxa that are in the taxa of interest
tissue_plot %>% mutate(avg = (Blue + Bud + Coloring + Flower + Green + Leaf + Stem)/7) %>% #getting an average abundance score across tissues for each taxon
  mutate(stdev = sqrt(((Blue-avg)^2 + (Bud-avg)^2 + (Coloring-avg)^2 + (Flower-avg)^2 + (Green-avg)^2 + (Leaf-avg)^2 + (Stem-avg)^2)/7))%>% #stdev of abundance scores for each taxon
  mutate(Blue = (Blue-avg)/stdev) %>% #getting affinity for each tissue by getting a value that demonstrates how many standard deviations away from the average the taxon is in each tissue
  mutate(Bud = (Bud-avg)/stdev) %>%
  mutate(Coloring = (Coloring-avg)/stdev) %>%
  mutate(Flower = (Flower-avg)/stdev) %>%
  mutate(Green = (Green-avg)/stdev) %>%
  mutate(Leaf = (Leaf-avg)/stdev) %>%
  mutate(Stem = (Stem-avg)/stdev) -> tissue_plot.test

coreplot_tissue<-ggplot(tissue_plot.test, aes(y=factor(Name), x=LDAscore))+ #generating plot of the LDA score for each taxon
  theme_classic()+
  geom_point(stat="identity", size=4)+
  labs(y="Species", x="LDAscore")+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0),
        panel.grid.major.y = element_line(color = 'grey', linetype=5),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(size = 16, face = "bold"))


tissue_plot.test %>% dplyr::select(Name, Blue, Bud, Green, Flower, Leaf, Coloring, Stem) %>%
  pivot_longer(!Name, names_to="Tissue")-> tissue_plot.side #prepares the data for doing a heatmap-esque presentation of those standard deviation numbers from above
plotside_tissue<-ggplot(tissue_plot.side,aes(x=Tissue, y=factor(Name), fill=value))+
  theme_classic()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_blank(),  axis.line = element_blank())+
  geom_tile(aes(fill = value, width = 0.7, height=0.7), linewidth=0.5, color="black")+
  scale_fill_gradient2(limits=c(min(plotdata3$value), max(plotdata3$value)),
                       breaks=c(min(plotdata3$value), 0, max(plotdata3$value)),
                       low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       label = function(x) sprintf("%.2f",x))+
  scale_x_discrete(position="top")+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=33, hjust=0.2, size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold", margin = margin(0,0,10,0)),
        legend.text = element_text(size = 14),
        legend.title.align = 0.5,
        plot.background = element_blank(),
        legend.position="right")+
  labs(fill=" Deviation \n from Mean")
library(ggpubr) #puts ggplots together. below is going to create a tiff with both portions of the figure which are then going to be edited briefly in paint to remove blank space
tiff(filename="Figures/Figure_5.tiff", width = 1500, height = 600, units="px")
ggarrange(coreplot_tissue, plotside_tissue, nrow=1, common.legend = FALSE, legend="right", align='hv', widths = c(1, 1))
dev.off()

#repeat above for fortnight

mbCourse<-PerformLefseAnal(mbCourse, 0.05, "fdr", 2.0, "Fortnight", "F", "NA", "Species")
mbCourse[["analSet"]][["lefse"]][["resTable"]]->fortnighttable

fortnight_plot<-fortnighttable[order(fortnighttable$LDAscore),]
fortnight_plot$Name=factor(row.names(fortnight_plot), levels=row.names(fortnight_plot))
fortnight_plot %>% filter(LDAscore >=2) %>% filter(FDR < 0.05) -> fnint #check what our valid fortnight LDA targets are
cbind(ofinterest, ofinterest %in% fnint$Name) #see which of the previous list of interest is in here and if we should change it at all
#we can see that a few members of the initial list aren't relevant, so lets swap out some 
ofinterestFN<-c("Epicoccum_nigrum","Aureobasidium_pullulans", "Filobasidium_magnum", "Sporobolomyces_roseus", "Taphrina_letifera", "Filobasidium_wieringae", "Peltaster_gemmifer","Filobasidium_floriforme","Colletotrichum_fioriniae","Bullera_alba", "Vishniacozyma_carnescens","Sclerotinia_sclerotiorum","Taphrina_inositophila","Sporobolomyces_ruberrimus","Botrytis_cinerea")
fortnight_plot %>% filter(LDAscore >=2) %>% filter(FDR < 0.05) %>% filter(Name %in% ofinterestFN) -> fortnight_plot
colnames(fortnight_plot)<-c("Pvalues", "FDR", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "LDAscore", "Name")
fortnight_plot.test = fortnight_plot %>% mutate(avg = (One + Two + Three + Four + Five + Six + Seven)/7) %>%
  mutate(stdev = sqrt(((Two-avg)^2 + (Three-avg)^2 + (Four-avg)^2 + (Five-avg)^2 + (Six-avg)^2 + (One-avg)^2 + (Seven - avg)^2)/7))%>%
  mutate(One = (One-avg)/stdev) %>%
  mutate(Two = (Two-avg)/stdev) %>%
  mutate(Three = (Three-avg)/stdev) %>%
  mutate(Four = (Four-avg)/stdev) %>%
  mutate(Five = (Five-avg)/stdev) %>%
  mutate(Six = (Six-avg)/stdev) %>%
  mutate(Seven = (Seven-avg)/stdev)

coreplot_fortnight<-ggplot(fortnight_plot.test, aes(y=factor(Name), x=LDAscore))+
  theme_classic()+
  geom_point(stat="identity", size=4)+
  labs(y="Species", x="LDAscore")+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0),
        panel.grid.major.y = element_line(color = 'grey', linetype=5),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(size = 16, face = "bold"))

fortnight_plot.test %>% dplyr::select(Name, One, Two, Three, Four,Five,Six, Seven) %>%
  pivot_longer(!Name, names_to="Fortnight")-> fortnight_plot.side

plotside_fortnight<-ggplot(fortnight_plot.side,aes(x=Fortnight, y=factor(Name), fill=value))+
  theme_classic()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_blank(),  axis.line = element_blank())+
  geom_tile(aes(fill = value, width = 0.7, height=0.7), linewidth=0.5, color="black")+
  scale_fill_gradient2(limits=c(min(fortnight_plot.side$value),max(fortnight_plot.side$value)),
                       low="darkblue",
                       mid="white",
                       high="darkred",
                       breaks=c(min(fortnight_plot.side$value),0,max(fortnight_plot.side$value)),
                       label = function(x) sprintf("%.2f",x))+
  scale_x_discrete(position="top", limits=c("One","Two","Three","Four","Five","Six", "Seven"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=33, hjust=0.2, size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold", margin = margin(0,0,10,0)),
        legend.text = element_text(size = 14),
        legend.title.align = 0.5,
        plot.background = element_blank(),
        legend.position="right")+
  labs(fill=" Deviation \n from Mean")

tiff(filename="Figures/Figure_6.tiff", width = 1500, height = 600, units="px")
ggarrange(coreplot_fortnight, plotside_fortnight, nrow=1, common.legend = FALSE, legend="right", align='hv', widths=c(1.2,1))
dev.off()

#as with the previous figure, this was adjusted in paint to remove blank space