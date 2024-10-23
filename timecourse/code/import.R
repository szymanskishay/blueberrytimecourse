library(phyloseq)
library(Biostrings)
library(data.table)
library(tidyverse)
library(microDecon)

##### 
options(scipen = 999) # sets the number of decimals displayed on screen 0
set.seed(0413) # setting the seed for reproducibility so analyses can be reproduced [0413]
#directory <- "/path/to/folder/timecourse" #set path to directory to be worked in
#setwd(directory) # set directory
fungi_otus<- read.delim("coredata/otutable_UPARSE_225bp.txt") #read in otu table from UPARSE

fungi_otus_reorder <- fungi_otus %>%
  relocate(all_of(c("Sample387","Sample388","Sample389","Sample390","Sample391","Sample392","Sample393")), .after = "X.OTU.ID")
#above line reorders "control" samples to be at the front for the purposes of microDecon
deconout<-decon(fungi_otus_reorder, numb.ind=383, numb.blanks=7, taxa=F) #filters out contaminate OTUs
#below code assumes a particular file structure
write.csv(deconout$reads.removed,file="Generated Tables/decon/reads_removed_counts.csv") #coutns reads removed per sample
write.csv(deconout$mean.per.group,file="Generated Tables/decon/mean_per_group.csv") #mean removed per OTU
write.csv(deconout$sum.per.group,file="Generated Tables/decon/sum_per_group.csv") #sum removed per OTU
write.csv(deconout$OTUs.removed,file="Generated Tables/decon/completely_removed_OTUs.csv") #OTUs completely eliminated

decon_tab<-deconout$decon.table
rownames(decon_tab)=decon_tab$X.OTU.ID
decon_otus<-decon_tab[,-c(1,2)]


fungi_otus_phy <-otu_table(decon_otus,
                           taxa_are_rows = TRUE) #read in decontaminated OTU table for phyloseq
fungi_metadata <-read.delim("coredata/map.txt",
                            row.names=1) #read in metadata

fungi_metadata_phy <-sample_data(fungi_metadata) #convert metadata for phyloseq
fungi_seqs<-readDNAStringSet("coredata/otus_225bp.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE) #DNA sequences
fungi_taxonomy<-read.delim("coredata/constax_taxonomy_225bp.txt",
                           header=TRUE, 
                           row.names=1) # read in taxonomy information to dataframe
fungi_taxonomy_phy <- tax_table(as.matrix(fungi_taxonomy)) # put into phyloseq readable format
course<-phyloseq(fungi_metadata_phy, fungi_otus_phy, fungi_taxonomy_phy, fungi_seqs) #assemble phyloseq
#df.mainseq<-as.data.frame(sample_data(fungi_metadata_phy))# Convert to data frame for exporting if you want
#these lines are more useful for identifying problematic samples (ie low reads)
#df_course <- as.data.frame(sample_data(course)) # Put sample_data into a easier to export data.frame
#df_course$LibrarySize_course <- sample_sums(course) #Add a column for Library Size
#df_course <- df_course[order(df_course$LibrarySize_course),] #sorts
#df_course$Index <- seq(nrow(df_course)) #Adds another layer of order
#write.csv(df_course, file = "Generated Tables/rank_sums_course.csv") #allows easier ID of samples with high/low reads
#the following code allows us to abolish bad samples. The above file shows a number of samples with less than 1000 reads after decontamination so we're cutting them
otu_table(course) <- subset(otu_table(course), select = -c(Sample100,
                                                           Sample146,
                                                           Sample27,
                                                           Sample29,
                                                           Sample2,
                                                           Sample328,
                                                           Sample5,
                                                           Sample103,
                                                           Sample101,
                                                           Sample1,
                                                           Sample200,
                                                           Sample293,
                                                           Sample6,
                                                           Sample292,
                                                           Sample4,
                                                           Sample83,
                                                           Sample308,
                                                           Sample306,
                                                           Sample288,
                                                           Sample307,
                                                           Sample275,
                                                           Sample199,
                                                           Sample212,
                                                           Sample202,
                                                           Sample309,
                                                           Sample344,
                                                           Sample97,
                                                           Sample136,
                                                           Sample49,
                                                           Sample40,
                                                           Sample277,
                                                           Sample41))

#samples for consideration to remove for strange/bad distribution of taxa suggesting poor sample
ps.course.naut<-course #backup of the phyloseq object at this point before further OTU removal
'%notin%'<- Negate('%in%') #useful for next step, also just in general! 
#below list is bunch of OTUs that are assigned to Viridiplantae or Metazoa or Malassezia
plant<-c("OTU_2578",
        "OTU_3509",
        "OTU_1484",
        "OTU_1896",
        "OTU_95",
        "OTU_2703",
        "OTU_2936",
        "OTU_3801",
        "OTU_2951",
        "OTU_3303",
        "OTU_2437",
        "OTU_3421",
        "OTU_2729",
        "OTU_3499",
        "OTU_953",
        "OTU_1366",
        "OTU_1474",
        "OTU_1863",
        "OTU_1873",
        "OTU_2161",
        "OTU_3038",
        "OTU_3135",
        "OTU_3538",
        "OTU_3623",
        "OTU_3638",
        "OTU_3789",
        "OTU_3854",
        "OTU_2891",
        "OTU_2209",
        "OTU_2535",
        "OTU_2615",
        "OTU_3012",
        "OTU_3126",
        "OTU_3289",
        "OTU_3293",
        "OTU_3399",
        "OTU_3404",
        "OTU_3526",
        "OTU_3634",
        "OTU_3679",
        "OTU_378",
        "OTU_2587",
        "OTU_3482",
        "OTU_1629",
        "OTU_2392",
        "OTU_300",
        "OTU_401",
        "OTU_462",
        "OTU_1054",
        "OTU_1093",
        "OTU_1151",
        "OTU_1698",
        "OTU_2093",
        "OTU_2242",
        "OTU_2300",
        "OTU_2303",
        "OTU_2494",
        "OTU_2628",
        "OTU_2850",
        "OTU_3414",
        "OTU_3443",
        "OTU_3517",
        "OTU_3598",
        "OTU_3769",
        "OTU_3797",
        "OTU_3847",
        "OTU_3860",
        "OTU_89",
        "OTU_170",
        "OTU_209",
        "OTU_233",
        "OTU_395",
        "OTU_455",
        "OTU_467",
        "OTU_484",
        "OTU_1423",
        "OTU_1832",
        "OTU_1916",
        "OTU_1184",
        "OTU_2200",
        "OTU_3204",
        "OTU_749",
        "OTU_490",
        "OTU_493",
        "OTU_824",
        "OTU_3369",
        "OTU_3378",
        "OTU_3699",
        "OTU_660",
        "OTU_582",
        "OTU_906",
        "OTU_908",
        "OTU_921",
        "OTU_1502",
        "OTU_1622",
        "OTU_3647",
        "OTU_985",
        "OTU_1166",
        "OTU_1307",
        "OTU_1413",
        "OTU_2252",
        "OTU_852") #OTU_89, 170, 233, 395, 455, 467, 484, 803, 852, 946, 1091, 2726, 2852, 209, 1715
plantframe<-cbind(plant,plant)
row.names(plantframe)<-plantframe[,1]
keep<-row.names(otu_table(course)) %notin% row.names(plantframe) #list of OTUs not in the above list
course_pruned<-prune_taxa(keep, course) #new phyloseq object with plants pruned
ps.course<-course_pruned #sets to name that will be used from now on
saveRDS(ps.course, file="pscourse.rds", compress=FALSE) #saves as an R object for easier re-use in other scripts!