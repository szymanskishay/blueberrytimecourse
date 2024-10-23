##Prioritizing core microbiome based on occupancy
##Pretty limited modifications from Shade, A and Stopnisek, N. Abundance-occupancy distributions to prioritize plant core microbiome membership. Curr Opin Microbiol; 49:50-8
library(vegan)
library(tidyverse)
#directory<-"path/to/timecourse"
#setwd(directory)
ps.course<-readRDS("pscourse.rds")
#The below script was executed three different times for the variables Fortnight, Tissue, SiteYear
#To repeat, run these but after running 'Tissue' , replace it with 'Fortnight' and 'SiteYear'
#After that, OTUs that were in all 3 output files were identified
#The final summary OTU table is available in the Generated Tables folder as "Core_All_Rank_Combined.csv"

mat.all <- as(t(otu_table(ps.course)), "matrix")
class(mat.all)
rowsums.all <- sort(rowSums(mat.all))
raremax.all <- min(rowSums(mat.all))
newmat.all<-rrarefy(mat.all, sample=raremax.all)

nReads=min(rowSums(newmat.all))                                                        # input dataset needs to be rarefied and the rarefaction depth included 
otu <- t(newmat.all)
map <- sample_data(ps.course)
map$sequence_name<-row.names(map)
otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('otu')

# Occupancy abundance plot:
ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill='white') +
  labs(x="log10(mean relative abundance)", y="Occupancy")

# Ranking OTUs based on their occupancy
# For calculating ranking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(sequence_name, abun, -otu) %>% #sequence name is a column designating to samples
  #I am deeply confused as to how this works because it works before referencing things needed further on
  left_join(map, by = 'sequence_name') %>% 
  group_by(otu, Tissue) %>% #sampling date is the conditional occupancy. so for tissue specific, it would be tissue here
  summarise(time_freq=sum(abun>0)/length(abun),            # frequency of detection between time points
            coreTime=ifelse(time_freq == 1, 1, 0)) %>%     # 1 only if occupancy 1 with specific time, 0 if not
  group_by(otu) %>%
  summarise(sumF=sum(time_freq),
            sumG=sum(coreTime),
            nS=length(Tissue),
            Index=(sumF+sumG)/nS)                 # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,                           # transmute creates a new data frame containing only the specified computations
            rank=Index) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the BC similarity
BCaddition <- NULL 

# calculating BC dissimilarity based on the 1st ranked OTU
otu_start=otu_ranked$otu[1]                   
start_matrix <- as.matrix(otu[otu_start,])
start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
#apply returns values parsed through a function. this takes combn as X, 2 is the margin (columns), and the function is the function applied
#combn generates all combinations of elements of X taken M at a time. Works with a supplied function. 
#function takes the sum of the absolute value of differences from element 1 and 2 in the start matrix and divides by 2xRarefying depth
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
#makes new names for rows based on the pairs in comparison?
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1  #replaces name of second column with 1?
BCaddition <- rbind(BCaddition,df_s)
# calculating BC dissimilarity based on addition of ranked OTUs from 2nd to 500th. 
# Can be set to the entire length of OTUs in the dataset, however it might take 
# some time if more than 5000 OTUs are included.
for(i in 2:500){                              
  otu_add=otu_ranked$otu[i]                       
  add_matrix <- as.matrix(otu[otu_add,])
  add_matrix <- t(add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
}
# calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])))
#Creating plot of Bray-Curtis similarity
ggplot(BC_ranked[1:500,], aes(x=factor(rank[1:500], levels=rank[1:500]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+14, y=.1, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=lastCall+3, y=.5, label=paste("Last 2% increase (",lastCall,")",sep=''), color="blue")

coreis.all<-otu_ranked[1:lastCall,]
write.csv(coreis.all, file="Generated Tables/corebyTissue.csv")
