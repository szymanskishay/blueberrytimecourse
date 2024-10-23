require(phyloseq)
require(vegan)
require(agricolae)
require(data.table)
require(ggplot2)
require(ggrepel)
require(ggpubr)
library(stringr)
#directory<-"/path/to/timecourse/"
#setwd(directory)
ps.course<-readRDS("pscourse.RDS")
#subsetting samples for down the line
PLP20<-subset_samples(ps.course, SiteYear%in%c("PLP20")) 
PLP21<-subset_samples(ps.course, SiteYear%in%c("PLP21"))
SWMREC<-subset_samples(ps.course, SiteYear%in%c("SWMREC"))

sd<-data.frame(sample_data(ps.course)) #Retrieving metadata from phyloseq
sd$shannon<-vegan::diversity(otu_table(ps.course), MARGIN=2, index="shannon") #shannon index for whole data set
sd$spec<-specnumber(otu_table(ps.course), MARGIN=2) #species number
sd$Pielou<-(sd$shannon/log(sd$spec)) #Pielou's evenness

shannon_anova_Tissue<-aov(shannon~Tissue, sd) #ANOVA relating Tissue and diversity index
pielou_anova_tissue<-aov(Pielou~Tissue, sd)

tukey_Shannon_tissue<-HSD.test(shannon_anova_Tissue,  "Tissue", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on tissue type
tukey_pielou_tissue<-HSD.test(pielou_anova_tissue, "Tissue", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.tissue<-cbind(tukey_Shannon_tissue$means[,2:9], tukey_Shannon_tissue$groups[order(row.names((tukey_Shannon_tissue$groups))),])
tk.pielou.tissue<-cbind(tukey_pielou_tissue$means[,2:9], tukey_pielou_tissue$groups[order(row.names((tukey_pielou_tissue$groups))),])

Shannon.tissue.plot<-ggplot(tk.Shannon.tissue, aes(x=row.names(tk.Shannon.tissue),
                                                   ymin=Min,
                                                   lower=Q25,
                                                   middle=Q50,
                                                   upper=Q75,
                                                   ymax=Max,
                                                   fill=row.names(tk.Shannon.tissue)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Shannon's Diversity")+
  labs(tag="A")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  guides(fill=guide_legend(title="Tissue"))

pielou.tissue.plot<-ggplot(tk.pielou.tissue, aes(x=row.names(tk.pielou.tissue),
                                                 ymin=Min,
                                                 lower=Q25,
                                                 middle=Q50,
                                                 upper=Q75,
                                                 ymax=Max,
                                                 fill=row.names(tk.pielou.tissue)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Pielou's Evenness")+
  labs(tag="B")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  guides(fill=guide_legend(title="Tissue"))

##repeat for Week
shannon_anova_Week<-aov(shannon~Week, sd) #ANOVA
pielou_anova_Week<-aov(Pielou~Week, sd)

tukey_Shannon_Week<-HSD.test(shannon_anova_Week,  "Week", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on Week type
tukey_pielou_Week<-HSD.test(pielou_anova_Week, "Week", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.Week<-cbind(tukey_Shannon_Week$means[,2:9], tukey_Shannon_Week$groups[order(row.names((tukey_Shannon_Week$groups))),])
tk.pielou.Week<-cbind(tukey_pielou_Week$means[,2:9], tukey_pielou_Week$groups[order(row.names((tukey_pielou_Week$groups))),])

real_weeks<-str_sort(row.names(tk.Shannon.Week), numeric=TRUE) #R assumes we want to see Week 1 next to Week 10 without this...

Shannon.Week.plot<-ggplot(tk.Shannon.Week, aes(x=row.names(tk.Shannon.Week),
                                               ymin=Min,
                                               lower=Q25,
                                               middle=Q50,
                                               upper=Q75,
                                               ymax=Max,
                                               fill=row.names(tk.Shannon.Week)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.1), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.1, max.overlaps = 50)+
  xlab("Week")+
  ylab("Shannon's Diversity")+
  labs(tag="C")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_weeks)+
  scale_x_discrete(limits=real_weeks)+
  guides(fill=guide_legend(title="Week"))
pielou.Week.plot<-ggplot(tk.pielou.Week, aes(x=row.names(tk.pielou.Week),
                                             ymin=Min,
                                             lower=Q25,
                                             middle=Q50,
                                             upper=Q75,
                                             ymax=Max,
                                             fill=row.names(tk.pielou.Week)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.1), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.1, max.overlaps = 50)+
  xlab("Week")+
  ylab("Pielou's Evenness")+
  labs(tag="D")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_weeks)+
  scale_x_discrete(limits=real_weeks)+
  guides(fill=guide_legend(title="Week"))

#tiff("Figures/Alpha/All/All-Tissues-Just-Shannon-Pielou.tiff", height = 300, width=600, units="px")
tissueall <-ggarrange(Shannon.tissue.plot, pielou.tissue.plot, common.legend=TRUE, legend = "right")
#dev.off()
#tiff("Figures/Alpha/All/All-Weeks-Just-Shannon-Pielou.tiff", height = 300, width=600, units="px")
weekall <-ggarrange(Shannon.Week.plot, pielou.Week.plot, common.legend = TRUE, legend="right")
#dev.off()
tiff("Figures/Figure_2.tiff", height=600, width=600, units='px')
ggarrange(tissueall, weekall, common.legend = FALSE, ncol = 1, nrow = 2)
dev.off()

#####
sd_PLP20<-data.frame(sample_data(PLP20)) #Retrieving metadata from phyloseq
sd_PLP20$shannon<-vegan::diversity(otu_table(PLP20), MARGIN=2, index="shannon") #shannon index for whole data set
sd_PLP20$spec<-specnumber(otu_table(PLP20), MARGIN=2)
sd_PLP20$Pielou<-(sd_PLP20$shannon/log(sd_PLP20$spec))

shannon_anova_Tissue_PLP20<-aov(shannon~Tissue, sd_PLP20) #asessing variance and the like
pielou_anova_tissue_PLP20<-aov(Pielou~Tissue, sd_PLP20)

tukey_Shannon_tissue_PLP20<-HSD.test(shannon_anova_Tissue_PLP20,  "Tissue", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on tissue type
tukey_pielou_tissue_PLP20<-HSD.test(pielou_anova_tissue_PLP20, "Tissue", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.tissue_PLP20<-cbind(tukey_Shannon_tissue_PLP20$means[,2:9], tukey_Shannon_tissue_PLP20$groups[order(row.names((tukey_Shannon_tissue_PLP20$groups))),])
tk.pielou.tissue_PLP20<-cbind(tukey_pielou_tissue_PLP20$means[,2:9], tukey_pielou_tissue_PLP20$groups[order(row.names((tukey_pielou_tissue_PLP20$groups))),])

Shannon.tissue_PLP20.plot<-ggplot(tk.Shannon.tissue_PLP20, aes(x=row.names(tk.Shannon.tissue_PLP20),
                                                               ymin=Min,
                                                               lower=Q25,
                                                               middle=Q50,
                                                               upper=Q75,
                                                               ymax=Max,
                                                               fill=row.names(tk.Shannon.tissue_PLP20)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Shannon's Diversity")+
  labs(tag="A")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  scale_y_continuous(limits=c(0,4.5))+
  guides(fill=guide_legend(title="Tissue"))
pielou.tissue_PLP20.plot<-ggplot(tk.pielou.tissue_PLP20, aes(x=row.names(tk.pielou.tissue_PLP20),
                                                             ymin=Min,
                                                             lower=Q25,
                                                             middle=Q50,
                                                             upper=Q75,
                                                             ymax=Max,
                                                             fill=row.names(tk.pielou.tissue_PLP20)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Pielou's Evenness")+
  labs(tag="B")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  scale_y_continuous(limits=c(0,1))+
  guides(fill=guide_legend(title="Tissue"))

shannon_anova_Week_PLP20<-aov(shannon~Week, sd_PLP20) #asessing variance and the like
pielou_anova_Week_PLP20<-aov(Pielou~Week, sd_PLP20)

tukey_Shannon_Week_PLP20<-HSD.test(shannon_anova_Week_PLP20,  "Week", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on Week_PLP20 type
tukey_pielou_Week_PLP20<-HSD.test(pielou_anova_Week_PLP20, "Week", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.Week_PLP20<-cbind(tukey_Shannon_Week_PLP20$means[,2:9], tukey_Shannon_Week_PLP20$groups[order(row.names((tukey_Shannon_Week_PLP20$groups))),])
tk.pielou.Week_PLP20<-cbind(tukey_pielou_Week_PLP20$means[,2:9], tukey_pielou_Week_PLP20$groups[order(row.names((tukey_pielou_Week_PLP20$groups))),])

real_Week_PLP20s<-str_sort(row.names(tk.Shannon.Week_PLP20), numeric=TRUE)

Shannon.Week_PLP20.plot<-ggplot(tk.Shannon.Week_PLP20, aes(x=row.names(tk.Shannon.Week_PLP20),
                                                           ymin=Min,
                                                           lower=Q25,
                                                           middle=Q50,
                                                           upper=Q75,
                                                           ymax=Max,
                                                           fill=row.names(tk.Shannon.Week_PLP20)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.1), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.1, max.overlaps = 50)+
  xlab("Week")+
  ylab("Shannon's Diversity")+
  labs(tag="A")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_Week_PLP20s)+
  scale_x_discrete(limits=real_Week_PLP20s)+
  scale_y_continuous(limits=c(0,4.5))+
  guides(fill=guide_legend(title="Week"))
pielou.Week_PLP20.plot<-ggplot(tk.pielou.Week_PLP20, aes(x=row.names(tk.pielou.Week_PLP20),
                                                         ymin=Min,
                                                         lower=Q25,
                                                         middle=Q50,
                                                         upper=Q75,
                                                         ymax=Max,
                                                         fill=row.names(tk.pielou.Week_PLP20)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.1), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.1, max.overlaps = 50)+
  xlab("Week")+
  ylab("Pielou's Evenness")+
  labs(tag="B")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_Week_PLP20s)+
  scale_x_discrete(limits=real_Week_PLP20s)+ scale_y_continuous(limits=c(0,1))+
  guides(fill=guide_legend(title="Week"))

#############
sd_PLP21<-data.frame(sample_data(PLP21)) #Retrieving metadata from phyloseq
sd_PLP21$shannon<-vegan::diversity(otu_table(PLP21), MARGIN=2, index="shannon") #shannon index for whole data set
sd_PLP21$spec<-specnumber(otu_table(PLP21), MARGIN=2)
sd_PLP21$Pielou<-(sd_PLP21$shannon/log(sd_PLP21$spec))

shannon_anova_Tissue_PLP21<-aov(shannon~Tissue, sd_PLP21) #asessing variance and the like
pielou_anova_tissue_PLP21<-aov(Pielou~Tissue, sd_PLP21)

tukey_Shannon_tissue_PLP21<-HSD.test(shannon_anova_Tissue_PLP21,  "Tissue", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on tissue type
tukey_pielou_tissue_PLP21<-HSD.test(pielou_anova_tissue_PLP21, "Tissue", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.tissue_PLP21<-cbind(tukey_Shannon_tissue_PLP21$means[,2:9], tukey_Shannon_tissue_PLP21$groups[order(row.names((tukey_Shannon_tissue_PLP21$groups))),])
tk.pielou.tissue_PLP21<-cbind(tukey_pielou_tissue_PLP21$means[,2:9], tukey_pielou_tissue_PLP21$groups[order(row.names((tukey_pielou_tissue_PLP21$groups))),])

Shannon.tissue_PLP21.plot<-ggplot(tk.Shannon.tissue_PLP21, aes(x=row.names(tk.Shannon.tissue_PLP21),
                                                               ymin=Min,
                                                               lower=Q25,
                                                               middle=Q50,
                                                               upper=Q75,
                                                               ymax=Max,
                                                               fill=row.names(tk.Shannon.tissue_PLP21)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Shannon's Diversity")+
  labs(tag="C")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  scale_y_continuous(limits=c(0,4.5))+
  guides(fill=guide_legend(title="Tissue"))
pielou.tissue_PLP21.plot<-ggplot(tk.pielou.tissue_PLP21, aes(x=row.names(tk.pielou.tissue_PLP21),
                                                             ymin=Min,
                                                             lower=Q25,
                                                             middle=Q50,
                                                             upper=Q75,
                                                             ymax=Max,
                                                             fill=row.names(tk.pielou.tissue_PLP21)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Pielou's Evenness")+
  labs(tag="D")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  scale_y_continuous(limits=c(0,1))+
  guides(fill=guide_legend(title="Tissue"))

shannon_anova_Week_PLP21<-aov(shannon~Week, sd_PLP21) #asessing variance and the like
pielou_anova_Week_PLP21<-aov(Pielou~Week, sd_PLP21)

tukey_Shannon_Week_PLP21<-HSD.test(shannon_anova_Week_PLP21,  "Week", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on Week_PLP21 type
tukey_pielou_Week_PLP21<-HSD.test(pielou_anova_Week_PLP21, "Week", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.Week_PLP21<-cbind(tukey_Shannon_Week_PLP21$means[,2:9], tukey_Shannon_Week_PLP21$groups[order(row.names((tukey_Shannon_Week_PLP21$groups))),])
tk.pielou.Week_PLP21<-cbind(tukey_pielou_Week_PLP21$means[,2:9], tukey_pielou_Week_PLP21$groups[order(row.names((tukey_pielou_Week_PLP21$groups))),])

real_Week_PLP21s<-str_sort(row.names(tk.Shannon.Week_PLP21), numeric=TRUE)

Shannon.Week_PLP21.plot<-ggplot(tk.Shannon.Week_PLP21, aes(x=row.names(tk.Shannon.Week_PLP21),
                                                           ymin=Min,
                                                           lower=Q25,
                                                           middle=Q50,
                                                           upper=Q75,
                                                           ymax=Max,
                                                           fill=row.names(tk.Shannon.Week_PLP21)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.2), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.1, max.overlaps = 50)+
  xlab("Week")+
  ylab("Shannon's Diversity")+
  labs(tag="C")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_Week_PLP21s)+
  scale_x_discrete(limits=real_Week_PLP21s)+
  scale_y_continuous(limits=c(0,4.5))+
  guides(fill=guide_legend(title="Week"))
pielou.Week_PLP21.plot<-ggplot(tk.pielou.Week_PLP21, aes(x=row.names(tk.pielou.Week_PLP21),
                                                         ymin=Min,
                                                         lower=Q25,
                                                         middle=Q50,
                                                         upper=Q75,
                                                         ymax=Max,
                                                         fill=row.names(tk.pielou.Week_PLP21)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.1), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.1, max.overlaps = 50)+
  xlab("Week")+
  ylab("Pielou's Evenness")+
  labs(tag="D")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_Week_PLP21s)+
  scale_x_discrete(limits=real_Week_PLP21s)+ scale_y_continuous(limits=c(0,1))+
  guides(fill=guide_legend(title="Week"))
#####
sd_SWMREC<-data.frame(sample_data(SWMREC)) #Retrieving metadata from phyloseq
sd_SWMREC$shannon<-vegan::diversity(otu_table(SWMREC), MARGIN=2, index="shannon") #shannon index for whole data set
sd_SWMREC$spec<-specnumber(otu_table(SWMREC), MARGIN=2)
sd_SWMREC$Pielou<-(sd_SWMREC$shannon/log(sd_SWMREC$spec))

shannon_anova_Tissue_SWMREC<-aov(shannon~Tissue, sd_SWMREC) #asessing variance and the like
pielou_anova_tissue_SWMREC<-aov(Pielou~Tissue, sd_SWMREC)

tukey_Shannon_tissue_SWMREC<-HSD.test(shannon_anova_Tissue_SWMREC,  "Tissue", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on tissue type
tukey_pielou_tissue_SWMREC<-HSD.test(pielou_anova_tissue_SWMREC, "Tissue", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.tissue_SWMREC<-cbind(tukey_Shannon_tissue_SWMREC$means[,2:9], tukey_Shannon_tissue_SWMREC$groups[order(row.names((tukey_Shannon_tissue_SWMREC$groups))),])
tk.pielou.tissue_SWMREC<-cbind(tukey_pielou_tissue_SWMREC$means[,2:9], tukey_pielou_tissue_SWMREC$groups[order(row.names((tukey_pielou_tissue_SWMREC$groups))),])

Shannon.tissue_SWMREC.plot<-ggplot(tk.Shannon.tissue_SWMREC, aes(x=row.names(tk.Shannon.tissue_SWMREC),
                                                                 ymin=Min,
                                                                 lower=Q25,
                                                                 middle=Q50,
                                                                 upper=Q75,
                                                                 ymax=Max,
                                                                 fill=row.names(tk.Shannon.tissue_SWMREC)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Shannon's Diversity")+
  labs(tag="E")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  scale_y_continuous(limits=c(0,4.5))+
  guides(fill=guide_legend(title="Tissue"))

pielou.tissue_SWMREC.plot<-ggplot(tk.pielou.tissue_SWMREC, aes(x=row.names(tk.pielou.tissue_SWMREC),
                                                               ymin=Min,
                                                               lower=Q25,
                                                               middle=Q50,
                                                               upper=Q75,
                                                               ymax=Max,
                                                               fill=row.names(tk.pielou.tissue_SWMREC)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  labs(tag="F", y="Pielou's")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"), axis.title.y=element_text(size=14, face="bold"), legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Blue"="Blue", "Green"="Green", "Coloring"="Red", "Bud"="Orange", "Flower"="Purple", "Stem"="Brown", "Leaf"="darkgreen"))+
  scale_y_continuous(limits=c(0,1))+
  ylab("Pielou's Evenness")+
  guides(fill=guide_legend(title="Tissue"))

shannon_anova_Week_SWMREC<-aov(shannon~Week, sd_SWMREC) #asessing variance and the like
pielou_anova_Week_SWMREC<-aov(Pielou~Week, sd_SWMREC)

tukey_Shannon_Week_SWMREC<-HSD.test(shannon_anova_Week_SWMREC,  "Week", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on Week_SWMREC type
tukey_pielou_Week_SWMREC<-HSD.test(pielou_anova_Week_SWMREC, "Week", group=TRUE, unbalanced = TRUE)

#make into new objects for easier plotting
tk.Shannon.Week_SWMREC<-cbind(tukey_Shannon_Week_SWMREC$means[,2:9], tukey_Shannon_Week_SWMREC$groups[order(row.names((tukey_Shannon_Week_SWMREC$groups))),])
tk.pielou.Week_SWMREC<-cbind(tukey_pielou_Week_SWMREC$means[,2:9], tukey_pielou_Week_SWMREC$groups[order(row.names((tukey_pielou_Week_SWMREC$groups))),])

real_Week_SWMRECs<-str_sort(row.names(tk.Shannon.Week_SWMREC), numeric=TRUE)

Shannon.Week_SWMREC.plot<-ggplot(tk.Shannon.Week_SWMREC, aes(x=row.names(tk.Shannon.Week_SWMREC),
                                                             ymin=Min,
                                                             lower=Q25,
                                                             middle=Q50,
                                                             upper=Q75,
                                                             ymax=Max,
                                                             fill=row.names(tk.Shannon.Week_SWMREC)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.1), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.1, max.overlaps = 50)+
  xlab("Week")+
  ylab("Shannon's Diversity")+
  labs(tag="E")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_Week_SWMRECs)+
  scale_x_discrete(limits=real_Week_SWMRECs)+
  scale_y_continuous(limits=c(0,4.5))+
  guides(fill=guide_legend(title="Week"))
pielou.Week_SWMREC.plot<-ggplot(tk.pielou.Week_SWMREC, aes(x=row.names(tk.pielou.Week_SWMREC),
                                                           ymin=Min,
                                                           lower=Q25,
                                                           middle=Q50,
                                                           upper=Q75,
                                                           ymax=Max,
                                                           fill=row.names(tk.pielou.Week_SWMREC)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text_repel(aes(y=Max, label=groups), direction="y", position = position_nudge_repel(y=0.07), min.segment.length = 1, segment.size = 0.1, ylim = c(-Inf, Inf), segment.linetype=5, box.padding = 0.0001, max.overlaps = 50)+
  xlab("Week")+
  ylab("Pielou's Evenness")+
  labs(tag="F")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_Week_SWMRECs)+
  scale_x_discrete(limits=real_Week_SWMRECs)+
  scale_y_continuous(limits=c(0,1))+
  guides(fill=guide_legend(title="Week"))
#the below figures were further modified in mspaint to add the site-year labels
tiff("Figures/Figure_3.tiff", width = 700, height = 900, units="px")
ggarrange(Shannon.tissue_PLP20.plot, pielou.tissue_PLP20.plot, 
          Shannon.tissue_PLP21.plot, pielou.tissue_PLP21.plot, 
          Shannon.tissue_SWMREC.plot, pielou.tissue_SWMREC.plot,  ncol=2, nrow=3, legend="right", common.legend = TRUE)
dev.off()
tiff("Figures/Figure_4.tiff", width = 700, height = 900, units="px")
ggarrange(Shannon.Week_PLP20.plot, pielou.Week_PLP20.plot, 
          Shannon.Week_PLP21.plot, pielou.Week_PLP21.plot, 
          Shannon.Week_SWMREC.plot, pielou.Week_SWMREC.plot,  ncol=2, nrow=3, legend="right", common.legend = FALSE)
dev.off()
