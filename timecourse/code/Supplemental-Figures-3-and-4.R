library(phyloseq)
library(ggplot2)
library(tidyverse)
ps.course<-readRDS("pscourse.RDS")

functionaltable<-read.csv(file="coredata/FunctionalGuilds.csv", row.names = 1)

tax_table(ps.course)<-as.matrix(functionaltable)
otu_table(ps.course) <- otu_table(ps.course)[which(rowSums(otu_table(ps.course)) >= 10),]### PCR Errors


ps.course.bar.Growth_form_template<- ps.course %>%
  tax_glom(taxrank = "Growth_form_template") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Growth_form_template)           # Sort data frame alphabetically by Genus

dat_ps.course.bar.Growth_form_template <- data.table(ps.course.bar.Growth_form_template)
dat_ps.course.out.Growth_form_template <- dat_ps.course.bar.Growth_form_template[(Abundance != 0)]
write.csv(dat_ps.course.out.Growth_form_template, "Generated Tables/Growth/All_Growth_form_template.csv") #allows examination of groupings

dat_ps.course.out.Growth_form_template %>%
  distinct(Growth_form_template) -> Growth_form_template_key
dat_ps.course.bar.Growth_form_template %>%
  mutate(Week= as.integer(Week)) %>%
  mutate(Growth_form_template = replace(Growth_form_template, Growth_form_template == "zoosporic-rhizomycelial_(chytrid-like)", "Zoosporic")) %>%
  mutate(Growth_form_template = replace(Growth_form_template, Growth_form_template == "thallus_photosynthetic", "Photosynthetic Thallus")) %>%
  mutate(Growth_form_template = replace(Growth_form_template, Growth_form_template == "dimorphic_yeast", "Dimorphic Yeast")) %>%
  mutate(Growth_form_template = replace(Growth_form_template, Growth_form_template == "yeast", "Yeast")) %>%
  mutate(Growth_form_template = replace(Growth_form_template, Growth_form_template == "filamentous_mycelium", "Filamentous")) %>%
  mutate(Growth_form_template = replace(Growth_form_template, Growth_form_template == "0", "#N/A"))


Tissue.Growth_form_template <- ggplot(dat_ps.course.bar.Growth_form_template, aes(x = Sample, y = Abundance, fill = Growth_form_template)) +
  theme_classic()+
  facet_wrap(~Tissue, strip.position="bottom", scales="free_x")+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#N/A"="#000000",
                             "0"="#ffffff",
                             "Photosynthetic Thallus"="#40f1c2",
                             "Filamentous"="#009736",
                             "Yeast"="#ffff00",
                             "Dimorphic Yeast"="#EE2A35",
                             "Zoosporic"="#a300ff"
  ))+
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 20, face = "bold")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=18)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) + theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 20, face = "bold")) +
  theme(legend.position="bottom")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, title = "Growth Form")) +
  ylab("Relative Abundance \n") +
  xlab("")+
  labs(title = "Relative Abundance of Growth Form by Tissue")


Fortnight.Growth_form_template <- ggplot(dat_ps.course.bar.Growth_form_template, aes(x = Sample, y = Abundance, fill = Growth_form_template)) +
  theme_classic()+
  facet_wrap(~Fortnight, strip.position="bottom", scales="free_x")+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#N/A"="#000000",
                             "0"="#ffffff",
                             "Photosynthetic Thallus"="#40f1c2",
                             "Filamentous"="#009736",
                             "Yeast"="#ffff00",
                             "Dimorphic Yeast"="#EE2A35",
                             "Zoosporic"="#a300ff"
  ))+
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 20, face = "bold")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=18)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")) + theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 20, face = "bold")) +
  theme(legend.position="bottom")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, title = "Growth Form")) +
  ylab("Relative Abundance \n") +
  xlab("")+
  labs(title = "Relative Abundance of Growth Form by Fortnight")

#figure out a good size for the tiff

#final figures modified slightly in MSPaint for better placement of the figure legends