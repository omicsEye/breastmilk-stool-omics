
library(reshape2)
library(phyloseq)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("phyloseq")
library(DT)
library(dplyr)
library(taxonomizr)
library(ggplot2)
library(readxl)
library(vegan)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
source("~/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/scripts/metaphlanToPhyloseq.R")
erie <- metaphlanToPhyloseq(infant_microbiome, infant_metadata_microbiome)
#otumat <- as.data.frame(mother_microbiome)
#taxmat = matrix(sample(letters, 7, replace = TRUE), nrow = nrow(otumat), ncol = 7)
#rownames(taxmat) <- rownames(otumat)
#colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#physeq <- phyloseq(otu_table(otumat, taxa_are_rows = T), tax_table(taxmat), sample_data(mother_metadata))
#plot_bar(erie, x="Sample", y="Abundance", fill = "Phylum")

# melt to long format (for ggploting) 
# prune out phyla below 2% in each sample
rank_level <- "Family"
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
erie_phylum <- erie %>%
  tax_glom(taxrank = rank_level) %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  dplyr::filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#410061"
)


# Plot 
abuddance_bar_plots <- vector('list', length(unique(colnames(infant_metadata_microbiome))))
names(abuddance_bar_plots) <- unique(colnames(infant_metadata_microbiome))
pdf(
  paste('analysis', '/', rank_level,'_infnat_abundance_barplot.pdf', sep = ''),
  width = 10,
  height = 10,
  onefile = TRUE
)

for (meta in unique(colnames(infant_metadata_microbiome))){
  tryCatch({
  abuddance_bar_plots[[meta]] <- ggplot(erie_phylum, aes(x = Sample, y = Abundance, fill = Family)) + 
    #facet_grid(Breast_milk_collected~.,scales="free") +
    facet_grid(.~get(meta),scales="free")+
    #facet_grid(Breast_milk_collected~.,scales="free") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = phylum_colors) +
    # scale_x_discrete(
    #   breaks = c("7/8", "8/4", "9/2", "10/6"),
    #   labels = c("Jul", "Aug", "Sep", "Oct"), 
    #   drop = FALSE
    # ) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) + 
    #
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance") +
    ggtitle(meta) + # Remove x axis title
    theme(
      axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks.x=element_blank()) #+ theme_omicsEye()#+theme(legend.position = "none")
  stdout <-
    capture.output(print(abuddance_bar_plots[[meta]]), type = "message")
  }, error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dev.off()

# Ordinate
erie_pcoa <- ordinate(
  physeq = erie, 
  method = "NMDS", 
  distance = "bray"
)

# Plot 
plot_donor_maternal <- plot_ordination(
  physeq = erie,
  ordination = erie_pcoa,
  color = "Gestational_age_week",
  #shape = "Gestational_age_week",
  title = "Gestational_age_week"
) + 
  #scale_color_manual(values = c("#a65628", "red", "#ffae19",
  #                              "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  #) +
  #geom_point(aes(color = Gestational_age_week, fill= Gestational_age_week), alpha = 0.7, size = 2, stroke = 1.5) +
  omicsArt::theme_omicsEye() + theme(legend.position = c(.55, .9))  #geom_point(colour = Gestational_age_week,size = 1.5)
plot_donor_maternal
ggsave(filename=paste('manuscript/figures/NMDS_plot_infnat_microbiome.png', sep=''), plot=plot_donor_maternal,
       width = 2.4, height = 2.2, units = "in", dpi = 350)

meta <-"Breast_milk_collected"  
infant_sample_collection <- ggplot(infant_metadata, aes(Gestational_age_week, External_ID, fill = Breast_milk_collected)) +
  geom_point( alpha = .75, shape = 21, size = 1, stroke = .1) +
  facet_grid(get(meta)~.,scales="free")+
  #scale_fill_manual(breaks = c(TRUE, FALSE),
  #                  values=c('white', 'maroon'))+
  xlab('Gestational age(day)') +  ylab('Infant stool') + #labs(title = "Infant sample collection" )+
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+ theme(legend.position = "none")
infant_sample_collection 
ggsave(filename=paste('manuscript/figures/infant_sample_collection.png', sep=''), plot=infant_sample_collection,
       width = 3.25, height = 2, units = "in", dpi = 350)

