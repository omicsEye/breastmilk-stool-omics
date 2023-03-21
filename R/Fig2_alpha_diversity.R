#library(ape)
library(vegan)
library(omicsArt)
library(Tweedieverse)

#run meatdata.R teh omics_process.R
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk")
# run omicsprocees.R script lines 1-140
microbiome_data <- as.data.frame(t(mother_microbiome))
#microbiome_data_metadata <- microbiome_data_metadata[microbiome_data_metadata$Visit =="a",]
intersect_samples <- intersect(rownames(microbiome_data), rownames(mother_metadata))
microbiome_data_metadata = mother_metadata[intersect_samples,]
microbiome_data = microbiome_data[intersect_samples,]

results_maternal_donor_alpha_divesity <- omicsArt::diversity_test(microbiome_data, microbiome_data_metadata[,c("Maternal vs. donor"), drop=F], output = NA)
ggsave(filename = 'manuscript/figures/fig3_diversity/results_maternal_donor_alpha_divesity.png', plot=results_maternal_donor_alpha_divesity$diversity_test_plots$`Maternal vs. donor`, width = 40, height = 40, units = "mm", dpi = 350)

microbiome_data_metadata_maternal <- microbiome_data_metadata[rownames(microbiome_data_metadata[microbiome_data_metadata$`Maternal vs. donor`=="Maternal", ]),]

maternal_microbiome <- microbiome_data[rownames(microbiome_data_metadata[microbiome_data_metadata$`Maternal vs. donor`=="Maternal", ]),]
microbiome_data_metadata_maternal <- microbiome_data_metadata_maternal[,colSums(is.na(microbiome_data_metadata_maternal))<nrow(microbiome_data_metadata_maternal)]
temp_filtered_metadata <- microbiome_data_metadata_maternal[, apply(microbiome_data_metadata_maternal, 2, entropy) > 0.5]
Excluded_metadata <- setdiff(colnames(metadata), colnames(temp_filtered_metadata))
logging::loginfo("Excluded metadata with entropy less or equal to 0.75: %s", Excluded_metadata)
microbiome_data_metadata_maternal <- temp_filtered_metadata
intersect_samples <- intersect(rownames(maternal_microbiome), rownames(microbiome_data_metadata_maternal))
microbiome_data_metadata_maternal = microbiome_data_metadata_maternal[intersect_samples,]
maternal_microbiome = microbiome_data[intersect_samples,]
microbiome_data_metadata_maternal <- microbiome_data_metadata_maternal[, ! colnames(microbiome_data_metadata_maternal) %in%
                                                                         c("External ID", "BMI category", "Maternal smoking")]



results <- omicsArt::diversity_test(maternal_microbiome, microbiome_data_metadata_maternal, output = NA)
maternal_overall_diversity_data <- results$alpha_diversity_data
maternal_overall_diversity_test <- results$alpha_diversity_test
maternal_overall_diversity_plots <- results$diversity_test_plots
maternal_overall_diversity_barplot <- results$overall_diversity_barplot
maternal_overall_diversity_barplot
microbiome_data_metadata_maternal <- microbiome_data_metadata_maternal[,maternal_overall_diversity_test$Metadata]
saveRDS(
  maternal_overall_diversity_plots,
  file = "analysis/diversity/maternal_overall_diversity_plots_gg_associations.RDS",
)

saveRDS(
  maternal_overall_diversity_barplot,
  file = "analysis/diversity/maternal_alpha_diversity_gg_associations.RDS",
)

pdf(
  paste('analysis/diversity', '/maternal_alpha_diversity.pdf', sep = ''),
  width = 2.4,
  height = 2,
  onefile = TRUE
)

for (meta in unique(colnames(mother_metadata))){
  tryCatch({
    stdout <-
      capture.output(print(maternal_overall_diversity_plots[[meta]]), type = "message")
  }, error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dev.off()
dev.off()
saveRDS(
  maternal_overall_diversity_barplot,
  file = "analysis/diversity/maternal_overall_diversity_barplot_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/maternal_overall_diversity_barplot.png', plot=maternal_overall_diversity_barplot,width = 2.5, height = 1.8, units = "in", dpi = 300)
ggsave(filename='analysis/diversity/maternal_overall_diversity_barplot.pdf', plot=maternal_overall_diversity_barplot, width = 2.5, height = 1.8, units = "in", dpi = 300)

###### infants
microbiome_data <- as.data.frame(t(infant_microbiome))


#tes alpha diversity of infant stool vs. Donor/Maternal
infant_results_maternal_donor_alpha_divesity <- omicsArt::diversity_test(microbiome_data[row.names(microbiome_data_metadata),], microbiome_data_metadata[,c("Maternal vs. donor"), drop=F], output = NA)

ggsave(filename = 'manuscript/figures/fig3_diversity/infant_results_maternal_donor_alpha_divesity.png', plot=infant_results_maternal_donor_alpha_divesity$diversity_test_plots$`Maternal vs. donor`, width = 40, height = 40, units = "mm", dpi = 350)
#infant_metadata2 <- infant_metadata[infant_metadata$Visit =="a",]
intersect_samples <- intersect(rownames(microbiome_data), rownames(infant_metadata))
infant_metadata = infant_metadata[intersect_samples,]
microbiome_data = microbiome_data[intersect_samples,]
microbiome_data_metadata <- infant_metadata[rownames(microbiome_data),]
infant_maternal_microbiome <- microbiome_data[rownames(microbiome_data_metadata[microbiome_data_metadata$`Maternal vs. donor`=="Maternal", ]),]
infant_microbiome_data_metadata_maternal <- microbiome_data_metadata[rownames(infant_maternal_microbiome),]

infant_microbiome_data_metadata_maternal <- infant_microbiome_data_metadata_maternal[,colSums(is.na(infant_microbiome_data_metadata_maternal))<nrow(microbiome_data_metadata_maternal)]
temp_filtered_metadata <- infant_microbiome_data_metadata_maternal[, apply(infant_microbiome_data_metadata_maternal, 2, entropy) > 0.5]
Excluded_metadata <- setdiff(colnames(metadata), colnames(temp_filtered_metadata))
logging::loginfo("Excluded metadata with entropy less or equal to 0.75: %s", Excluded_metadata)
infant_microbiome_data_metadata_maternal <- temp_filtered_metadata
intersect_samples <- intersect(rownames(infant_maternal_microbiome), rownames(infant_microbiome_data_metadata_maternal))
infant_maternal_microbiome = microbiome_data[intersect_samples,]
maternal_microbiome = microbiome_data[intersect_samples,]
infant_microbiome_data_metadata_maternal <- infant_microbiome_data_metadata_maternal[intersect_samples,]



infant_results <- omicsArt::diversity_test(infant_maternal_microbiome, infant_microbiome_data_metadata_maternal)
infant_alpha_diversity_data <- infant_results$alpha_diversity_data
infant_alpha_diversity_test <- infant_results$alpha_diversity_test
infant_alpha_diversity_plots <- infant_results$diversity_test_plots
infant_overall_diversity_barplot <- infant_results$overall_diversity_barplot
infant_microbiome_data_metadata_maternal <- infant_microbiome_data_metadata_maternal[,infant_alpha_diversity_test$Metadata]

saveRDS(
  infant_overall_diversity_barplot,
  file = "analysis/infant_overall_diversity_barplot_gg_associations.RDS",
)
saveRDS(
  infant_alpha_diversity_plots,
  file = "analysis/infant_overall_diversity_plots_gg_associations.RDS",
)
pdf(
  paste('analysis/diversity', '/infant_maternal_alpaha_diversity.pdf', sep = ''),
  width = 2.4,
  height = 2.25,
  onefile = TRUE
)

for (meta in unique(colnames(infant_microbiome_data_metadata_maternal))){
  tryCatch({
    stdout <-
      capture.output(print(infant_alpha_diversity_plots[[meta]]), type = "message")
  }, error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dev.off()
dev.off()
saveRDS(
  infant_overall_diversity_barplot,
  file = "analysis/diversity/infant_maternal_overall_diversity_barplot_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/infant_maternal_overall_diversity_barplot.png', plot=infant_overall_diversity_barplot,width = 2.5, height = 1.8, units = "in", dpi = 300)
ggsave(filename='analysis/diversity/infant_maternal_overall_diversity_barplot.pdf', plot=infant_overall_diversity_barplot,width = 2.5, height = 1.8, units = "in", dpi = 300)

intersect_samples <- intersect(rownames(infant_maternal_microbiome), rownames(microbiome_data_metadata_maternal))
microbiome_data_metadata_maternal2 <- microbiome_data_metadata_maternal[intersect_samples,]
infant_maternal_microbiome <- infant_maternal_microbiome[intersect_samples,]

infant_microbiomre_maternal_metadata_results <- omicsArt::diversity_test(infant_maternal_microbiome, microbiome_data_metadata_maternal2, order = F)
infantMicrobiome_maternalMetadata_alpha_diversity_data <- infant_microbiomre_maternal_metadata_results$alpha_diversity_data
infantMicrobiome_maternalMetadata_alpha_diversity_test <- infant_microbiomre_maternal_metadata_results$alpha_diversity_test
infantMicrobiome_maternalMetadata_alpha_diversity_plots <- infant_microbiomre_maternal_metadata_results$diversity_test_plots
infantMicrobiome_maternalMetadata_alpha_diversity_barplot <- infant_microbiomre_maternal_metadata_results$overall_diversity_barplot

saveRDS(
  infantMicrobiome_maternalMetadata_alpha_diversity_barplot,
  file = "analysis/infantMicrobiome_maternalMetadata_alpha_diversity_barplot_gg_associations.RDS",
)
saveRDS(
  infantMicrobiome_maternalMetadata_alpha_diversity_plots,
  file = "analysis/infantMicrobiome_maternalMetadata_alpha_diversity_plots_gg_associations.RDS",
)

pdf(
  paste('analysis/diversity', '/infantMicrobiome_maternalMetadata_alpha_diversity.pdf', sep = ''),
  width = 2.4,
  height = 2.25,
  onefile = TRUE
)

for (meta in unique(colnames(microbiome_data_metadata_maternal))){
  tryCatch({
    stdout <-
      capture.output(print(infantMicrobiome_maternalMetadata_alpha_diversity_plots[[meta]]), type = "message")
  }, error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dev.off()
dev.off()
saveRDS(
  infantMicrobiome_maternalMetadata_alpha_diversity_barplot,
  file = "analysis/diversity/infantMicrobiome_maternalMetadata_alpha_diversity_barplot_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/infantMicrobiome_maternalMetadata_alpha_diversity_barplot.png', plot=infantMicrobiome_maternalMetadata_alpha_diversity_barplot,width = 2.5, height = 1.8, units = "in", dpi = 300)
ggsave(filename='analysis/diversity/infantMicrobiome_maternalMetadata_alpha_diversity_barplot.pdf', plot=infantMicrobiome_maternalMetadata_alpha_diversity_barplot,width = 2.5, height = 1.8, units = "in", dpi = 300)

microbiome_data <- as.data.frame(t(mother_microbiome))
intersect_samples <- intersect(rownames(maternal_microbiome), rownames(infant_metadata))
microbiome_data_metadata_maternal2 <- microbiome_data_metadata_maternal[intersect_samples,]
infant_maternal_microbiome <- infant_maternal_microbiome[intersect_samples,]
maternal_microbiome <- microbiome_data[rownames(microbiome_data_metadata[microbiome_data_metadata$`Maternal vs. donor`=="Maternal", ]),]
microbiome_data <- maternal_microbiome[intersect_samples,]
maternalMicrobiome_InfantMetadata <- infant_metadata[intersect_samples,colnames(infant_microbiome_data_metadata_maternal)]
maternalMicrobiome_InfantMetadata_results <- omicsArt::diversity_test(microbiome_data, maternalMicrobiome_InfantMetadata, order = F)
maternalMicrobiome_InfantMetadata_alpha_diversity_data <- maternalMicrobiome_InfantMetadata_results$alpha_diversity_data
maternalMicrobiome_InfantMetadata_alpha_diversity_test <- maternalMicrobiome_InfantMetadata_results$alpha_diversity_test
maternalMicrobiome_InfantMetadata_alpha_diversity_plots <- maternalMicrobiome_InfantMetadata_results$diversity_test_plots
maternalMicrobiome_InfantMetadata_alpha_diversity_barplot <- maternalMicrobiome_InfantMetadata_results$overall_diversity_barplot

saveRDS(
  maternalMicrobiome_InfantMetadata_alpha_diversity_barplot,
  file = "analysis/diversity/maternalMicrobiome_InfantMetadata_alpha_diversity_barplot_gg_associations.RDS",
)
saveRDS(
  maternalMicrobiome_InfantMetadata_alpha_diversity_plots,
  file = "analysis/maternalMicrobiome_InfantMetadata_alpha_diversity_plots_gg_associations.RDS",
)

pdf(
  paste('analysis/diversity', '/maternalMicrobiome_infantMetadata_alpaha_diversity.pdf', sep = ''),
  width = 2.4,
  height = 2.25,
  onefile = TRUE
)

for (meta in unique(colnames(infant_microbiome_data_metadata_maternal))){
  tryCatch({
    stdout <-
      capture.output(print(maternalMicrobiome_InfantMetadata_alpha_diversity_plots[[meta]]), type = "message")
  }, error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dev.off()
dev.off()
saveRDS(
  maternalMicrobiome_InfantMetadata_alpha_diversity_barplot,
  file = "analysis/diversity/maternalMicrobiome_InfantMetadata_alpha_diversity_barplot_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/maternalMicrobiome_InfantMetadata_alpha_diversity_barplot.png', plot=maternalMicrobiome_InfantMetadata_alpha_diversity_barplot,width = 2.5, height = 1.8, units = "in", dpi = 300)
ggsave(filename='analysis/diversity/maternalMicrobiome_InfantMetadata_alpha_diversity_barplot.pdf', plot=maternalMicrobiome_InfantMetadata_alpha_diversity_barplot,width = 2.5, height = 1.8, units = "in", dpi = 300)

dev.off()
dev.off()
fig2 <- ggdraw() +
  draw_plot(maternal_overall_diversity_barplot + coord_flip(),
            x = 0, y = 0, width = .22, height = 1) +
  draw_plot(infantMicrobiome_maternalMetadata_alpha_diversity_barplot + coord_flip() + 
              theme(axis.text.y=element_blank(),
                    axis.line.y = element_blank(),
                    axis.ticks.y=element_blank()),
              x = 0.21, y = 0, width = .08, height = 1)  +
  draw_plot(maternal_overall_diversity_plots$BMI + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), 
    x = .3, y = .5, width = .1, height = .5)+
  draw_plot(maternal_overall_diversity_plots$GDM+ theme(
    #axis.title.x = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .3, y = 0, width = .1, height = .5)+
  draw_plot(infant_overall_diversity_barplot + coord_flip(),
            x = .4, y = 0, width = .22, height = 1) +
  draw_plot(maternalMicrobiome_InfantMetadata_alpha_diversity_barplot + coord_flip() + 
              theme(axis.text.y=element_blank(),
                    axis.line.y = element_blank(),
                    axis.ticks.y=element_blank()),
            x = 0.61, y = 0, width = .08, height = 1) +
  draw_plot(maternalMicrobiome_InfantMetadata_alpha_diversity_plots$Sex + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .7, y = 0, width = .15, height = .5)+
  draw_plot(results_maternal_donor_alpha_divesity$diversity_test_plots$`Maternal vs. donor` + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .7, y = .5, width = .15, height = .5) +
  draw_plot(infant_alpha_diversity_plots$Visit + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 5)), x = .85, y = 0, width = .15, height = .5)+
  draw_plot(maternalMicrobiome_InfantMetadata_alpha_diversity_plots$`Fortifer in feed` + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 5)), x = .85, y = .5, width = .15, height = .5)+
  draw_plot_label((label = c("a", "b", "c", "d", "e", "f", "g", "h")),
                  size = 8,x = c(0, .3, .3, .4, .7, .7, .85, .85), y = c(1, 1, 0.5, 1, 1, .5, 1, .5))
fig2
ggsave(filename = 'manuscript/figures/fig3_diversity/fig2.pdf', plot=fig2, width = 183, height = 60, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/figures/fig3_diversity/fig2.png', plot=fig2, width = 183, height = 60, units = "mm", dpi = 350)



