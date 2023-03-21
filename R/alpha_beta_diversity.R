#library(ape)
library(vegan)
library(omicsArt)
library(Tweedieverse)

# run omicsprocees.R script lines 1-140
microbiome_data <- as.data.frame(t(mother_microbiome))
microbiome_data_metadata <- mother_metadata[rownames(microbiome_data),]


microbiome_data_metadata$Breast_milk_collected <- infant_metadata[gsub("_m","", rownames(microbiome_data_metadata)), "Breast_milk_collected"]
microbiome_data_metadata$Gestational_age_at_birth_week <- infant_metadata[rownames(microbiome_data_metadata), "Gestational_age_at_birth_week"]
microbiome_data_metadata$Day_of_life_at_time_of_breast_milk_sample <- infant_metadata[rownames(microbiome_data_metadata), "Day_of_life_at_time_of_breast_milk_sample"]

results_maternal_donor_alpha_divesity <- omicsArt::diversity_test(microbiome_data, microbiome_data_metadata[,c("Breast_milk_collected"), drop=F], output = NA)

#veg_dist <- as.matrix(vegdist(microbiome_data, method="bray"))
#write.table(veg_dist, 'data/mother_microbiome_bray_dist.txt',
#            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

#diversity_data$Visit <- microbiome_data_metadata[rownames(microbiome_data),"Visit"] #"Country of living from birth" Visit  "Breast_milk_collected"
#diversity_data$Visit <- microbiome_data_metadata[rownames(microbiome_data),"Day_of_life_at_time_of_breast_milk_sample"]
microbiome_data_metadata_maternal <- microbiome_data_metadata[rownames(microbiome_data_metadata[microbiome_data_metadata$Breast_milk_collected=="Maternal", ]),]

maternal_microbiome <- microbiome_data[rownames(microbiome_data_metadata[microbiome_data_metadata$Breast_milk_collected=="Maternal", ]),]
microbiome_data_metadata_maternal <- microbiome_data_metadata_maternal[,colSums(is.na(microbiome_data_metadata_maternal))<nrow(microbiome_data_metadata_maternal)]
temp_filtered_metadata <- microbiome_data_metadata_maternal[, apply(microbiome_data_metadata_maternal, 2, entropy) > 0.5]
Excluded_metadata <- setdiff(colnames(metadata), colnames(temp_filtered_metadata))
logging::loginfo("Excluded metadata with entropy less or equal to 0.75: %s", Excluded_metadata)
microbiome_data_metadata_maternal <- temp_filtered_metadata
microbiome_data_metadata_maternal <- microbiome_data_metadata_maternal[rownames(microbiome_data),]
microbiome_data_metadata_maternal <- microbiome_data_metadata_maternal[, ! colnames(microbiome_data_metadata_maternal) %in% 
                                                                         c("External_ID", "BMI_category", "Maternal_smoking")]



results <- omicsArt::diversity_test(microbiome_data, microbiome_data_metadata_maternal, output = NA)
alpha_diversity_data <- results$alpha_diversity_data
alpha_diversity_test <- results$alpha_diversity_test
alpha_diversity_plots <- results$diversity_test_plots
overall_diversity_barplot <- results$overall_diversity_barplot


saveRDS(
  alpha_diversity_plots,
  file = "Analysis/diversity/alpha_diversity_maternal_gg_associations.RDS",
)

pdf(
  paste('analysis/diversity', '/mathernal_alpaha_diversity.pdf', sep = ''),
  width = 2.4,
  height = 2,
  onefile = TRUE
)

for (meta in unique(colnames(mother_metadata))){
  tryCatch({
    stdout <-
      capture.output(print(alpha_diversity_plots[[meta]]), type = "message")
  }, error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dev.off()
dev.off()
saveRDS(
  overall_diversity_barplot,
  file = "Analysis/diversity/overall_diversity_barplot_maternal_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/overall_diversity_barplot_maternal.png', plot=overall_diversity_barplot, width = 3, height = 1.8, units = "in", dpi = 300)

###### infants
microbiome_data <- as.data.frame(t(infant_microbiome))
microbiome_data_metadata <- infant_metadata[rownames(microbiome_data),]

#microbiome_data_metadata$Breast_milk_collected <- infant_metadata[rownames(microbiome_data_metadata), "Breast_milk_collected"]
#microbiome_data_metadata$Gestational_age_at_birth_week <- infant_metadata[rownames(microbiome_data_metadata), "Gestational_age_at_birth_week"]
#microbiome_data_metadata$Day_of_life_at_time_of_breast_milk_sample <- infant_metadata[rownames(microbiome_data_metadata), "Day_of_life_at_time_of_breast_milk_sample"]
#microbiome_data_metadata$Breast_milk_collected[microbiome_data_metadata$Breast_milk_collected=="1"] <- "Maternal breast milk"
#microbiome_data_metadata$Breast_milk_collected[microbiome_data_metadata$Breast_milk_collected=="2"] <- "Donor breast milk"
veg_dist <- as.matrix(vegdist(microbiome_data, method="bray"))



#infant_diversity_data$Visit <- microbiome_data_metadata[rownames(microbiome_data),"Visit"] #"Country of living from birth" Visit  "Breast_milk_collected"
#diversity_data$Visit <- microbiome_data_metadata[rownames(microbiome_data),"Day_of_life_at_time_of_breast_milk_sample"]
infant_maternal_microbiome <- microbiome_data[rownames(microbiome_data_metadata[microbiome_data_metadata$Breast_milk_collected=="Maternal", ]),]
infant_microbiome_data_metadata_maternal <- microbiome_data_metadata[rownames(infant_maternal_microbiome),]

infant_microbiome_data_metadata_maternal <- infant_microbiome_data_metadata_maternal[,colSums(is.na(infant_microbiome_data_metadata_maternal))<nrow(microbiome_data_metadata_maternal)]
temp_filtered_metadata <- infant_microbiome_data_metadata_maternal[, apply(infant_microbiome_data_metadata_maternal, 2, entropy) > 0.5]
Excluded_metadata <- setdiff(colnames(metadata), colnames(temp_filtered_metadata))
logging::loginfo("Excluded metadata with entropy less or equal to 0.75: %s", Excluded_metadata)
infant_microbiome_data_metadata_maternal <- temp_filtered_metadata
infant_microbiome_data_metadata_maternal <- infant_microbiome_data_metadata_maternal[rownames(infant_maternal_microbiome),]



infant_results <- omicsArt::diversity_test(infant_maternal_microbiome, infant_microbiome_data_metadata_maternal)
infant_alpha_diversity_data <- infant_results$alpha_diversity_data
infant_alpha_diversity_test <- infant_results$alpha_diversity_test
infant_alpha_diversity_plots <- infant_results$diversity_test_plots
infant_overall_diversity_barplot <- infant_results$overall_diversity_barplot

pdf(
  paste('analysis/diversity', '/infnat_mathernal_alpaha_diversity.pdf', sep = ''),
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
  file = "Analysis/diversity/infant_overall_diversity_barplot_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/infant_overall_diversity_barplot.png', plot=infant_overall_diversity_barplot, width = 3, height = 1.8, units = "in", dpi = 300)


infant_microbiomre_maternal_metadata_results <- omicsArt::diversity_test(infant_maternal_microbiome, microbiome_data_metadata_maternal)
infantMicrobiome_maternalMetadata_alpha_diversity_data <- infant_microbiomre_maternal_metadata_results$alpha_diversity_data
infantMicrobiome_maternalMetadata_alpha_diversity_test <- infant_microbiomre_maternal_metadata_results$alpha_diversity_test
infantMicrobiome_maternalMetadata_alpha_diversity_plots <- infant_microbiomre_maternal_metadata_results$diversity_test_plots
infantMicrobiome_maternalMetadata_alpha_diversity_barplot <- infant_microbiomre_maternal_metadata_results$overall_diversity_barplot



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
  file = "Analysis/diversity/infantMicrobiome_maternalMetadata_alpha_diversity_barplot_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/infantMicrobiome_maternalMetadata_alpha_diversity_barplot.png', plot=infantMicrobiome_maternalMetadata_alpha_diversity_barplot, width = 3, height = 1.8, units = "in", dpi = 300)

microbiome_data <- microbiome_data[rownames(infant_microbiome_data_metadata_maternal),]
maternalMicrobiome_InfantMetadata_results <- omicsArt::diversity_test(microbiome_data, infant_microbiome_data_metadata_maternal)
maternalMicrobiome_InfantMetadata_alpha_diversity_data <- maternalMicrobiome_InfantMetadata_results$alpha_diversity_data
maternalMicrobiome_InfantMetadata_alpha_diversity_test <- maternalMicrobiome_InfantMetadata_results$alpha_diversity_test
maternalMicrobiome_InfantMetadata_alpha_diversity_plots <- maternalMicrobiome_InfantMetadata_results$diversity_test_plots
maternalMicrobiome_InfantMetadata_alpha_diversity_barplot <- maternalMicrobiome_InfantMetadata_results$overall_diversity_barplot



pdf(
  paste('analysis/diversity', '/maternalMicrobiome_InfnatMetadata_alpaha_diversity.pdf', sep = ''),
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
  file = "Analysis/diversity/maternalMicrobiome_InfantMetadata_alpha_diversity_barplot_gg_associations.RDS",
)
ggsave(filename='analysis/diversity/maternalMicrobiome_InfantMetadata_alpha_diversity_barplot.png', plot=maternalMicrobiome_InfantMetadata_alpha_diversity_barplot, width = 3, height = 1.8, units = "in", dpi = 300)
