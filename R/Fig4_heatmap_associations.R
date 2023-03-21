# heatmap of Tweedieverse results
setwd("~/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk")
colnames(infant_metadata) <- gsub(" ", "_", colnames(infant_metadata))

time_based <- c( "Visit", "Gestational_age_week", "Gestational_age_day", "Gestational_age_category", "Gestational_by_weight_percentage", "Gestational_by_weight_category",
                 "Sex", "Days_in_NICU", "Breast_milk_collected", "Diet_at_sample", "All_medications_infant", "Recent_follow_up_date")

not_time_based <- c("Birth_weight", "Birth_weight_category", "Breast_milk_sample_date", "Age_at_breast_milk_sample", "Sample_fortified_with", "Stool_sample_date", "Age_stool_sample", "Method_of_feeding", 
                    "Last_antibiotic_received", "Gastric_acid_suppresion", "GI_related_diagnosis",  "GI_surgery_prior", "Metabolic_Condition", "Cystic_Fibrosis",   "Underwent_sepsis", "Diagnosed_with_Sepsis", 
                    "Infection_due_resistant",  "MRSA_swab_result", "VRE_swab_result", "Problems_at_discharge",  "Other_medical_conditions")

infant_Metabolite_results <- NULL
for (meta in colnames(infant_metadata)){
  #if (meta %in% time_based) next
  tryCatch({
    temp_res <- read.delim(
      paste0('analysis/Tweedieverse/Infant/Metabolite/Tweedieverse_meatbolites_', meta, '/all_results.tsv',sep = ""),
      #paste('/Users/rah/Box/HMH-CDI COVID Data/Analysis/maaslin2_metadata_vs_microbiome/maaslin2_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      #row.names = NA
    )
    infant_Metabolite_results <- rbind(infant_Metabolite_results, temp_res)
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dim(infant_Metabolite_results)
write.table( infant_Metabolite_results,"Manuscript/figures/fig#_associations/Tweedieverse_meatbolites_infant_results.tsv",
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
for (meta in colnames(infant_metadata)){
  #if (meta %in% not_time_based) next
  tryCatch({
    temp_res <- read.delim(
      paste0('analysis/Tweedieverse/Infant/Metabolite/Tweedieverse_meatbolites_visit_a_', meta, '/all_results.tsv',sep = ""),
      #paste('/Users/rah/Box/HMH-CDI COVID Data/Analysis/maaslin2_metadata_vs_microbiome/maaslin2_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      #row.names = NA
    )
    infant_Metabolite_results <- rbind(infant_Metabolite_results, temp_res)
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dim(infant_Metabolite_results)
#[1] 924  12

df1 = infant_Metabolite_results[infant_Metabolite_results$qval<=0.1,]
df1 <- df1[order(df1$qval,decreasing = F),]
dim(df1)
#[1] 260  12
infant_Metabolite_associations_heatmap <- Tweedieverse::Tweedieverse_heatmap(df1,
  title = "-log(q-value)",
  cell_value = 'qval',
  data_label = 'data',
  metadata_label = 'metadata',
  border_color = 'grey93',
  color = colorRampPalette(c("darkblue", "grey90", "darkred")),
  col_rotate = 90,
  first_n = 50,
  write_to = NA
)

##use Tweedieverse_heatmap in viz.R Tweedieverse
saveRDS(infant_Metabolite_associations_heatmap, file = paste('manuscript/figures/fig#_associations', "/infnat_metabolites_gg_heatmap_associations.RDS", sep = ""))
ggsave(filename=paste('manuscript/figures/fig#_associations', "/infant_metabolites_gg_heatmap_associations.png", sep = ""), plot=associations_heatmap$gtable,
       width = 5, height = 7, units = "in", dpi = 350)

infant_Microbiome_results <- NULL
for (meta in colnames(infant_metadata)){
  #if (meta %in% time_based) next
  tryCatch({
    temp_res <- read.delim(
      paste0('analysis/Tweedieverse/Infant/Microbiome/Tweedieverse_infant_microbiome_', meta, '/all_results.tsv',sep = ""),
      #paste('/Users/rah/Box/HMH-CDI COVID Data/Analysis/maaslin2_metadata_vs_microbiome/maaslin2_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      #row.names = NA
    )
    infant_Microbiome_results <- rbind(infant_Microbiome_results, temp_res)
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dim(infant_Microbiome_results)
write.table( infant_Microbiome_results,"Manuscript/figures/fig#_associations/Tweedieverse_microbiome_infant_results.tsv",
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
for (meta in colnames(infant_metadata)){
  if (meta %in% not_time_based) next
  tryCatch({
    temp_res <- read.delim(
      paste0('analysis/Tweedieverse/Infant/Microbiome/Tweedieverse_infant_microbiome_visit_a_', meta, '/all_results.tsv',sep = ""),
      #paste('/Users/rah/Box/HMH-CDI COVID Data/Analysis/maaslin2_metadata_vs_microbiome/maaslin2_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      #row.names = NA
    )
    infant_Microbiome_results <- rbind(infant_Microbiome_results, temp_res)
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dim(infant_Microbiome_results)
#[1] 924  12

df2 <- infant_Microbiome_results[infant_Microbiome_results$qval<=0.1,]
df2 <- df[order(df2$qval,decreasing = F),]
dim(df2)
#[1] 260  12
infant_Microbiome_associations_heatmap <- Tweedieverse::Tweedieverse_heatmap(df2,
                                                                             title = "-log(q-value)",
                                                                             cell_value = 'qval',
                                                                             data_label = 'data',
                                                                             metadata_label = 'metadata',
                                                                             border_color = 'grey93',
                                                                             color = colorRampPalette(c("darkblue", "grey90", "darkred")),
                                                                             col_rotate = 90,
                                                                             first_n = 50,
                                                                             write_to = NA
)

##use Tweedieverse_heatmap in viz.R Tweedieverse
saveRDS(infant_Microbiome_associations_heatmap, file = paste('manuscript/figures/fig#_associations', "/infnat_microbiome_gg_heatmap_associations.RDS", sep = ""))
ggsave(filename=paste('manuscript/figures/fig#_associations', "/infant_microbiome_gg_heatmap_associations.png", sep = ""), plot=infant_Microbiome_associations_heatmap$gtable,
       width = 5, height = 5, units = "in", dpi = 350)


cell_value = "Q.value"
data_label = 'Data'
metadata_label = 'Metadata'
border_color = "grey93"
color = colorRampPalette(c("blue", "whitesmoke", "red"))(51) #whitesmoke
# read MaAsLin output
gaps <- c(0, 0)
mat <- NA
Omics_color_order = list(
  "Omics" = c(
    "Metabolite" = "blue",
    "Microbiome" = "red"
  )
)

mat <- convert_maaslin_output2matrix(df1)

gaps[1] <- length(colnames(mat))
method_order = data.frame(ID = factor(rep(c("Microbiome"), each = length(colnames(mat)))))

mat2 <- convert_maaslin_output2matrix(df2)
if (!is.na(mat2) && dim(mat2)[1] > 0) {
  #mat2 <- mat2[, colnames(mat2)%in%names(bugs_to_show_assoc[["Buccal mucosa"]])]
  print ("Metabolite")
  #colnames(mat2) <- bugs_to_show_assoc[["Buccal mucosa"]]
  #mat <- cbind(mat, mat2)
  mat3 <- combine2df(mat, mat2)
  gaps[2] <- length(colnames(mat3))
  method_order = rbind(method_order, data.frame(ID = factor(rep(
    c("Metabolite"), each = length(colnames(mat2))
  ))))
}
  
colnames(method_order)[1] <- "Omics"
labels_col <- colnames(mat3)
colnames(mat3) <- rownames(method_order)
mylegend = TRUE
my_show_rownames = TRUE
mat3[is.na(mat3)] <- 0
temp <- mat3
plot_result <-
pheatmap(
  temp,
  cellwidth = 5,
  cellheight = 5,
  # changed to 3
  main = title,
  annotation_col = method_order["Omics"],
  #annotation_colors = Omics_color_order["Omics"],
  annotation_legend = T,
  fontsize = 6,
  kmeans_k = NA,
  #border=TRUE,
  labels_col = labels_col,
  show_rownames = T,
  show_colnames = T,
  scale = "none",
  #clustering_method = "complete",
  cluster_rows = FALSE,
  cluster_cols = F,
  #clustering_distance_rows = "euclidean",
  #clustering_distance_cols = "euclidean",
  legend = TRUE,
  border_color = border_color,
  color = color,
  treeheight_row = 0,
  treeheight_col = 0,
  gaps_col = gaps,
  display_numbers = matrix(ifelse(temp > 0.0, "+", ifelse(temp <= 0.0, "-", "")),  nrow(temp))
)
pheatmap(temp)
dev.off()
ggsave(
  filename = 'associations_overview.pdf',
  plot = plot_result$gtable,
  width = 12,
  limitsize = FALSE,
  height = 3,
  units = "in",
  dpi = 300
)
return (plot_result)


m_time_based <- c("Visit", "Maternal_education", "Gravidity", "Parity_after_delivery", "Maternal_smoking", "Alcohol_during_pregnancy", "Hypertensive_disorder", "GBS_positive", 
                  "PROM", "Chorioamnionitis", "Supplements_during_pregnancy")

m_not_time_based <- c("Age", "Maritial_status", "Weight_1Year_Prior", "Height", "BMI", "BMI_category", "Ethnicity", "Race", "Country_lived", "Household_income", "Method_of_delivery", "Breast_milk_collected", 
                      "Antibiotic_use_prenatally","Antibiotic_during_delivery")

breastmilk_Microbiome_results <- NULL
for (meta in colnames(mother_metadata)){
  if (meta %in% m_time_based) next
  tryCatch({
    temp_res <- read.delim(
      paste0('analysis/Tweedieverse/Breastmilk/Tweedieverse_mother_microbiome_visit_a_', meta, '/all_results.tsv',sep = ""),
      #paste('/Users/rah/Box/HMH-CDI COVID Data/Analysis/maaslin2_metadata_vs_microbiome/maaslin2_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      #row.names = NA
    )
    breastmilk_Microbiome_results <- rbind(breastmilk_Microbiome_results, temp_res)
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dim(breastmilk_Microbiome_results)

for (meta in colnames(mother_metadata)){
  if (meta %in% m_not_time_based) next
  tryCatch({
    temp_res <- read.delim(
      paste0('analysis/Tweedieverse/Breastmilk/Tweedieverse_mother_microbiome_', meta, '/all_results.tsv',sep = ""),
      #paste('/Users/rah/Box/HMH-CDI COVID Data/Analysis/maaslin2_metadata_vs_microbiome/maaslin2_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      #row.names = NA
    )
    breastmilk_Microbiome_results <- rbind(breastmilk_Microbiome_results, temp_res)
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
dim(breastmilk_Microbiome_results)
#[1] 924  12

df = breastmilk_Microbiome_results[breastmilk_Microbiome_results$qval<=0.1,]
df <- df[order(df$qval,decreasing = F),]
dim(df)
#[1] 260  12
breastmilk_Microbiome_associations_heatmap <- Tweedieverse::Tweedieverse_heatmap(df,
                                                                             title = "-log(q-value)",
                                                                             cell_value = 'qval',
                                                                             data_label = 'data',
                                                                             metadata_label = 'metadata',
                                                                             border_color = 'grey93',
                                                                             color = colorRampPalette(c("darkblue", "grey90", "darkred")),
                                                                             col_rotate = 90,
                                                                             first_n = 50,
                                                                             write_to = NA
)

##use Tweedieverse_heatmap in viz.R Tweedieverse
saveRDS(breastmilk_Microbiome_associations_heatmap, file = paste('manuscript/figures/fig#_associations', "/breastmilk_microbiome_gg_heatmap_associations.RDS", sep = ""))
ggsave(filename=paste('manuscript/figures/fig#_associations', "/breastmilk_microbiome_gg_heatmap_associations.png", sep = ""), plot=breastmilk_Microbiome_associations_heatmap$gtable,
       width = 5, height = 5, units = "in", dpi = 350)



# ## read association
association_Azithromycin <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/hackensack/analysis/Tweedieverse/Tweedieverse_Azithromycin/figures/Azithromycin_gg_associations.RDS")
association_Hypertension <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/hackensack/analysis/Tweedieverse/Tweedieverse_Hypertension/figures/Hypertension_gg_associations.RDS")
association_Diabetes <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/hackensack/analysis/Tweedieverse/Tweedieverse_Diabetes/figures/Diabetes_gg_associations.RDS")
association_acetaminophen <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/hackensack/analysis/Tweedieverse/Tweedieverse_acetaminophen/figures/acetaminophen_gg_associations.RDS")

# do plots
fig_microbiome_metadata_association <- ggdraw() +
  draw_plot(associations_heatmap$gtable,
            x = 0, y = .0, width = .8, height = 1)+ 
  draw_plot(association_Azithromycin[[1]] + scale_y_log10() , x = .82 , y = 0, width = .2, height = .25) +
  draw_plot(association_Hypertension[[1]] + scale_y_log10() + theme(
    axis.title.x = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 5),
    axis.text.y = element_text(size = 4)), x = .82, y = .25, width = .2, height = .25) +
  draw_plot(association_Diabetes[[1]]+ scale_y_log10() + theme(
    axis.title.x = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 5),
    axis.text.y = element_text(size = 4)), x = .82, y = 0.5, width = .2, height = .25) +
  draw_plot(association_acetaminophen[[1]]+ scale_y_log10() + theme(
    axis.title.x = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 5),
    axis.text.y = element_text(size = 4)), x = .82, y = 0.75, width = .2, height = .25) +

  draw_plot_label((label = c("a", "b", "c", "d", "e")),
                  size = 7,x = c(0, .82, .82, .82, 0.82), y = c(1, 1, .75, .5, .25))
fig_microbiome_metadata_association

ggsave(filename = paste0('manuscript/figures/fig#_associations', "/Fig5.png", sep = ""), plot=fig_microbiome_metadata_association, width = 183, height = 120, units = "mm", dpi = 350)
