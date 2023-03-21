
library(Tweedieverse)

colnames(infant_metadata) <- gsub(" ", "_", colnames(infant_metadata))



infant_metabolite [is.na(infant_metabolite)] <- 0
# for (meta in colnames(infant_metadata)){
#   tryCatch({
#     
#     Tweedieverse::Tweedieverse(infant_metabolite,
#                                infant_metadata[infant_metadata$Visit=='a', meta , drop=F],
#                                paste('analysis/Tweedieverse/Infant/Metabolite/Tweedieverse_meatbolites_visit_a_', meta, sep = ""),
#                                max_significance = 0.25,
#                                plot_heatmap = T,
#                                plot_scatter = T,
#                                standardize = F,
#                                fixed_effects = meta
#     )
#     
#   },
#   error = function(e) {
#     print(meta)
#     print(paste('error:', e))
#   })
# }

for (meta in colnames(infant_metadata)){
  tryCatch({
    Tweedieverse::Tweedieverse(infant_metabolite,
                               infant_metadata[, meta, drop=F],
                               paste('analysis/Tweedieverse/Infant/Metabolite/Tweedieverse_metabolites', meta , sep = ""),
                               max_significance = 0.25,
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F,
                               fixed_effects = meta
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}

for (meta in colnames(infant_metadata)){
  tryCatch({
    #if (meta=="Date_of_breast_milk_sample") next
    Tweedieverse::Tweedieverse(infant_microbiome,
                               infant_metadata[, meta, drop=F],
                               paste('analysis/Tweedieverse/Infant/Microbiome/Tweedieverse_infant_microbiome_', meta, sep = ""),
                               max_significance = 0.1,
                               base_model = "CPLM",
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F,
                               fixed_effects = meta
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}



for (meta in colnames(infant_metadata)){
  tryCatch({
    #if (meta=="Date_of_breast_milk_sample") next
    Maaslin2::Maaslin2(infant_microbiome,
                       infant_metadata[, meta, drop=F],
                       paste('analysis/Tweedieverse/Infant/Microbiome/Maaslin2_infant_microbiome_', meta, sep = ""),
                       min_abundance = 0.0,
                       min_prevalence = 0.0,
                       max_significance = .1,
                       standardize = T,
                       transform = 'Log',
                       analysis_method = 'LM'
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
for (meta in colnames(infant_metadata)){
  tryCatch({
    Maaslin2::Maaslin2(infant_metabolite,
                       infant_metadata[, meta, drop=F],
                       paste('analysis/Tweedieverse/Infant/Metabolite/Maaslin2_metabolites', meta , sep = ""),
                       min_abundance = 0.0,
                       min_prevalence = 0.0,
                       max_significance = .1,
                       standardize = T,
                       transform = 'Log',
                       analysis_method = 'LM'
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}


#visit a only
# for (meta in colnames(infant_metadata)){
#   tryCatch({
#     #if (meta=="Date_of_breast_milk_sample") next
#     Tweedieverse::Tweedieverse(infant_microbiome,
#                                infant_metadata[infant_metadata$Visit=='a', meta, drop=F],
#                                paste('analysis/Tweedieverse/Infant/Microbiome/Tweedieverse_infant_microbiome_visit_a_', meta, sep = ""),
#                                max_significance = 0.1,
#                                base_model = "CPLM",
#                                plot_heatmap = T,
#                                plot_scatter = T,
#                                standardize = F,
#                                fixed_effects = meta
#     )
#   },
#   error = function(e) {
#     print(meta)
#     print(paste('error:', e))
#   })
# }
#meta <- "Country_of_living_from_birth"
colnames(mother_metadata) <- gsub(" ", "_", colnames(mother_metadata))

for (meta in colnames(mother_metadata)){
  if (meta=="Race")
    Tweedieverse::Tweedieverse(mother_microbiome,
                               mother_metadata[, meta , drop=F],
                               paste('analysis/Tweedieverse/Breastmilk/Tweedieverse_mother_microbiome_', meta, sep = ""),
                               max_significance = 0.25,
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F,
                               fixed_effects = meta,
                               reference = c('Race,White_caucasion')
    )
  else
  tryCatch({
    Tweedieverse::Tweedieverse(mother_microbiome,
                               mother_metadata[, meta, drop=F],
                               paste('analysis/Tweedieverse/Breastmilk/Tweedieverse_mother_microbiome_', meta, sep = ""),
                               abd_threshold = 0,
                               prev_threshold = 0.0,
                               max_significance = 0.25,
                               base_model = "CPLM",
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F,
                               fixed_effects = meta
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}


for (meta in colnames(mother_metadata)){
  if (meta=="Race")
    Maaslin2::Maaslin2(mother_microbiome,
                       mother_metadata[, meta , drop=F],
                       paste('analysis/Tweedieverse/Breastmilk/Maaslin2_mother_microbiome_', meta, sep = ""),
                       min_abundance = 0.0,
                       min_prevalence = 0.0,
                       max_significance = .1,
                       standardize = T,
                       transform = 'Log',
                       analysis_method = 'LM',
                       reference = c('Race,White_caucasion')
    )
  else
    tryCatch({
      Maaslin2::Maaslin2(mother_microbiome,
                         mother_metadata[, meta, drop=F],
                         paste('analysis/Tweedieverse/Breastmilk/Maaslin2_mother_microbiome_', meta, sep = ""),
                         min_abundance = 0.0,
                         min_prevalence = 0.0,
                         max_significance = .1,
                         standardize = T,
                         transform = 'Log',
                         analysis_method = 'LM'
      )
    },
    error = function(e) {
      print(meta)
      print(paste('error:', e))
    })
}

# Visit a
# 
# for (meta in colnames(mother_metadata)){
#   if (meta=="Race")
#     Tweedieverse::Tweedieverse(mother_microbiome,
#                                mother_metadata[mother_metadata$Visit=='a', meta , drop=F],
#                                paste('analysis/Tweedieverse/Breastmilk/Tweedieverse_mother_microbiome_visit_a_', meta, sep = ""),
#                                max_significance = 0.25,
#                                plot_heatmap = T,
#                                plot_scatter = T,
#                                standardize = F,
#                                fixed_effects = meta,
#                                reference = c('Race,White_caucasion')
#     )
#   else
#     tryCatch({
#       Tweedieverse::Tweedieverse(mother_microbiome,
#                                  mother_metadata[mother_metadata$Visit=='a', meta, drop=F],
#                                  paste('analysis/Tweedieverse/Breastmilk/Tweedieverse_mother_microbiome_visit_a_', meta, sep = ""),
#                                  abd_threshold = 0,
#                                  prev_threshold = 0.0,
#                                  max_significance = 0.25,
#                                  base_model = "CPLM",
#                                  plot_heatmap = T,
#                                  plot_scatter = T,
#                                  standardize = F,
#                                  fixed_effects = meta
#       )
#     },
#     error = function(e) {
#       print(meta)
#       print(paste('error:', e))
#     })
# }

complete_samples <- intersect(rownames(mother_metadata), rownames(infant_metadata))
mother_metadata2 <- cbind(mother_metadata[complete_samples,], infant_metadata[complete_samples,])
mother_metadata2 <- mother_metadata2[mother_metadata2$Breast_milk_collected=="Maternal",]
mother_microbiome2 <- mother_microbiome[ ,colnames(mother_microbiome) %in% complete_samples]
intersect(rownames(mother_metadata2), colnames(mother_microbiome))
meta <- "Country_of_living_from_birth"

# for (meta in colnames(mother_metadata2)){
#   tryCatch({
#     Maaslin2(mother_microbiome,
#              mother_metadata[mother_metadata$Visit=='a', meta, drop=F],
#              paste('analysis/Maaslin2/Breastmilk/Maaslin2_maternal_microbiome_visit_a_',meta, sep =""),
#              min_abundance = 0.0,
#              min_prevalence = 0.0,
#              max_significance = .25,
#              standardize = T,
#              transform = 'Log',
#              analysis_method = 'LM',
#              #normalization = 'AST',
#              #transform = 'NONE',
#              #fixed_effects = c("Visit", "Race"),
#              heatmap_first_n = 50,  #"Has_subj_ever_had_a_yeast_infection"),
#     )
#     
#   },
#   error = function(e) {
#     print(meta)
#     print(paste('error:', e))
#   })
# }


##Cross analysis
for (meta in colnames(infant_metadata)){
  tryCatch({
    Maaslin2::Maaslin2(mother_microbiome,
             infant_metadata[, meta, drop=F],
             paste('analysis/Tweedieverse/Cross/Breastmilk/Maaslin2_maternal_microbiome_infnat_',meta, sep =""),
             min_abundance = 0.0,
             min_prevalence = 0.0,
             max_significance = .1,
             standardize = T,
             transform = 'Log',
             analysis_method = 'LM',
             #normalization = 'AST',
             #transform = 'NONE',
             #fixed_effects = c("Visit", "Race"),
             heatmap_first_n = 50,  #"Has_subj_ever_had_a_yeast_infection"),
    )
    
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
for (meta in colnames(infant_metadata)){
  tryCatch({
    Tweedieverse::Tweedieverse(mother_microbiome,
                               infant_metadata[, meta, drop=F],
                               paste('analysis/Tweedieverse/Cross/Breastmilk/Tweedieverse_mother_microbiome_infant_', meta, sep = ""),
                               max_significance = 0.1,
                               base_model = "CPLM",
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
for (meta in colnames(infant_metadata)){
  tryCatch({
    Tweedieverse::Tweedieverse(mother_microbiome,
                               infant_metadata[, meta, drop=F],
                               paste('analysis/Tweedieverse/Cross/Breastmilk/Tweedieverse_mother_microbiome_infant_', meta, sep = ""),
                               max_significance = 0.25,
                               base_model = "CPLM",
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F)
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}

for (meta in colnames(mother_metadata)){
  tryCatch({
    Maaslin2::Maaslin2(infant_microbiome,
             mother_metadata[, meta, drop=F],
             paste('analysis/Tweedieverse/Cross/Stool/Maaslin2_infant_microbiome_mother_',meta, sep =""),
             min_abundance = 0.0,
             min_prevalence = 0.0,
             max_significance = .1,
             standardize = T,
             transform = 'Log',
             analysis_method = 'LM',
             #normalization = 'AST',
             #transform = 'NONE',
             #fixed_effects = c("Visit", "Race"),
             heatmap_first_n = 50,  #"Has_subj_ever_had_a_yeast_infection"),
    )

  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
for (meta in colnames(mother_metadata)){
  tryCatch({
    Tweedieverse::Tweedieverse(infant_microbiome,
                               mother_metadata[, meta, drop=F],
                               paste('analysis/Tweedieverse/Stool/Breastmilk/Tweedieverse_infant_microbiome_mother_', meta, sep = ""),
                               max_significance = 0.25,
                               base_model = "CPLM",
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F)
    
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}


# for (meta in colnames(mother_metadata)){
#   tryCatch({
#     Maaslin2(infant_metabolite,
#              mother_metadata[, meta, drop=F],
#              paste('analysis/Maaslin2/Cross/Stool/Maaslin2_infant_metabolite_mother_',meta, sep =""),
#              min_abundance = 0.0,
#              min_prevalence = 0.0,
#              max_significance = .25,
#              standardize = T,
#              transform = 'Log',
#              analysis_method = 'LM',
#              #normalization = 'AST',
#              #transform = 'NONE',
#              #fixed_effects = c("Visit", "Race"),
#              heatmap_first_n = 50,  #"Has_subj_ever_had_a_yeast_infection"),
#     )
#     
#   },
#   error = function(e) {
#     print(meta)
#     print(paste('error:', e))
#   })
# }
for (meta in colnames(mother_metadata)){
  tryCatch({
    Tweedieverse::Tweedieverse(infant_metabolite,
                               infant_metadata[, meta, drop=F],
                               paste('analysis/Tweedieverse/Cross/Stool/Tweedieverse_infant_metabolite_mother_', meta, sep = ""),
                               max_significance = 0.25,
                               base_model = "CPLM",
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F)
    
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
