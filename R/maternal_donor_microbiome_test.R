library(ape)



library(vegan)
moms <- as.data.frame(t(mother_microbiome))
moms_metadata <- mother_metadata[rownames(moms),]

moms_metadata$Breast_milk_collected <- infant_metadata[rownames(moms_metadata), "Breast_milk_collected"]
moms_metadata$Gestational_age_at_birth_week <- infant_metadata[rownames(moms_metadata), "Gestational_age_at_birth_week"]
moms_metadata$Breast_milk_collected[moms_metadata$Breast_milk_collected=="1"] <- "Maternal breast milk"
moms_metadata$Breast_milk_collected[moms_metadata$Breast_milk_collected=="2"] <- "Donor breast milk"
moms_metadata$Breast_milk_collected[moms_metadata$Breast_milk_collected=="2"] <- "Donor breast milk"
veg_dist <- as.matrix(vegdist(moms, method="bray"))
species_1 <- moms / 100.00 #max(species)
species_2 <- species_1[,colSums(species_1>0.01) > nrow(species_1)*.01]
# 9- write the  in you computer as a tab-delimited file
fakepcl <- list(meta=as.data.frame(moms_metadata), x=as.matrix((species_2)),
                ns=dim(species_2)[1], nf=dim(species_2)[2])

meta <- "Breast_milk_collected"
 tryCatch({
  Tweedieverse::Tweedieverse(moms,
                             moms_metadata[, meta, drop=F],
                             paste('analysis/Tweedieverse_moms_microbiome_', meta, sep = ""),
                             max_significance = 0.62,
                             base_model = "CPLM",
                             prev_threshold = 0.0,
                             abd_threshold = 0.0,
                             link = "log",
                             plot_heatmap = T,
                             plot_scatter = T,
                             standardize = F
  )
},
error = function(e) {
  print(meta)
  print(paste('error:', e))
})

tryCatch({
  Maaslin2(moms,
           moms_metadata[, c("Breast_milk_collected"), drop=F],
           paste('analysis/Maaslin2_moms_microbiome_', meta, sep = ""),
           min_abundance = 0.0,
           min_prevalence = 0.0,
           max_significance = .37,
           standardize = T,
           #transform = 'Log',
           analysis_method = 'LM',
           #normalization = 'TSS',
           #transform = 'NONE',
           #fixed_effects = c("Visit", "Race"),
           heatmap_first_n = 50,  #"Has_subj_ever_had_a_yeast_infection"),
  )

},
error = function(e) {
  print(meta)
  print(paste('error:', e))
})

