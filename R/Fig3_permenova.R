library(vegan)
library(permute)
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")

#infant
infant_microbiome <- as.data.frame(t(infant_microbiome))
infant_metadata2 = infant_metadata#[infant_metadata$Visit=="a",]
infant_metadata2 <- infant_metadata2[, apply(infant_metadata2, 2, Tweedieverse::entropy) > 0.5]
intersect_samples <- intersect(rownames(infant_microbiome2), rownames(infant_metadata2))
infant_metadata2 = infant_metadata[intersect_samples,]
infant_microbiome2 = infant_microbiome2[intersect_samples,]


M_perm <- matrix(NA, nrow=length(infant_metadata2), ncol=3)

rownames(M_perm) <- colnames(infant_metadata2)
colnames(M_perm) <- c("Breast milk microbiome", "Infant stool microbiome", "Infant stool metabolite")
Nperms <- 4999
R_perm = P_perm = M_perm

#infant_metabolite [is.na(infant_metabolite)] <- 0
#infant_microbiome2 [is.na(infant_microbiome2)] <- 0
#infant_microbiome2 <- infant_microbiome2[ , colSums(is.na(infant_microbiome2)) < 40]
#infant_metabolite <- infant_metabolite[ , colSums(is.na(infant_metabolite)) < 40]
omic <- "Infant stool microbiome"
data <- infant_microbiome2
for (meta in colnames(infant_metadata2)){
  #if (meta =="Visit time point") next
  tryCatch({
    ad <- adonis2(data ~ . , data= infant_metadata2[, meta, drop=F] , permutations=4999, mehod='euclidean')
    temp_R <- ad$R2 #aov.tab$R2
    #names(temp_R) <- rownames(ad$aov.tab)
    #n <- length(temp_R)
    temp_P <- ad$`Pr(>F)` # aov.tab$`Pr(>F)`
    P_perm[meta, omic] <- temp_P[1]
    R_perm[meta, omic] <- temp_R[1]

  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
infant_metadata2 = infant_metadata#[infant_metadata$Visit=="a",]
infant_metadata2 <- infant_metadata2[, apply(infant_metadata2, 2, Tweedieverse::entropy) > 0.5]
intersect_samples <- intersect(rownames(infant_metabolite), rownames(infant_metadata2))
infant_metadata2 = infant_metadata2[intersect_samples,]
infant_metabolite = infant_metabolite[intersect_samples,]
omic <- "Infant stool metabolite"
data <- infant_metabolite
data[is.na(data)]<-0
for (meta in colnames(infant_metadata2)){
  #if (meta =="Visit") next
  tryCatch({
    ad <- adonis2(data ~ . , data= infant_metadata2[, meta, drop=F] , permutations=4999, mehod='euclidean')
    temp_R <- ad$R2  #aov.tab$R2
    names(temp_R) <- rownames(ad$aov.tab)
    n <- length(temp_R)
    temp_P <- ad$`Pr(>F)` #aov.tab$`Pr(>F)`
    P_perm[meta, omic] <- temp_P[1]
    R_perm[meta, omic] <- temp_R[1]

  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}

#####################Mother###############################
mother_microbiome_data <- as.data.frame(t(mother_microbiome))
infant_metadata2 = infant_metadata#[infant_metadata$Visit=="a",]
infant_metadata2 <- infant_metadata2[, apply(infant_metadata2, 2, Tweedieverse::entropy) > 0.5]
intersect_samples <- intersect(rownames(mother_microbiome_data), rownames(infant_metadata2))
infant_metadata2 = infant_metadata2[intersect_samples,]
mother_microbiome_data = mother_microbiome_data[intersect_samples,]
omic <- "Breast milk microbiome"
data <- mother_microbiome_data

for (meta in colnames(infant_metadata2)){
  #if (meta =="Visit") next
  tryCatch({
    ad <- adonis2(data ~ . , data= infant_metadata2[, meta, drop=F] , permutations=4999, mehod='euclidean')
    temp_R <- ad$R2
    #names(temp_R_mother) <- rownames(ad$aov.tab)
    #n <- length(temp_R_mother)
    temp_P <- ad$`Pr(>F)`
    P_perm[meta, omic] <- temp_P[1]
    R_perm[meta, omic] <- temp_R[1]

  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}

R_perm[is.na(R_perm)] <- 0.001
P_perm[is.na(P_perm)] <- 0.99


omnibus_heatmap_Infant <- omicsArt::ps2heatmap(t(R_perm), t(P_perm), FDR = F)
ggsave(filename=paste('manuscript/figures/fig3_PERMANOVA/omnibus_heatmap_Infant_euclidean_FDR_Va.pdf', sep=''), plot=omnibus_heatmap_Infant,
       width = 8, height = .75, units = "in", dpi = 350)

write.table( P_perm,"manuscript/figures/fig3_PERMANOVA/omnibus_heatmap_Infant_P.txt",
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

write.table( R_perm,"manuscript/figures/fig3_PERMANOVA/omnibus_heatmap_Infant_R.txt",
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

##########################Maternal Metadata #######################


M_perm <- matrix(NA, nrow=length(mother_metadata), ncol=3)

rownames(M_perm) <- colnames(mother_metadata)
colnames(M_perm) <- c("Breast milk microbiome", "Infant stool microbiome", "Infant stool metabolite")
R_perm = P_perm = M_perm

#infant
infant_microbiome2 <- infant_microbiome
mother_metadata2 = mother_metadata#[mother_metadata$Visit=="a",]

intersect_samples <- intersect(rownames(infant_microbiome2), rownames(mother_metadata2))
mother_metadata2 = mother_metadata2[intersect_samples,]
infant_microbiome2 = infant_microbiome2[intersect_samples,]

omic <- "Infant stool microbiome"
data <- infant_microbiome2
for (meta in colnames(mother_metadata2)){
  #if (meta =="Visit") next
  tryCatch({
    mother_metadata_complete  <- mother_metadata2[rowSums(is.na(mother_metadata2[, meta, drop=F])) == 0,]
    data_complete <- data[row.names(mother_metadata_complete),]
    ad <- adonis2(data_complete ~ . , data= mother_metadata_complete[, meta, drop=F] , permutations=4999, mehod='euclidean')
    temp_R <- ad$R2 #aov.tab$R2
    #names(temp_R) <- rownames(ad$aov.tab)
    #n <- length(temp_R)
    temp_P <- ad$`Pr(>F)` # aov.tab$`Pr(>F)`
    P_perm[meta, omic] <- temp_P[1]
    R_perm[meta, omic] <- temp_R[1]

  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}

intersect_samples <- intersect(rownames(infant_metabolite), rownames(mother_metadata))
mother_metadata2 = mother_metadata[intersect_samples,]
infant_metabolite = infant_metabolite[intersect_samples,]
omic <- "Infant stool metabolite"
data <- infant_metabolite
data[is.na(data)]<-0
for (meta in colnames(mother_metadata2)){
  #if (meta =="Visit") next
  tryCatch({
    mother_metadata_complete  <- mother_metadata2[rowSums(is.na(mother_metadata2[, meta, drop=F])) == 0,]
    data_complete <- data[row.names(mother_metadata_complete),]
    ad <- adonis2(data_complete ~ . , data= mother_metadata_complete[, meta, drop=F] , permutations=4999, mehod='euclidean')
    temp_R <- ad$R2  #aov.tab$R2
    names(temp_R) <- rownames(ad$aov.tab)
    n <- length(temp_R)
    temp_P <- ad$`Pr(>F)` #aov.tab$`Pr(>F)`
    P_perm[meta, omic] <- temp_P[1]
    R_perm[meta, omic] <- temp_R[1]

  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}

# Mother
mother_microbiome_data <- as.data.frame(t(mother_microbiome))
#mother_metadata2 = mother_metadata[mother_metadata$Visit=="a",]
intersect_samples <- intersect(rownames(mother_microbiome_data), rownames(mother_metadata2))
mother_metadata2 = mother_metadata[intersect_samples,]
mother_microbiome_data = mother_microbiome_data[intersect_samples,]
omic <- "Breast milk microbiome"
data <- mother_microbiome_data
data[is.na(data)]<-0

for (meta in colnames(mother_metadata2)){
  if (meta =="Visit") next
  tryCatch({
    mother_metadata_complete  <- mother_metadata2[rowSums(is.na(mother_metadata2[, meta, drop=F])) == 0,]
    data_complete <- data[row.names(mother_metadata_complete),]
    ad <- adonis2(data_complete ~ . , data= mother_metadata_complete[, meta, drop=F] , permutations=4999, mehod='euclidean')
    temp_R <- ad$R2
    #names(temp_R_mother) <- rownames(ad$aov.tab)
    n <- length(temp_R)
    temp_P <- ad$`Pr(>F)`
    P_perm[meta, omic] <- temp_P[1]
    R_perm[meta, omic] <- temp_R[1]

  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}

omnibus_heatmap_maternal <- omicsArt::ps2heatmap(t(R_perm), t(P_perm), FDR = T)
ggsave(filename=paste('manuscript/figures/fig3_PERMANOVA/omnibus_heatmap_maternal_euclidean_FDR.pdf', sep=''), plot=omnibus_heatmap_maternal,
       width = 8, height = .75, units = "in", dpi = 350)
write.table( P_perm,"manuscript/figures/fig3_PERMANOVA/omnibus_heatmap_maternal_P.txt",
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

write.table( R_perm,"manuscript/figures/fig3_PERMANOVA/omnibus_heatmap_maternal_R.txt",
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

###############################################################################
##############combine figures ##############################
fig_permenova <-   plot_grid(omnibus_heatmap_Infant,
                             omnibus_heatmap_maternal,
                             #rel_widths = c(1, 1, 1,1),
                             labels = c('a, infnat clinical information', 'b, maternal clinical information'),
                             label_size = 7, ncol = 1, nrow = 2)

fig_permenova
ggsave(filename=paste('manuscript/figures/fig3_PERMANOVA/fig3_omnibus_heatmap_euclidean_FDR.pdf', sep=''), plot=fig_permenova,
       width = 7.2, height = 3.75, units = "in", dpi = 350)

ggsave(filename=paste('manuscript/figures/fig3_PERMANOVA/fig3_omnibus_heatmap_euclidean_FDR.png', sep=''), plot=fig_permenova,
       width = 7.2, height = 3.75, units = "in", dpi = 350)
##############################




