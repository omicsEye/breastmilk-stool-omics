
bSpecies_05 <- breastmilk_species[,colSums(breastmilk_species>0.0) > nrow(breastmilk_species)*.25]
dim(bSpecies_05)

sSpecies_05 <- stool_species[,colSums(stool_species>0.0) > nrow(stool_species)*.25]
dim(sSpecies_05)
species_75 <- species[,colSums(species>0.0) > nrow(species)*.25]
dim(species_75)
species_filtered <- species[,union(colnames(sSpecies_05), union(colnames(bSpecies_05), colnames(species_75)))]

dim(species_filtered)

fakepcl <- list(meta=mom_or_infant, x=as.matrix(species_filtered),
                ns=dim(species_filtered)[2], nf=dim(species_filtered)[1])
speices_heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= T, show_colnames = F, show_rownames = T,treeheight_row = 0, treeheight_col= 5 )
#ggsave(filename='analysis/species_heatmap_05.png', plot=heat_plot, width = 12, height = 8, units = "in", dpi = 300)
#ggsave(filename='analysis/species_heatmap_05.pdf', plot=heat_plot, width = 12, height = 8, units = "in", dpi = 300)

fakepcl <- list(meta=mom_or_infant, x=as.matrix(species),
                ns=dim(species)[2], nf=dim(species)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= T, show_colnames = F, show_rownames = T,treeheight_row = 0, treeheight_col= 5)
ggsave(filename='analysis/species_heatmap.png', plot=heat_plot, width = 11, height = 18, units = "in", dpi = 300)
ggsave(filename='analysis/species_heatmap.png', plot=heat_plot, width = 11, height = 18, units = "in", dpi = 300)
ord <- omicsArt:::pcl.pcoa(fakepcl)
ord_plot <- omicsArt:::pcl.ordplot(fakepcl, ord, colour="Sample")
ggsave(filename='analysis/species_pcoa.png', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)
ggsave(filename='analysis/species_pcoa.pdf', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)
ord <- omicsArt:::pcl.tsne(fakepcl)
ord_plot <- omicsArt:::pcl.ordplot(fakepcl, ord, colour="Sample")
pcoa_theme <-   theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      plot.title =element_text(face="bold.italic"), 
                      axis.ticks.y=element_blank(),
                      axis.text.y=element_blank())
ord_plot <- ord_plot+  guides(fill = guide_legend(title = "", keywidth=0.1 ,keyheight=0.1, default.unit="inch",  override.aes = list(size=4)))+ pcoa_theme
ggsave(filename='analysis/species_tsne.png', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)
ggsave(filename='analysis/species_tsne.pdf', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)

fakepcl <- list(meta=mom_or_infant, x=as.matrix(genera),
                ns=dim(genera)[2], nf=dim(genera)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= T, show_colnames = F, show_rownames = T,treeheight_row = 0, treeheight_col= 5)
ggsave(filename='analysis/genera.png', plot=heat_plot, width = 12, height = 12, units = "in", dpi = 300)
ggsave(filename='analysis/genera.pdf', plot=heat_plot, width = 12, height = 12, units = "in", dpi = 300)

genera_05 <- genera / 100.00 #max(species)
genera_05 <- genera[,colSums(genera_05>0.0) > nrow(genera_05)*.05]
fakepcl <- list(meta=mom_or_infant, x=as.matrix(genera_05),
                ns=dim(genera_05)[2], nf=dim(genera_05)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= T, show_colnames = F, show_rownames = T,treeheight_row = 0, treeheight_col= 5)
ggsave(filename='analysis/genera_05.png', plot=heat_plot, width = 14, height = 8, units = "in", dpi = 300)
ggsave(filename='analysis/genera_05.pdf', plot=heat_plot, width = 14, height = 8, units = "in", dpi = 300)

fakepcl <- list(meta=mom_or_infant_BGC, x=as.matrix(BGCs),
                ns=dim(BGCs)[2], nf=dim(BGCs)[1])
ord <- omicsArt:::pcl.tsne(fakepcl)
ord_plot <- omicsArt:::pcl.ordplot(fakepcl, ord, colour="Sample")
ord_plot <- ord_plot+  guides(fill = guide_legend(title = "", keywidth=0.1 ,keyheight=0.1, default.unit="inch",  override.aes = list(size=4)))+ pcoa_theme
ggsave(filename='analysis/BGC_tsne.png', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)
ggsave(filename='analysis/BGC_tsne.pdf', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)


bBGCs_05 <- breastmilk_BGC[,colSums(breastmilk_BGC>0.0) > nrow(breastmilk_BGC)*.85]

sBGCs_05 <- stool_BGC[,colSums(stool_BGC>0.0) > nrow(stool_BGC)*.85]

BGCs_75 <- BGCs[,colSums(BGCs>0.0) > nrow(BGCs)*.85]
BGCs_05 <- BGCs[,union(colnames(sBGCs_05), union(colnames(bBGCs_05), colnames(BGCs_75)))]
fakepcl <- list(meta=mom_or_infant_BGC, x=as.matrix(BGCs_05),
                ns=dim(BGCs)[2], nf=dim(BGCs)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 10, meta= "Sample", show_colnames = F, show_rownames = T,treeheight_row = 0, treeheight_col= 5 )
ggsave(filename='analysis/BGCs_heatmap_85.png', plot=heat_plot, width = 16, height = 6, units = "in", dpi = 300)
ggsave(filename='analysis/BGCs_heatmap_85.pdf', plot=heat_plot, width = 16, height = 6, units = "in", dpi = 300)

## species  overall plot  ######

# Adonis/PERMANOVA omnibus association with outcomes wrt demographic variables
library(vegan)
sort(diversity(infant_metabolite, index = 'shannon'))
spa <- specaccum(infant_metabolite)
plot(spa)
mod <- decorana(infant_metabolite)
plot(mod)
dist.metabolites<-vegdist(infant_metabolite, mehod='bray')
D <-  cor(t(infant_metabolite), method = 'spearman')
library(pheatmap)
pheatmap(dist.metabolites)
pheatmap(D)

infant_metadata <- infant_metadata[ , colSums(is.na(infant_metadata)) == 0]
fakepcl <- list(meta=as.data.frame(infant_metadata), x=as.matrix(infant_metabolite),
                ns=dim(infant_metabolite)[1], nf=dim(infant_metabolite)[2])

# heatmap plot
heat_plot <- pcl.heatmap(fakepcl, sqrtspace = T, gamma = 5, meta = c( "Gestational_age_at_birth_day",
                                                                      "Birth_weight",
                                                                      "Gender_of_child",
                                                                      "Infant_diet_at_age_of_breastmilk_sample_obtained"))
heat_plot
ggsave(filename=paste('figures/heatmap_metabolites_infant.png', sep=''), plot=heat_plot,
       width = 15, height = 6, units = "in", dpi = 350)

ord <- omicsArt:::pcl.pcoa(fakepcl)
ord_plot <- omicsArt:::pcl.ordplot(fakepcl, ord, colour="Method_of_feeding")
ggsave(filename='analysis/infant_metabolites_pcoa.png', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)
ggsave(filename='analysis/infant_metabolites_pcoa.pdf', plot=ord_plot, width = 4, height = 3, units = "in", dpi = 300)
#pcl.tsne(fakepcl)
#
# wdbc.pr <- prcomp(infant_metabolite, center = F, scale = F, )
# summary(wdbc.pr)
# screeplot(wdbc.pr, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
# abline(h = 1, col="red", lty=5)
# legend("topright", legend=c("Eigenvalue = 1"),
#        col=c("red"), lty=5, cex=0.6)
# cumpro <- cumsum(wdbc.pr$sdev^2 / sum(wdbc.pr$sdev^2))
# plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
# abline(v = 6, col="blue", lty=5)
# abline(h = 0.88759, col="blue", lty=5)
# legend("topleft", legend=c("Cut-off @ PC6"),
#        col=c("blue"), lty=5, cex=0.6)
# library("factoextra")
# pdf("figures/infant_metabolites.pdf", 6, 4)
# fig <- fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 21,
#                     pointsize = 2,
#                     #fill.ind = infant_metadata$Gender_of_child,
#                     col.ind = "black",
#                     palette = "jco",
#                     addEllipses = TRUE,
#                     label = "var",
#                     col.var = "black",
#                     repel = TRUE,
#                     legend.title = "Gender of child") +
#   ggtitle("2D PCA-plot from infant metabolit dataset") +
#   theme(plot.title = element_text(hjust = 0.5))
# print(fig)
# dev.off()

metadata_all <- readxl::read_excel("data/external_metabolon_gw_matchup.xlsx")
metadata_all <- as.data.frame(metadata_all)
rownames(metadata_all) <- metadata_all$`Event ID (barcode scanned)`

rownames(metadata) <- metadata$`CLIENT IDENTIFIER`

diff_bar <- diff_bar_plot(stats_table, threshold = threshold, method = method, orderby = 'logFC', y_label ='')
ggsave(filename=paste(output_path,'/','diff_bar_CUSTOM_ATTRIBUTE_2_MO_1.pdf', sep=''), plot=diff_bar,
       width = 80, height = 60, units = "mm", dpi = 350)
volcano_p <- volcano_plot(stats_table, threshold = threshold, method = method)
ggsave(filename=paste(output_path,'/','volcano_CUSTOM_ATTRIBUTE_2_MO_1.pdf', sep=''), plot=volcano_p,
       width = 100, height = 85, units = "mm", dpi = 300)



###### infant microbial species vs metadata

fakepcl <- list(meta=mother_metadata, x=as.matrix(infant_microbiome2),
                ns=dim(infant_microbiome2)[1], nf=dim(infant_microbiome2)[2])
heat_plot <- pcl.heatmap(fakepcl, sqrtspace = T, gamma = 3, meta = c("Birth_weight", "Gender_of_child", "Metabolic_Condition", "Method_of_feeding_at_time_of_breatmilk_sample"))
# "Active_problem_list_at_time_of_discharge"
ordplots(infant_microbiome2, mother_metadata, "analysis/species_pcoa", method = "pcoa")


pcoa_plots <- ordplots(metabolites, metadata, output = 'analysis/meatbolite/', outputname = 'pcoa', method = 'pcoa')


fakepcl <- list(meta=mother_metadata_microbiome, x=as.matrix(t(mother_microbiome)),
                ns=dim(mother_microbiome)[1], nf=dim(mother_microbiome)[2])
heat_plot <- pcl.heatmap(fakepcl, sqrtspace = T, gamma = 3, meta = c("Method of delivery", "BMI", "Visit","Parity (after delivery)", "Race","GDM","PROM","Maritial Status") )
pcoa_plots <- ordplots(t(mother_microbiome), mother_metadata_microbiome, output = 'analysis/', outputname = 'mother_microbiome_pcoa', method = 'pcoa')
tsne_plots <- ordplots(metabolites, metadata, output = 'analysis/meatbolite', outputname = 'tsne_scaled', method = 'tsne')

for (meta in colnames(mother_metadata_microbiome)){
  tryCatch({
    Tweedieverse::Tweedieverse(t(mother_microbiome),
                               mother_metadata_microbiome[, meta, drop=F],
                               paste('analysis/Tweedieverse_mother_microbiome', meta, sep ="_"),
                               base_model = "CPLM",
                               prev_threshold = 0.01,
                               max_significance = .1,
                               standardize = T,
                               
                               #fixed_effects = c("Visit_type", "Race"),
                               plot_scatter =  T,  #"Has_subj_ever_had_a_yeast_infection"),
                               plot_heatmap = T
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}
###############################################################################

metadata_all <- readxl::read_excel("data/external_metabolon_gw_matchup.xlsx")
metadata_all <- as.data.frame(metadata_all)
rownames(metadata_all) <- metadata_all$`Event ID (barcode scanned)`

rownames(metadata) <- metadata$`CLIENT IDENTIFIER`

diff_bar <- diff_bar_plot(stats_table, threshold = threshold, method = method, orderby = 'logFC', y_label ='')
ggsave(filename=paste(output_path,'/','diff_bar_CUSTOM_ATTRIBUTE_2_MO_1.pdf', sep=''), plot=diff_bar,
       width = 80, height = 60, units = "mm", dpi = 350)
volcano_p <- volcano_plot(stats_table, threshold = threshold, method = method)
ggsave(filename=paste(output_path,'/','volcano_CUSTOM_ATTRIBUTE_2_MO_1.pdf', sep=''), plot=volcano_p,
       width = 100, height = 85, units = "mm", dpi = 300)



###### infant microbial species vs metadata

fakepcl <- list(meta=infant_metadata, x=as.matrix(infant_microbiome2),
                ns=dim(infant_microbiome2)[1], nf=dim(infant_microbiome2)[2])
heat_plot <- pcl.heatmap(fakepcl, sqrtspace = T, gamma = 3, meta = c("Birth_weight", "Gender_of_child", "Metabolic_Condition", "Method_of_feeding_at_time_of_breatmilk_sample"))
# "Active_problem_list_at_time_of_discharge"
ordplots(infant_microbiome2, infant_metadata, "analysis/species_pcoa", method = "pcoa")


pcoa_plots <- ordplots(metabolites, metadata, output = 'analysis/meatbolite/', outputname = 'pcoa', method = 'pcoa')


fakepcl <- list(meta=mother_metadata_microbiome, x=as.matrix(t(mother_microbiome)),
                ns=dim(mother_microbiome)[1], nf=dim(mother_microbiome)[2])
heat_plot <- pcl.heatmap(fakepcl, sqrtspace = T, gamma = 3, meta = c("Method of delivery", "BMI", "Visit","Parity (after delivery)", "Race","GDM","PROM","Maritial Status") )
pcoa_plots <- ordplots(t(mother_microbiome), mother_metadata_microbiome, output = 'analysis/', outputname = 'mother_microbiome_pcoa', method = 'pcoa')
tsne_plots <- ordplots(metabolites, metadata, output = 'analysis/meatbolite', outputname = 'tsne_scaled', method = 'tsne')

for (meta in colnames(mother_metadata_microbiome)){
  tryCatch({
    Tweedieverse::Tweedieverse(t(mother_microbiome),
                               mother_metadata_microbiome[, meta, drop=F],
                               paste('analysis/Tweedieverse_mother_microbiome', meta, sep ="_"),
                               base_model = "CPLM",
                               prev_threshold = 0.01,
                               max_significance = .1,
                               standardize = T,
                               
                               #fixed_effects = c("Visit_type", "Race"),
                               plot_scatter =  T,  #"Has_subj_ever_had_a_yeast_infection"),
                               plot_heatmap = T
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
}


