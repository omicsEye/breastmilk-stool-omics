#library(vegan)
library (corrplot)
library(Hmisc)
library(omicsPattern)
library(psych)

setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/btest_results/infant_metabolite_mother_microbiome_btest_pair/")

X <- read.table('pvalues_table.txt', header = TRUE,
                row.names = 1,   sep = "\t", fill = FALSE, comment.char = "" , check.names = FALSE)
pheatmap::pheatmap(X)
method = 'spearman'
x_corr <- rcorr(as.matrix(t(X)), type=method)
D <- as.matrix(vegdist(t(X), method="bray")) #as.matrix(cor2dist(x_corr$r))
#col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE","#4393C3", "#2166AC", "#053061"))(200))
dev.off()
pdf(height=5, width=6, file = 'corr_X.pdf')
cor_plot <- corrplot(x_corr$r, type = "lower", order = "hclust",
         tl.col = "black", tl.srt = 60,
         tl.cex = .5,
         #col = col,
         addrect = 5,
         p.mat = x_corr$P,
         pch.cex = .7,
         sig.level = 0.01,
         insig = "label_sig"
         ) #diag = FALSE
dev.off()
p_plot <- pheatmap::pheatmap(x_corr$r)
#p_plot2 <- p_plot+omicsPattern::theme_nature()
#pheatmap::pheatmap(X)
ggsave(filename='manuscript/Figures/corr_heat_proteomics.pdf', plot=p_plot, width = 450, height = 400, units = "mm", dpi = 300) #width = 183, height = 126.8,
dev.off()


##### metabolites correlation plot
number_of_sig_to_keep <- 15
all_assoc <- NULL
for (time_point_base in c('Week 1', 'Month 3')){
  maaslin_wm_time_side <- read.table(
    paste('analysis/metabolomics/P/', time_point_base,'/', 'significant_results.tsv', sep = ""),
    sep = '\t',
    header = TRUE,
    fill = FALSE,
    comment.char = "" ,
    check.names = FALSE,
    #row.names = 1
  )
  number_of_sig_to_keep <- min(number_of_sig_to_keep, dim(maaslin_wm_time_side)[1])
  all_assoc <- rbind(all_assoc, maaslin_wm_time_side[1:number_of_sig_to_keep ,])
}

#}

length(unique(all_assoc$feature))
assoc <- unique(all_assoc$feature)[1:number_of_sig_to_keep]
# X is the dataframe metabolites (rows) and samples (columns)
source('~/Documents/omicsEye/omicsPattern/R/utils.R')
loaded_data <- load_data(input = '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/RYGB/data/Metabolites_netome_format.xlsx',
                                      type = 'all', sheet = 1, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
features_info <- loaded_data$feature_metadata
metadata <- loaded_data$sample_metadata

X2 <- data [, colnames(data) %in% all_assoc$feature]
method = 'spearman'
#cor.test(X$Albumin.Base, X$Age_Est.FU3)
x_corr <- rcorr(as.matrix(X2), type=method)
write.table(x_corr$r,
            "data/metabolite_corr.txt",
            sep = "\t", eol = "\n", col.names = NA, row.names = T)
dev.off()
pdf(height = 5.5, width = 6, file = 'manuscript/Figures/corr_metabolites.pdf')
cor_plot <- corrplot(x_corr$r, type = "upper", order = "hclust",
                     tl.col = "black", tl.srt = 60,
                     tl.cex = .5,
                     col = col,
                     addrect = 5,
                     p.mat = x_corr$P,
                     pch.cex = .7,
                     sig.level = 0.01,
                     insig = "label_sig"
) #diag = FALSE
dev.off()
p_plot <- pheatmap::pheatmap(x_corr$r)
#p_plot2 <- p_plot+omicsPattern::theme_nature()
#pheatmap::pheatmap(X)
ggsave(filename='manuscript/Figures/corr_heat_metabolites.pdf', plot=p_plot, width = 220, height = 200, units = "mm", dpi = 300) #width = 183, height = 126.8,
dev.off()

common_sample <- intersect(row.names(X), row.names(X2))
X <- X[common_sample,]
X2 <- X2[common_sample,]
x_y_corr <- rcorr(as.matrix(X2),as.matrix(X), type=method)
x_y_corr_sub <- NULL
x_y_corr_sub$r <- x_y_corr$r[colnames(X2), colnames(X)]
x_y_corr_sub$P <- x_y_corr$P[colnames(X2), colnames(X)]
x_y_corr_sub$n <- x_y_corr$n[colnames(X2), colnames(X)]
write.table(x_y_corr_sub$r,
            "data/x_y_corr.txt",
            sep = "\t", eol = "\n", col.names = NA, row.names = T)
dev.off()
pdf(height = 5, width = 8, file = 'manuscript/Figures/corr_m2.pdf')
cor_plot <- corrplot(x_y_corr_sub$r,
                     tl.col = "black", tl.srt = 45,
                     tl.cex = .65,
                     col = col,
                     addrect = 5,
                     p.mat = x_y_corr_sub$P,
                     pch.cex = .9,
                     sig.level = 0.01,
                     insig = "label_sig"
) #diag = FALSE

dev.off()


max_value <- ceiling(max(x_y_corr_sub$r))
min_value <- ceiling(min(x_y_corr_sub$r))
range_value <- max(c(abs(max_value),abs(min_value)))
breaks <- seq(range_value, range_value, by = 0.1)
color = colorRampPalette(c("darkblue", "grey90", "darkred"))
p <- pheatmap::pheatmap(
  x_y_corr_sub$r,
  cellwidth = 5,
  cellheight = 5,
  # changed to 3
  #main = title,
  fontsize = 6,
  kmeans_k = NA,
  border = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  legend = TRUE,
  border_color = 'grey93',
  color = color(range_value*100),
  #breaks = breaks,
  treeheight_row = 0,
  treeheight_col = 0,
  display_numbers = matrix(ifelse(
    x_y_corr_sub$P< 0.05, "*",  ""), nrow(x_y_corr_sub$r))
)
ggsave(filename='manuscript/Figures/corr_heat_m2.pdf', plot=p, width = 210, height = 120, units = "mm", dpi = 300) #width = 183, height = 126.8,

dev.off()


