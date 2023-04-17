library(ggplot2)
library(ggridges)
library(cowplot)
## read association
box_association <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Tweedieverse/Cross/Breastmilk/Tweedieverse_mother_microbiome_infant_Breast_milk_collected/figures/Breast_milk_collected_gg_associations.RDS")
## do plots

paird_microbiome <- read.delim(
  '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/btest_results/infant_microbiome_mother_microbiome_btest_pair/X_Y.tsv', #'data/pathab.tsv', #
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)



paird_microbiome <- subset(paird_microbiome, paird_microbiome$Feature_1==paird_microbiome$Feature_2)
paird_microbiome$gwfill[paird_microbiome$Correlation > 0.0 & paird_microbiome$pval < 0.05 ] <- "#002654"
paird_microbiome$gwfill[paird_microbiome$Correlation < 0.0  & paird_microbiome$pval < 0.05] <- "#E5D19D"
paird_microbiome$gwfill[is.na(paird_microbiome$gwfill)] <- "gray"
paird_microbiome <- paird_microbiome[order(paird_microbiome$Correlation),]
paird_microbiome <- within(paird_microbiome,
                           Feature_1 <- factor(Feature_1,
                                      levels=paird_microbiome$Feature_1))
corr_bar_plot <- ggplot2::ggplot(data=paird_microbiome, aes(x= Feature_1, y=Correlation )) +
  ggplot2::geom_bar(stat="identity", fill = paird_microbiome$gwfill, alpha = 0.5, size=0.1) +
  xlab("Common species in infant stool and breast milk")+ylab("Spearman Coef")+
  omicsArt::theme_omicsEye()+
  theme(axis.title = element_text(face="bold"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))

corr_bar_plot

fig2 <- ggdraw() +
  draw_plot(speices_heat_plot$gtable,
            x = 0.02, y = .6, width = .99, height = .4) +
  draw_plot(corr_bar_plot,
            x = 0.02, y = 0, width = .5, height = .6) +
  draw_plot(box_association[[1]] + scale_y_log10() + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 7, face ="bold"),
    axis.text.y = element_text(size = 5)), x = .54, y = 0.01, width = .23, height = .5) +
  draw_plot(box_association[[2]] + scale_y_log10() + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 5),
    axis.title.y = element_text(size = 7, face ="bold"),
    axis.text.y = element_text(size = 5)), x = .77, y = 0, width = .23, height = .5)+
  draw_plot_label((label = c("a", "b", "c", "d")),
                  size = 9,x = c(0, 0, .53, .76), y = c(1, 0.6, 0.6, .6))

# fig2 <- ggdraw() +
#   draw_plot(speices_heat_plot$gtable,
#             x = 0.02, y = .6, width = .99, height = .4) +
#   draw_plot(corr_bar_plot,
#             x = 0.02, y = 0, width = 1, height = .6) +
#     draw_plot_label((label = c("a", "b")),
#                   size = 9,x = c(0, 0), y = c(1, 0.6))


fig2

setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")

ggsave(filename = 'manuscript/figures/fig_overviwe_microbiome/fig2_v3.pdf', plot=fig2, width = 4.8, height = 3, units = "in", dpi = 350)
ggsave(filename = 'manuscript/figures/fig_overviwe_microbiome/fig2_v3.png', plot=fig2, width = 4.8, height = 3, units = "in", dpi = 350)
