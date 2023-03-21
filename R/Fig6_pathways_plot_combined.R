library(ggplot2)
library(cowplot)


setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")

BMI_box_association <- readRDS("analysis/deepath/deepath_output_BMI/figures/gg_enrichment_rank.RDS")
Donor_Maternal_box_association <- readRDS("analysis/deepath/deepath_output_Breast_milk_collected/figures/gg_enrichment_rank.RDS")

## do plots
a_theme <- theme(
  plot.title = element_text(size =7, face = 'bold'),
  axis.title.x = element_text(size = 6),
  axis.text.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  axis.text.y = element_text(size = 5))

fig6 <-   plot_grid(BMI_box_association$`GO:0006099|BP|03|tricarboxylic acid cycle` +a_theme,
                    BMI_box_association$`GO:0009231|BP|04|riboflavin biosynthetic process`,
                    BMI_box_association$`GO:0019323|BP|06|pentose catabolic process`,
                    BMI_box_association$`GO:0046080|BP|07|dUTP metabolic process`,
                    Donor_Maternal_box_association$`GO:0044718|BP|04|siderophore transmembrane transport`+a_theme,
                    Donor_Maternal_box_association$`GO:0042732|BP|06|D-xylose metabolic process`,
                    Donor_Maternal_box_association$`GO:0019509|BP|07|L-methionine biosynthetic process from methylthioadenosine`+a_theme,
                    Donor_Maternal_box_association$`GO:2001059|BP|05|D-tagatose 6-phosphate catabolic process`,
          #rel_widths = c(1, 1, 1,1),
          labels = c('a, BMI',"","","","b, Maternal vs. donor"),
          label_size = 7, ncol = 4, nrow = 2)

ggsave(filename = 'manuscript/figures/figure_6_pathway/fig6_pathways_rank.pdf', plot=fig6, width = 350, height = 150, units = "mm", dpi = 350)
ggsave(filename = 'manuscript/figures/figure_6_pathway/fig6_pathways_rank.png', plot=fig6, width = 350, height = 150, units = "mm", dpi = 350)

