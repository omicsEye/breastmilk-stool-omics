library(omicsArt)
library(ggplot2)
library(ggrepel)
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")


path_to_file <- "data/INOV-0102-19VW+ CDT_21AUG2020_Cleanded_10_15_2020.xlsx"
output_path <- paste('analysis/meatbolite', sep = '')
dir.create(file.path(output_path), showWarnings = FALSE)

##############################
#AO and AV label missing as well as data
loaded_data <- load_data(input=path_to_file, type = 'all', sheet = 3, ID = 'Metabolite')
data <- loaded_data$data
data <- numeric_dataframe(data)
metadata <- loaded_data$sample_metadata
rownames(metadata) <- metadata$`SUBJECT ID TIME POINT`
rownames(data) <- rownames(metadata)
metabilotes_info <- loaded_data$feature_metadata
stats_table <- read.delim(
  'analysis/Tweedieverse/Infant/Metabolite/Tweedieverse_metabolitesBirth_weight/all_results.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  #row.names = 1
)

stats_table$RT <- metabilotes_info[match(stats_table$feature, metabilotes_info$Metabolite), "RT"]

#stats_table <- stats_table[stats_table$pval < 0.05,]

p <- omicsArt::lollipop_plot(stats_table, x= "RT", y= "coef",  x_label = "Retention time", y_label = "Effect size (Gestational age)") 
p <- p +
  geom_text_repel(
    size = 1,
    force = .5,
    fontface =  "italic")#,
#box.padding = unit(.01, "lines"),
#point.padding = unit(0.01, "lines")
#)
ggsave(filename = 'manuscript/figures/fig#_metabolites/coef_RT_metabolites_annotated_Gestational_age.png', plot = p,
       height = 2,  width = 7.2 ,units = "in", dpi = 350)
ggsave(filename = 'manuscript/figures/fig#_metabolites/coef_RT_metabolites_annotated_Gestational_age.pdf', plot = p,
       height = 2,  width = 7.2 ,units = "in", dpi = 350)


