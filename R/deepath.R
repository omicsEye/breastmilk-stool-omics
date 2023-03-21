
# details at https://github.com/omicsEye/deepath

# Install devtools
install.packages('devtools')
library(devtools)

# Install the Bioconuctor dependencies
install.packages('BiocManager'); library('BiocManager');


# Install the CRAN dependencies:
install.packages(c('future', 'downloader', 'reader', 'backports', 'gsEasy','pscl','pbapply','car','nlme','dplyr','vegan','chemometrics','ggplot2','pheatmap','cplm','hash','logging','data.table'), repos='http://cran.r-project.org')

# Install deepath library from GitHub
devtools::install_github('omicsEye/deepath', force = TRUE)

# load deepath library
library(deepath)

setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")

# read effect size values from Tweediverse output 
Effect_size_result <- read.delim(
  "analysis/Maaslin2/Breastmilk/Maaslin2_mother_genes_Breast_milk_collected/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
score_data_filtered <- Effect_size_result[Effect_size_result$metadata=="Breast_milk_collected" & Effect_size_result$value=="Maternal" ,]
row.names(score_data_filtered)=score_data_filtered$feature


pathway_col = "Pathway"
feature_col = "Feature_ID"

mapper_file <- read.delim('/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/deepath/GO_UNIREF90_MAP.tsv', 
                          sep = '\t',
                          header = T,
                          fill = F,
                          comment.char = "" ,
                          check.names = F,
                          #row.names = NA
)

deepath_result <- deepath::deepath(
  input_data = score_data_filtered,
  output = "analysis/deepath_output_Breast_milk_collected",
  score_col = 'coef',
  pval_threshold = 0.05,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = mapper_file,
  method = "ks",
  min_member = 2,
  pathway_col = pathway_col,
  feature_col = feature_col)






Effect_size_result_BMI <- read.delim(
  "analysis/Maaslin2/Breastmilk/Maaslin2_mother_genes_BMI/all_results.tsv",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
score_data_filtered <- Effect_size_result_BMI[Effect_size_result_BMI$metadata=="BMI"  ,]
row.names(score_data_filtered)=score_data_filtered$feature

deepath_result <- deepath::deepath(
  input_data = score_data_filtered,
  output = "analysis/deepath/deepath_output_BMI",
  score_col = 'coef',
  pval_threshold = 0.05,
  fdr_threshold = NA,
  Pathway.Subject = NA,#'Metabolic',
  do_plot = TRUE,
  mapper_file = mapper_file,
  method = "ks",
  min_member = 2,
  pathway_col = pathway_col,
  feature_col = feature_col)
