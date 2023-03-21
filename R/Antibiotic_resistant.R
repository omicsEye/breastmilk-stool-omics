
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")


AR_genes <- read.delim(
  "data/AR_genes/uniref-yourlist_M20220115F248CABF64506F29A91F8037F07B67D13CBFFCY.tab",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)
all_R_genes <- read.delim(
  "data/AR_genes/resistance_genes_UniPortKB_2_UniRef90.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = NA
)

infant_genes <- read.delim( 'data/infant_genes_mother_genes.txt',
                           sep = '\t',
                           header = T,
                           fill = F,
                           comment.char = "" ,
                           check.names = F,
                           row.names = 1
)
infant_genes <- as.data.frame(t(infant_genes))

mother_genes <- read.delim( 'data/mother_genes_infant_genes.txt',
                           sep = '\t',
                           header = T,
                           fill = F,
                           comment.char = "" ,
                           check.names = F,
                           row.names = 1
)
mother_genes <- as.data.frame(t(mother_genes))


AR_genes_study_f <- intersect(AR_genes$`Cluster ID`, colnames(infant_genes))
# [1] "UniRef90_P31121" "UniRef90_Q8ZPP2" "UniRef90_P52477" "UniRef90_Q9I3Y3" "UniRef90_P56976" "UniRef90_O31855"
# [7] "UniRef90_Q9HTR0" "UniRef90_Q8DPQ6" "UniRef90_P25396" "UniRef90_P52067" "UniRef90_G3XD21" "UniRef90_A7ZKF6"
# [13] "UniRef90_Q7AFA1"

AR_genes_study_m <-  intersect(AR_genes$`Cluster ID`, colnames(mother_genes))
# [1] "UniRef90_P31121" "UniRef90_Q8ZPP2" "UniRef90_P52477" "UniRef90_Q9I3Y3" "UniRef90_P56976" "UniRef90_O31855"
# [7] "UniRef90_Q9HTR0" "UniRef90_Q8DPQ6" "UniRef90_P25396" "UniRef90_P52067" "UniRef90_G3XD21" "UniRef90_A7ZKF6"
# [13] "UniRef90_Q7AFA1"

AR_genes_study <- union(AR_genes_study_f,AR_genes_study_m)

resistant_genes_study_f <- intersect(all_R_genes$`Cluster ID`, colnames(infant_genes))
# [1] "UniRef90_Q8XFG0"     "UniRef90_P52002"     "UniRef90_Q5KJE4"     "UniRef90_P0CW83"     "UniRef90_O34777"    
# [6] "UniRef90_O35018"     "UniRef90_Q8Y9K8"     "UniRef90_Q8DPQ6"     "UniRef90_A0A0P0F8B5" "UniRef90_Q8ZPP2"    
# [11] "UniRef90_Q8ZNQ0"     "UniRef90_O06967"     "UniRef90_O07550"     "UniRef90_O07549"     "UniRef90_P52003"    
# [16] "UniRef90_O31855"     "UniRef90_P0AAA2"     "UniRef90_Q8X7J5"     "UniRef90_P52067"     "UniRef90_P31121"    
# [21] "UniRef90_P15177"     "UniRef90_Q71YH0"     "UniRef90_Q4K7U5"     "UniRef90_A0A069QB48" "UniRef90_P76352"    
# [26] "UniRef90_P0AET6"     "UniRef90_P75685"     "UniRef90_P56976"     "UniRef90_A7ZKF6"     "UniRef90_P0AB42"    
# [31] "UniRef90_P32157"     "UniRef90_P76573"     "UniRef90_P37621"     "UniRef90_Q7AFA1"     "UniRef90_P75863"    
# [36] "UniRef90_P75821"     "UniRef90_P77549"     "UniRef90_P25396"     "UniRef90_Q8ZKY1"     "UniRef90_P96712"    
# [41] "UniRef90_G3XD21"     "UniRef90_Q9K015"     "UniRef90_P52477"     "UniRef90_Q8DR23"     "UniRef90_Q9HTR0"    
# [46] "UniRef90_Q9I3Y3"     "UniRef90_A0A0C6FA12" "UniRef90_A0A0C6ETG9"

resistant_genes_study_m <- intersect(all_R_genes$`Cluster ID`, colnames(mother_genes))
# [1] "UniRef90_Q8XFG0"     "UniRef90_P52002"     "UniRef90_Q5KJE4"     "UniRef90_P0CW83"     "UniRef90_O34777"    
# [6] "UniRef90_O35018"     "UniRef90_Q8Y9K8"     "UniRef90_Q8DPQ6"     "UniRef90_A0A0P0F8B5" "UniRef90_Q8ZPP2"    
# [11] "UniRef90_Q8ZNQ0"     "UniRef90_O06967"     "UniRef90_O07550"     "UniRef90_O07549"     "UniRef90_P52003"    
# [16] "UniRef90_O31855"     "UniRef90_P0AAA2"     "UniRef90_Q8X7J5"     "UniRef90_P52067"     "UniRef90_P31121"    
# [21] "UniRef90_P15177"     "UniRef90_Q71YH0"     "UniRef90_Q4K7U5"     "UniRef90_A0A069QB48" "UniRef90_P76352"    
# [26] "UniRef90_P0AET6"     "UniRef90_P75685"     "UniRef90_P56976"     "UniRef90_A7ZKF6"     "UniRef90_P0AB42"    
# [31] "UniRef90_P32157"     "UniRef90_P76573"     "UniRef90_P37621"     "UniRef90_Q7AFA1"     "UniRef90_P75863"    
# [36] "UniRef90_P75821"     "UniRef90_P77549"     "UniRef90_P25396"     "UniRef90_Q8ZKY1"     "UniRef90_P96712"    
# [41] "UniRef90_G3XD21"     "UniRef90_Q9K015"     "UniRef90_P52477"     "UniRef90_Q8DR23"     "UniRef90_Q9HTR0"    
# [46] "UniRef90_Q9I3Y3"     "UniRef90_A0A0C6FA12" "UniRef90_A0A0C6ETG9"

resistant_genes_study <- union(resistant_genes_study_f,resistant_genes_study_m)

mapper_file <- read.delim('/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/omePath/GO_UNIREF90_MAP.tsv', 
                          sep = '\t',
                          header = T,
                          fill = F,
                          comment.char = "" ,
                          check.names = F,
                          #row.names = NA
)
intersect(AR_genes$`Cluster ID`, mapper_file$Feature_ID)


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

omePath_result <- omePath::omePath(
  input_data = score_data_filtered,
  output = "analysis/omePath/omePath_output_BMI",
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

