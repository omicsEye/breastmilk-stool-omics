

setwd('/Users/rah/Library/CloudStorage/Box-Box/GW Genomics Core Projects/Projects Currently being Analyzed/Projects for Pay/INOVA_Stool_and_Breastmilk_0074 (Ali)/Analysis/Inova QC')
qc_data <- read.delim(
  'INOVA_OVERALL_RESULT2023-01-29.tsv', #
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)

wilcox.test(qc_data[qc_data$Source=='m', "Viral"] , qc_data[qc_data$Source=='f', "Viral"], paired=F)
wilcox.test(qc_data[qc_data$Source=='m', "Bacterial"] , qc_data[qc_data$Source=='f', "Bacterial"], paired=F)
wilcox.test(qc_data[qc_data$Source=='m', "Fungi"] , qc_data[qc_data$Source=='f', "Fungi"], paired=F)
wilcox.test(qc_data[qc_data$Source=='m', "Archaea"] , qc_data[qc_data$Source=='f', "Archaea"], paired=F)
