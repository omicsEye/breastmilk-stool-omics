


source('~/Documents/omicsEye/omicsArt/R/pcl_utils.R')

species <- read.delim(
  '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/data/merged_species.tsv', #'/Users/rah/Box/HMH-CDI COVID Data/Analysis/HackensackMetaData_v2.txt',
  sep = '\t',
  header = TRUE,
  fill = T,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

species <- species[grepl("\\|s__[^\\|]+$", rownames(species)),]
rownames(species) <- gsub("_", " ", gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "", rownames(species)))


## species pcoa plot
fakepcl <- list(meta=NA, x=as.matrix(t(species)),
                ns=dim(species)[1], nf=dim(species)[2])
heat_plot <- pcl.heatmap(fakepcl, sqrtspace = T, gamma = 3)

ggsave(filename='/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/species.pdf', plot=heat_plot, 
       width = 5, height = 4, units = "in", dpi = 300)
