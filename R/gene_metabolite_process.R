
#require(gdata)
library(gdata)
library(pheatmap)
library(vegan)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(Hmisc)
library(Maaslin2)
library(Tweedieverse)
library(omicsArt)
library(stringr)
library(readxl)

setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")

output_path <- paste('analysis/genes', sep = '')
dir.create(file.path(output_path), showWarnings = FALSE)

####Read metabolites ##########################
#AO and AV label missing as well as data
path_to_file <- "data/INOV-0102-19VW+ CDT_21AUG2020_Cleanded_10_15_2020.xlsx"
loaded_data <- load_data(input=path_to_file, type = 'all', sheet = 2, ID = 'Metabolite')
metabolites <- loaded_data$data
metabolites <- numeric_dataframe(metabolites)
metadata <- loaded_data$sample_metadata
rownames(metadata) <- metadata$`SUBJECT ID TIME POINT`
rownames(metabolites) <-  loaded_data$sample_metadata$`SUBJECT ID TIME POINT`

metadata_file <- read_xlsx("data/EXTERNAL_Breast_Milk_Data_Final_4.16.19.xlsx", skip =2)
metadata_file <- metadata_file[-1, ]

 ###RUN Metadata.R

##### metaphlan species profiles ###
genes <- read.delim(
  'data/gene_list.tsv', #'data/pathab.tsv', #
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
#clean sample anmes and species names
cols_temp <- str_split(colnames(genes), "-")
#paste(sapply(cols_temp, "[[", 2), sapply(cols_temp, "[[", 3), sep = "_")

colnames(genes) <-  paste(sapply(cols_temp, "[[", 2), sapply(cols_temp, "[[", 3), sep = "_") #sapply(cols_temp, "[[", 2)#


# id_mapper <- read_xlsx("/Users/rah/Library/CloudStorage/Box-Box/GW Genomics Core Projects/Projects Currently being Analyzed/INOVA_Stool_and_Breastmilk_0074/data/metadata/external_metabolon_gw_matchup.xlsx")
  #"~/Box/GW Genomics Core Projects/Projects Currently being Analyzed/INOVA_Stool_and_Breastmilk_0074/
  #data/metadata/external_metabolon_gw_matchup.xlsx")
#mapped_sample_ids <- as.list.data.frame( id_mapper[match(sapply(cols_temp, "[[", 2), id_mapper$`GW Sample ID`), "MGX_MBX_Map"] )
#colnames(microbiome) <- mapped_sample_ids$MGX_MBX_Map

mom_or_infant <- genes
mom_or_infant[1,] <- sapply(cols_temp, "[[", 3)
genes <- genes[, !is.na(colnames(genes))]

mom_or_infant <- mom_or_infant[!is.na(colnames(mom_or_infant))]
#colnames(mom_or_infant) <- gsub("\\.1","", colnames(mom_or_infant))
#colnames(microbiome) <- colnames(mom_or_infant)
#mom_or_infant <- as.data.frame(t(mom_or_infant))
mom_or_infant <- mom_or_infant[1,, drop=F]
row.names(mom_or_infant)[1] <- "Sample"
mom_or_infant[mom_or_infant=="f"] <- "Stool"
mom_or_infant[mom_or_infant=="m"] <- "Breast milk"
mom_or_infant <- as.data.frame(t(mom_or_infant))
###visualization
mom_or_infant <- mom_or_infant[mom_or_infant$Sample %in% c("Stool", "Breast milk"),, drop=F]
# Antibiotic resistant genes
AR_genes <- genes[AR_genes_study,]
AR_genes <- as.data.frame(t(AR_genes))
AR_genes <- AR_genes[rownames(mom_or_infant),]
fakepcl <- list(meta=mom_or_infant, x=as.matrix(AR_genes),
                ns=dim(AR_genes)[2], nf=dim(AR_genes)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= c("Sample"), show_colnames = F, show_rownames = T,treeheight_row = 5, treeheight_col= 5 )
ggsave(filename='analysis/AR_genes.png', plot=heat_plot, width = 18, height = 3.1, units = "in", dpi = 300)
ggsave(filename='analysis/AR_genes.pdf', plot=heat_plot, width = 18, height = 3.1, units = "in", dpi = 300)

# Antibiotic resistant genes
R_genes <- genes[resistant_genes_study,]
R_genes <- as.data.frame(t(R_genes))
R_genes <- R_genes[rownames(mom_or_infant),]
fakepcl <- list(meta=mom_or_infant, x=as.matrix(R_genes),
                ns=dim(AR_genes)[2], nf=dim(AR_genes)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 2, meta= c("Sample"), show_colnames = T, show_rownames = T,treeheight_row = 5, treeheight_col= 5 )
ggsave(filename='analysis/R_genes.png', plot=heat_plot, width = 18, height = 9, units = "in", dpi = 300)
ggsave(filename='analysis/R_genes.pdf', plot=heat_plot, width = 18, height = 9, units = "in", dpi = 300)


###### data pairing #########
# infant microbial species and samples information
infant_genes_samples <- intersect(colnames(genes[,colnames(mom_or_infant[,mom_or_infant[1,]=="Infant",drop=F])]), rownames(infant_metadata))
infant_metadata_genes <- infant_metadata[infant_genes_samples,]
infant_genes <- genes[,infant_genes_samples]
rownames(infant_metadata) <- gsub("_f","", rownames(infant_metadata))
rownames(infant_metadata) <- paste0(infant_metadata$External_ID, infant_metadata$Visit, spe="")
mother_genes_samples <- intersect(colnames(genes[,colnames(mom_or_infant[,mom_or_infant[1,]=="Mother",drop=F])]), rownames(mother_metadata))
mother_metadata_genes <- mother_metadata[mother_genes_samples,]
mother_genes <- genes[,mother_genes_samples]




infant_genes <- infant_genes[!grepl("\\|",rownames(infant_genes)),]
colnames(infant_genes) <- gsub("_f","", colnames(infant_genes))

# mother microbial species and samples information
mother_genes <- mother_genes[!grepl("\\|",rownames(mother_genes)),]
rownames(mother_metadata) <- gsub("_m","", rownames(mother_metadata))
colnames(mother_genes) <- gsub("_m","", colnames(mother_genes))

write.table(infant_genes, 'data/infant_genes.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(mother_genes, 'data/mother_genes.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

# pair infant's genes with mother genes
infant_genes_mother_genes_samples <- intersect(colnames(infant_genes), colnames(mother_genes))
infant_genes_mother_genes <- infant_genes[,infant_genes_mother_genes_samples]
mother_genes_infant_genes <- mother_genes[,infant_genes_mother_genes_samples]
write.table(infant_genes_mother_genes, 'data/infant_genes_mother_genes.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(mother_genes_infant_genes, 'data/mother_genes_infant_genes.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
mother_metadata$Breast_milk_collected <- infant_metadata[rownames(mother_metadata), "Breast_milk_collected"]

Maaslin2(mother_genes,
         mother_metadata[, "Breast_milk_collected", drop=F],
         paste('analysis/Maaslin2/Breastmilk/Maaslin2_mother_genes_',"Breast_milk_collected", sep =""),
         min_abundance = 0.0,
         min_prevalence = 0.1,
         max_significance = .1,
         standardize = T,
         transform = 'Log',
         analysis_method = 'LM',
         #normalization = 'AST',
         #transform = 'NONE',
         #fixed_effects = c("Visit", "Race"),
         heatmap_first_n = 50,  #"Has_subj_ever_had_a_yeast_infection"),
)
Maaslin2(mother_genes,
         mother_metadata[, "BMI", drop=F],
         paste('analysis/Maaslin2/Breastmilk/Maaslin2_mother_genes_',"BMI", sep =""),
         min_abundance = 0.0,
         min_prevalence = 0.1,
         max_significance = .1,
         standardize = T,
         transform = 'Log',
         analysis_method = 'LM',
         #normalization = 'AST',
         #transform = 'NONE',
         #fixed_effects = c("Visit", "Race"),
         heatmap_first_n = 50,  #"Has_subj_ever_had_a_yeast_infection"),
)

