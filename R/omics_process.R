
#require(gdata)
library(gdata)
library(pheatmap)
library(vegan)
#library(corrplot)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
#library(Hmisc)
library(omicsArt)
library(stringr)
library(readxl)
library(massSight)
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")

#path_to_file <- "data/INOV-0101-19INTR Client Data Table (Frozen)-Raw Data_Cleaned_Ali.xlsx"
output_path <- paste('analysis/metabolite', sep = '')
dir.create(file.path(output_path), showWarnings = FALSE)

do_write = F

####Read metabolites ##########################
#AO and AV label missing as well as data
path_to_file <- "data/INOV-0102-19VW+ CDT_21AUG2020_Cleanded_10_15_2020.xlsx"
loaded_data <- load_data(input=path_to_file, type = 'all', sheet = 2, id = 'Metabolite')
metabolites <- loaded_data$data
metabolites <- numeric_dataframe(metabolites)
metadata <- loaded_data$sample_metadata
rownames(metadata) <- metadata$`SUBJECT ID TIME POINT`
rownames(metabolites) <-  loaded_data$sample_metadata$`SUBJECT ID TIME POINT`

metadata_file <- read_xlsx("data/EXTERNAL_Breast_Milk_Data_Final_4.16.19.xlsx", skip =2)
metadata_file <- metadata_file[-1, ]


temp_filtered_metadata <- mother_metadata[, apply(mother_metadata, 2, entropy) > 0.5]
Excluded_metadata <- setdiff(colnames(mother_metadata), colnames(temp_filtered_metadata))
#logging::loginfo("Excluded metadata with entropy less or equal to 0.75: %s", Excluded_metadata)
mother_metadata <- temp_filtered_metadata

#Maritial_status, "External_ID","Maternal_smoking",

##### metaphlan species profiles ###
microbiome <- read.delim(
  'data/species_abundance_table.txt', #'data/pathab.tsv', #
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
#clean sample anmes and species names
cols_temp <- str_split(colnames(microbiome), "-")
#paste(sapply(cols_temp, "[[", 2), sapply(cols_temp, "[[", 3), sep = "_")

colnames(microbiome) <-  paste(sapply(cols_temp, "[[", 2), sapply(cols_temp, "[[", 3), sep = "_") #sapply(cols_temp, "[[", 2)#

# pair mother's microbiome sample with infant's microbiome and metabolites
# pair mothers' metagenomics samples' ids to infants' samples (metabolomics and metagenomics) ids

id_mapper <- read_xlsx("/Users/rah/Library/CloudStorage/Box-Box/GW Genomics Core Projects/Projects Currently being Analyzed/Projects for Pay/INOVA_Stool_and_Breastmilk_0074 (Ali)/data/metadata/external_metabolon_gw_matchup.xlsx")
  #"~/Box/GW Genomics Core Projects/Projects Currently being Analyzed/INOVA_Stool_and_Breastmilk_0074/
  #data/metadata/external_metabolon_gw_matchup.xlsx")
#mapped_sample_ids <- as.list.data.frame( id_mapper[match(sapply(cols_temp, "[[", 2), id_mapper$`GW Sample ID`), "MGX_MBX_Map"] )
#colnames(microbiome) <- mapped_sample_ids$MGX_MBX_Map

mom_or_infant <- microbiome
mom_or_infant[1,] <- sapply(cols_temp, "[[", 3)
microbiome <- microbiome[, !is.na(colnames(microbiome))]

mom_or_infant <- mom_or_infant[!is.na(colnames(mom_or_infant))]
#colnames(mom_or_infant) <- gsub("\\.1","", colnames(mom_or_infant))
#colnames(microbiome) <- colnames(mom_or_infant)
#mom_or_infant <- as.data.frame(t(mom_or_infant))
mom_or_infant <- mom_or_infant[1,, drop=F]
row.names(mom_or_infant)[1] <- "Sample"
mom_or_infant[mom_or_infant=="f"] <- "Infant stool"
mom_or_infant[mom_or_infant=="m"] <- "Breast milk"
#mom_or_infant[!mom_or_infant%in% c("Infant stool", "Breast milk")] <- "Control"


species <- microbiome
species <- species[grepl("\\|s__[^\\|]+$", rownames(species)),]
rownames(species) <- gsub("_", " ", gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "", rownames(species)))
colnames(species) <- paste(sapply(cols_temp, "[[", 2), sapply(cols_temp, "[[", 3), sep = "_")
genera <- microbiome
genera <- genera[grepl("\\|g__[^\\|]+$", rownames(genera)),]
rownames(genera) <- gsub("_", " ", gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "", rownames(genera)))
colnames(genera) <- paste(sapply(cols_temp, "[[", 2), sapply(cols_temp, "[[", 3), sep = "_")



#colnames(mom_or_infant) <- paste(sapply(cols_temp, "[[", 2), sapply(cols_temp, "[[", 3), sep = "_")
# transpose to have row as samples and columns as microbial species

#mom_or_infant = as.data.frame(t(mom_or_infant))

#species <- species[!rownames(species) %in% c("negative", "Zymo"), ]
#species <- species[, !grepl("negative", colnames(species))]
#species <- species[, !grepl("Zymo", colnames(species))]
#genera <- genera[!rownames(genera) %in% c("negative", "Zymo"), ]
#genera <- genera[,!grepl("negative", colnames(genera))]
#genera <- genera[,!grepl("Zymo", colnames(genera))]

#mom_or_infant <- mom_or_infant[!rownames(mom_or_infant) %in% c("negative", "Zymo"),, drop=F]
#mom_or_infant <- mom_or_infant[,!grepl("negative", colnames(mom_or_infant)), drop=F]
#mom_or_infant <- mom_or_infant[,!grepl("Zymo", colnames(mom_or_infant)), drop=F]


mom_or_infant <-  as.data.frame(t(mom_or_infant))
mom_or_infant$Sample[!mom_or_infant$Sample %in% c("Infant stool", "Breast milk")] <- "Control"
species <- as.data.frame(t(species))
rownames(species) <- gsub("X", "", rownames(species))
genera <- as.data.frame(t(genera))
rownames(genera) <- gsub("X", "", rownames(genera))

##### explanatory plots ###

ad <- adonis2(species ~ . , data= mom_or_infant , permutations=4999, mehod='euclidean')
ad

############## Filter species for
breastmilk_species <- species[grepl("_m", row.names(species)),]
stool_species <- species[grepl("_f", row.names(species)),]

rownames(infant_metadata) <- paste0(rownames(infant_metadata),"_f")

# infant microbial species and samples information
infant_microbiome_samples <- intersect(colnames(microbiome[,rownames(mom_or_infant[mom_or_infant$Sample=="Infant stool",,drop=F])]), rownames(infant_metadata))
infant_metadata_microbiome <- infant_metadata[infant_microbiome_samples,]
infant_microbiome <- microbiome[,infant_microbiome_samples]
rownames(infant_metadata) <- gsub("_f","", rownames(infant_metadata))
rownames(infant_metadata) <- paste0(infant_metadata$`External ID`, infant_metadata$`Visit time point`, spe="")
mother_microbiome_samples <- intersect(colnames(microbiome[,rownames(mom_or_infant[mom_or_infant$Sample=="Breast milk",,drop=F])]), rownames(mother_metadata))
mother_metadata_microbiome <- mother_metadata[mother_microbiome_samples,]
mother_microbiome <- microbiome[,mother_microbiome_samples]




infant_microbiome <- infant_microbiome[grepl("\\|s__[^\\|]+$", rownames(infant_microbiome)),]
rownames(infant_microbiome) <- gsub("_", " ", gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "", rownames(infant_microbiome)))
colnames(infant_microbiome) <- gsub("_f","", colnames(infant_microbiome))
# infant metabolites and samples information
infant_metabolites_samples <- intersect(rownames(metabolites), rownames(infant_metadata))
infant_metadata_metabolites = infant_metadata[infant_metabolites_samples,]
infant_metabolite = metabolites[infant_metabolites_samples,]

# mother microbial species and samples information
mother_microbiome <- mother_microbiome[grepl("\\|s__[^\\|]+$", rownames(mother_microbiome)),]
rownames(mother_microbiome) <- gsub("_", " ", gsub(
  "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "", rownames(mother_microbiome)))
rownames(mother_metadata) <- gsub("_m","", rownames(mother_metadata))
colnames(mother_microbiome) <- gsub("_m","", colnames(mother_microbiome))



#read BGCs

BGCs <- read.delim(
  'data/BGC_Abundance_Across_Samples.tsv', #'data/pathab.tsv', #
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
BGCs <- as.data.frame(t(BGCs))
rownames(BGCs) <- gsub("X","", rownames(BGCs))

mom_or_infant_BGC <- BGCs[,1, drop=F]
colnames(mom_or_infant_BGC)[1] <- "Sample"
mom_or_infant_BGC[grepl("_f", rownames(BGCs)),"Sample"] <- "Infant stool"
mom_or_infant_BGC[grepl("_m", rownames(BGCs)), "Sample"] <- "Breast milk"
mom_or_infant_BGC[!mom_or_infant_BGC$Sample %in% c("Infant stool", "Breast milk") , "Sample"] <- "Control"

breastmilk_BGC <- BGCs[grepl("_m", row.names(BGCs)),]
stool_BGC <- BGCs[grepl("_f", row.names(BGCs)),]

if (do_write){

  write.table(mother_microbiome, 'data/breastmilk_microbiome.tsv',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  write.table(infant_microbiome, 'data/infant_microbiome.tsv',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  write.table(infant_metabolite, 'data/infant_metabolite.tsv',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  write.table(breastmilk_BGC, 'data/breastmilk_BGC.tsv',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  write.table(stool_BGC, 'data/stool_BGC.tsv',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  write.table(stool_BGC, 'data/stool_BGC.tsv',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  write.table(BGCs, 'data/BGCs.tsv',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
}


