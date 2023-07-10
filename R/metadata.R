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
######## Infant information#####
do_write = F
infant_metadata <- read.delim(
  #'~/Box/GW Genomics Core Projects/Projects Currently being Analyzed/INOVA_Stool_and_Breastmilk_0074(Ali)/data/metadata/infant_metadata.txt',
  '/Users/rah/Library/CloudStorage/Box-Box/GW Genomics Core Projects/Projects Currently being Analyzed/Projects for Pay/INOVA_Stool_and_Breastmilk_0074 (Ali)/data/metadata/infant_metadata.txt',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
#infant_metadata$`Person ID` <- gsub("a|b|c|d", "", rownames(infant_metadata))

#temp_filtered_metadata <- infant_metadata[, apply(infant_metadata, 2, omicsArt::entropy) > 0.5]
#Excluded_metadata <- setdiff(colnames(infant_metadata), colnames(temp_filtered_metadata))
#logging::loginfo("Excluded metadata with entropy less or equal to 0.75: %s", Excluded_metadata)
#infant_metadata <- temp_filtered_metadata
infant_metadata$`Gestational age` <- infant_metadata$`Gestational age week` * 7 + infant_metadata$`Gestational age day`

# time_based (all vitis) <- c:
#   Visit (Visit time point), Gestational_age_week (Gestational age (weeks))", "Gestational_age_day", "Gestational_age_category (Gestational age category)", "Gestational_by_weight_percentage (Weight percentile by gestational age)", "Gestational_by_weight_category (Weight category by gestational age)",  "Sex (Sex)", "Days_in_NICU", "Breast_milk_collected", "Diet_at_sample (Diet) ", "All_medications_infant", "Recent_follow_up_date"), "Birth_weight (Birth weight)", "Birth_weight_category (Birth weight category)",   "Sample_fortified_with (Fortifer in feed)",  "Age_stool_sample (Age at sample)",  Method_of_feeding (Feeding method)",      "Last_antibiotic_received (Time since last antibiotic)",          "Gastric_acid_suppresion (Gastric acid suppression)",  "Underwent_sepsis (Underwent sepsis rule out)", "Diagnosed_with_Sepsis (Diagnosed with sepsis)"
#
# not_time_based (visit a) <- c "Breast_milk_sample_date", "Age_at_breast_milk_sample", "Stool_sample_date",
# "GI_related_diagnosis",  "GI_surgery_prior", "Metabolic_Condition", "Cystic_Fibrosis",
# "Infection_due_resistant",  "MRSA_swab_result", "VRE_swab_result", "Problems_at_discharge",  "Other_medical_conditions")
#
# For breastmilk:
#
#   m_time_based: (all vitis) <- c("Visit (Visit time point)", "Maternal_education", "Gravidity (Gravidity)", "Parity_after_delivery (Parity after delivery)", "Maternal_smoking (Maternal smoking)", "Alcohol_during_pregnancy (Alcohol during pregnancy)", "Hypertensive_disorder (Hypertensive disorder)", "GBS_positive", "PROM (Premature rupture of membranes)", "Chorioamnionitis (Chorioamnionitis)", "Supplements_during_pregnancy"),  "Age (Maternal Age)",  "Method_of_delivery (Method of Delivery)",              "Antibiotic_use_prenatally (Prenatal antibiotic use)","Antibiotic_during_delivery (Peripartum antibiotic use)"  "BMI( BMI), "BMI_category (BMI category)", "Ethnicity (Ethnicity)", "Race (Race)",
#
# m_not_time_based (visit a) <- c( "Maritial_status", "Weight_1Year_Prior", "Height", "Country_lived", "Household_income", "Breast_milk_collected",
#           )

infant_metadata <- infant_metadata[, ! colnames(infant_metadata) %in% c("BMI category", "Maternal smoking",
                                                                        "Other medical conditions",
                                                                        "Gestational age week",
                                                                        "Gestational age day",
                                                                        "Gestational age category",
                                                                        "Birth weight category",
                                                                        "Weight category by gestational age",
                                                                        "Stool sample date",
                                                                        "Breast milk sample date",
                                                                        "Other medical conditions",
                                                                        "All medications infant",
                                                                        "Age at sample",
                                                                        "Problems at discharge",
                                                                        "Recent follow up date",
                                                                        "VRE swab result",
                                                                        "Days in NICU",
                                                                        "Breast milk sample date",
                                                                        "Age at breast milk sample",
                                                                        "Stool sample date",
                                                                        "GI related diagnosis",
                                                                        "GI surgery prior",
                                                                        "Metabolic Condition",
                                                                        "Cystic Fibrosis",
                                                                        "Infection due resistant",
                                                                        "MRSA swab result", "VRE swab result",
                                                                        "Problems at discharge",  "Other medical conditions",
                                                                        "Diet"

)]

#colnames(infant_metadata) <- gsub("_", " ", colnames(infant_metadata))
infant_metadata2 <- infant_metadata[, ! colnames(infant_metadata) %in% c("External ID")]

if (do_write){
  result_infant <- omicsArt::metadataCorr(infant_metadata2, entropy_threshold = 0.5, p_threshold = 0.01) #[infant_metadata$Visit=='a',]
  ggsave(filename='Manuscript/figures/SFig1/SFig1a.png', plot=result_infant$pval_hetamap, width = 7.2, height = 6, units = "in", dpi = 300)
  write.table( result_infant$P_perm,"Manuscript/figures/SFig1/SFig1a_infant_stool_metadataCorrelation_Pvalue.txt",
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  #
  write.table(infant_metadata, '/Users/rah/Library/CloudStorage/Box-Box/GW Genomics Core Projects/Projects Currently being Analyzed/Projects for Pay/INOVA_Stool_and_Breastmilk_0074 (Ali)/data/metadata/infant_metadata_processed.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
}

######## Mothers information#####

mother_metadata <- read.delim(
  #'~/Box/GW Genomics Core Projects/Projects Currently being Analyzed/INOVA_Stool_and_Breastmilk_0074/data/metadata/mother_metadata.txt',
  '/Users/rah/Library/CloudStorage/Box-Box/GW Genomics Core Projects/Projects Currently being Analyzed/Projects for Pay/INOVA_Stool_and_Breastmilk_0074 (Ali)/data/metadata/mother_metadata.txt',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
#mother_metadata$`Person ID` <- gsub("a|b|c|d", "", rownames(mother_metadata))

mother_metadata$`Gestational age` <- infant_metadata[match(mother_metadata$`External ID`, infant_metadata$`External ID`), "Gestational age"]
mother_metadata$`Maternal vs. donor` <- infant_metadata[match(mother_metadata$`External ID`, infant_metadata$`External ID`), "Maternal vs. donor"]
mother_metadata$Gestational_age_at_birth_week <- infant_metadata[match(mother_metadata$External_ID, infant_metadata$External_ID), "Gestational_age_at_birth_week"]
mother_metadata$Day_of_life_at_time_of_breast_milk_sample <- infant_metadata[match(mother_metadata$External_ID, infant_metadata$External_ID), "Day_of_life_at_time_of_breast_milk_sample"]

#mother_metadata$Day_of_life_at_time_of_breast_milk_sample <- infant_metadata[rownames(mother_metadata), "Day_of_life_at_time_of_breast_milk_sample"]
#temp_filtered_metadata <- infant_metadata[, apply(infant_metadata, 2, omicsArt::entropy) > 0.5]
#temp_filtered_metadata <- mother_metadata[, apply(mother_metadata, 2, omicsArt::entropy) > 0.5]
#Excluded_metadata <- setdiff(colnames(mother_metadata), colnames(temp_filtered_metadata))
#logging::loginfo("Excluded metadata with entropy less or equal to 0.75: %s", Excluded_metadata)
#mother_metadata <- temp_filtered_metadata

mother_metadata <- mother_metadata[, ! colnames(mother_metadata) %in% c("Supplements during pregnancy",
                                                                        "Maternal education",
                                                                        "Maritial status",
                                                                        "Meconium stain"
)]
colnames(mother_metadata) <- gsub("_", " ", colnames(mother_metadata))

#Maritial_status, "External_ID","Maternal_smoking",
mother_metadata2 <- mother_metadata[mother_metadata$`Maternal vs. donor`=="Maternal", !colnames(mother_metadata) %in% c("External ID")]
result_breastmilk <- omicsArt::metadataCorr(mother_metadata2, entropy_threshold = 0.5, p_threshold = 0.01)

if (do_write){
  ggsave(filename='Manuscript/figures/SFig1/SFig1b.png', plot=result_breastmilk$pval_hetamap, width = 7.2, height = 6, units = "in", dpi = 300)
  write.table( result_breastmilk$P_perm,"Manuscript/figures/SFig1/SFig1a_breastmilk_metadataCorrelation_Pvalue.txt",
               sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
  write.table(mother_metadata, '/Users/rah/Library/CloudStorage/Box-Box/GW Genomics Core Projects/Projects Currently being Analyzed/Projects for Pay/INOVA_Stool_and_Breastmilk_0074 (Ali)/data/metadata/mother_mother_processed.txt',
              sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

  SFig1 <- plot_grid(result_infant$pval_hetamap$gtable,
                     result_breastmilk$pval_hetamap$gtable,
                    ncol = 2, nrow = 1, align = 'h', labels = c("a, infant metadata", "b, maternal-donor metadata"),
                    hjust=-0.1, label_size=14)
  SFig1
  ggsave(filename='Manuscript/figures/SFig1/SFig1_Visit_a.png',
          plot=SFig1, width = 12, height =5., units = "in", dpi = 350)
  ggsave(filename='Manuscript/figures/SFig1/SFig1_Visit_a.pdf',
         plot=SFig1, width = 12, height =5., units = "in", dpi = 350)
}

mother_metadata_a <- mother_metadata
mother_metadata_b <- mother_metadata

mother_metadata_b$`Visit time point` <- "b"
rownames(mother_metadata_b)  <- gsub("a", "b", rownames(mother_metadata_b))
mother_metadata_c <- mother_metadata
mother_metadata_c$`Visit time point` <- "c"
rownames(mother_metadata_c)  <- gsub("a", "c", rownames(mother_metadata_c))
mother_metadata_d <- mother_metadata
mother_metadata_d$`Visit time point` <- "d"
rownames(mother_metadata_d)  <- gsub("a", "d", rownames(mother_metadata_d))
mother_metadata <- rbind(mother_metadata_a, mother_metadata_b, mother_metadata_c, mother_metadata_d)
row.names(mother_metadata)<- paste0(row.names(mother_metadata), "_", "m", spe="")


