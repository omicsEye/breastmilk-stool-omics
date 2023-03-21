library(ggplot2)
library(cowplot)

# Gestational age (column C “in weeks” and column E “in groups” in infant tab of metadata spreadsheet. If no differences are found you may want to look at the first sample “a” for these infants only where the biggest differences may be found)
# Mode of delivery, vaginal delivery versus CS (column U from maternal tab of metadata spreadsheet. I recommend combining all C-sections together i.e. look at vaginal delivery versus C section (both scheduled and unscheduled combined).
# Mother’s country of birth (column J from maternal tab of metadata spreadsheet)
# Race of Mother (column I from maternal tab of metadata spreadsheet)
# Ethnicity of Mother (column H from maternal tab of metadata spreadsheet
# BMI category of Mother (column G from maternal tab of metadata spreadsheet)
# Maternal antibiotic use during pregnancy (column Z from maternal tab of metadata spreadsheet)
# Maternal antibiotic use during delivery ( column AB from maternal tab of metadata spreadsheet)




Gestational_age_category_associations<- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Tweedieverse/Cross/Breastmilk/Tweedieverse_mother_microbiome_infant_Gestational_age_day/figures/Gestational_age_day_gg_associations.RDS")
Birth_weight_category_associations<- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Birth_weight_category/figures/Birth_weight_category_gg_associations.RDS")
Visit_associations<- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Visit/figures/Visit_gg_associations.RDS")
Method_of_delivery_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Method_of_delivery/figures/Method_of_delivery_gg_associations.RDS")
Country_of_living_from_birth_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Country_of_living_from_birth/figures/Country_of_living_from_birth_gg_associations.RDS")
Rac_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Race/figures/Race_gg_associations.RDS")
Ethnicity_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Ethnicity/figures/Ethnicity_gg_associations.RDS")
BMI_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_BMI_category/figures/BMI_category_gg_associations.RDS")
Antibiotic_use_during_delivery_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Antibiotic_use_during_delivery/figures/Antibiotic_use_during_delivery_gg_associations.RDS")
PROM_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_PROM/figures/PROM_gg_associations.RDS")
If_sample_fortified_with_what_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_If_sample_fortified_with_what/figures/If_sample_fortified_with_what_gg_associations.RDS")
VRE_swab_result_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_VRE_swab_result/figures/VRE_swab_result_gg_associations.RDS")
Other_medical_conditions_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Other_medical_conditions/figures/Other_medical_conditions_gg_associations.RDS")
Breast_milk_collected_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Breast_milk_collected/figures/Breast_milk_collected_gg_associations.RDS")
Alcohol_use_during_pregnancy_associations <- readRDS("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/analysis/Maaslin2_maternal_microbiome_Alcohol_use_during_pregnancy/figures/Alcohol_use_during_pregnancy_gg_associations.RDS")

## do plots
a_theme <- theme(
  plot.title = element_text(size =7, face = 'bold'),
  axis.title.x = element_text(size = 6),
  axis.text.x = element_text(size = 6),
  axis.title.y = element_text(size = 6),
  axis.text.y = element_text(size = 5))

fig2 <-   plot_grid(Gestational_age_category_associations[[1]] + a_theme,
                    Gestational_age_category_associations[[3]] + a_theme , #C-reactive protein boxplox
                    Method_of_delivery_associations[[1]] + a_theme, #glucose vs. uridine scatterplot
                    #Country_of_living_from_birth_associations[[1]] + a_theme, #glucose vs. citrulline scatterplot
                    #Method_of_delivery_association[[2]] + a_theme, #C-reactive protein vs. Kynurenine
                    Rac_associations[[1]] + a_theme,
                    #Rac_associations[[2]] + a_theme,
                    BMI_associations[[1]] + a_theme,
                    #Antibiotic_use_during_delivery_associations[[3]] + a_theme,
                    #PROM_associations[[1]] + a_theme
                    VRE_swab_result_associations[[1]] + a_theme,
                    VRE_swab_result_associations[[2]] + a_theme ,
                    #Other_medical_conditions_associations[[2]]+ a_theme,
                    Alcohol_use_during_pregnancy_associations[[1]]+ a_theme,
          #rel_widths = c(1, 1, 1,1),
          #labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', "i"),
          label_size = 7, ncol = 4, nrow = 2)
fig2
fig2_2 <- plot_grid(fig2, alpha_diversity_plot,
                    rel_heights = c(4,3),
                    label_size = 7, ncol = 1, nrow = 2)
ggsave(filename = '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/figures/fig2_2_associations.pdf', plot=fig2_2, width = 183, height = 200, units = "mm", dpi = 350)
ggsave(filename = '/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/figures/fig2_2_associations.png', plot=fig2_2, width = 183, height = 200, units = "mm", dpi = 350)



draw_plot(infant_alpha_diversity_plots$Visit + theme(
  axis.title.x = element_text(size = 7),
  axis.text.x = element_text(size = 7),
  axis.title.y = element_text(size = 7),
  axis.text.y = element_text(size = 5)), x = .28, y = 0, width = .4, height = .45) +
  draw_plot(alpha_diversity_plots$Visit + theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 5)), x = .68, y = 0, width = .32, height = .45)
