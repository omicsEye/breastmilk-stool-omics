
library(xlsx)
setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")
mother_metadata_a <- mother_metadata[mother_metadata$Visit=='a',]
cont_cols <- colnames(mother_metadata_a)[sapply(mother_metadata_a, function(x) length(unique(x)))>10]
cat_cols <- setdiff(colnames(mother_metadata_a),cont_cols)

cont_cols <- cont_cols[cont_cols!="External_ID"]
#create summary data for continuous variables
res_cont <- data.frame(do.call(rbind,lapply(mother_metadata_a[,cont_cols],summary)))

#add column with variable name (unwise to store in rownames)



res_cat <- do.call(rbind,lapply(cat_cols, function(x){
  res <- data.frame(table(mother_metadata_a[,x],useNA="always")) #added it to deal with missings, can be changes
  res$variable <- x
  colnames(res)[1] <- "Category"
  res
}
))
sum_table <- as.data.frame(summary(mother_metadata_a))

write.xlsx2(res_cont, 'manuscript/Tables/STable1_V2.xlsx', sheetName = "summary_breastmilk_continuous_vars",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(res_cat, 'manuscript/Tables/STable1_V2.xlsx', sheetName = "summary_breatmilk_categorical_vars",
            col.names = TRUE, row.names = TRUE, append = T)





infant_metadata_a <- infant_metadata[infant_metadata$Visit=='a',]
cont_cols <- colnames(infant_metadata)[sapply(infant_metadata, function(x) length(unique(x)))>10]
cat_cols <- setdiff(colnames(infant_metadata),cont_cols)

cont_cols <- cont_cols[cont_cols!="External_ID"]
#create summary data for continuous variables
res_cont <- data.frame(do.call(rbind,lapply(infant_metadata[,cont_cols],summary)))

#add column with variable name (unwise to store in rownames)



res_cat <- do.call(rbind,lapply(cat_cols, function(x){
  res <- data.frame(table(infant_metadata[,x],useNA="always")) #added it to deal with missings, can be changes
  res$variable <- x
  colnames(res)[1] <- "Category"
  res
}
))
sum_table <- as.data.frame(summary(mother_metadata_a))


write.xlsx2(res_cont, 'manuscript/Tables/STable1_V2.xlsx', sheetName = "summary_stool_continuous_vars",
            col.names = TRUE, row.names = TRUE, append = T)
write.xlsx2(res_cat, 'manuscript/Tables/STable1_V2.xlsx', sheetName = "summary_stool_categorical_vars",
            col.names = TRUE, row.names = TRUE, append = T)
write.table(summary(infant_metadata), 'data/summary_infant_stool.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

