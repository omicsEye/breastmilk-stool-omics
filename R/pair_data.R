

###### data pairing #########

# pair infants' microbiome and metabolites
mapped_infant_mgx_mbx <- as.list(id_mapper[match(colnames(infant_microbiome), id_mapper$`GW Sample ID`), "MGX_MBX_Map"])$MGX_MBX_Map
metabolites_infant <- metabolites[mapped_infant_mgx_mbx,]
colnames(metabolites_infant) <- as.list(id_mapper[match(rownames(metabolites_infant), id_mapper$MGX_MBX_Map), "GW Sample ID"])$`GW Sample ID`

infant_microbiome_metabolite_samples <- intersect(colnames(infant_microbiome), rownames(metabolites))
infant_microbiome_metabolite <- infant_microbiome[,infant_microbiome_metabolite_samples]
infant_metabolite_microbiome <- infant_metabolite[infant_microbiome_metabolite_samples,]
infant_metabolite_microbiome <- as.data.frame(t(infant_metabolite_microbiome))
infant_metabolite_microbiome <- infant_metabolite_microbiome[rowSums(!is.na(infant_metabolite_microbiome)) > 10,]
write.table(infant_microbiome_metabolite, 'data/infant_microbiome_metabolite.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(infant_metabolite_microbiome, 'data/infant_metabolite_microbiome.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

#Pair BGCs Mother - Infant
rownames(breastmilk_BGC) <- gsub("_m", "", rownames(breastmilk_BGC))
rownames(stool_BGC) <- gsub("_f", "", rownames(stool_BGC))

BGCs_samples <- intersect(row.names(breastmilk_BGC), rownames(stool_BGC))
stool_BGC2 <- stool_BGC[BGCs_samples,]
breastmilk_BGC2 <- breastmilk_BGC[BGCs_samples,]
write.table(t(stool_BGC2), 'data/stool_BGC.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(t(breastmilk_BGC2), 'data/breastmilk_BGC.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

#Pair BGCs Mother - Infant
common_samples <- intersect(row.names(breastmilk_BGC), rownames(infant_metabolite))
stool_BGC_metabolite <- stool_BGC[BGCs_samples,]
infant_metabolite_BGC <- infant_metabolite[BGCs_samples,]
write.table(t(stool_BGC_metabolite), 'data/stool_BGC_metabolite.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(t(infant_metabolite_BGC), 'data/infant_metabolite_BGC.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)



# pair infant's microbiome with mother microbiome
infant_microbiome_mother_microbiome_samples <- intersect(colnames(infant_microbiome), colnames(mother_microbiome))
infant_microbiome_mother_microbiome <- infant_microbiome[,infant_microbiome_mother_microbiome_samples]
mother_microbiome_infant_microbiome <- mother_microbiome[,infant_microbiome_mother_microbiome_samples]
write.table(infant_microbiome_mother_microbiome, 'data/infant_microbiome_mother_microbiome.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(mother_microbiome_infant_microbiome, 'data/mother_microbiome_infant_microbiome.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

# pair infant's metabolite with mothers' microbiome
mapped_mother_mgx_mbx <- as.list(id_mapper[match(colnames(mother_microbiome), id_mapper$`GW Sample ID`), "MGX_MBX_Map"])$MGX_MBX_Map
metabolites_with_mother_microbiome <- metabolites[mapped_mother_mgx_mbx,]
rownames(metabolites_with_mother_microbiome) <- as.list(id_mapper[match(rownames(metabolites_with_mother_microbiome), id_mapper$MGX_MBX_Map), "GW Sample ID"])$`GW Sample ID`
#colnames(mother_microbiome_withmbx_id) <- mapped_mother_mgx_mbx$MGX_MBX_Map
#infant_microbiome_withmbx_id <- as.data.frame(t(infant_microbiome_withmbx_id))
infant_metabolite_mother_microbiome_samples <- intersect(colnames(mother_microbiome), rownames(metabolites_with_mother_microbiome))


infant_metabolite_mother_microbiome <- metabolites[infant_metabolite_mother_microbiome_samples,]
mother_microbiome_infant_metabolite <- mother_microbiome[,infant_metabolite_mother_microbiome_samples]
infant_metabolite_mother_microbiome <- as.data.frame(t(infant_metabolite_mother_microbiome))
infant_metabolite_mother_microbiome <- infant_metabolite_mother_microbiome[ rowSums(!is.na(infant_metabolite_mother_microbiome)) > 10,]
write.table(infant_metabolite_mother_microbiome, 'data/infant_metabolite_mother_microbiome.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(mother_microbiome_infant_metabolite, 'data/mother_microbiome_infant_metabolite.txt',
            sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)





