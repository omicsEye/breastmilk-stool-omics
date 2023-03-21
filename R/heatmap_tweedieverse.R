# heatmap of Tweedieverse results

all_results <- NULL
meta <- "Breast_milk_collected"
for (meta in colnames(mother_metadata)){
  tryCatch({
    temp_res <- read.delim(
      paste('analysis/Tweedieverse_mother_microbiome_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      row.names = NA
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
  all_results <- rbind(all_results, temp_res)
}

df = all_results[all_results$qval<=0.25,]
title = NA
cell_value = 'qval'
data_label = 'data'
metadata_label = 'metadata'
border_color = 'grey93'
color = colorRampPalette(c("darkblue", "grey90", "darkred"))
col_rotate = 90
first_n = 75
write_to = NA
ggsave(filename=paste('manuscript/figures/Tweedieverse_mother_microbiome.png', sep=''), plot=p,
       width = 5, height = 5, units = "in", dpi = 350)



all_results <- NULL
for (meta in colnames(infant_metadata)){
  tryCatch({
    temp_res <- read.delim(
      paste('analysis/Tweedieverse_infant_microbiome_', meta, '/all_results.tsv',sep = ""),
      sep = '\t',
      header = TRUE,
      fill = FALSE,
      comment.char = "" ,
      check.names = FALSE,
      #row.names = NA
    )
  },
  error = function(e) {
    print(meta)
    print(paste('error:', e))
  })
  all_results <- rbind(all_results, temp_res)
}
df = all_results[all_results$qval<=0.25,]
title = NA
cell_value = 'qval'
data_label = 'data'
metadata_label = 'metadata'
border_color = 'grey93'
color = colorRampPalette(c("darkblue", "grey90", "darkred"))
col_rotate = 90
first_n = 50
write_to = NA
ggsave(filename=paste('manuscript/figures/Tweedieverse_infant_microbiome.png', sep=''), plot=p,
       width = 4.5, height = 4.5, units = "in", dpi = 350)

