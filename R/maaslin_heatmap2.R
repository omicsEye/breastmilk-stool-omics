maaslin_heatmap <- function(maaslin_output, title = "", cell_value = "Q.value", data_label = 'Data', metadata_label = 'Metadata',
                            border_color = "grey93", color = colorRampPalette(c("blue","grey90", "red"))(50)) {#)
  library(ggplot2)
  library(pheatmap)
  # read MaAsLin output
  df <- read.table( maaslin_output,
                    header = TRUE, sep = "\t", fill = TRUE, comment.char = "" , check.names = FALSE)
  metadata <- df$Variable
  data <- df$Feature
  value <- NA
  # values to use for coloring the heatmap
  if (cell_value == "P.value"){
    value <- -log(df$P.value)*sign(df$Coefficient)
    value <- pmax(-2, pmin(2, value))
  }else if(cell_value == "Q.value"){
    value <- -log(df$Q.value)*sign(df$Coefficient)
    value <- pmax(-2, pmin(2, value))
  }else if(cell_value == "Coefficient"){
    value <- df$Coefficient
  }
  n <- length(unique(metadata))
  m <- length(unique(data))
  a = matrix(0, nrow=n, ncol=m)
  for (i in 1:length(metadata)){
    #if (abs(a[as.numeric(metadata[i]), as.numeric(metadata[i])]) > abs(value[i]))
    #  next
    a[as.numeric(metadata)[i], as.numeric(data)[i]] <- value[i]
  }
  rownames(a) <- levels(metadata)
  colnames(a) <- levels(data)
  
  #colnames(a) <-  sapply(colnames(a), pcl.sub )
  
  p <- pheatmap(a, cellwidth = 7, cellheight = 7,   # changed to 3
                main = title,
                fontsize = 6,
                kmeans_k = NA,
                border=TRUE,
                show_rownames = T, show_colnames = T,
                scale="none",
                #clustering_method = "complete",
                cluster_rows = FALSE, cluster_cols = TRUE,
                clustering_distance_rows = "euclidean", 
                clustering_distance_cols = "euclidean",
                legend=TRUE,
                border_color = border_color,
                color = color,
                treeheight_row=0,
                treeheight_col=0,
                display_numbers = matrix(ifelse(a > 0.0, "+", ifelse(a < 0.0, "|", "")),  nrow(a)))
  return(p)
} # srb
maaslin_output <- '~/Downloads/ali.txt'
srb_ddr <- maaslin_heatmap ('~/Downloads/ali.txt', title = "", cell_value = "Q.value", data_label = 'Data', 
                            metadata_label = 'Metadata', border_color = "grey93", color = colorRampPalette(c("blue","grey90", "red"))(500))

ggsave(filename='~/Documents/01_srb_ddr.pdf', #'~/Dropbox (Huttenhower Lab)/hutlab/long/sulfur/figures/bugs/01_srb_ddr.pdf', 
       plot=srb_ddr$gtable, width = 135, height = 100, units = "mm", dpi = 350)

dev.off()