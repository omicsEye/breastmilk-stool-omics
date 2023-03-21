#read in required packages
require(readxl)
require(tidyverse)

#set the working directory from which the files will be read from
setwd("/Users/rah/Box/INOVA_Stool_and_Breastmilk_0074/Analysis/Metagenomics/Pathoscope_MS/pathoscope_output_BM/fungal/")

#create a list of the files from your target directory
file_list <- list.files(path="/Users/rah/Box/INOVA_Stool_and_Breastmilk_0074/Analysis/Metagenomics/Pathoscope_MS/pathoscope_output_BM/fungal/")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- data.frame( )

#had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  temp_data <-read.delim(
    file_list[i],
    sep = '\t',
    header = T,
    fill = F,
    comment.char = "" ,
    check.names = F,
    row.names = 1,
    skip = 1
  )
 if (i==1){
   dataset <- temp_data[,"Final Guess", drop = F]
   colnames(dataset)[i] <- file_list[i]
 }
 else{
   dataset[rownames(temp_data), file_list[i]] <- temp_data[,"Final Guess", drop = F]
   colnames(dataset)[i] <- file_list[i]
 }
}

### clean smaple snames to reflect names in the study files ####
x <- sapply(colnames(dataset), function(x) strsplit(x, "-")[[1]], USE.NAMES=FALSE)
colnames(dataset) <- x[2,]


#######
colnames(dataset)[colnames(dataset) %in% id_mapper$`GW Sample ID`]
barcodes <- id_mapper[match(colnames(dataset),id_mapper$`GW Sample ID`), "Event ID (barcode scanned)"]
colnames(dataset) <-  unlist(as.list(barcodes))

setwd("/Users/rah/Dropbox/Ali-Docs/Research_docs/Projects/INOVA_Breastmilk/")

write.table( dataset,'data/BM_fungal.txt',
             sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("compbiomed/animalcules")

library(animalcules)
