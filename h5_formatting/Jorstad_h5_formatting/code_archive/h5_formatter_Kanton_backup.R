suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))

biological_replicates <- c("HSB106", "HSB628", "HSB189", "HSB340")

technical_replicates <- c(1:4)

data_path    <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Kanton/raw_data"
h5_data_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Kanton/"

meta_data    <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Kanton/metadata/metadata_file.xlsx")

set1         <- c("CHF 3", "HC1", "MI1")
set2         <- c("CHD2", "HE1", "MG2")
set3         <- c("HE1", "MH1", "BA1")


for (set in list(set1, set2, set3){
	print(set)
}
#current_data_matrix <- Read10X(paste0(data_path), gene.column = 1,  cell.column = 1)

#write10xCounts(paste0(h5_data_path, "Kanton_adult.h5"  ),
#               current_data_matrix, barcodes = colnames(current_data_matrix),  gene.id = rownames(current_data_matrix),
#               type = "HDF5")

