suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))
suppressWarnings(suppressMessages(library(readxl)))


technical_replicates <- c(1:4)

data_path    <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Kanton/raw_data"
h5_data_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Kanton/Individual/"

meta_data    <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Kanton/metadata/metadata_file_withSpecies_symbol.xlsx")
data_matrix  <- Read10X(paste0(data_path), gene.column = 1,  cell.column = 1)

set1         <- c("3", "HC1", "MI1")
set2         <- c("CHD2", "HG1", "MG2")
set3         <- c("HE1", "MH1", "BA1")

sets         <- list(set1, set2, set3)

for (set in c(set1, set2, set3)){
	print(set)
	barcodes_in_set <- paste0(meta_data[meta_data$Individual %in% set, ]$Barcode, 
				  meta_data[meta_data$Individual %in% set, ]$Symbol)

	counts_in_set   <- data_matrix[, barcodes_in_set]
	write10xCounts(paste0(h5_data_path, "set_", set, ".h5"  ),
                       counts_in_set, barcodes = barcodes_in_set,  gene.id = rownames(data_matrix),
                       type = "HDF5")
	
}
#data_matrix <- Read10X(paste0(data_path), gene.column = 1,  cell.column = 1)

#write10xCounts(paste0(h5_data_path, "Kanton_adult.h5"  ),
#               data_matrix, barcodes = colnames(data_matrix),  gene.id = rownames(data_matrix),
#               type = "HDF5")

