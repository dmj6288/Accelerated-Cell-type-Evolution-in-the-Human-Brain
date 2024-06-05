suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))
suppressWarnings(suppressMessages(library(readxl)))


technical_replicates <- c(1:4)

data_path    <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Jorstad/"
h5_data_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Jorstad/"

meta_data    <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Jorstad/Human_metadata.xlsx")
data_matrix  <- readRDS(paste0(data_path, "human_SCT_UMI_expression_matrix.RDS"))

data_matrix  <- data_matrix@assays$RNA@counts

donors       <- unique(meta_data$donor)

for (set in donors){

	print(set)
	barcodes_in_set <- meta_data[meta_data$donor %in% set, ]$sample_id

	counts_in_set   <- data_matrix[, barcodes_in_set]
	write10xCounts(paste0(h5_data_path, "set_", set, ".h5"  ),
                       counts_in_set, barcodes = barcodes_in_set,  gene.id = rownames(data_matrix),
                       type = "HDF5")
	
}
#data_matrix <- Read10X(paste0(data_path), gene.column = 1,  cell.column = 1)

#write10xCounts(paste0(h5_data_path, "Kanton_adult.h5"  ),
#               data_matrix, barcodes = colnames(data_matrix),  gene.id = rownames(data_matrix),
#               type = "HDF5")

