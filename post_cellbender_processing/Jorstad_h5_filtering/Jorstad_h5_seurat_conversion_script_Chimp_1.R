suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path  <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN015/Jorstad/"
Jorstad_seurat      <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Jorstad/potentially_processed_data/chimp_SCT_UMI_expression_matrix.RDS")

Jorstad_count_matrix<-Jorstad_seurat@assays$RNA@counts

metadata_file       <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Jorstad/potentially_processed_data/Chimp_metadata.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

Chimp_donors        <- unique(metadata_file$donor)

for (replicate in Chimp_donors[c(1:4)]){

	replicate_barcodes <- metadata_file[metadata_file$donor %in% replicate, ]$sample_id

	Chimp_data         <- Jorstad_count_matrix[, replicate_barcodes]

	Chimp_seurat       <- CreateSeuratObject(Chimp_data)

	Chimp_seurat_filt  <- comprehensive_filtering(Chimp_seurat, paste0("Jorstad_MTG_Chimp_", replicate))

	saveRDS(Chimp_seurat_filt, paste0(seurat_object_path, "Jorstad_", replicate, ".rds"))
}
