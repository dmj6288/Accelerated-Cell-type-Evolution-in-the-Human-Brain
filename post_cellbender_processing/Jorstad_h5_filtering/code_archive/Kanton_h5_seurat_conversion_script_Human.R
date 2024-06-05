suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path  <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN015/Jorstad/"
Jorstad_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Jorstad/human_counts.h5")

metadata_file       <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Jorstad/potentially_processed_data/Human_metadata.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

Human_donors        <- unique(metadata_file$donor)

print(head(metadata_file))
print(Human_barcodes)

for (replicate in Human_donors){

	replicate_barcodes <- metadata_file[metadata_file$donor %in% replicate, ]$sample_id

	Human_data         <- Jorstad_count_matrix[, replicate_barcodes]

	Human_seurat       <- CreateSeuratObject(Human_data)

	Human_seurat@meta.data$orig.ident <- metadata_file[metadata_file$Species == "human",]$Individual

	Human_seurat_filt  <- comprehensive_filtering(Human_seurat, paste0("Jorstad_MTG_Human_", donor))

	saveRDS(Human_seurat_filt, paste0(seurat_object_path, "Jorstad_", replicate, ".rds"))
}
