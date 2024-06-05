suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_biological/RN006/Kanton/"
Kanton_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Kanton/Kanton_adult.h5")
metadata_file      <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Kanton/metadata_file.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

Human_barcodes     <- paste0(metadata_file[metadata_file$Species == "human",]$Barcode, "H")

Human_Individual      <- metadata_file[metadata_file$Species == "human",]$Individual
Human_segment      <- metadata_file[metadata_file$Species == "human",]$Segment

print(head(metadata_file))
print(Human_barcodes)

Human_data         <- Kanton_count_matrix[, Human_barcodes]
colnames(Human_data) <- paste(colnames(Human_data), Human_Individual, sep = "_")
colnames(Human_data) <- paste(colnames(Human_data), Human_segment, sep = "_")

Human_seurat       <- CreateSeuratObject(Human_data)

Human_seurat@meta.data$orig.ident <- metadata_file[metadata_file$Species == "human",]$Individual

Human_seurat_filt  <- comprehensive_filtering(Human_seurat, "Kanton_PFC_Human")

saveRDS(Human_seurat_filt, paste0(seurat_object_path, "Kanton_Human.rds"))
