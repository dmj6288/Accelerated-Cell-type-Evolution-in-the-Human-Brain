suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_biological/RN006/Kanton/"
Kanton_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Kanton/Kanton_adult.h5")

metadata_file      <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Kanton/metadata_file.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

Macaque_barcodes     <- paste0(metadata_file[metadata_file$Species == "macaque",]$Barcode, "M")
Macaque_Individual      <- metadata_file[metadata_file$Species == "macaque",]$Individual
Macaque_segment      <- metadata_file[metadata_file$Species == "macaque",]$Segment

print(head(metadata_file))
print(Macaque_barcodes)

Macaque_data         <- Kanton_count_matrix[, Macaque_barcodes]
colnames(Macaque_data) <- paste(colnames(Macaque_data), Macaque_Individual, sep = "_")
colnames(Macaque_data) <- paste(colnames(Macaque_data), Macaque_segment, sep = "_")


Macaque_seurat       <- CreateSeuratObject(Macaque_data)

Macaque_seurat@meta.data$orig.ident <- metadata_file[metadata_file$Species == "macaque",]$Individual

Macaque_seurat_filt  <- comprehensive_filtering(Macaque_seurat, "Kanton_PFC_Macaque")

saveRDS(Macaque_seurat_filt, paste0(seurat_object_path, "Kanton_Macaque.rds"))
