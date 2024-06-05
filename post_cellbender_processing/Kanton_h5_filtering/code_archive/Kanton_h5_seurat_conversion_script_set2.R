suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN011/Kanton/"
Kanton_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Kanton/set_2.h5")

meta_data           <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Kanton/metadata/metadata_file_withSpecies_symbol.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

set2              <- c("CHD2", "HE1", "MG2")

set2_barcodes     <- paste0(meta_data[meta_data$Individual %in% set2, ]$Barcode,
                            meta_data[meta_data$Individual %in% set2, ]$Symbol)

set2_data         <- Kanton_count_matrix[, set2_barcodes]
set2_seurat       <- CreateSeuratObject(set2_data)

set2_seurat@meta.data$orig.ident <- meta_data[meta_data$Individual %in% set2,]$Individual

set2_seurat_filt  <- comprehensive_filtering(set2_seurat, "Kanton_PFC_set2")

saveRDS(set2_seurat_filt, paste0(seurat_object_path, "Kanton_set2.rds"))
