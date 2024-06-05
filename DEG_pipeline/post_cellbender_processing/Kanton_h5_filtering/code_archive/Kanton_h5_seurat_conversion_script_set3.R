suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN011/Kanton/"
Kanton_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Kanton/set_3.h5")

meta_data           <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data_processing_folder/Kanton/metadata/metadata_file_withSpecies_symbol.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

set3              <- c("HE1", "MH1")

set3_barcodes     <- paste0(meta_data[meta_data$Individual %in% set3, ]$Barcode,
                            meta_data[meta_data$Individual %in% set3, ]$Symbol)

set3_data         <- Kanton_count_matrix[, set3_barcodes]
set3_seurat       <- CreateSeuratObject(set3_data)

set3_seurat@meta.data$orig.ident <- meta_data[meta_data$Individual %in% set3, ]$Individual

set3_seurat_filt  <- comprehensive_filtering(set3_seurat, "Kanton_PFC_set3")

saveRDS(set3_seurat_filt, paste0(seurat_object_path, "Kanton_set3.rds"))
