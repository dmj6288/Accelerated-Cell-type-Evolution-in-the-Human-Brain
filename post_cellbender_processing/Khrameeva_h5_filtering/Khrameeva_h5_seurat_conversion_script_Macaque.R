suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_biological/RN009/Khrameeva/"

Khrameeva_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Khrameeva/Macaque.h5")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

Macaque_data         <- Khrameeva_count_matrix
colnames(Macaque_data) <- paste0(colnames(Macaque_data), "M")

Macaque_seurat       <- CreateSeuratObject(Macaque_data)

Macaque_seurat_filt  <- comprehensive_filtering(Macaque_seurat, "Khrameeva_ACC_Macaque")

saveRDS(Macaque_seurat_filt, paste0(seurat_object_path, "Khrameeva_Macaque.rds"))
