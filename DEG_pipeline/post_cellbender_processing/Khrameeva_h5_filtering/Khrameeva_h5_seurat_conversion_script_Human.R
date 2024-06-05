suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_biological/RN009/Khrameeva/"

Khrameeva_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Khrameeva/Human.h5")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

Human_data         <- Khrameeva_count_matrix
colnames(Human_data) <- paste0(colnames(Human_data), "H")

Human_seurat       <- CreateSeuratObject(Human_data)

Human_seurat_filt  <- comprehensive_filtering(Human_seurat, "Khrameeva_ACC_Human")

saveRDS(Human_seurat_filt, paste0(seurat_object_path, "Khrameeva_Human.rds"))
