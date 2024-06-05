suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(hash)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(stringr)))


support_functions <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

source(paste(support_functions,"comprehensive_filtering_Evo_TS.R", sep =""))

#Chimp_data <- Read10X("raw_data/Chimp/",   gene.column = 1,  cell.column = 1)
Chimp_data <- Read10X_h5("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_data/Ma_Chimp_filtered.h5", use.names = TRUE, unique.features = TRUE)

colnames(Chimp_data) <- paste0(colnames(Chimp_data), "C")

Chimp_seurat <- CreateSeuratObject(Chimp_data)

Chimp_seurat <- comprehensive_filtering(Chimp_seurat, "Shaojie_chimp")
print("Filtering Done")
saveRDS(Chimp_seurat, "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Ma/Shaojie_Chimp.rds")
rm(Chimp_seurat)


#Human_data <- Read10X("raw_data/Human/",   gene.column = 1,  cell.column = 1)
Human_data <- Read10X_h5("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_data/Ma_Human_filtered.h5", use.names = TRUE, unique.features = TRUE)

colnames(Human_data) <- paste0(colnames(Human_data), "C")

Human_seurat <- CreateSeuratObject(Human_data)

Human_seurat <- comprehensive_filtering(Human_seurat, "Shaojie_Human")
print("Filtering Done")
saveRDS(Human_seurat, "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Ma/Shaojie_Human.rds")
rm(Human_seurat)

#Macaque_data <- Read10X("raw_data/Macaque/",   gene.column = 1,  cell.column = 1)
Macaque_data <- Read10X_h5("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_data/Ma_Macaque_filtered.h5", use.names = TRUE, unique.features = TRUE)

colnames(Macaque_data) <- paste0(colnames(Macaque_data), "C")

Macaque_seurat <- CreateSeuratObject(Macaque_data)

Macaque_seurat <- comprehensive_filtering(Macaque_seurat, "Shaojie_Macaque")
print("Filtering Done")
saveRDS(Macaque_seurat, "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Ma/Shaojie_Macaque.rds")
rm(Macaque_seurat)
