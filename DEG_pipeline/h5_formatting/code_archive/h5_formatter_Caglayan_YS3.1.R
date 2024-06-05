suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))

method <- "gprofiler"
cellbender_ready_data <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Caglayan"
Cagleyan_count_matrix <- Read10X("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Caglayan/",   gene.column = 2,  cell.column = 2)
metadata_file         <- read.table("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Caglayan/metadata.txt", sep = "\t")

human_barcodes     <- rownames(metadata_file[metadata_file$Species == "human",])
chimp_barcodes     <- rownames(metadata_file[metadata_file$Species == "chimp",])
macaque_barcodes   <- rownames(metadata_file[metadata_file$Species == "macaque",])

Human_data <- Cagleyan_count_matrix[, human_barcodes]
Chimp_data <- Cagleyan_count_matrix[, chimp_barcodes]
Macak_data <- Cagleyan_count_matrix[, macaque_barcodes]

write10xCounts(paste0(cellbender_ready_data, "_Human.h5"  ), Human_data, barcodes = colnames(Human_data),  gene.id = rownames(Human_data), type = "HDF5")
write10xCounts(paste0(cellbender_ready_data, "_Chimp.h5"  ), Chimp_data, barcodes = colnames(Chimp_data),  gene.id = rownames(Chimp_data), type = "HDF5")
write10xCounts(paste0(cellbender_ready_data, "_Macaque.h5"), Macak_data, barcodes = colnames(Macak_data),  gene.id = rownames(Macak_data), type = "HDF5")

