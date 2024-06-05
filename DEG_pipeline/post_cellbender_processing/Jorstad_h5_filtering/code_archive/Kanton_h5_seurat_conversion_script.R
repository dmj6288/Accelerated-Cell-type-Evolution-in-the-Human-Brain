suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_001/"
Kanton_count_matrix<- Read10X("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Kanton", gene.column = 1,  cell.column = 1)
metadata_file      <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Kanton/metadata_file.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_001.R")
human_protein_coding_genes        <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/rds_files/human_PC_genes.rds")
method <- "gprofiler"

human_barcodes     <- paste0(metadata_file[metadata_file$Species == "human",]$Barcode, "H")
chimp_barcodes     <- paste0(metadata_file[metadata_file$Species == "chimp",]$Barcode, "C")
macaque_barcodes   <- paste0(metadata_file[metadata_file$Species == "macaque",]$Barcode, "M")

Human_data <- Kanton_count_matrix[, human_barcodes]
Chimp_data <- Kanton_count_matrix[, chimp_barcodes]
Macak_data <- Kanton_count_matrix[, macaque_barcodes]

common_gene_set <- Reduce(intersect, list(rownames(Human_data),
                                          rownames(Chimp_data),
                                          rownames(Macak_data)))

