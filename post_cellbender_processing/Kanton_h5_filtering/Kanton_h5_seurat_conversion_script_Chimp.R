suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN022/Kanton/"
Kanton_count_matrix<- Read10X_h5("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Kanton/Kanton_adult.h5")

metadata_file      <- read_excel("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Kanton/metadata_file.xlsx")

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_004.R")

Chimp_barcodes     <- paste0(metadata_file[metadata_file$Species == "chimp",]$Barcode, "C")

Chimp_Individual   <- metadata_file[metadata_file$Species == "chimp",]$Individual
Chimp_segment      <- metadata_file[metadata_file$Species == "chimp",]$Segment

print(head(metadata_file))
#print(Chimp_barcodes)

Chimp_data         <- Kanton_count_matrix[, Chimp_barcodes]
colnames(Chimp_data) <- paste(colnames(Chimp_data), Chimp_Individual, sep = "_")
Chimp_data         <- Chimp_data[, !duplicated(colnames(Chimp_data))]

#colnames(Chimp_data) <- paste(colnames(Chimp_data), Chimp_segment, sep = "_")

Chimp_seurat       <- CreateSeuratObject(counts = Chimp_data)
print(str(Chimp_seurat))

#Chimp_seurat@meta.data$orig.ident <- metadata_file[metadata_file$Species == "chimp",]$Individual

Chimp_seurat_filt  <- comprehensive_filtering(Chimp_seurat, "Kanton_PFC_Chimp")

saveRDS(Chimp_seurat_filt, paste0(seurat_object_path, "Kanton_Chimp.rds"))
