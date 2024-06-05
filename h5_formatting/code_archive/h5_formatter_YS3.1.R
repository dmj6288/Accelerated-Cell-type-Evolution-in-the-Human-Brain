suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))

cellbender_ready_data <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Khrameeva"

Human_data <- Read10X("../raw_data/Khrameeva/Human/",   gene.column = 2,  cell.column = 2)
Chimp_data <- Read10X("../raw_data/Khrameeva/Chimp/",   gene.column = 2,  cell.column = 2)
Macak_data <- Read10X("../raw_data/Khrameeva/Macaque/",   gene.column = 2,  cell.column = 2)

print(dim(Human_data))
print(dim(Chimp_data))
print(dim(Macak_data))

colnames(Human_data) <- paste0(colnames(Human_data), "H")
colnames(Chimp_data) <- paste0(colnames(Chimp_data), "C")
colnames(Macak_data) <- paste0(colnames(Macak_data), "M")

write10xCounts(paste0(cellbender_ready_data, "_Human.h5"  ), Human_data, barcodes = colnames(Human_data),  gene.id = rownames(Human_data), type = "HDF5")
write10xCounts(paste0(cellbender_ready_data, "_Chimp.h5"  ), Chimp_data, barcodes = colnames(Chimp_data),  gene.id = rownames(Chimp_data), type = "HDF5")
write10xCounts(paste0(cellbender_ready_data, "_Macaque.h5"), Macak_data, barcodes = colnames(Macak_data),  gene.id = rownames(Macak_data), type = "HDF5")
