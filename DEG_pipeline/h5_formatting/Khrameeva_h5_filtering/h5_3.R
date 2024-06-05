suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))

data_path    <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/aligned_10x_data/Khrameeva/Macaque/"
h5_data_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Khrameeva/"

mtx      <- paste0(data_path, "matrix.mtx.gz")
cells    <- paste0(data_path, "barcodes.tsv.gz")
features <- paste0(data_path, "features.tsv.gz")

#counts   <- ReadMtx(mtx = mtx, cells = cells, features = features)

counts <- Read10X(data_path) 
print(str(counts)) 
write10xCounts(paste0(h5_data_path, "Macaque.h5"),
               counts, barcodes = colnames(counts),  gene.id = rownames(counts),
               type = "HDF5")
