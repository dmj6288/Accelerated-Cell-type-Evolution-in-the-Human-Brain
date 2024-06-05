suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))

biological_replicates <- c("RMB161", "RMB196", "RMB295", "RMB307")

technical_replicates <- c(1:4)

data_path    <- "../aligned_10x_data/Ma/data/cellranger_counts/"
h5_data_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Ma/"

for (biological_replicate in biological_replicates){
        print(biological_replicate)
        for (technical_replicate in technical_replicates){
                print(technical_replicate)

                current_data_matrix <- Read10X(paste0(data_path, biological_replicate, "_", technical_replicate, "/raw_feature_bc_matrix/"))
                write10xCounts(paste0(h5_data_path, biological_replicate, "_", technical_replicate, ".h5"  ),
                               current_data_matrix, barcodes = colnames(current_data_matrix),  gene.id = rownames(current_data_matrix),
                               type = "HDF5")

        }
}
