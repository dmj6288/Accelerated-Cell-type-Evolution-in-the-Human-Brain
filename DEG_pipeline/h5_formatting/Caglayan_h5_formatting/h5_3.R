suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DropletUtils)))

biological_replicates <- c("m1", "m2", "m3", "m4")

data_path    <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/aligned_10x_data/Caglayan"

h5_data_path <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/h5_data/Caglayan/"

for (biological_replicate in biological_replicates){

	current_data_matrix <- Read10X(paste0(data_path, "/03_seurat_qc_annotate_", biological_replicate), gene.column = 1,  cell.column = 1)
	print(head(current_data_matrix))

	write10xCounts(paste0(h5_data_path, biological_replicate, ".h5"  ), 
                               current_data_matrix, barcodes = colnames(current_data_matrix),  gene.id = rownames(current_data_matrix), 
			       type = "HDF5")

}


