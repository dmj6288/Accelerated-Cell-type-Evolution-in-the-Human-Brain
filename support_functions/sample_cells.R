library(Seurat)

# Function to sample cells from each cell type
sample_cells <- function(seurat_object, cells_per_type = 100) {
	# Find unique cell types
	cell_types <- unique(Idents(seurat_object))
  
	# Initialize a list to store sampled cells
	sampled_cells <- list()
  
	# Loop through each cell type and sample cells
	for(cell_type in cell_types) {

		cells_in_type <- WhichCells(seurat_object, ident = cell_type)
		if(length(cells_in_type) > cells_per_type) {
			sampled_cells[[cell_type]] <- sample(cells_in_type, cells_per_type)
		} else {
			sampled_cells[[cell_type]] <- cells_in_type
		}
	}

  
	# Combine all sampled cell barcodes
	sampled_cells_all <- unlist(sampled_cells)
  
	# Subset the Seurat object to only include the sampled cells
	sampled_seurat_object <- subset(seurat_object, cells = sampled_cells_all)
  
	return(sampled_seurat_object)
}

