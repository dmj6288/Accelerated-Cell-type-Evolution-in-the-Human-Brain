suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

seurat_folder_path           <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN015/Jorstad/"
seurat_files                 <- list.files(seurat_folder_path)
human_protein_coding_genes   <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/rds_files/human_PC_genes.rds")
output_path                  <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_final/RN015/"

all_files                    <- list.files(seurat_folder_path)

species_seurat_object_list <- list()

for (sample in all_files){
	species_seurat_object_list[[sample]] <- readRDS(paste0(seurat_folder_path, sample))
}

species_seurat <- merge(species_seurat_object_list[[1]], y = c(species_seurat_object_list[[2]],
                                                               species_seurat_object_list[[3]]), project = "Jorstad")
	
print(str(species_seurat))
Jorstad_counts                  <- species_seurat@assays$RNA@counts
protein_coding_common_gene_set <- intersect(human_protein_coding_genes, rownames(Jorstad_counts))

print(protein_coding_common_gene_set)
combined_matrix <- Jorstad_counts[protein_coding_common_gene_set, ]

combined_matrix <- combined_matrix[, !duplicated(colnames(combined_matrix))]

print(colnames(combined_matrix))
print(dim(combined_matrix))

saveRDS(CreateSeuratObject(combined_matrix), paste0(output_path, "Jorstad_seurat_object.rds"))
