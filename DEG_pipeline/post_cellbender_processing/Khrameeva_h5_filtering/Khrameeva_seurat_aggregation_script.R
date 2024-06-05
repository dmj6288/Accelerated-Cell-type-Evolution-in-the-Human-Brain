suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

seurat_folder_path           <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_biological/RN009/Khrameeva/"
seurat_files                 <- list.files(seurat_folder_path)
human_protein_coding_genes   <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/rds_files/human_PC_genes.rds")
output_path                  <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_final/RN009/"

species_seurat_object_list <- list()
for (sample in seurat_files){
	species_seurat_object_list[[sample]] <- readRDS(paste0(seurat_folder_path, sample))
}

species_seurat <- merge(species_seurat_object_list[[1]], y = c(species_seurat_object_list[[2]],
                                                               species_seurat_object_list[[3]]), project = "Khrameeva")
	
print(str(species_seurat))
Khrameeva_counts                  <- species_seurat@assays$RNA@counts
protein_coding_common_gene_set <- intersect(human_protein_coding_genes, rownames(Khrameeva_counts))

print(protein_coding_common_gene_set)
combined_matrix <- Khrameeva_counts[protein_coding_common_gene_set, ]

print(dim(combined_matrix))

saveRDS(CreateSeuratObject(combined_matrix), paste0(output_path, "Khrameeva_seurat_object.rds"))
