suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

seurat_object_path                <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN018/Khrameeva/"
output_path                       <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_final/RN018/"

filtered_human_seurat             <- readRDS(paste0(seurat_object_path, "Khrameeva_Human.rds"))
filtered_chimp_seurat             <- readRDS(paste0(seurat_object_path, "Khrameeva_Chimp.rds"))
filtered_macak_seurat             <- readRDS(paste0(seurat_object_path, "Khrameeva_Macaque.rds"))

human_protein_coding_genes        <- readRDS("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/rds_files/human_PC_genes.rds")
method                            <- "gprofiler"

human_counts                      <- filtered_human_seurat@assays$RNA@counts
chimp_counts                      <- filtered_chimp_seurat@assays$RNA@counts
macak_counts                      <- filtered_macak_seurat@assays$RNA@counts

#colnames(human_counts)            <- paste0(colnames(human_counts), "H")
#colnames(chimp_counts)            <- paste0(colnames(chimp_counts), "C")
#colnames(macak_counts)            <- paste0(colnames(macak_counts), "M")


gene_df_chimp_human               <- orthogene::convert_orthologs(gene_df = chimp_counts,
                                                    gene_input = "rownames",
                                                    gene_output = "rownames",
                                                    input_species = "chimp",
                                                    output_species = "human",
                                                    non121_strategy = "drop_both_species",
                                                    method = method)


gene_df_macak_human               <- orthogene::convert_orthologs(gene_df = macak_counts,
                                                    gene_input = "rownames",
                                                    gene_output = "rownames",
                                                    input_species = "macaque",
                                                    output_species = "human",
                                                    non121_strategy = "drop_both_species",
                                                    method = method)

common_gene_set                   <- Reduce(intersect, list(rownames(human_counts),
                                          rownames(gene_df_chimp_human),
                                          rownames(gene_df_macak_human)))

protein_coding_common_gene_set    <- intersect(human_protein_coding_genes, common_gene_set)

print(protein_coding_common_gene_set)

combined_matrix <- cbind(human_counts[protein_coding_common_gene_set, ],
                         gene_df_chimp_human[protein_coding_common_gene_set, ],
                         gene_df_macak_human[protein_coding_common_gene_set, ])

print(head(combined_matrix))
print(dim(combined_matrix))

saveRDS(CreateSeuratObject(combined_matrix), paste0(output_path, "Khrameeva_seurat_object.rds"))
