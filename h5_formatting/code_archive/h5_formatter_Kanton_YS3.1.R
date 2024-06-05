suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))

seurat_object_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_custom_doublet/"
Kanton_count_matrix<- Read10X("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Kanton", gene.column = 1,  cell.column = 1)
metadata_file      <- read_excel("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Kanton/metadata_file.xlsx")

source("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_Evo_TS.R")
human_protein_coding_genes        <- readRDS("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/rds_files/human_PC_genes.rds")
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

Human_seurat <- CreateSeuratObject(Human_data)
Chimp_seurat <- CreateSeuratObject(Chimp_data)
Macak_seurat <- CreateSeuratObject(Macak_data)

filtered_human_seurat <- comprehensive_filtering(Human_seurat, "Kanton_PFC")
filtered_chimp_seurat <- comprehensive_filtering(Chimp_seurat, "Kanton_PFC")
filtered_macak_seurat <- comprehensive_filtering(Macak_seurat, "Kanton_PFC")

human_counts <- filtered_human_seurat@assays$RNA@counts
chimp_counts <- filtered_chimp_seurat@assays$RNA@counts
macak_counts <- filtered_macak_seurat@assays$RNA@counts

gene_df_chimp_human <- orthogene::convert_orthologs(gene_df = chimp_counts,
                                                    gene_input = "rownames",
                                                    gene_output = "rownames",
                                                    input_species = "chimp",
                                                    output_species = "human",
                                                    non121_strategy = "drop_both_species",
                                                    method = method)


gene_df_macak_human <- orthogene::convert_orthologs(gene_df = macak_counts,
                                                    gene_input = "rownames",
                                                    gene_output = "rownames",
                                                    input_species = "macaque",
                                                    output_species = "human",
                                                    non121_strategy = "drop_both_species",
                                                    method = method)

gene_df_macak_human <- orthogene::convert_orthologs(gene_df = macak_counts,
                                                    gene_input = "rownames",
                                                    gene_output = "rownames",
                                                    input_species = "macaque",
                                                    output_species = "human",
                                                    non121_strategy = "drop_both_species",
                                                    method = method)


common_gene_set <- Reduce(intersect, list(rownames(human_counts),
                                          rownames(gene_df_chimp_human),
                                          rownames(gene_df_macak_human)))
#                                         rownames(filtered_bonobo_seurat)))

protein_coding_common_gene_set <- intersect(human_protein_coding_genes, common_gene_set)

print(protein_coding_common_gene_set)

combined_matrix <- cbind(human_counts[protein_coding_common_gene_set, ],
                         gene_df_chimp_human[protein_coding_common_gene_set, ],
                         gene_df_macak_human[protein_coding_common_gene_set, ])
#                        filtered_bonobo_seurat@assays$RNA@counts[common_gene_set, ])

print(dim(combined_matrix))

saveRDS(CreateSeuratObject(combined_matrix), paste0(seurat_object_path, "Kanton_seurat_object.rds"))
