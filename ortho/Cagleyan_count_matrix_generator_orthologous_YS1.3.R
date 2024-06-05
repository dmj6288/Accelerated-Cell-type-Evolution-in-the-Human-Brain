suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(orthogene)))
suppressWarnings(suppressMessages(library(biomaRt)))

method <- "gprofiler"

Cagleyan_human_seurat   <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN018/Caglayan/h1.rds")
Cagleyan_chimp_seurat   <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN018/Caglayan/c1.rds")
Cagleyan_macaque_seurat <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN018/Caglayan/m1.rds")

Caglayan_human_count_matrix   <- Cagleyan_human_seurat@assays$RNA@counts
Caglayan_chimp_count_matrix   <- Cagleyan_chimp_seurat@assays$RNA@counts
Caglayan_macaque_count_matrix <- Cagleyan_macaque_seurat@assays$RNA@counts

colnames(Caglayan_human_count_matrix)   <- paste(tstrsplit(colnames(Caglayan_human_count_matrix), "_")[[2]],   tstrsplit(colnames(Caglayan_human_count_matrix), "_")[[1]],sep = "_")
colnames(Caglayan_chimp_count_matrix)   <- paste(tstrsplit(colnames(Caglayan_chimp_count_matrix), "_")[[2]],   tstrsplit(colnames(Caglayan_chimp_count_matrix), "_")[[1]],sep = "_")
colnames(Caglayan_macaque_count_matrix) <- paste(tstrsplit(colnames(Caglayan_macaque_count_matrix), "_")[[2]], tstrsplit(colnames(Caglayan_macaque_count_matrix), "_")[[1]],sep = "_")

seurat_object_path    <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_final/RN018/"
metadata_file         <- read.table("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/raw_data/Caglayan/metadata.txt", sep = "\t")

human_protein_coding_genes        <- readRDS("~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/rds_files/human_PC_genes.rds")

human_barcodes     <- intersect(rownames(metadata_file[metadata_file$Species == "human",]),   colnames(Caglayan_human_count_matrix))
chimp_barcodes     <- intersect(rownames(metadata_file[metadata_file$Species == "chimp",]),   colnames(Caglayan_chimp_count_matrix))
macaque_barcodes   <- intersect(rownames(metadata_file[metadata_file$Species == "macaque",]), colnames(Caglayan_macaque_count_matrix))

print(human_barcodes)
print(Caglayan_human_count_matrix)

human_counts <- Caglayan_human_count_matrix[, human_barcodes]
chimp_counts <- Caglayan_chimp_count_matrix[, chimp_barcodes]
macaque_counts <- Caglayan_macaque_count_matrix[, macaque_barcodes]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

human_genes  <- rownames(human_counts)
G_list_human <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=human_genes,mart= mart)
G_list_human <- G_list_human[G_list_human$hgnc_symbol != "",]
human_counts <- human_counts[G_list_human[G_list_human$hgnc_symbol != "",]$ensembl_gene_id, ]
rownames(human_counts) <- G_list_human[G_list_human$hgnc_symbol != "",]$hgnc_symbol

chimp_genes  <- rownames(chimp_counts)
G_list_chimp <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=chimp_genes,mart= mart)
G_list_chimp <- G_list_chimp[G_list_chimp$hgnc_symbol != "",]
chimp_counts <- chimp_counts[G_list_chimp[G_list_chimp$hgnc_symbol != "",]$ensembl_gene_id, ]
rownames(chimp_counts) <- G_list_chimp[G_list_chimp$hgnc_symbol != "",]$hgnc_symbol

macaque_genes  <- rownames(macaque_counts)
G_list_macaque <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=macaque_genes,mart= mart)
G_list_macaque <- G_list_macaque[G_list_macaque$hgnc_symbol != "",]
macaque_counts <- macaque_counts[G_list_macaque[G_list_macaque$hgnc_symbol != "",]$ensembl_gene_id, ]
rownames(macaque_counts) <- G_list_macaque[G_list_macaque$hgnc_symbol != "",]$hgnc_symbol

print(dim(human_counts))
print(dim(chimp_counts))
print(dim(macaque_counts))

gene_df_chimp_human <- orthogene::convert_orthologs(gene_df = chimp_counts,
                                                    gene_input = "rownames",
                                                    gene_output = "rownames",
#						    standardise_genes = TRUE,
                                                    input_species = "chimp",
                                                    output_species = "human",
                                                    non121_strategy = "drop_both_species",
                                                    method = method)


gene_df_macak_human <- orthogene::convert_orthologs(gene_df = macaque_counts,
                                                    gene_input = "rownames",
                                                    gene_output = "rownames",
#						    standardise_genes = TRUE,
                                                    input_species = "macaque",
                                                    output_species = "human",
                                                    non121_strategy = "drop_both_species",
                                                    method = method)

common_gene_set <- Reduce(intersect, list(rownames(human_counts),
                                          rownames(gene_df_chimp_human),
                                          rownames(gene_df_macak_human)))

protein_coding_common_gene_set <- intersect(human_protein_coding_genes, common_gene_set)

#print(common_gene_set)

combined_matrix <- cbind(human_counts[protein_coding_common_gene_set, ],
                         gene_df_chimp_human[protein_coding_common_gene_set, ],
                         gene_df_macak_human[protein_coding_common_gene_set, ])

print(dim(combined_matrix))

protein_coding_common_gene_set <- intersect(human_protein_coding_genes, common_gene_set)

saveRDS(CreateSeuratObject(combined_matrix), paste0(seurat_object_path, "Cagleyan_seurat_object.rds"))
