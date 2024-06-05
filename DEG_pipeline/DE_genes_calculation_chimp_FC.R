AsuppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(this.path)))

study_results     <- list.files("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_020/")

path              <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_020/"
output_path       <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/original_020/"

coldata_path      <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/coldata/coldata_020/"

cell_types        <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

FDR           <- 0.05
FC            <- 0.3

#Inhibit_#human-vs-Inhibit_#macak_compare_B_vs_Aboot#_79.csv
for (FDR in c(0.05, 0.1, 0.2)){
for (study in study_results[c(5)]){
	human_chimp_results            <- paste0(path, study, "/", study, "_df_human_chimp.rds_sce.rdscounts_ls.rds/")
	chimp_macak_results            <- paste0(path, study, "/", study, "_df_chimp_macak.rds_sce.rdscounts_ls.rds/")
	human_macak_results            <- paste0(path, study, "/", study, "_df_human_macak.rds_sce.rdscounts_ls.rds/")

	print(human_chimp_results)

	human_chimp_results_considered <- list.files(human_chimp_results)
	chimp_macak_results_considered <- list.files(chimp_macak_results)
	human_macak_results_considered <- list.files(human_macak_results)
	
	print(head(human_chimp_results_considered))

	up_regulated_human_gene_number   <- list()
        down_regulated_human_gene_number <- list()
	up_regulated_human_genes         <- list()
        down_regulated_human_genes       <- list()
	up_regulated_human_genes_table   <- list()
        down_regulated_human_genes_table <- list()
	considered_genes_table           <- list()

#       CUSTOM CELL TYPE
	
	cell_types <- unique(tstrsplit(human_chimp_results_considered, "_#")[[1]])

	for (cell_type in cell_types){

		cell_type_human_chimp_results_considered <- list.files(human_chimp_results)[grepl(cell_type, human_chimp_results_considered, fixed = TRUE)]
		cell_type_chimp_macak_results_considered <- list.files(chimp_macak_results)[grepl(cell_type, chimp_macak_results_considered, fixed = TRUE)]
		cell_type_human_macak_results_considered <- list.files(human_macak_results)[grepl(cell_type, human_macak_results_considered, fixed = TRUE)]
		
		print(cell_type_human_chimp_results_considered)

		cell_type_human_chimp_results_considered <- read_csv(paste0(human_chimp_results, cell_type_human_chimp_results_considered[grepl("all_genes", cell_type_human_chimp_results_considered, fixed = TRUE)]))
                cell_type_chimp_macak_results_considered <- read_csv(paste0(chimp_macak_results, cell_type_chimp_macak_results_considered[grepl("all_genes", cell_type_chimp_macak_results_considered, fixed = TRUE)]))
                cell_type_human_macak_results_considered <- read_csv(paste0(human_macak_results, cell_type_human_macak_results_considered[grepl("all_genes", cell_type_human_macak_results_considered, fixed = TRUE)]))
		
		print(head(cell_type_human_chimp_results_considered))
			
		up_raw_chimp_human_genes                        <- na.omit(cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$padj        < FDR &
                                                                                                                      cell_type_human_chimp_results_considered$log2FoldChange > FC,]$gene)

                up_high_expressing_chimp_macak_genes            <- na.omit(cell_type_chimp_macak_results_considered[cell_type_chimp_macak_results_considered$padj        < FDR &
                                                                                                                      cell_type_chimp_macak_results_considered$log2FoldChange > FC,]$gene)

                down_raw_chimp_human_genes                      <- na.omit(cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$padj        < FDR &
                                                                                                                      cell_type_human_chimp_results_considered$log2FoldChange < -FC,]$gene)

                down_high_expressing_chimp_macak_genes          <- na.omit(cell_type_chimp_macak_results_considered[cell_type_chimp_macak_results_considered$padj        < FDR &
                                                                                                                      cell_type_chimp_macak_results_considered$log2FoldChange < -FC,]$gene)

                potentially_interfering_human_macak_genes       <- na.omit(cell_type_human_macak_results_considered[cell_type_human_macak_results_considered$padj   > 0.1,]$gene)
		
		
                up_Berto_human_genes                            <- Reduce(intersect, list(up_raw_chimp_human_genes,
                                                                                                  up_high_expressing_chimp_macak_genes,
                                                                                                  potentially_interfering_human_macak_genes))

                down_Berto_human_genes                          <- Reduce(intersect, list(down_raw_chimp_human_genes,
                                                                                                  down_high_expressing_chimp_macak_genes,
                                                                                                  potentially_interfering_human_macak_genes))	
		#print(cell_type_human_chimp_results_considered)

		up_gene_res_table                               <- cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$gene %in% up_Berto_human_genes  , ]
		down_gene_res_table                             <- cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$gene %in% down_Berto_human_genes, ]
		
		print("up_gene_res_table")
		print(up_gene_res_table)

		up_regulated_human_gene_number[[cell_type]]    <- length(up_Berto_human_genes)
                down_regulated_human_gene_number[[cell_type]]  <- length(down_Berto_human_genes)
                up_regulated_human_genes[[cell_type]]          <- up_Berto_human_genes
                down_regulated_human_genes[[cell_type]]        <- down_Berto_human_genes

		up_regulated_human_genes_table[[cell_type]]    <- up_gene_res_table
	        down_regulated_human_genes_table[[cell_type]]  <- down_gene_res_table


        }
        saveRDS(up_regulated_human_gene_number,   paste0(output_path, study, "/", study, FDR, "_original_CHIMP_URGene_numbers.rds"))
        saveRDS(down_regulated_human_gene_number, paste0(output_path, study, "/", study, FDR, "_original_CHIMP_DRGene_numbers.rds"))
        saveRDS(up_regulated_human_genes,         paste0(output_path, study, "/", study, FDR, "_original_CHIMP_URGenes.rds"))
        saveRDS(down_regulated_human_genes,       paste0(output_path, study, "/", study, FDR, "_original_CHIMP_DRGenes.rds"))
	saveRDS(up_regulated_human_genes_table,   paste0(output_path, study, "/", study, FDR, "_original_CHIMP_URGenes_restable.rds"))
        saveRDS(down_regulated_human_genes_table, paste0(output_path, study, "/", study, FDR, "_original_CHIMP_DRGenes_restable.rds"))
}
}
