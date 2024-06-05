suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(this.path)))

study_results <- list.files("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_015BT/")

path          <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_015BT/"
output_path   <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/bootstrapping_results_015BT/"

cell_types    <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

FDR           <- 0.05
FC            <- 0.3

#Inhibit_#human-vs-Inhibit_#macak_compare_B_vs_Aboot#_79.csv

for (study in study_results[4]){
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
	considered_genes                 <- list()

	for (cell_type in cell_types){

		cell_type_human_chimp_results_considered <- list.files(human_chimp_results)[grepl(cell_type, human_chimp_results_considered, fixed = TRUE)][[1]][1]
		cell_type_chimp_macak_results_considered <- list.files(chimp_macak_results)[grepl(cell_type, chimp_macak_results_considered, fixed = TRUE)][[1]][1]
		cell_type_human_macak_results_considered <- list.files(human_macak_results)[grepl(cell_type, human_macak_results_considered, fixed = TRUE)][[1]][1]
		
		print(cell_type_human_chimp_results_considered)
	
		cell_type_human_chimp_results_base <- strsplit(cell_type_human_chimp_results_considered, "100")[[1]][1]
                cell_type_chimp_macak_results_base <- strsplit(cell_type_chimp_macak_results_considered, "100")[[1]][1]
                cell_type_human_macak_results_base <- strsplit(cell_type_human_macak_results_considered, "100")[[1]][1]

		print(cell_type_human_chimp_results_base)

		for (index in seq(6100, 7000, by=100)){

			cell_type_human_chimp_results_current <- paste0(cell_type_human_chimp_results_base, index)
			cell_type_chimp_macak_results_current <- paste0(cell_type_chimp_macak_results_base, index)
			cell_type_human_macak_results_current <- paste0(cell_type_human_macak_results_base, index)

			cell_type_human_chimp_results_res_tables <- readRDS(paste0(human_chimp_results, cell_type_human_chimp_results_current, ".rds"))

			cell_type_chimp_macak_results_res_tables <- readRDS(paste0(chimp_macak_results, cell_type_chimp_macak_results_current, ".rds"))

			cell_type_human_macak_results_res_tables <- readRDS(paste0(human_macak_results, cell_type_human_macak_results_current, ".rds"))

			print(paste0(human_chimp_results, cell_type_human_chimp_results_current, ".rds"))
			print(head(c((index):(index+400))))

			for (list_index in c((index-99):(index))){
				print(list_index)
				print(head(cell_type_human_chimp_results_res_tables[[list_index]]))
				
				up_raw_chimp_human_genes                        <- na.omit(cell_type_human_chimp_results_res_tables[[list_index]][cell_type_human_chimp_results_res_tables[[list_index]]$padj        < FDR &
                                                                                                                                          cell_type_human_chimp_results_res_tables[[list_index]]$log2FoldChange > FC,]$gene)

				up_high_expressing_chimp_macak_genes            <- na.omit(cell_type_chimp_macak_results_res_tables[[list_index]][cell_type_chimp_macak_results_res_tables[[list_index]]$padj        < FDR &
                                                                                                                                          cell_type_chimp_macak_results_res_tables[[list_index]]$log2FoldChange > FC,]$gene)

				down_raw_chimp_human_genes                      <- na.omit(cell_type_human_chimp_results_res_tables[[list_index]][cell_type_human_chimp_results_res_tables[[list_index]]$padj        < FDR &
                                                                                                                                          cell_type_human_chimp_results_res_tables[[list_index]]$log2FoldChange < -FC,]$gene)

				down_high_expressing_chimp_macak_genes          <- na.omit(cell_type_chimp_macak_results_res_tables[[list_index]][cell_type_chimp_macak_results_res_tables[[list_index]]$padj        < FDR &
                                                                                                                                          cell_type_chimp_macak_results_res_tables[[list_index]]$log2FoldChange < -FC,]$gene)

				potentially_interfering_human_macak_genes       <- na.omit(cell_type_human_macak_results_res_tables[[list_index]][cell_type_human_macak_results_res_tables[[list_index]]$padj   > 0.1,]$gene)

				up_Berto_human_genes                            <- Reduce(intersect, list(up_raw_chimp_human_genes,
                                                                                                  up_high_expressing_chimp_macak_genes,
                                                                                                  potentially_interfering_human_macak_genes))

				down_Berto_human_genes                          <- Reduce(intersect, list(down_raw_chimp_human_genes,
                                                                                                  down_high_expressing_chimp_macak_genes,
                                                                                                  potentially_interfering_human_macak_genes))
				up_regulated_human_gene_number[[list_index]]    <- length(up_Berto_human_genes)
				down_regulated_human_gene_number[[list_index]]  <- length(down_Berto_human_genes)
			}
		}
	paste0(output_path, study, "/", cell_type, "/", study, "_", cell_type, "_", index, "_boot_CHIMP_URGene_numbers_cellbender.rds")

        saveRDS(up_regulated_human_gene_number,   paste0(output_path, study, "/", cell_type, "/",  study, "_", cell_type, "_", index, "_boot_CHIMP_URGene_numbers_cellbender.rds"))
        saveRDS(down_regulated_human_gene_number, paste0(output_path, study, "/", cell_type, "/",  study, "_", cell_type, "_", index, "_boot_CHIMP_DRGene_numbers_cellbender.rds"))
	
	}
}
