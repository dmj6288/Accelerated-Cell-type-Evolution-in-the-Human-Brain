suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(this.path)))

study_results     <- list.files("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_018/")

path              <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_018/"
output_path       <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/original_018/Repeat_1/restables_all_comparisons/"

cell_types        <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

FDR           <- 0.05
FC            <- 0.3

#Inhibit_#human-vs-Inhibit_#macak_compare_B_vs_Aboot#_79.csv
for (FDR in c(0.05, 0.1, 0.2)){
for (study in study_results[c(1, 5)]){
	human_chimp_results            <- paste0(path, study, "/", study, "_df_human_chimp.rds_sce.rdscounts_ls.rds/")
	chimp_macak_results            <- paste0(path, study, "/", study, "_df_chimp_macak.rds_sce.rdscounts_ls.rds/")
	human_macak_results            <- paste0(path, study, "/", study, "_df_human_macak.rds_sce.rdscounts_ls.rds/")

	print(human_chimp_results)

	human_chimp_results_considered <- list.files(human_chimp_results)
	chimp_macak_results_considered <- list.files(chimp_macak_results)
	human_macak_results_considered <- list.files(human_macak_results)
	
	print(head(human_chimp_results_considered))

	cell_type_human_chimp_results_considered_all_cells <- c()
	cell_type_chimp_macak_results_considered_all_cells <- c()
	cell_type_human_macak_results_considered_all_cells <- c()	

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

		cell_type_human_chimp_results_considered$Cell <- cell_type
		cell_type_chimp_macak_results_considered$Cell <- cell_type
		cell_type_human_macak_results_considered$Cell <- cell_type

		cell_type_human_chimp_results_considered_all_cells <- rbind(cell_type_human_chimp_results_considered_all_cells, cell_type_human_chimp_results_considered)
		cell_type_chimp_macak_results_considered_all_cells <- rbind(cell_type_chimp_macak_results_considered_all_cells, cell_type_chimp_macak_results_considered)
		cell_type_human_macak_results_considered_all_cells <- rbind(cell_type_human_macak_results_considered_all_cells, cell_type_human_macak_results_considered)
		


        }

	saveRDS(cell_type_human_chimp_results_considered_all_cells, paste0(output_path, study, "/", study, FDR, "_original_HC_DESEQ2_results_table.rds"))
	saveRDS(cell_type_human_macak_results_considered_all_cells, paste0(output_path, study, "/", study, FDR, "_original_HM_DESEQ2_results_table.rds"))
	saveRDS(cell_type_chimp_macak_results_considered_all_cells, paste0(output_path, study, "/", study, FDR, "_original_CM_DESEQ2_results_table.rds"))

}
}
