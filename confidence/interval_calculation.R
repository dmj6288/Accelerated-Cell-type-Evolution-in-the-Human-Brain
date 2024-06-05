suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(this.path)))

specific_gene_numbers_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/bootstrap_species_specific_genes_number/"

cell_types                 <- list.files(specific_gene_numbers_path)
print(cell_types)

regulations                 <- c("UR", "DR")

for (regulation in regulations){
	for (cell_type in cell_types){

		cell_type_files           <- list.files(paste0(specific_gene_numbers_path, cell_type))
		regulated_cell_type_files <- cell_type_files[grepl(regulation, cell_type_files, fixed = TRUE)]
		full_list_of_gene_numbers <- c()

		for (regulated_cell_type_file in regulated_cell_type_files){

			current_gene_list <- readRDS(paste0(specific_gene_numbers_path, cell_type, "/", regulated_cell_type_file))
			full_list_of_gene_numbers <- c(full_list_of_gene_numbers, current_gene_list)
			full_list_of_gene_numbers[sapply(full_list_of_gene_numbers, is.null)] <- NULL
			saveRDS(full_list_of_gene_numbers, paste0("interval_results/", cell_type, "_", regulation, ".rds"))
		}
	}
}
