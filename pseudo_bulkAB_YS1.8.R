#suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(S4Vectors)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(apeglm)))
suppressMessages(suppressWarnings(library(ashr)))
suppressMessages(suppressWarnings(library(png)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(comprehenr)))
suppressMessages(suppressWarnings(library(hash)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(tibble)))

support_functions_path     <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

pseudo_bulk_ready_path     <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbds/pseudo_bulk_ready_data_021/"
results_path               <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_021/"
#bootstrapped_gene_selection_matrix_list <- readRDS("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/rds_files/bootstrapped_gene_selection_matrix_list.rds")


source(paste(support_functions_path, "get_dds_resultsAvsB_age_sex_018.R", sep = ""))

target_folders               <- list.files(pseudo_bulk_ready_path)

print(target_folders)
for (folder in target_folders[c(5)]){
	
	target_files               <- list.files(paste0(pseudo_bulk_ready_path, folder))
	counts_ls_files            <- sort(target_files[grepl("counts_ls", target_files)])
	metadata_ls_files          <- sort(target_files[grepl("metadata_ls", target_files)])

	pseudo_bulk_ready_data_to_study_map <- hash(counts_ls_files, metadata_ls_files)
	print(pseudo_bulk_ready_data_to_study_map)	

	for (counts_ls in keys(pseudo_bulk_ready_data_to_study_map)){
		print(counts_ls)
	
		species <- substr(strsplit(counts_ls, ".rds")[[1]], start = 3, stop = nchar(strsplit(counts_ls, ".rds")[[1]]))
	
		counts_ls_file    <- readRDS(paste(pseudo_bulk_ready_path, folder, "/", counts_ls, sep = ""))
		metadata_ls_file <- readRDS(paste(pseudo_bulk_ready_path, folder, "/", pseudo_bulk_ready_data_to_study_map[[counts_ls]], sep = ""))

		get_dds_resultsAvsB(padj_cutoff = 0.05, folder, 
                                    counts_ls = counts_ls_file, metadata_ls = metadata_ls_file, 
                                    count_ls  = counts_ls,      path = results_path)
	}
}
