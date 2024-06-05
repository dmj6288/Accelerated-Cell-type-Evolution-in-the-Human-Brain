suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))


study_results <- list.files("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_bootAB/")

path          <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_bootAB/Ma/"
output_path   <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/bootstrap_species_specific_genes_number/"

cell_types    <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

FDR           <- 0.05
FC            <- 0.3
current_path               <- this.path()
bootstrap_Step             <- as.numeric(strsplit(strsplit(a, "_")[[1]][8], ".R")[[1]])

print(bootstrap_Step)

step_size                  <- 500
#Inhibit_#human-vs-Inhibit_#macak_compare_B_vs_Aboot#_79.csv

for (study in study_results[4]){
	human_chimp_results            <- paste0(path, study, "/", study, "_df_human_chimp.rds_sce.rdscounts_ls.rds/")
	chimp_macak_results            <- paste0(path, study, "/", study, "_df_chimp_macak.rds_sce.rdscounts_ls.rds/")
	human_macak_results            <- paste0(path, study, "/", study, "_df_human_macak.rds_sce.rdscounts_ls.rds/")

	human_chimp_results_considered <- list.files(human_chimp_results)
	chimp_macak_results_considered <- list.files(chimp_macak_results)
	human_macak_results_considered <- list.files(human_macak_results)

	for (cell_type in cell_types){

		cell_type_human_chimp_results_considered <- list.files(human_chimp_results)[grepl(cell_type, human_chimp_results_considered, fixed = TRUE)][1]
		cell_type_chimp_macak_results_considered <- list.files(chimp_macak_results)[grepl(cell_type, chimp_macak_results_considered, fixed = TRUE)][1]
		cell_type_human_macak_results_considered <- list.files(human_macak_results)[grepl(cell_type, human_macak_results_considered, fixed = TRUE)][1]

		cell_type_human_chimp_results_considered_base <- strsplit(cell_type_human_chimp_results_considered, "#_")[1]
		cell_type_chimp_macak_results_considered_base <- strsplit(cell_type_human_chimp_results_considered, "#_")[1]
		cell_type_human_macak_results_considered_base <- strsplit(cell_type_human_chimp_results_considered, "#_")[1]
		
		up_regulated_human_genes                 <- hash()
		down_regulated_human_genes               <- hash()		

		for (bootstrap_sample in c(((bootstrap_Step*step_size)+1):((bootstrap_Step+1)*step_size))){

			boot_cell_type_human_chimp_results_considered <-  read_csv(paste0(cell_type_human_chimp_results_considered_base, "#_", bootstrap_sample, ".csv"))
	                boot_cell_type_chimp_macak_results_considered <-  read_csv(paste0(cell_type_chimp_macak_results_considered_base, "#_", bootstrap_sample, ".csv"))
		        boot_cell_type_human_macak_results_considered <-  read_csv(paste0(cell_type_human_macak_results_considered_base, "#_", bootstrap_sample, ".csv"))
			
			up_raw_chimp_human_genes                        <- boot_cell_type_human_chimp_results_considered[boot_cell_type_human_chimp_results_considered$padj        < FDR &
                                                                                                                      boot_cell_type_human_chimp_results_considered$log2FoldChange < -FC,]$gene

			up_high_expressing_chimp_macak_genes            <- boot_cell_type_chimp_macak_results_considered[boot_cell_type_chimp_macak_results_considered$padj        < FDR &
                                                                                                                      boot_cell_type_chimp_macak_results_considered$log2FoldChange < -FC,]$gene

			down_raw_chimp_human_genes                      <- boot_cell_type_human_chimp_results_considered[boot_cell_type_human_chimp_results_considered$padj        < FDR &
                                                                                                                      boot_cell_type_human_chimp_results_considered$log2FoldChange > FC,]$gene 

			down_high_expressing_chimp_macak_genes          <- boot_cell_type_chimp_macak_results_considered[boot_cell_type_chimp_macak_results_considered$padj        < FDR &
                                                                                                                      boot_cell_type_chimp_macak_results_considered$log2FoldChange > FC,]$gene

			potentially_interfering_human_macak_genes       <- boot_cell_type_human_macak_results_considered[boot_cell_type_human_macak_results_considered$padj   > 0.1,]$gene
			
			up_Berto_human_genes                            <- Reduce(intersect, list(up_raw_chimp_human_genes,
                                                                                                  up_high_expressing_chimp_macak_genes,
                                                                                                  potentially_interfering_human_macak_genes))

			down_Berto_human_genes                          <- Reduce(intersect, list(down_raw_chimp_human_genes,
                                                                                                  down_high_expressing_chimp_macak_genes,
                                                                                                  potentially_interfering_human_macak_genes))

			up_regulated_human_gene_number[[bootstrap_sample]]    <- length(up_Berto_human_genes)
			down_regulated_human_gene_number[[bootstrap_sample]]  <- length(down_Berto_human_genes)
		}
		saveRDS(up_regulated_human_gene_number, paste0(cell_type, "bootstrap_URGene_numbers.rds"))
		saveRDS(down_regulated_human_gene_number, paste0(cell_type, "bootstrap_DRGene_numbers.rds"))
	}
}
