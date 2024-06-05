suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(this.path)))


path          <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/bootstrapping_results_015BT"

output_path   <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/consolidated_lists_015BT"

study_results <- list.dirs(path, full.names = FALSE, recursive = FALSE)
print(study_results)
#output_path   <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/bootstrapping_results_015BT/"

cell_types    <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

FDR           <- 0.05
FC            <- 0.3

#Inhibit_#human-vs-Inhibit_#macak_compare_B_vs_Aboot#_79.csv

up_regulated_human_study_wise_list   <- list()
down_regulated_human_study_wise_list <- list()
up_regulated_chimp_study_wise_list   <- list()
down_regulated_chimp_study_wise_list <- list()

considered_genes_study_wise_list     <- list()

for (study in study_results){

	up_regulated_human_gene_number   <- list()
        down_regulated_human_gene_number <- list()
	up_regulated_chimp_gene_number   <- list()
        down_regulated_chimp_gene_number <- list()

	considered_genes                 <- list()

	for (cell_type in cell_types){
		
		print(cell_type)
		cell_type_path <- paste(path, study,  cell_type, sep = "/")
		all_gene_files <- list.files(cell_type_path)

		#print(all_gene_files)
		
		chimp_DR_files   <- all_gene_files[grepl("CHIMP_DR",   all_gene_files, fixed = TRUE)]
		chimp_UR_files   <- all_gene_files[grepl("CHIMP_UR",   all_gene_files, fixed = TRUE)]
		human_DR_files   <- all_gene_files[grepl("HUMAN_DR",   all_gene_files, fixed = TRUE)]
		human_UR_files   <- all_gene_files[grepl("HUMAN_UR",   all_gene_files, fixed = TRUE)]
		considered_files <- all_gene_files[grepl("considered", all_gene_files, fixed = TRUE)]
		
		#print(human_UR_files)
		
		current_up_regulated_human_gene_number   <- c()
		current_down_regulated_human_gene_number <- c()
		current_up_regulated_chimp_gene_number   <- c()
                current_down_regulated_chimp_gene_number <- c()

	        current_considered_genes                 <- c()

		for (step in seq(1000, 10000, by = 1000)){
			
			current_up_regulated_human_gene_number   <- Filter(Negate(is.null), c(current_up_regulated_human_gene_number,   readRDS(paste(cell_type_path, human_UR_files[tstrsplit(human_UR_files, "_")[[3]] %in% as.character(step)], sep = "/"))))
			current_down_regulated_human_gene_number <- Filter(Negate(is.null), c(current_down_regulated_human_gene_number, readRDS(paste(cell_type_path, human_DR_files[tstrsplit(human_DR_files, "_")[[3]] %in% as.character(step)], sep = "/"))))
			current_up_regulated_chimp_gene_number   <- Filter(Negate(is.null), c(current_up_regulated_chimp_gene_number,   readRDS(paste(cell_type_path, chimp_UR_files[tstrsplit(chimp_UR_files, "_")[[3]] %in% as.character(step)], sep = "/"))))
			current_down_regulated_chimp_gene_number <- Filter(Negate(is.null), c(current_down_regulated_chimp_gene_number, readRDS(paste(cell_type_path, chimp_DR_files[tstrsplit(chimp_DR_files, "_")[[3]] %in% as.character(step)], sep = "/"))))

			current_considered_genes                 <- Filter(Negate(is.null), c(current_considered_genes, readRDS(paste(cell_type_path, considered_files[tstrsplit(considered_files, "_")[[3]] %in% as.character(step)], sep = "/"))))
			
			
		}
		up_regulated_human_gene_number[[cell_type]]   <- current_up_regulated_human_gene_number
		down_regulated_human_gene_number[[cell_type]] <- current_down_regulated_human_gene_number
		up_regulated_chimp_gene_number[[cell_type]]   <- current_up_regulated_chimp_gene_number
		down_regulated_chimp_gene_number[[cell_type]] <- current_down_regulated_chimp_gene_number

		considered_genes[[cell_type]]                 <- current_considered_genes
	}

	up_regulated_human_study_wise_list[[study]]   <- up_regulated_human_gene_number
	down_regulated_human_study_wise_list[[study]] <- down_regulated_human_gene_number
	up_regulated_chimp_study_wise_list[[study]]   <- up_regulated_chimp_gene_number
	down_regulated_chimp_study_wise_list[[study]] <- down_regulated_chimp_gene_number

	considered_genes_study_wise_list[[study]]     <- considered_genes	
}
print(up_regulated_human_study_wise_list)
print(paste(output_path, "consolidated_human_UR.rds", sep ="/"))

saveRDS(up_regulated_human_study_wise_list,   paste(output_path, "consolidated_human_UR.rds", sep = "/"))
saveRDS(down_regulated_human_study_wise_list, paste(output_path, "consolidated_human_DR.rds", sep = "/"))
saveRDS(up_regulated_chimp_study_wise_list,   paste(output_path, "consolidated_chimp_UR.rds", sep = "/"))
saveRDS(down_regulated_chimp_study_wise_list, paste(output_path, "consolidated_chimp_DR.rds", sep = "/"))
saveRDS(considered_genes_study_wise_list,     paste(output_path, "consolidated_considered.rds", sep = "/"))
