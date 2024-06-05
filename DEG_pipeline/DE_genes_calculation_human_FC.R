suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(this.path)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(openxlsx)))

study_results     <- list.files("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_018/")

path              <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbres/pseudo_bulk_analysis_results_018/"
output_path       <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/original_018/Repeat_1/"

coldata_path      <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/coldata/coldata_018/"

cell_types        <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")


FDR               <- 0.05
FC                <- 0.3
coldata_folders   <- list.files(coldata_path)
#Inhibit_#human-vs-Inhibit_#macak_compare_B_vs_Aboot#_79.csv

for (FDR in c(0.05)){
for (study in study_results[c(1)]){
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
	considered_human_genes_table     <- list()

	considered_genes_human_chimp     <- list()
	considered_genes_human_macaque   <- list()
	considered_genes_chimp_macaque   <- list()

	three_species_considered_genes          <- list()
	three_species_considered_genes_number   <- list()

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
		
		print(cell_type_human_chimp_results_considered)		
			
		up_raw_human_chimp_genes                        <- na.omit(cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$padj        < FDR &
                                                                                                                      cell_type_human_chimp_results_considered$log2FoldChange < -FC,]$gene)
		print(length(up_raw_human_chimp_genes))
		
		up_high_expressing_human_macak_genes            <- na.omit(cell_type_human_macak_results_considered[cell_type_human_macak_results_considered$padj        < FDR &
                                                                                                                      cell_type_human_macak_results_considered$log2FoldChange >  FC,]$gene)
		print(length(up_high_expressing_human_macak_genes))		
		down_raw_human_chimp_genes                      <- na.omit(cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$padj        < FDR &
                                                                                                                      cell_type_human_chimp_results_considered$log2FoldChange > FC,]$gene)
		
		down_high_expressing_human_macak_genes          <- na.omit(cell_type_human_macak_results_considered[cell_type_human_macak_results_considered$padj        < FDR &
                                                                                                                      cell_type_human_macak_results_considered$log2FoldChange < -FC,]$gene)
		
		potentially_interfering_chimp_macak_genes       <- na.omit(cell_type_chimp_macak_results_considered[cell_type_chimp_macak_results_considered$padj   > 0.1,]$gene)

		print(length(potentially_interfering_chimp_macak_genes))

		up_Berto_human_genes                            <- Reduce(intersect, list(up_raw_human_chimp_genes,
                                                                                                  up_high_expressing_human_macak_genes,
                                                                                                  potentially_interfering_chimp_macak_genes))
		
		down_Berto_human_genes                          <- Reduce(intersect, list(down_raw_human_chimp_genes,
                                                                                                  down_high_expressing_human_macak_genes,
                                                                                                  potentially_interfering_chimp_macak_genes))
			
		up_gene_res_table                               <- cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$gene %in% up_Berto_human_genes  , ]
                down_gene_res_table                             <- cell_type_human_chimp_results_considered[cell_type_human_chimp_results_considered$gene %in% down_Berto_human_genes, ]		
		considered_res_table                            <- cell_type_human_chimp_results_considered
		
		up_regulated_human_gene_number[[cell_type]]    <- length(up_Berto_human_genes)
		down_regulated_human_gene_number[[cell_type]]  <- length(down_Berto_human_genes)
		up_regulated_human_genes[[cell_type]]          <- up_Berto_human_genes
		down_regulated_human_genes[[cell_type]]        <- down_Berto_human_genes

		up_regulated_human_genes_table[[cell_type]]    <- up_gene_res_table
                down_regulated_human_genes_table[[cell_type]]  <- down_gene_res_table
		considered_human_genes_table[[cell_type]]      <- considered_res_table

#		print(cell_type)
#		print(up_Berto_human_genes)
#		print(down_Berto_human_genes)
		
		considered_genes_human_chimp[[cell_type]]      <- na.omit(cell_type_human_chimp_results_considered)
		considered_genes_human_macaque[[cell_type]]    <- na.omit(cell_type_human_macak_results_considered)
		considered_genes_chimp_macaque[[cell_type]]    <- na.omit(cell_type_chimp_macak_results_considered)
		
		three_species_considered_genes_number[[cell_type]]  <- length(Reduce(intersect, list(  considered_genes_human_chimp[[cell_type]]$gene,
										                considered_genes_human_macaque[[cell_type]]$gene,
										                considered_genes_chimp_macaque[[cell_type]]$gene)))
		
		three_species_considered_genes[[cell_type]]         <- Reduce(intersect, list(  considered_genes_human_chimp[[cell_type]]$gene,
                                                                                                considered_genes_human_macaque[[cell_type]]$gene,
                                                                                                considered_genes_chimp_macaque[[cell_type]]$gene))
		
	}
	print(up_regulated_human_gene_number)
	print(down_regulated_human_gene_number)
	
	print(paste0(output_path, study, "/", study, "_original_HUMAN_URGene_numbers.rds"))

	saveRDS(up_regulated_human_gene_number,          paste0(output_path, study, "/", study, FDR, "_original_HUMAN_URGene_numbers.rds"))
        saveRDS(down_regulated_human_gene_number,        paste0(output_path, study, "/", study, FDR, "_original_HUMAN_DRGene_numbers.rds"))
	saveRDS(up_regulated_human_genes,                paste0(output_path, study, "/", study, FDR, "_original_HUMAN_URGenes.rds"))
        saveRDS(down_regulated_human_genes,              paste0(output_path, study, "/", study, FDR, "_original_HUMAN_DRGenes.rds"))
	saveRDS(up_regulated_human_genes_table,          paste0(output_path, study, "/", study, FDR, "_original_HUMAN_URGenes_restable.rds"))
        saveRDS(down_regulated_human_genes_table,        paste0(output_path, study, "/", study, FDR, "_original_HUMAN_DRGenes_restable.rds"))
	saveRDS(considered_human_genes_table,            paste0(output_path, study, "/", study, FDR, "_original_considered_Genes_restable.rds"))

	saveRDS(considered_genes_human_chimp,            paste0(output_path, study, "/", study, FDR, "_human_chimp_considered_Genes.rds"))
	saveRDS(considered_genes_human_macaque,          paste0(output_path, study, "/", study, FDR, "_human_macaque_considered_Genes.rds"))
	saveRDS(considered_genes_human_chimp,            paste0(output_path, study, "/", study, FDR, "_chimp_macaque_considered_Genes.rds"))

	saveRDS(three_species_considered_genes,          paste0(output_path, study, "/", study, FDR, "_three_species_considered_Genes.rds"))
	saveRDS(three_species_considered_genes_number,   paste0(output_path, study, "/", study, FDR, "_three_species_considered_Gene_Number.rds"))

}
}
for(study in coldata_folders[c(5)]){
  
	coldata_files          <- list.files(paste0(coldata_path, study))
	chimp_macak_coldata    <- readRDS(paste0(coldata_path, study, "/", coldata_files[1]))
	human_chimp_coldata    <- readRDS(paste0(coldata_path, study, "/", coldata_files[2]))
  
	chimp_macak_count      <- count(chimp_macak_coldata$cluster_id)
	human_chimp_count      <- count(human_chimp_coldata$cluster_id)
  
	cells                  <- unique(tstrsplit(chimp_macak_count$x, "_#")[[1]])
	
	print(chimp_macak_count)	
  
	counts                 <- data.frame(round(100*chimp_macak_count[grepl("chimp", chimp_macak_count$x), "freq"]/sum(chimp_macak_count[grepl("chimp", chimp_macak_count$x), "freq"]), 2),
                                             round(100*chimp_macak_count[grepl("macak", chimp_macak_count$x),"freq"]/sum(chimp_macak_count[grepl("macak", chimp_macak_count$x), "freq"]), 2),
                                             round(100*human_chimp_count[grepl("human", human_chimp_count$x), "freq"]/sum(human_chimp_count[grepl("human", human_chimp_count$x), "freq"]), 2))
  
  
	colnames(counts)       <- c("chimp", "macaque", "human")
  
  
	counts <- rbind(counts, c(sum(chimp_macak_count[grepl("chimp", chimp_macak_count$x), "freq"]), 
                            sum(chimp_macak_count[grepl("macak", chimp_macak_count$x), "freq"]), 
                            sum(human_chimp_count[grepl("human", human_chimp_count$x), "freq"])))
    
  
	rownames(counts)       <- c(cells, "TOTAL")
	write.xlsx(as.data.frame(counts), paste0(output_path, study, "/", study, "_MajorCellTypes.xlsx"), sheetName=study, append = TRUE, rowNames = TRUE)
}

