suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))

input_path        <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/original_018/"
output_path       <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Functional_Analysis/Functional_Analysis_018/"

cell_types        <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

study_results     <- list.files(input_path)

chimp_combined_regulated_GO            <- list()
human_combined_regulated_GO            <- list()

for (FDR1 in c(0.1)){
for (FDR2 in c(0.05, 0.1, 0.2)){
for (ont in c("MF", "BP", "CC")){

        chimp_combined_regulated_GO_ont            <- list()
        human_combined_regulated_GO_ont            <- list()

for (study in study_results[c(1, 5)]){
   
	chimp_down_regulated          <- readRDS(paste0(input_path, study, "/", study, FDR1, "_original_CHIMP_DRGenes.rds"))
	human_down_regulated          <- readRDS(paste0(input_path, study, "/", study, FDR1, "_original_HUMAN_DRGenes.rds"))
	chimp_up_regulated            <- readRDS(paste0(input_path, study, "/", study, FDR1, "_original_CHIMP_URGenes.rds"))
	human_up_regulated            <- readRDS(paste0(input_path, study, "/", study, FDR1, "_original_HUMAN_URGenes.rds"))
	
	considered_Genes              <- readRDS(paste0(input_path, study, "/", study, FDR1, "_three_species_considered_Genes.rds"))

	#considered_Genes_cell_type_entrez       <- bitr(considered_Genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")		
  
	chimp_combined_regulated_GO_cell_type            <- list()
	human_combined_regulated_GO_cell_type            <- list()
  
	for (cell_type in names(chimp_down_regulated)){
		
		print(length(considered_Genes[[cell_type]]))

		chimp_down_regulated_cell_type          <- chimp_down_regulated[[cell_type]]
                human_down_regulated_cell_type          <- human_down_regulated[[cell_type]] 
		chimp_up_regulated_cell_type            <- chimp_up_regulated[[cell_type]]
		human_up_regulated_cell_type            <- human_up_regulated[[cell_type]]

		human_combined_regulated_cell_type            <- c(human_up_regulated_cell_type, human_down_regulated_cell_type)
		chimp_combined_regulated_cell_type            <- c(chimp_up_regulated_cell_type, chimp_down_regulated_cell_type)
        
		chimp_combined_regulated_cell_type_entrez     <- bitr(chimp_combined_regulated_cell_type,   fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
		human_combined_regulated_cell_type_entrez     <- bitr(human_combined_regulated_cell_type,   fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    
		considered_Genes_cell_type_entrez       <- bitr(c(considered_Genes[[cell_type]], "EEIG1", "EEIG2"), 
								fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

		print(length(considered_Genes_cell_type_entrez$ENTREZID))
      
		chimp_combined_regulated_GO_cell_type[[cell_type]]            <- enrichGO(gene          = chimp_combined_regulated_cell_type_entrez$ENTREZID,
										    universe      = considered_Genes_cell_type_entrez$ENTREZID,
                                                                                    OrgDb         = org.Hs.eg.db,
                                                                                    ont           = ont,
                                                                                    pAdjustMethod = "BH",
                                                                                    pvalueCutoff  = FDR2,
                                                                                    qvalueCutoff  = 0.05,
                                                                                    readable      = TRUE)

		human_combined_regulated_GO_cell_type[[cell_type]]            <- enrichGO(gene          = human_combined_regulated_cell_type_entrez$ENTREZID,
										    universe      = considered_Genes_cell_type_entrez$ENTREZID,
                                                                                    OrgDb         = org.Hs.eg.db,
                                                                                    ont           = ont,
                                                                                    pAdjustMethod = "BH",
                                                                                    pvalueCutoff  = FDR2,
                                                                                    qvalueCutoff  = 0.05,
                                                                                    readable      = TRUE)

		}

		chimp_combined_regulated_GO_ont[[study]]            <- chimp_combined_regulated_GO_cell_type
		human_combined_regulated_GO_ont[[study]]            <- human_combined_regulated_GO_cell_type

	}

        chimp_combined_regulated_GO[[ont]]            <- chimp_combined_regulated_GO_ont
        human_combined_regulated_GO[[ont]]            <- human_combined_regulated_GO_ont

}

saveRDS(chimp_combined_regulated_GO,   paste0(output_path, FDR1, FDR2, "chimp_combined_regulated_GO.rds"))
saveRDS(human_combined_regulated_GO,   paste0(output_path, FDR1, FDR2, "human_combined_regulated_GO.rds"))
}
}
