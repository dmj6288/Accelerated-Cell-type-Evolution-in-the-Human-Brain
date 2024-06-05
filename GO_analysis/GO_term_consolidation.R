suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(writexl)))
suppressWarnings(suppressMessages(library(clusterProfiler)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))

input_path        <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Gene_Lists/cellbender_processed_data/original_018/"
GO_term_path      <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Functional_Analysis/Functional_Analysis_018/"

output_path       <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Functional_Analysis/Functional_Analysis_018/Excel/"

cell_types        <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

study_results     <- list.files(input_path)

#0.050.05human_up_regulated_GO.rds
for (FDR1 in c(0.05, 0.1, 0.2)){
	for (FDR2 in c(0.05, 0.1, 0.2)){
		
		print(paste0(GO_term_path, FDR1, FDR2, "chimp_down_regulated_GO.rds"))
		chimp_down_regulated_GO      <- readRDS(paste0(GO_term_path, FDR1, FDR2, "chimp_down_regulated_GO.rds"))
		human_down_regulated_GO      <- readRDS(paste0(GO_term_path, FDR1, FDR2, "human_down_regulated_GO.rds"))   
		chimp_up_regulated_GO        <- readRDS(paste0(GO_term_path, FDR1, FDR2, "chimp_up_regulated_GO.rds"))
		human_up_regulated_GO        <- readRDS(paste0(GO_term_path, FDR1, FDR2, "human_up_regulated_GO.rds"))

#		print(chimp_down_regulated_GO)
		for (cell_type in cell_types){
			#print(ont)
			chimp_down_summary  <- DataFrame()
                        human_down_summary  <- DataFrame()
                        chimp_up_summary    <- DataFrame()
                        human_up_summary    <- DataFrame()

			for (study in study_results[c(1, 4)]){
				print(study)
				print(colnames(chimp_up_summary))

				#print(chimp_down_regulated_GO[[ont]][[study]])
				for (ont in c("MF", "BP", "CC")){
					print(cell_type)
					#print(head(chimp_down_regulated_GO[[ont]][[study]][[cell_type]]@result, 10))					
					
					chimp_down <- head(chimp_down_regulated_GO[[ont]][[study]][[cell_type]]@result, 10)
					human_down <- head(human_down_regulated_GO[[ont]][[study]][[cell_type]]@result, 10)
					chimp_up   <- head(chimp_up_regulated_GO[[ont]][[study]][[cell_type]]@result, 10)
                                        human_up   <- head(human_up_regulated_GO[[ont]][[study]][[cell_type]]@result, 10) 
					
					chimp_down$Ont  <- ont
					human_down$Ont  <- ont
					chimp_up$Ont    <- ont
					human_up$Ont    <- ont

					chimp_down$Study  <- study
                                        human_down$Study  <- study
                                        chimp_up$Study    <- study
                                        human_up$Study    <- study

					chimp_down_summary  <- rbind(chimp_down_summary, chimp_down)
                                        human_down_summary  <- rbind(human_down_summary, human_down)
                                        chimp_up_summary    <- rbind(chimp_up_summary,   chimp_up)
                                        human_up_summary    <- rbind(human_up_summary,   human_up)

					#print(colnames(chimp_up_summary))

				}

			}
			write_xlsx(data.frame(chimp_down_summary), paste0(output_path, FDR1, FDR2, cell_type, "chimp_down_regulated_GO.xlsx"))
                        write_xlsx(data.frame(human_down_summary), paste0(output_path, FDR1, FDR2, cell_type, "human_down_regulated_GO.xlsx"))
                        write_xlsx(data.frame(chimp_up_summary),   paste0(output_path, FDR1, FDR2, cell_type, "chimp_up_regulated_GO.xlsx"))
                        write_xlsx(data.frame(human_up_summary),   paste0(output_path, FDR1, FDR2, cell_type, "human_up_regulated_GO.xlsx"))
		}
	}
}
