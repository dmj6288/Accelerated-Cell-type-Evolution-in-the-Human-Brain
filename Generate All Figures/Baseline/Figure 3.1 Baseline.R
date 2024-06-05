rm(list=ls())
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(rstatix)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggprism)))


regulation         <- c("UR", "DR")
studies            <- c("Caglayan", "Ma")
region_list        <- list("PCC", "dlPFC")
names(region_list) <- studies
species_type       <- c("HUMAN", "CHIMP")

folder_name        <- "../../Cloud Data/major/"
FDR                <- 0.05
considered_list    <- list()
for (study in studies){
  
  considered_list[[study]] <- readRDS(paste0(folder_name, study, "/", study, FDR, "_three_species_considered_Genes.rds"))
}

mini_cell_type_list <- c( "Astro",     "Excite",    "Inhibit",   "Microglia", "Oligo",     "OPC")

cell_plot_order <- c("Excite PCC",  "Inhibit PCC", "Astro PCC", "Microglia PCC", "Oligo PCC", "OPC PCC",
                     "Excite dlPFC", "Inhibit dlPFC", "Astro dlPFC", "Microglia dlPFC","Oligo dlPFC","OPC dlPFC")

barplot_data <- readRDS("../../Visualization_data/Figure 3/data_stats_results_all_considered_genes_newdNdS.rds")
stats_data   <- readRDS("../../Visualization_data/Figure 3/stats_data_pvalues_all_considered_genes_newdNdS.rds")

