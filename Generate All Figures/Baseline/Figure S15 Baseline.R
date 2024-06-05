rm(list=ls())
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(hash)))
suppressMessages(suppressWarnings(library(biomaRt)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(tidyr)))

regulation         <- c("UR", "DR")
studies            <- c("Caglayan", "Ma")
region_list        <- list("PCC", "dlPFC")
names(region_list) <- studies
species_type       <- c("HUMAN", "CHIMP")

folder_name        <- "../../Cloud Data/"
FDR                <- 0.05
considered_list    <- list()
for (study in studies){
  
  considered_list[[study]] <- readRDS(paste0(folder_name, study, "/", study, FDR, "_three_species_considered_Genes.rds"))
}


cell_plot_order <- c("Excite PCC",  "Inhibit PCC", "Astro PCC", "Microglia PCC", "Oligo PCC", "OPC PCC",
                     "Excite dlPFC", "Inhibit dlPFC", "Astro dlPFC", "Microglia dlPFC","Oligo dlPFC","OPC dlPFC")
