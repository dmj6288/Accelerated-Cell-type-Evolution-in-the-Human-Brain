rm(list=ls())
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(rstatix)))
suppressMessages(suppressWarnings(library(rlang)))


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

