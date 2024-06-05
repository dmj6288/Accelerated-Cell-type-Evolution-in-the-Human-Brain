rm(list=ls())

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(readxl)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressMessages(suppressWarnings(library(rstatix)))
suppressMessages(suppressWarnings(library(ggprism)))
suppressMessages(suppressWarnings(library(ggVennDiagram)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(VennDiagram)))

  
# USE PACKAGE - chromoMAP

library(parallel)
# Example using mclapply to parallelize operations across cores

studies            <- c("Caglayan", "Ma")
species            <- c("HUMAN", "CHIMP")
regions            <- c("PCC", "dlPFC")
names(regions)     <- studies
regulation         <- c("UR", "DR")
folder_name        <- "../Gene_List_results/"
FDR                <- 0.05

cell_names_of_interest_Excite  <- c("L2-3 IT", "L3-5 IT-1", "L3-5 IT-2", "L3-5 IT-3")

# "L2-3 IT L3-5 IT-1 L3-5 IT-2 L3-5 IT-3"
cell_names_of_interest_Inhibit <- c("PVALB", "SST", "VIP")

cell_names_of_interest         <- list(cell_names_of_interest_Excite,
                                       cell_names_of_interest_Inhibit)

names(cell_names_of_interest)  <- c("Excite", "Inhibit")

subtype_folder_name   <- "../../Cloud Data/subtype/"
majortype_folder_name <- "../../Cloud Data/major/"
