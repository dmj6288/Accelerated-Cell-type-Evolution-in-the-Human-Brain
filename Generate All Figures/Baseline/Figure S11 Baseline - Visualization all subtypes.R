# Initialize helper variables
rm(list = ls())
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(readxl)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(writexl)))
suppressMessages(suppressWarnings(library(rstatix)))
suppressMessages(suppressWarnings(library(ggprism)))
suppressMessages(suppressWarnings(library(ggVennDiagram)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(VennDiagram)))
suppressMessages(suppressWarnings(library(ggrepel)))

studies             <- c("Caglayan", "Ma")
cellnames           <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")
cellfullnames       <- c("Astrocyte", "Excitatory Neuron", "Inhibitory Neuron", "Microglia", "Oligodendrocyte", "OPC")
cellfullnames_factor   <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte",  "Microglia", "Oligodendrocyte", "OPC")

names(cellfullnames)    <- cellnames

species            <- c("HUMAN", "CHIMP")
regions            <- c("PCC", "dlPFC")
names(regions)     <- studies
regulation         <- c("UR", "DR")
folder_name        <- "../Gene_List_results/"
FDR                <- 0.05

combination_lengths    <- list(c(2:4), c(2:3))
cellnames_of_interest  <- c("Excite", "Inhibit")

names(combination_lengths) <- cellnames_of_interest

study <- "Ma"


