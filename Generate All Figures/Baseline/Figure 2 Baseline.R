
rm(list = ls())

suppressMessages(suppressWarnings(library(ggplot2)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressMessages(suppressWarnings(library(ggprism)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(openxlsx)))

studies             <- c("Caglayan", "Ma")
cellnames           <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")
cellfullnames       <- c("Astrocyte", "Excitatory Neuron", "Inhibitory Neuron", "Microglia", "Oligodendrocyte", "OPC")
cellfullnames_factor   <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte",  "Microglia", "Oligodendrocyte", "OPC")

names(cellfullnames)    <- cellnames

species            <- c("human", "chimp")
regions            <- c("PCC", "dlPFC")
names(regions)     <- studies
regulation         <- c("UR", "DR")
folder_name        <- "../../Cloud Data/major/"
figure_path        <- "../../Figures/Figure 2/"
FDR                <- 0.05

colour_list <- list(list("#a3ebe7", "#54c5e0"), list("#eb934b", "#d95734"))
names(colour_list) <- species