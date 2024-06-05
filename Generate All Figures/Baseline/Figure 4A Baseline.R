suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(readxl)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(grid)))
suppressMessages(suppressWarnings(library(stringr)))
suppressWarnings(suppressMessages(library(dplyr)))



studies             <- c("Caglayan", "Ma")
cellnames           <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")
cellfullnames       <- c("Astrocyte", "Excitatory Neuron", "Inhibitory Neuron", "Microglia", "Oligodendrocyte", "OPC")

cells_of_interest_excite    <- c("L2-3 IT", "L3-5 IT-1", "L3-5 IT-2", "L3-5 IT-3")
cells_of_interest_inhibit   <- c("PVALB", "SST", "VIP")

cells_of_interest_list      <- list(cells_of_interest_inhibit, cells_of_interest_excite)
names(cells_of_interest_list) <- c("Inhibit", "Excite")

regions             <- c("PCC", "dlPFC")

names(regions)      <- studies
names(cellfullnames)    <- cellnames

geom_bar_width = 0.6
geom_bar_alpha_prop    = 0.8
geom_bar_alpha_overlap = 0.2
geom_errorbar_width = 0.15
geom_errorbar_dodge = 0.6
save_flag <- 0