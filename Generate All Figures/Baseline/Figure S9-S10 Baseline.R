rm(list = ls())

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(ggpubr)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(rstatix)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))

studies             <- c("Caglayan", "Ma")
cellnames           <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")
cellfullnames       <- c("Astrocyte", "Excitatory Neuron", "Inhibitory Neuron", "Microglia", "Oligodendrocyte", "OPC")

names(cellfullnames)    <- cellnames

checkpoints_path    <- "../../Visualization_data/Figure S9-S10/"
figure_output_path  <- "../../Figures/Figure S9-S10/"


FDR <- 5
gene_count_value <- 10

regulation.labs <- c("chimp", "human combined up/down regulated genes")
names(regulation.labs) <- c("chimp", "human")

FDR_value <- 0.05