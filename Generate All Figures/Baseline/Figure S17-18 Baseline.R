rm(list=ls())

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(readxl)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(grid)))
suppressMessages(suppressWarnings(library(stringr)))
suppressWarnings(suppressMessages(library(dplyr)))

geom_bar_width = 0.6
geom_bar_alpha_prop    = 0.8
geom_bar_alpha_overlap = 0.2

studies             <- c("Caglayan", "Ma")
cellnames           <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")
cellfullnames       <- c("Astrocyte", "Excitatory Neuron", "Inhibitory Neuron", "Microglia", "Oligodendrocyte", "OPC")

                       
regions             <- c("PCC", "dlPFC")

names(regions)      <- studies
names(cellfullnames)    <- cellnames

major_folder_name         <- "../../Cloud Data/major/"
subtype_folder_name       <- "../../Cloud Data/subtype/"

