rm(list = ls())

suppressWarnings(suppressMessages(library(readr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(grid)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(xlsx)))

folder_name                    <- "../../Cloud Data/major/"
study_results                  <- list.files(folder_name)[c(1, 6)]

FDR                 <-  0.05
studies             <- c("Caglayan", "Ma")
cell_types          <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

cellfullnames       <- c("Astrocyte", "Excitatory Neuron", "Inhibitory Neuron", "Microglia", "Oligodendrocyte", "OPC")

regions             <- c("PCC", "dlPFC")

names(regions)      <- studies


bootstrap_results_path <- "../../Cloud Data/major/consolidated_lists_015BT/"
checkpoints_path       <- "../../Visualization_data/Figure 1/checkpoints/"
gene_table_path        <- "../../Visualization_data/Figure 1/"
figure_output_path     <- "../../Figures/Figure 1/"

cell_full_names_hash   <- readRDS("../../common_support_files/cell_full_names.rds")


geom_bar_width = 0.6
geom_bar_alpha = 0.8
geom_errorbar_width = 0.15
geom_errorbar_dodge = 0.6

