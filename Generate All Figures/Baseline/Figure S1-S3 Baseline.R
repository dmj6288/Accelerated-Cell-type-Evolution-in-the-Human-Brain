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
suppressWarnings(suppressMessages(library(writexl)))
suppressWarnings(suppressMessages(library(xlsx)))
suppressWarnings(suppressMessages(library(ggpattern)))


FDR <- 0.05

folder_name                    <- "../../../Visualization - cellbender processed data - 018/Gene_List_results/"
cell_maps                      <- list.files("../../../Visualization - cellbender processed data - 018/cellname_maps/")
cell_counts_path               <- "../../Visualization_data/Figure S1-S3/checkpoints/"
gene_table_path                <- "../../Visualization_data/Figure S1-S3/"
cell_counts_files              <- list.files("../../Visualization_data/Figure S1-S3/checkpoints/")
cell_full_names_hash           <- readRDS("../../../Visualization - cellbender processed data - 019/rds_files/cell_full_names.rds")
figure_output_path             <- "../../Figures/Figure S1-S3/"

cell_types <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")

study_results                        <- list.files(folder_name)[c(2, 3, 4)]
study_to_gene_number_hash            <- hash(study_results, c(13088, 11141, 9420))
name_to_region_hash                  <- hash(study_results, c("MTG", "dlPFC", "ACC"))

study_to_human_replicate_number_hash <- hash(study_results, c(7, 3, 3))
study_to_chimp_replicate_number_hash <- hash(study_results, c(7, 2, 3))


geom_bar_width = 0.6
geom_bar_alpha = 0.8

