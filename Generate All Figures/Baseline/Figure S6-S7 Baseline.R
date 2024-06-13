rm(list = ls())

suppressMessages(suppressWarnings(library(ComplexHeatmap)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggplotify)))
suppressMessages(suppressWarnings(library(UpSetR)))

powerset <- function(x) {
  sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = F))
  unlist(sets, recursive = F)
}

folder_name     <- "../../Cloud Data/"
FDR             <- 0.05
study_results   <- list.files(folder_name)

cell_full_names <- c("Astrocyte", "Excitatory Neuron", "Inhibitory Neuron",  "Microglia", "Oligodendrocyte", "OPC")

checkpoints_path         <- "../../Visualization_data/Figure S8-S9/checkpoints/"

cell_order <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Oligodendrocyte", "Microglia", "OPC")

