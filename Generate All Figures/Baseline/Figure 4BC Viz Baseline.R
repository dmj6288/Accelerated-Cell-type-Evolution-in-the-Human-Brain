rm(list = ls())

suppressMessages(suppressWarnings(library(ComplexHeatmap)))
suppressMessages(suppressWarnings(library(ggplot2)))
library(ggplotify) 

library(UpSetR)

powerset <- function(x) {
  sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = F))
  unlist(sets, recursive = F)
}