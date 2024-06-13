rm(list = ls())

suppressWarnings(suppressMessages(library(VennDiagram)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggplotify)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggVennDiagram)))
suppressWarnings(suppressMessages(library(gmp)))
suppressWarnings(suppressMessages(library(hash)))
suppressWarnings(suppressMessages(library(gridExtra)))

enrich_pvalue <- function(N, A, B, k)
{
  m <- A + k
  n <- B + k
  i <- k:min(m,n)
  
  as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}    

folder_name   <- "../../Cloud Data/major/"
FDR           <- 0.05

study_results           <- list.files(folder_name)[c(1, 6)]

traditional_cell_types   <- c("Astro", "Excite", "Inhibit", "Microglia", "Oligo", "OPC")
cell_full_names_hash     <- readRDS("../../common_support_files/cell_full_names.rds")
checkpoints_path         <- "../../Visualization_data/Figure S2-S5/checkpoints/"
figure_output_path       <- "../../Figures/Figure S2-5/"
supplementary_table_path <- "../../Supplementary Tables/"

study_to_gene_number_hash <- hash(study_results, c(13791, 13316))
study_to_human_replicate_number_hash <- hash(study_results, c(4, 4))
study_to_chimp_replicate_number_hash <- hash(study_results, c(4, 4))

name_to_region_hash <- hash(study_results, c("PCC", "dlPFC"))



