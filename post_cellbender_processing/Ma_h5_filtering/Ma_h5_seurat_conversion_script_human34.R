suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

h5_folder_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Ma/"

h5_files       <- list.files(h5_folder_path)

samples        <- unique(tstrsplit(h5_files, "_")[[1]])

output_path    <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN003/Ma/"

source("/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

human_samples   <- samples[grepl("HSB", samples, fixed = TRUE)]
chimp_samples   <- samples[grepl("PTB", samples, fixed = TRUE)]
macaque_samples <- samples[grepl("RMB", samples, fixed = TRUE)]

for (sample in human_samples){

        for (technical_replicate in c(3:4)){
                print(paste0(h5_folder_path, sample, "_", technical_replicate, ".h5"))
                current_h5_file <- CreateSeuratObject(Read10X_h5(paste0(h5_folder_path, sample, "_", technical_replicate, ".h5")))
                print(str(current_h5_file))

                filtered_seurat <- comprehensive_filtering(current_h5_file, paste0("Ma_dlPFC", sample, "_", technical_replicate))
                saveRDS(filtered_seurat, paste0(output_path, sample, "_", technical_replicate, ".rds"))

        }
}
