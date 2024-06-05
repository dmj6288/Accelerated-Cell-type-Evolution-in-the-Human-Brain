suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(stringi)))

seurat_folder_path           <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN015/Jorstad/"
output_path                  <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN015/"

all_files                    <- list.files(seurat_folder_path)

Human_samples                <- all_files[stri_detect_fixed(all_files, "Jorstad_H")]

species_seurat_object_list <- list()

for (sample in Human_samples){
        species_seurat_object_list[[sample]] <- readRDS(paste0(seurat_folder_path, sample))
}

species_seurat <- merge(species_seurat_object_list[[1]], y = c(species_seurat_object_list[[2]],
                                                               species_seurat_object_list[[3]],
							       species_seurat_object_list[[4]],
                                                               species_seurat_object_list[[5]],
                                                               species_seurat_object_list[[6]],
                                                               species_seurat_object_list[[7]]), project = "Jorstad_Human")

saveRDS(species_seurat, paste0(output_path, "Jorstad_Human.rds"))
