suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

h5_folder_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/cellbender_finished_RAW_data/Caglayan/"

h5_files       <- list.files(h5_folder_path)

output_path    <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN011/Caglayan/"

source("~/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/comprehensive_filtering_003.R")

for (technical_replicate in c("m1", "m2")){

	print(paste0(h5_folder_path, technical_replicate, ".h5"))
        current_h5_file <- CreateSeuratObject(Read10X_h5(paste0(h5_folder_path, technical_replicate, ".h5")))
        print(str(current_h5_file))

        filtered_seurat <- comprehensive_filtering(current_h5_file, paste0("Caglayan_PCC", technical_replicate))
        saveRDS(filtered_seurat, paste0(output_path, technical_replicate, ".rds"))

}

