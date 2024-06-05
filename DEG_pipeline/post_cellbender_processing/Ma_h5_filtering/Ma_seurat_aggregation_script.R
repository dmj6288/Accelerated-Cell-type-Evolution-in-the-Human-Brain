suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

seurat_folder_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN003/Ma/"
seurat_files       <- list.files(seurat_folder_path)

print(seurat_files)
samples        <- unique(tstrsplit(seurat_files, "_")[[1]])

output_path    <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_biological/RN003/Ma/"

human_samples   <- samples[grepl("HSB", samples, fixed = TRUE)]
chimp_samples   <- samples[grepl("PTB", samples, fixed = TRUE)]
macaque_samples <- samples[grepl("RMB", samples, fixed = TRUE)]

all_samples     <- list(human_samples, 
			chimp_samples,
			macaque_samples)


for (samples in all_samples){

	species_seurat_object_list <- list()
	for (sample in samples){
		
		print(sample)
		seurat_object_list <- list()
		for (technical_replicate in c(1:4)){

			current_seurat_file <- readRDS(paste0(seurat_folder_path, sample, "_", technical_replicate, ".rds"))
			#print(str(current_seurat_file))
			seurat_object_list[[technical_replicate]] <- current_seurat_file
		}
		
		biological_replicate <- merge(seurat_object_list[[1]], y = c(seurat_object_list[[2]],
									     seurat_object_list[[3]],
                                                                             seurat_object_list[[4]]), add.cell.ids = c("1", "2", "3", "4"), project = sample)
		print(str(biological_replicate))

		species_seurat_object_list[[sample]] <- biological_replicate
		#saveRDS(biological_replicate, paste0(output_path, sample, ".rds"))
	}
	species_seurat <- merge(species_seurat_object_list[[1]], y = c(species_seurat_object_list[[2]],
                                                                       species_seurat_object_list[[3]],
                                                                       species_seurat_object_list[[4]]), add.cell.ids = samples, project = sample)
	
	print(str(species_seurat))
	saveRDS(species_seurat, paste0(output_path, samples[1], ".rds"))
}
