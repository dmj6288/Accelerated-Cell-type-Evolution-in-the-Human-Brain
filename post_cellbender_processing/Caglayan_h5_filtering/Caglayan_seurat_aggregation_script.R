suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))

seurat_folder_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_technical/RN011/Caglayan/"
seurat_files       <- list.files(seurat_folder_path)

#print(seurat_files)
samples        <- unique(tstrsplit(seurat_files, ".rds")[[1]])

output_path    <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_species/RN011/Caglayan/"

human_samples   <- samples[grepl("h", samples, fixed = TRUE)]
chimp_samples   <- samples[grepl("c", samples, fixed = TRUE)]
macaque_samples <- samples[grepl("m", samples, fixed = TRUE)]

all_samples     <- list(human_samples, 
			chimp_samples,
			macaque_samples)


for (samples in all_samples){

	species_seurat_object_list <- list()
#	print(samples)
	for (sample in samples){
		
		#print(sample)
		#for (technical_replicate in c(1:4)){
		print(paste0(seurat_folder_path, sample, ".rds"))
		current_seurat_object <- readRDS(paste0(seurat_folder_path, sample, ".rds"))
		print(str(current_seurat_object))
		species_seurat_object_list[[sample]] <- current_seurat_object
#		print(str(species_seurat_object_list[[sample]]))		
	}
	species_seurat <- merge(species_seurat_object_list[[samples[1]]], y = c(species_seurat_object_list[[samples[2]]],
                                                                                species_seurat_object_list[[samples[3]]],
                                                                                species_seurat_object_list[[samples[4]]]), add.cell.ids = samples, project = sample)
	
	print(str(species_seurat))
	saveRDS(species_seurat, paste0(output_path, samples[1], ".rds"))
}
