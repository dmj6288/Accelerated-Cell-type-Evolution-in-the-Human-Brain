rm(list = ls())
suppressMessages(suppressWarnings(library(stringr)))
print("STRINGR IMPORTED")
suppressMessages(suppressWarnings(library(Seurat)))
#suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(hash)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

support_functions <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

metadata_path                <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/mdata/mdata_011/"
seurat_objects_path          <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_final/RN018/"
sce_path                     <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/sces/sce_folder_021/"

metadata_to_seurat_hash      <- hash(list.files(metadata_path), 
                                     list.files(seurat_objects_path))

print(metadata_to_seurat_hash)

filtering <- 1

for (study_metadata in c(list.files(metadata_path))[c(5)]){
	
	filtered_seurat_object <- readRDS(paste0(seurat_objects_path, metadata_to_seurat_hash[[study_metadata]]))
	
	print(str(filtered_seurat_object))
	
	for (metadata_file in list.files(paste0(metadata_path, study_metadata))){

		print(list.files(paste0(metadata_path, study_metadata)))
		print(paste(metadata_path, study_metadata, "/", metadata_file, sep =""))
		df                           <- readRDS(paste(metadata_path, study_metadata, "/", metadata_file, sep =""))

		filtered_seurat_object_with_identified_cells <- subset(filtered_seurat_object, cells = rownames(df))
		df                           <- df[colnames(filtered_seurat_object_with_identified_cells@assays$RNA@counts), ]

		print(identical(colnames(filtered_seurat_object_with_identified_cells), rownames(df)))
		
		filtered_seurat_object_with_identified_cells@meta.data$orig.ident <- df$sample_id
		filtered_seurat_object_with_identified_cells@meta.data$human_age  <- df$human_age
		filtered_seurat_object_with_identified_cells@meta.data$sex        <- df$sex

		filtered_seurat_object_with_identified_cells <- SetIdent(filtered_seurat_object_with_identified_cells, value = df$trad)
	
		# print(str(filtered_seurat_object_with_identified_cells))

		# FILTERING THE DATA
		
		counts                       <- filtered_seurat_object_with_identified_cells@assays$RNA@counts
		metadata                     <- filtered_seurat_object_with_identified_cells@meta.data

	
		# print("------FILTERED SEURAT STRUCTURE------")
		# print(str(filtered_seurat_object_with_identified_cells))


		metadata$cluster_id          <- df[colnames(counts),"trad"]
		metadata$sample_id           <- df[colnames(counts),"sample_id"]
	
		sce                          <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
	
		# print("------SCE STRUCTURE------")
		# print(str(sce))
		
		print(paste(sce_path, study_metadata, "_", metadata_file, "_sce.rds", sep =""))
		saveRDS(sce, paste(sce_path, study_metadata, "/", study_metadata, "_", metadata_file, "_sce.rds", sep =""))
	}
}
