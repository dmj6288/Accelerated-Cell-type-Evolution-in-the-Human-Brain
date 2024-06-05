rm(list = ls())
suppressMessages(suppressWarnings(library(stringr)))
print("STRINGR IMPORTED")
suppressMessages(suppressWarnings(library(Seurat)))
#suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(hash)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(writexl)))

support_functions <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

metadata_path                <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/mdata/mdata_010/"
seurat_objects_path          <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/filtered_seurat/filtered_seurat_objects_final/RN018/"
sce_path                     <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/sces/sce_folder_019/"
sample_wise_cell_count_path  <- "~/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Evolution_Paper/how_any_cells/cell_count_by_sample/"

metadata_to_seurat_hash      <- hash(list.files(metadata_path), 
                                     list.files(seurat_objects_path))

print(metadata_to_seurat_hash)

filtering <- 1

for (study_metadata in c(list.files(metadata_path))[c(1)]){
	
	filtered_seurat_object <- readRDS(paste0(seurat_objects_path, metadata_to_seurat_hash[[study_metadata]]))
	
	print(str(filtered_seurat_object))
	study_wise_cell_counts_per_sample <- c()
	
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

		filtered_seurat_orig.ident <- filtered_seurat_object_with_identified_cells@meta.data$orig.ident
		cell_type_by_sample <- table(filtered_seurat_object_with_identified_cells@meta.data$orig.ident, filtered_seurat_object_with_identified_cells@active.ident)

		df_counts <- as.data.frame.matrix(cell_type_by_sample)		
		df_counts$row_names <- rownames(df_counts)
		
		cell_types <- unique(gsub("_#.*", "", names(df_counts)))
		df_aggregated <- data.frame(row = 1:nrow(df_counts))		

		for(cell_type in cell_types) {
			# Identify columns belonging to the current cell type
			cols <- grep(paste0("^", cell_type, "_#"), names(df_counts), value = TRUE)

			# Combine these columns by row-wise sum (or other appropriate aggregation)
			combined_column <- rowSums(df_counts[, cols], na.rm = TRUE)

			# Add the combined column to the aggregated data frame
			df_aggregated[[cell_type]] <- combined_column
			print(cell_type)
			print(combined_column)
		}

		rownames(df_aggregated) <- df_counts$row_names
		df_aggregated$Row_names <- df_counts$row_names

		study_wise_cell_counts_per_sample <- rbind(study_wise_cell_counts_per_sample, df_aggregated)
		
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
	
	#	print(paste(sce_path, study_metadata, "_", metadata_file, "_sce.rds", sep =""))
	#	saveRDS(sce, paste(sce_path, study_metadata, "/", study_metadata, "_", metadata_file, "_sce.rds", sep =""))
	}

	study_wise_cell_counts_per_sample <- study_wise_cell_counts_per_sample[!duplicated(study_wise_cell_counts_per_sample$Row_names), ]
	print(study_wise_cell_counts_per_sample)

	write_xlsx(study_wise_cell_counts_per_sample, paste0(sample_wise_cell_count_path, study_metadata, "sub_celltype_by_sample.xlsx"))
}
