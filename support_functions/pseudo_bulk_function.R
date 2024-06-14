# Single-cell RNA-seq analysis - Pseudobulk DE analysis with DESeq2

#suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(S4Vectors)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(apeglm)))
suppressMessages(suppressWarnings(library(png)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(comprehenr)))


pseudo_bulk_sce_data <- function(i){
	support_functions_path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/support_functions/"

	# INPUT

	sce_path               <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/sces/sce_folder_019/"

	# OUTPUT

	pseudo_bulk_ready_data <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/pbds/pseudo_bulk_ready_data_019/"
	coldata_path           <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/coldata/coldata_019/"
	aggmatrix_path         <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/Aggregate_matrices/aggregate_019/"

	source(paste(support_functions_path, "aggregate_function_Matrix.Utils.R", sep = ""))

	#source(paste(support_functions_path, "cell_type_replicability_check_Emrev2.R", sep = ""))

	sce_folders_to_analyze   <- list.files(sce_path)

	for (sce_folder in sce_folders_to_analyze[c(1, 5)]){

		sce_files_to_analyze <- list.files(paste0(sce_path, sce_folder))
		print(sce_files_to_analyze)

		for (sce_file in sce_files_to_analyze){

			print(sce_file)
			sce <- readRDS(paste(sce_path, sce_folder, "/", sce_file, sep =""))
			saveRDS(colData(sce), paste0(coldata_path, "/", sce_folder, "/", sce_file))
			print("----READ COMPLETE----")
	
			groups <- colData(sce)[, c("cluster_id", "sample_id")]
			print(groups)
	
			groups$both <- as.factor(paste(groups$cluster_id, groups$sample_id, sep ="_"))
	
			# ----------- BEFORE REPLICABILITY CHECK ----------- #

			print("GROUP LEVELS, BOTH, CLUSTER, SAMPLE")
			print((levels(groups$both)))
			print((levels(groups$cluster_id)))	
			print((levels(groups$sample_id)))	

		        cluster_names <- levels(colData(sce)$cluster_id)
			print(cluster_names)
	        	#cluster_names <- chartr(" ", "_", cluster_names)

	        	# Extract unique names of samples (= levels of sample_id factor variable)
		        sample_names <- levels(colData(sce)$sample_id)
		#        print(sample_names)

			groups <- colData(sce)[, c("cluster_id", "sample_id")]
		#        print(groups$cluster_id)	

			# Aggregate across cluster-sample groups
			# transposing row/columns to have cell_ids as row names matching those of groups
			aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")

			aggr_counts <- t(aggr_counts)


			print("#--- AGGR COUNT HEAD ---#")
			print(head(aggr_counts))

		#	print(colnames(aggr_counts))
			counts_ls <- list()
	
		#	cluster_combinations <- readRDS(paste0("rds_files/",clust_combo_hash[[sce_file]]))
			cells   <- unique(tstrsplit(levels(groups$cluster_id), "_#")[[1]])
			species <- unique(tstrsplit(levels(groups$cluster_id), "_#")[[2]])

			cluster_combinations <- rbind(paste(cells, species[1], sep = "_#"),
						      paste(cells, species[2], sep = "_#"))

			print("CLUST_COMB")
			print(cluster_combinations)	

			for (i in 1:dim(cluster_combinations)[2]) {
  
				print(i)
				## Extract indexes of columns in the global matrix that match a given cluster
				print(cluster_combinations[1,i])
				print(cluster_combinations[2,i])
				column_id1 <- which(substr(colnames(aggr_counts), start = 1, stop = nchar(cluster_combinations[1,i])) == cluster_combinations[1,i])
				column_id2 <- which(substr(colnames(aggr_counts), start = 1, stop = nchar(cluster_combinations[2,i])) == cluster_combinations[2,i])
		                print(column_id1)
				print(column_id2)

				## Store corresponding sub-matrix as one element of a list

				print("# ------------ STORING AGGCOUNTS ------------ #")

				#filter         <- rowSums(aggr_counts >= 10) >= 6
				counts_ls[[i]] <- aggr_counts[, c(column_id1, column_id2)]

				#filter         <- rowSums(counts_ls[[i]] >= 10) >= 6

				#counts_ls[[i]] <- counts_ls[[i]][filter, ]
			
				names(counts_ls)[i] <- paste(cluster_combinations[1, i], cluster_combinations[2, i], sep ="-vs-")
		#		print(names(counts_ls)[i])
		#		print(head(counts_ls[[i]]))

			}
			counts_ls <- counts_ls[lapply(counts_ls, length) > 0]

			# Reminder: explore structure of metadata
			head(colData(sce))

			# Extract sample-level variables
			# metadata <- colData(sce) %>% as.data.frame() %>% dplyr::select(cluster_id, sample_id, lib_batch, human_age, sex)
			metadata <- colData(sce) %>% as.data.frame() %>% dplyr::select(cluster_id, sample_id, human_age, sex)
			# print("#---metadata---#")
			# print(metadata)
		
	
			metadata <- metadata[!duplicated(metadata), ]

			# Number of cells per sample and cluster
			t <- table(colData(sce)$sample_id, colData(sce)$cluster_id)
			print("#---- T ----#")
			print(t)
			# Creating metadata list

			## Initiate empty list
			metadata_ls <- list()
			for (i in 1:length(counts_ls)) {
  
				## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  
				print(i)
				print(colnames(counts_ls[[i]]))
				df <- data.frame(cluster_sample_id = sort(colnames(counts_ls[[i]])))
				print(df)		
	
				## Use tstrsplit() to separate cluster (cell type) and sample IDs
				df$cluster_id <- tstrsplit(df$cluster_sample_id, "_#")[[1]]
				df$sample_id  <- tstrsplit(df$cluster_sample_id, "_#")[[2]]

				cell_species <- unique(paste(tstrsplit(df$cluster_sample_id, "_#")[[1]],
		                                             tstrsplit(tstrsplit(df$cluster_sample_id, "_#")[[2]], "_")[[1]], sep = "_#"))
				print(cell_species)
				covariate_df  <- unique(metadata[metadata$cluster_id %in% cell_species, ])
				covariate_df  <- covariate_df[order(covariate_df$cluster_id),]
				print(covariate_df)

				#df$lib_batch  <- covariate_df$lib_batch
				df$human_age  <- scale(covariate_df$human_age, center = TRUE, scale = TRUE)[, 1]
				# covariate_df$human_age
				df$sex        <- covariate_df$sex
#				df$Segment    <- covariate_df$Segment
			
				species       <- unique(tstrsplit(df$sample_id, "_")[[1]])
				print(species)
				num1          <- length(df$sample_id[tstrsplit(df$sample_id, "_")[[1]] %in% species[1]])
				num2          <- length(df$sample_id[tstrsplit(df$sample_id, "_")[[1]] %in% species[2]])
				print(species[1] %in% tstrsplit(df$sample_id, "_")[[1]])
				print(df)
				df$compare <- factor(c(to_vec(for(i in 1:num1) "B"), to_vec(for(i in 1:num2) "A")))
				print(df)
  
				print("#------------ COMPARE COMPLETE ------------#")
				## Retrieve cell count information for this cluster from global cell count table
				idx <- pmatch(unique(substr(df$cluster_sample_id, start = 1, stop = nchar(df$cluster_sample_id) - 3)), colnames(t))
				print("# ---- IDX ---- #")
				print(idx)

				cell_counts <- t[, idx]

				print(cell_counts)

				df <- plyr::join(df, metadata, by = intersect(names(df), names(metadata)))
  
				print(df)
				## Update rownames of metadata to match colnames of count matrix, as needed later for DE
				rownames(df) <- df$cluster_sample_id
  
				## Store complete metadata for cluster i in list
				metadata_ls[[i]] <- df
				names(metadata_ls)[i] <- paste(cluster_combinations[1, i], cluster_combinations[2, i], sep ="-vs-")
			}
		
			saveRDS(counts_ls,   paste(pseudo_bulk_ready_data, sce_folder, "/", sce_file, "counts_ls.rds" ,sep = ""))
			saveRDS(metadata_ls, paste(pseudo_bulk_ready_data, sce_folder, "/", sce_file, "metadata_ls.rds",sep = ""))
			print("# ---- SAVED ----#")
		}
	}
}
