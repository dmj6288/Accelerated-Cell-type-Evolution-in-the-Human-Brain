get_dds_resultsAvsB <- function(padj_cutoff = 0.05, study, counts_ls, metadata_ls, count_ls, path, bootstrapped_gene_selection_matrix_list, bootstrap_Step, step_size) {

        results_path <- paste(path, study, "/", count_ls, "/", sep = "")
	print(results_path)
	print(paste0(results_path, "bootstrap_results_", metadata_ls[1], ".rds"))

	if (!dir.exists(results_path)) { dir.create(results_path) }
	#print(bootstrap_Step)
	#print(step_size)
	print(c((bootstrap_Step*step_size):((bootstrap_Step+1)*step_size))[1:10])
	bootstrap_results <- c()
				
	for (i in c(1:length(metadata_ls))){
		clustx <- names(counts_ls)[i]
		print(clustx) # useful for debugging
		# print(head(bootstrapped_gene_selection_matrix_list[[i]]))
		# Extract counts matrix and metadata for cluster x
		idx <- which(names(counts_ls) == clustx)
    
		print(idx)
		    
		clust_mat <- counts_ls[[idx]]

		#cluster_counts <- clust_mat[1:max(which(rowSums(clust_mat != 0) > 0)), ]

		#print(any(is.na(cluster_counts)))

		#cluster_counts <- clust_mat[which(rowSums(clust_mat) != 0),]
	    
		#cluster_counts <- as.matrix(clust_df)
		cluster_metadata <- metadata_ls[[idx]]

#		if ( all(colnames(cluster_counts) != rownames(cluster_metadata))) {
#			print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
#		}
		
		for (bootstrap_sample in c(((bootstrap_Step*step_size)+1):((bootstrap_Step+1)*step_size))){
			#print(bootstrap_sample)
			#cluster_counts <- clust_mat[which(rowSums(clust_mat) != 0),]
			cluster_counts_strapped <- clust_mat[bootstrapped_gene_selection_matrix_list[bootstrap_sample, ], ]
			cluster_counts_strapped <- cluster_counts_strapped[which(rowSums(cluster_counts_strapped) != 0),]
			#print(head(cluster_counts_strapped))
			dds <- DESeqDataSetFromMatrix(cluster_counts_strapped, 
		                                      colData = cluster_metadata, 
                                                      design = ~human_age + sex + compare)
			
			#print("BOOTSTRAP SAMPLE : ")
			#print(bootstrap_sample)

			print("# -------- DESeqDataSetFromMatrix COMPLETE -------- #")
			# Transform counts for data visualization
			rld <- rlog(dds, blind = TRUE)
    
    
			# Generate QC plots
    
			rld_mat <- assay(rld)
			rld_cor <- cor(rld_mat)
			
			keep <- rowSums(counts(dds) >= 10) >= 4
	                dds <- dds[keep,]
    
			## Plot and save heatmap
    
			# Run DESeq2 differential expression analysis
			dds <- DESeq(dds)
			#print("# -------- DESEQ COMPLETE -------- #")
                ## Plot dispersion estimates
			#png(paste0(results_path, clustx, bootstrap_sample, "_dispersion_plot.png"),
                        #   height = 5, width = 6, units = "in", res = 300)
			#plotDispEsts(dds)
			#dev.off()
			print("# -------- DESEQ COMPLETE -------- #")
		        contrast  <- "compare_B_vs_A"

		
		        res <- results(dds, name = contrast, alpha = 0.05)
		
		        # print(res)
				
		        print("# -------- RES COMPLETE -------- #")
    
		        res <- lfcShrink(dds, coef = contrast, res = res, type = c("ashr"))
    
			res_tbl <- res %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble()
    
			#write.csv(res_tbl, paste0(results_path, clustx, "_", contrast, "boot#_", bootstrap_sample, ".csv"),
			#	  quote = FALSE, row.names = FALSE)
			#print(res_tbl)
			#print(paste0(results_path, "bootstrap_results_", clustx, ".rds"))
			bootstrap_results[[bootstrap_sample]] <- res_tbl


			if (bootstrap_sample %% 5 == 0){
				print(bootstrap_sample)
			}
		}
	
	print(paste0(results_path, "bootstrap_results_", clustx, bootstrap_sample, ".rds"))
	saveRDS(bootstrap_results, paste0(results_path, "bootstrap_results_", clustx, bootstrap_sample, ".rds"))
	}
}
