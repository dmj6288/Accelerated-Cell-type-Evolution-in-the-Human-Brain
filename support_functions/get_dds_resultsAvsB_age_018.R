get_dds_resultsAvsB <- function(padj_cutoff = 0.05, study, counts_ls, metadata_ls, count_ls, path) {

        results_path <- paste(path, study, "/", count_ls, "/", sep = "")
	if (!dir.exists(results_path)) { dir.create(results_path) }
			
	for (i in c(1:length(metadata_ls))){
		clustx <- names(counts_ls)[i]
		print(clustx) # useful for debugging
		# Extract counts matrix and metadata for cluster x
		idx <- which(names(counts_ls) == clustx)
    
		print(idx)
    
		clust_mat <- counts_ls[[idx]]
		#clust_df <- clust_df[!apply(clust_df, 1, function(x) all(x == 0)), ]
		print("# -------- CLUST MAT -------- #")
		print(head(clust_mat, 20))

		#cluster_counts <- clust_mat[1:max(which(rowSums(clust_mat != 0) > 0)), ]

		cluster_counts <- clust_mat[which(rowSums(clust_mat) != 0),]
		print(cluster_counts[which(rowSums(cluster_counts) == 0),])
		print(any(is.na(cluster_counts)))		
    
		#cluster_counts <- as.matrix(clust_df)
		cluster_metadata <- metadata_ls[[idx]]

		print("# -------- CLUST METADATA -------- #")
		print(cluster_metadata)        
#		print(dim(cluster_counts))
#		print(head(cluster_counts))
    
		# Print error message if sample names do not match
		if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
			print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
		}
    
#		dds <- DESeqDataSetFromMatrix(cluster_counts,
#                                              colData = cluster_metadata, 
#                                              design = ~human_age + sex + compare)

		dds <- DESeqDataSetFromMatrix(cluster_counts,
                                              colData = cluster_metadata,
                                              design = ~human_age + compare)

		print("# -------- DESeqDataSetFromMatrix COMPLETE -------- #")
		# Transform counts for data visualization
		rld <- rlog(dds, blind = TRUE)
    
    
		# Generate QC plots
    
		## Plot and save PCA plot
		DESeq2::plotPCA(rld, intgroup = "cluster_id")
    
		if (!dir.exists(results_path)) { dir.create(results_path) }
		ggsave(paste0(results_path, clustx, "_specific_PCAplot.png"))
    
    		## Extract rlog matrix from the object and compute pairwise correlation values
		rld_mat <- assay(rld)
		rld_cor <- cor(rld_mat)
    
		## Plot and save heatmap
		png(paste0(results_path, clustx, "_specific_heatmap.png"),
			   height = 6, width = 7.5, units = "in", res = 300)
		pheatmap(rld_cor, annotation = cluster_metadata[, c("cluster_id"), drop = FALSE])
		dev.off()
		
#		keep <- rowSums(counts(dds) >= 10) >= 4
#		dds <- dds[keep,]    
    
		# Run DESeq2 differential expression analysis
		dds <- DESeq(dds)
		print("# -------- DESEQ COMPLETE -------- #")
		## Plot dispersion estimates
		png(paste0(results_path, clustx, "_dispersion_plot.png"),
			   height = 5, width = 6, units = "in", res = 300)
		plotDispEsts(dds)
		dev.off()
    
		## Output and shrink results of Wald test for contrast A vs B
    
		#contrast <- paste(c("group_id", unique(cluster_metadata$cluster_id)[2], "vs", unique(cluster_metadata$cluster_id)[1]), collapse = "_")
		contrast  <- "compare_B_vs_A"

		print("# -------- RESNAMES DDS -------- #")
		print(resultsNames(dds))
		# print(dds)
		
		res <- results(dds, name = contrast, alpha = 0.05)
		
		# print(res)
				
		print("# -------- RES COMPLETE -------- #")
    
		res <- lfcShrink(dds, coef = contrast, res = res, type = c("ashr"))
		print(res)

		## Turn the results object into a tibble for use with tidyverse functions
    
		res_tbl <- res %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble()
    
		write.csv(res_tbl, paste0(results_path, clustx, "_", contrast, "_all_genes.csv"),
			  quote = FALSE, row.names = FALSE)
    
		## Subset the significant results
		sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj)
    
		write.csv(sig_res, paste0(results_path, clustx, "_", contrast, "_signif_genes.csv"),
                          quote = FALSE, row.names = FALSE)
    
    
		# Generate results visualization plots
    
		## Extract normalized counts from dds object
		normalized_counts <- counts(dds, normalized = TRUE)
    
		## Extract top 20 DEG from resLFC (make sure to order by padj)
		top20_sig_genes <- sig_res %>% dplyr::arrange(padj) %>% dplyr::pull(gene) %>% head(n = 20)
    
		## Extract matching normalized count values from matrix 
		top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
    
    		## Convert wide matrix to long data frame for ggplot2
		top20_sig_df <- data.frame(top20_sig_counts)
		top20_sig_df$gene <- rownames(top20_sig_counts)
    
		#top20_sig_df <- melt(setDT(top20_sig_df), id.vars = c("gene"), variable.name = "cluster_sample_id") %>% data.frame()
    
		## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
    
		#top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
		#top20_sig_df$cluster_sample_id <- gsub("\\  ", "+ ", top20_sig_df$cluster_sample_id)
    
    		## Join counts data frame with metadata
    
		#top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)), by = "cluster_sample_id")
    
		## Generate plot
		#ggplot(top20_sig_df, aes(y = value, x = cluster_id)) + geom_jitter(height = 0, width = 0.15)        +
		#scale_y_continuous(trans = 'log10')                  + ylab("log10 of normalized expression level") +
		#xlab("condition")                                    + ggtitle("Top 20 Significant DE Genes")       +
		#theme(plot.title = element_text(hjust = 0.5))        + facet_wrap(~ gene)
    
		#ggsave(paste0(results_path, clustx, "_", contrast, "_top20_DE_genes.png"))
		#return(cluster_counts) 
	}
}
