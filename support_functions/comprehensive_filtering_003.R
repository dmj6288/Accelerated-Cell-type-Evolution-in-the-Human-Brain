comprehensive_filtering <- function(seurat_object, study_and_region_name) {
      
      # PIPELINE SOURCE : https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
      #path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/filtered/data_quality_figures/"
	path <- "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/Evolution/data_quality_figures_active_ident/"
	seurat_object <- PercentageFeatureSet(seurat_object, "^MT-", col.name = "percent_mito")
      
      # Column sum of indicated gene name pattern divided by total sum of features
	      
	seurat_object <- PercentageFeatureSet(seurat_object, "^RP[SL]", col.name = "percent_ribo")
      
	seurat_object <- PercentageFeatureSet(seurat_object, "^HB[^(P)]", col.name = "percent_hb")
      
	seurat_object <- PercentageFeatureSet(seurat_object, "PECAM1|PF4", col.name = "percent_plat")
      
	# REMOVE GENES WHICH ARE NOT EXPRESSED IN LESS THAN OR EQUAL TO 3 CELLS
	
	selected_c <- WhichCells(seurat_object, expression = nFeature_RNA > 200)
	selected_f <- rownames(seurat_object)[Matrix::rowSums(seurat_object) > 3]
	
	print(head(selected_f))      
	data.filt <- subset(seurat_object, features = selected_f, cells = selected_c)
      
	print(head(data.filt))

	par(mar = c(4, 8, 2, 1))
	C <- data.filt@assays$RNA@counts
	C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
	most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
      
	print("Boxplot params")
	#print(C)
	#print(most_expressed)

	most_expressed_genes <- rownames(C)[most_expressed]
	print(most_expressed_genes)
	
	C_max <- t(as.matrix(C[most_expressed_genes, ]))

	png(paste(path,study_and_region_name, "high_expression_genes.png", sep = ""))
	boxplot(C_max, cex = 0.1, las = 1, xlab = "% total count per cell",
		col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
	dev.off()

	print("By Detection")
	print(dim(data.filt))
      
	selected_mito <- WhichCells(data.filt, expression = percent_mito < 5)
	#selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 0.025)
      
	# and subset the object to only keep those cells
	data.filt <- subset(data.filt, cells = selected_mito)
	#data.filt <- subset(data.filt, cells = selected_ribo)
      
	print("By Mito/Ribo genes")
	print(dim(data.filt))
      
	feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
      
	VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
		NoLegend()
	ggsave(paste(path,study_and_region_name, "gene_percent_metrics.png", sep = ""))
      
	# FILTER BY GENES
      
	print(dim(data.filt))
      
	# Filter MALAT1
	data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]
      
	# Filter Mitocondrial
	#data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
      
	# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
	# <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]
      
	# Filter Hemoglobin gene (optional if that is a problem on your data)
	#data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]
      
	print(dim(data.filt))
      
	genes.file = "/home/local/ADS/dennis00/Yi_Lab/Lab_WorkDir/Dennis/pseudo_single_pipeline/sex_chromosome_gene_table/results/genes.table.csv"
	genes.table = read.csv(genes.file)
	genes.table <- genes.table[genes.table$external_gene_name %in% rownames(data.filt),]
      
	chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
      
	chrY.gene = chrY.gene[chrY.gene != ""]
      
	chrY.gene = chrY.gene[chrY.gene %in% rownames(data.filt)]
      
	print(chrY.gene)

#      data.filt$pct_chrY = colSums(data.filt@assays$RNA@counts[chrY.gene, ])/colSums(data.filt@assays$RNA@counts)
      
#      FeatureScatter(data.filt, feature2 = "pct_chrY")
      
#      VlnPlot(data.filt, group.by = "orig.ident", features = c("pct_chrY"))      

#      ggsave(paste(path, study_and_region_name, "XY.png", sep = ""))
      
	data.filt = NormalizeData(data.filt)
      
      
#	data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes,
#                                    s.features = cc.genes$s.genes)
      
#	VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
#              ncol = 4, pt.size = 0.1)
      
#	ggsave(paste(path, study_and_region_name, "division_cycle.png", sep = ""))
      
	suppressMessages(require(DoubletFinder))
      
	data.filt = FindVariableFeatures(data.filt, verbose = F)
	data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"),
                            verbose = F)
	data.filt = RunPCA(data.filt, verbose = F, npcs = 20)

	# get significant PCs
	stdv         <- data.filt[["pca"]]@stdev
	sum.stdv     <- sum(data.filt[["pca"]]@stdev)
	percent.stdv <- (stdv / sum.stdv) * 100
	cumulative   <- cumsum(percent.stdv)
	co1          <- which(cumulative > 90 & percent.stdv < 5)[1]
	co2          <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                                    percent.stdv[2:length(percent.stdv)]) > 0.1), 
                                    decreasing = T)[1] + 1
	min.pc       <- min(co1, co2)
	min.pc

	data.filt = RunUMAP(data.filt, dims = 1:min.pc, verbose = F)

	sweep.res <- paramSweep_v3(data.filt) 
	sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
	bcmvn <- find.pK(sweep.stats)

	# Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions
	bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric),]
	pK_value <- bcmvn_max$pK
	pK_value <- as.numeric(levels(pK_value))[pK_value]

	png(paste(path,study_and_region_name, "param_sweep.png", sep = ""))
	barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
	dev.off()

	nExp <- round(ncol(data.filt) * 0.04)  # expect 4% doublets

	data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = pK_value, nExp = nExp, PCs = 1:min.pc)
	DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
      
      
	cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                         DimPlot(data.filt, group.by = DF.name) + NoAxes())
      
	ggsave(paste(path, study_and_region_name, "doublets.png", sep = ""))
      
	VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
      
	ggsave(paste(path, study_and_region_name, "doublet_violin.png", sep =""))
	data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
	print(dim(data.filt))
      
	return(data.filt)
}
