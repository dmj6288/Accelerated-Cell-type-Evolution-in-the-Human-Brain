library(Seurat)
library(hash)
library(readr)

Bakken_metadata_file <- read.table("author_source/metadata.txt", sep = "\t")


traditional_cell_types <- c("Excite", "Astrocyte", "Inhibit", "Microglia", "OPC", "Oligo")
emre_broadannot        <- c("Excitatory", "Astrocyte", "Inhibitory",  "Microglia"  ,"OPC", "MOL")


Bakken_metadata_file <- Bakken_metadata_file[Bakken_metadata_file$broadAnnot %in% emre_broadannot,]

local_to_standard_name_hash <- hash(emre_broadannot, traditional_cell_types)

modified_cluster_names <- c()
for (cell_label in Bakken_metadata_file$broadAnnot){
  modified_cluster_names <- c(modified_cluster_names, 
                              local_to_standard_name_hash[[cell_label]])
}

df <- data.frame(as.factor(paste(Bakken_metadata_file$CellType, Bakken_metadata_file$Species, sep = "_#")),
                 as.factor(Bakken_metadata_file$Sample), 
                 as.factor(Bakken_metadata_file$lib_batch),
                 Bakken_metadata_file$human_age,
                 as.factor(Bakken_metadata_file$sex))



# df <- data.frame(as.factor(Bakken_metadata_file$CellType),
#                  as.factor(Bakken_metadata_file$Sample))

colnames(df) <- c("trad", "sample_id", "lib_batch", "human_age", "sex")
rownames(df) <- rownames(Bakken_metadata_file)

df_human_chimp <- df[df$sample_id %in% c("h1", "h2", "h3", "h4", "c1", "c2", "c3", "c4"), ]
df_chimp_macak <- df[df$sample_id %in% c("c1", "c2", "c3", "c4", "m1", "m2", "m3", "m4"), ]
df_human_macak <- df[df$sample_id %in% c("h1", "h2", "h3", "h4", "m1", "m2", "m3", "m4"), ]


df_human_chimp$trad <- droplevels(df_human_chimp$trad)
df_chimp_macak$trad <- droplevels(df_chimp_macak$trad)
df_human_macak$trad <- droplevels(df_human_macak$trad)


df_human_chimp$sample_id <- droplevels(df_human_chimp$sample_id)
df_chimp_macak$sample_id <- droplevels(df_chimp_macak$sample_id)
df_human_macak$sample_id <- droplevels(df_human_macak$sample_id)


saveRDS(df_human_chimp, "pairwise_metadata_files/df_human_chimp.rds")
saveRDS(df_chimp_macak, "pairwise_metadata_files/df_chimp_macak.rds")
saveRDS(df_human_macak, "pairwise_metadata_files/df_human_macak.rds")


