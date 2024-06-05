library(Seurat)
library(hash)
library(readr)

Bakken_metadata_file <- read.table("author_source/metadata.txt", sep = "\t")


emre_broadannot        <- c("Excitatory", "Astrocyte", "Inhibitory",  "Microglia"  ,"OPC", "MOL")

broad_annot_hash              <- readRDS("Caglayan_cell_name_hash.rds")

Bakken_metadata_file <- Bakken_metadata_file[Bakken_metadata_file$broadAnnot %in% emre_broadannot,]


modified_cluster_names <- c()
for (cell_label in Bakken_metadata_file$CellType){
  modified_cluster_names <- c(modified_cluster_names, 
                              broad_annot_hash[[cell_label]])
}

df <- data.frame(as.factor(paste(modified_cluster_names, Bakken_metadata_file$Species, sep = "_#")),
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


saveRDS(df_human_chimp, "pairwise_metadata_major_cell_types/df_human_chimp.rds")
saveRDS(df_chimp_macak, "pairwise_metadata_major_cell_types/df_chimp_macak.rds")
saveRDS(df_human_macak, "pairwise_metadata_major_cell_types/df_human_macak.rds")


