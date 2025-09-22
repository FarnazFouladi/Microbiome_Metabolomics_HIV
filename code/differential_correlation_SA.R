source("code/load.R")
source("code/differential_correlations_functions_v1.R")



g1_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g1.rds"))
g2_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g2.rds"))
g3_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g3.rds"))
g4_bias_corrected <- read_rds(file.path("Outputs/Processed_data", "processed_g4.rds"))


output_dir <- file.path("Outputs", "Differential_Correlations_SA","Tables")
dir.create(output_dir, showWarnings = F, recursive = T)


# Gut species
for (i in 2:4){
  
  g_bias_corrected <- get(paste0("g",i,"_bias_corrected"))
  
  
  df.g1 <- t(g1_bias_corrected[[1]])
  df.g <- t(g_bias_corrected[[1]])
  common_features <- intersect(colnames(df.g1), colnames(df.g))
  common_features <- common_features[!grepl("s__sp.", common_features)]
  df.g1 <- df.g1[, common_features]
  df.g <- df.g[, common_features]
  species_res <- disco_score(data1 = list(df.g1), data2 = list(df.g),
                             output_dir = output_dir,prefix = paste0("g",i,"_gut_species"))
}


# Gut GO
for (i in 2:4){
  
  g_bias_corrected <- get(paste0("g",i,"_bias_corrected"))
  
  
  df.g1 <- t(g1_bias_corrected[[2]])
  df.g <- t(g_bias_corrected[[2]])
  common_features <- intersect(colnames(df.g1), colnames(df.g))
  df.g1 <- df.g1[, common_features]
  df.g <- df.g[, common_features]
  go_res <- disco_score(data1 = list(df.g1), data2 = list(df.g),
                             output_dir = output_dir,prefix = paste0("g",i,"_gut_go"))
}



# Gut Go
go_annotation <- readRDS("data/processed_data/gut/GO_phy.rds")
go_annotation <- tax_table(go_annotation) %>% as.data.frame() %>% rownames_to_column("t2") %>%
  dplyr::select(c("t2","Description","Ontology"))

go.nc <- t(nc_bias_corrected[[2]])
go.sc <- t(sc_bias_corrected[[2]])
common_features <- intersect(colnames(go.nc), colnames(go.sc))
go.nc <- go.nc[, common_features]
go.sc <- go.sc[, common_features]
go_res <- disco_score(
  data1 = list(df.nc, go.nc), data2 = list(df.sc, go.sc), 
  output_dir = output_dir,prefix = "gut_go"
)
# Add annotation to significant correlations
go_res_sig <- go_res[[2]] %>%
  left_join(go_annotation) 
write.csv(go_res_sig, file.path(output_dir, "gut_go_cor_diff_sig.csv"), quote = T, row.names = F)

# Gut Metabolites - Positive mode
mtb_ann <- read.table(file.path("data/metabolomics_annotation/metabolites_annotation.txt"), sep = "\t", comment.char = "", quote = "", check.names = F, header = T)
mtb_ann <- mtb_ann %>%
  dplyr::rename(t2 = accession, Metabolites = Compound_name) %>%
  dplyr::select(-c("name"))

meta.p.nc <- t(nc_bias_corrected[[3]])
meta.p.sc <- t(sc_bias_corrected[[3]])
common_features <- intersect(colnames(meta.p.nc), colnames(meta.p.sc))
meta.p.nc <- meta.p.nc[, common_features]
meta.p.sc <- meta.p.sc[, common_features]
meta.p_res <- disco_score(
  data1 = list(df.nc, meta.p.nc), data2 = list(df.sc, meta.p.sc),
  output_dir = output_dir,prefix = "gut_meta.p"
)
meta.p_res_sig <- meta.p_res[[2]] %>%
  left_join(mtb_ann)
write.csv(meta.p_res_sig, file.path(output_dir, "gut_meta.p_cor_diff_sig.csv"), quote = T, row.names = F)


# Gut Metabolites - Negative mode
meta.n.nc <- t(nc_bias_corrected[[4]])
meta.n.sc <- t(sc_bias_corrected[[4]])
common_features <- intersect(colnames(meta.n.nc), colnames(meta.n.sc))
meta.n.nc <- meta.n.nc[, common_features]
meta.n.sc <- meta.n.sc[, common_features]
meta.n_res <- disco_score(
  data1 = list(df.nc, meta.n.nc), data2 = list(df.sc, meta.n.sc),
  output_dir = output_dir,prefix = "gut_meta.n"
)
meta.n_res_sig <- meta.n_res[[2]] %>%
  left_join(mtb_ann) 
write.csv(meta.n_res_sig, file.path(output_dir, "gut_meta.n_cor_diff_sig.csv"), quote = T, row.names = F)

# Plasma Metabolites - Positive mode
meta.plasma.p.nc <- t(nc_bias_corrected[[5]])
meta.plasma.p.sc <- t(sc_bias_corrected[[5]])
common_features <- intersect(colnames(meta.plasma.p.nc), colnames(meta.plasma.p.sc))
meta.plasma.p.nc <- meta.plasma.p.nc[, common_features]
meta.plasma.p.sc <- meta.plasma.p.sc[, common_features]
meta.plasma.p_res <- disco_score(
  data1 = list(df.nc, meta.plasma.p.nc), data2 = list(df.sc, meta.plasma.p.sc),
  output_dir = output_dir,prefix = "gut_meta.plasma.p"
)
meta.plasma.p_res_sig <- meta.plasma.p_res[[2]] %>%
  left_join(mtb_ann) 
write.csv(meta.plasma.p_res_sig , file.path(output_dir, "gut_meta.plasma.p_cor_diff_sig.csv"), quote = T, row.names = F)

# Plasma Metabolites - Negative mode
meta.plasma.n.nc <- t(nc_bias_corrected[[6]])
meta.plasma.n.sc <- t(sc_bias_corrected[[6]])
common_features <- intersect(colnames(meta.plasma.n.nc), colnames(meta.plasma.n.sc))
meta.plasma.n.nc <- meta.plasma.n.nc[, common_features]
meta.plasma.n.sc <- meta.plasma.n.sc[, common_features]
meta.plasma.n_res <- disco_score(
  data1 = list(df.nc, meta.plasma.n.nc), data2 = list(df.sc, meta.plasma.n.sc),
  output_dir = output_dir,prefix = "gut_meta.plasma.n"
)
meta.plasma.n_res_sig <- meta.plasma.n_res[[2]] %>%
  left_join(mtb_ann) 
write.csv(meta.plasma.n_res_sig, file.path(output_dir, "gut_meta.plasma.n_cor_diff_sig.csv"), quote = T, row.names = F)

# Oral species
df.nc.oral <- t(nc_bias_corrected[[7]])
df.sc.oral <- t(sc_bias_corrected[[7]])
common_features <- intersect(colnames(df.nc.oral), colnames(df.sc.oral))
common_features <- common_features[!grepl("s__sp.", common_features)]
df.nc.oral <- df.nc.oral[, common_features]
df.sc.oral <- df.sc.oral[, common_features]
oral_species_res <- disco_score(
  data1 = list(df.nc.oral),
  data2 = list(df.sc.oral),
  output_dir = output_dir,prefix = "oral_species"
)

# Oral Go
go_annotation <- readRDS("data/processed_data/oral/GO_phy.rds")
go_annotation <- tax_table(go_annotation) %>% as.data.frame() %>% rownames_to_column("t2") %>%
  dplyr::select(c("t2","Description","Ontology"))

go.nc.oral <- t(nc_bias_corrected[[8]])
go.sc.oral <- t(sc_bias_corrected[[8]])
common_features <- intersect(colnames(go.nc.oral), colnames(go.sc.oral))
go.nc.oral <- go.nc.oral[, common_features]
go.sc.oral <- go.sc.oral[, common_features]
go_res_oral <- disco_score(
  data1 = list(df.nc.oral, go.nc.oral), 
  data2 = list(df.sc.oral, go.sc.oral),
  output_dir = output_dir,prefix = "oral_go"
)
go_res_oral_sig <- go_res_oral[[2]] %>%
  left_join(go_annotation) 

write.csv(go_res_oral_sig, file.path(output_dir, "oral_go_cor_diff_sig.csv"), quote = T, row.names = F)

# Oral Metabolites - Positive mode
meta.p.nc.oral <- t(nc_bias_corrected[[9]])
meta.p.sc.oral <- t(sc_bias_corrected[[9]])
common_features <- intersect(colnames(meta.p.nc.oral), colnames(meta.p.sc.oral))
meta.p.nc.oral <- meta.p.nc.oral[, common_features]
meta.p.sc.oral <- meta.p.sc.oral[, common_features]
meta.p_oral <- disco_score(
  data1 = list(df.nc.oral, meta.p.nc.oral), 
  data2 = list(df.sc.oral, meta.p.sc.oral),
  output_dir = output_dir,prefix = "oral_meta.p"
)
meta.p_oral_sig <- meta.p_oral[[2]] %>%
  left_join(mtb_ann) 

write.csv(meta.p_oral_sig, file.path(output_dir, "oral_meta.p_cor_diff_sig.csv"), quote = T, row.names = F)

# Gut Metabolites - Negative mode
meta.n.nc.oral <- t(nc_bias_corrected[[10]])
meta.n.sc.oral <- t(sc_bias_corrected[[10]])
common_features <- intersect(colnames(meta.n.nc.oral), colnames(meta.n.sc.oral))
meta.n.nc.oral <- meta.n.nc.oral[, common_features]
meta.n.sc.oral <- meta.n.sc.oral[, common_features]
meta.n_oral <- disco_score(
  data1 = list(df.nc.oral, meta.n.nc.oral), 
  data2 = list(df.sc.oral, meta.n.sc.oral),
  output_dir = output_dir,prefix = "oral_meta.n"
)
meta.n_oral_sig <- meta.n_oral[[2]] %>%
  left_join(mtb_ann) 
write.csv( meta.n_oral_sig, file.path(output_dir, "oral_meta.n_cor_diff_sig.csv"), quote = T, row.names = F)

# Plasma Metabolites - Positive mode
meta.plasma.p_res.oral <- disco_score(
  data1 = list(df.nc.oral, meta.plasma.p.nc), 
  data2 = list(df.sc.oral, meta.plasma.p.sc),
  output_dir = output_dir,prefix = "oral_meta.plasma.p"
)
meta.plasma.p_res.oral_sig <- meta.plasma.p_res.oral[[2]] %>%
  left_join(mtb_ann) 

write.csv(meta.plasma.p_res.oral_sig, file.path(output_dir, "oral_meta.plasma.p_cor_diff_sig.csv"), quote = T, row.names = F)

# Plasma Metabolites - Negative mode
meta.plasma.n_res.oral <- disco_score(
  data1 = list(df.nc.oral, meta.plasma.n.nc), 
  data2 = list(df.sc.oral, meta.plasma.n.sc),
  output_dir = output_dir,prefix = "oral_meta.plasma.n"
)
meta.plasma.n_res.oral_sig <- meta.plasma.n_res.oral[[2]] %>%
  left_join(mtb_ann) 

write.csv(meta.plasma.n_res.oral_sig, file.path(output_dir, "oral_meta.plasma.n_cor_diff_sig.csv"), quote = T, row.names = F)

# Oral and gut species
colnames(df.nc) <- paste0("gut_", colnames(df.nc))
colnames(df.sc) <- paste0("gut_", colnames(df.sc))

gut.oral <- disco_score(
  data1 = list(df.nc.oral, df.nc), data2 = list(df.sc.oral, df.sc),
  output_dir = output_dir,prefix = "gut.oral"
)
gut.oral_sig1 <- gut.oral[[2]] %>% mutate(oral_genus = unlist(purrr::map(stringi::stri_split_fixed(t1,"_"),3))) %>%
  mutate(gut_genus = unlist(purrr::map(stringi::stri_split_fixed(t2,"_"),4)))%>%
  group_by(oral_genus) %>% mutate(oral_genus_n = n()) %>%
  group_by(gut_genus) %>% mutate(gut_genus_n = n())
  
  
