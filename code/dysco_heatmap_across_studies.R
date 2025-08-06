Fulcher_neg_species <- read.csv("Outputs/External_Dataset/Fulcher_2022/Negative/Fulcher_2022_Negative_gut_species_dysco.csv", check.names = F)
Fulcher_neg_go <- read.csv("Outputs/External_Dataset/Fulcher_2022/Negative/Fulcher_2022_Negative_gut_go_dysco.csv", check.names = F)
Fulcher_pos_species <- read.csv("Outputs/External_Dataset/Fulcher_2022/Positive/Fulcher_2022_Positive_gut_species_dysco.csv", check.names = F)
Fulcher_pos_go <- read.csv("Outputs/External_Dataset/Fulcher_2022/Positive/Fulcher_2022_Positive_gut_go_dysco.csv", check.names = F)
Garcia_species <- read.csv("Outputs/External_Dataset/Garcia_2024/Garcia_2024_gut_species_dysco.csv", check.names = F)
Garcia_go <- read.csv("Outputs/External_Dataset/Garcia_2024/Garcia_2024_gut_go_dysco.csv", check.names = F)
macs_species <- read.csv("Outputs/Differential_Correlations/Tables/gut_species_dysco.csv", check.names = F)
macs_go <- read.csv("Outputs/Differential_Correlations/Tables/gut_go_dysco.csv", check.names = F)
boston_msm <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/boston/MSM", paste0("Rocafort_boston_MSM", "_gut_species_dysco.csv")), check.names = F)
boston_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/boston/heterosexual", paste0("Rocafort_boston_heterosexual", "_gut_species_dysco.csv")), check.names = F)
bostwana_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/botswana/heterosexual", paste0("Rocafort_botswana_heterosexual", "_gut_species_dysco.csv")), check.names = F)
uganda_hetr <- read.csv(file.path("Outputs/External_Dataset/Rocafort_2024/uganda_2/heterosexual", paste0("Rocafort_uganda_2_heterosexual", "_gut_species_dysco.csv")), check.names = F)
no_ART <- read.csv(file.path("Outputs/External_Dataset/Armstrong_2018/MSM_HIV+ART-", paste0("Armstrong_MSM_HIV+ART-", "_gut_species_dysco.csv")), check.names = F)
with_ART <- read.csv(file.path("Outputs/External_Dataset/Armstrong_2018/MSM_HIV+ART", paste0("Armstrong_MSM_HIV+ART", "_gut_species_dysco.csv")), check.names = F)

datasets <- list(
  macs_species, macs_go,
  Fulcher_neg_species, Fulcher_neg_go, Fulcher_pos_species, Fulcher_pos_go, 
  Garcia_species, Garcia_go, boston_msm, boston_hetr, bostwana_hetr, uganda_hetr, no_ART, with_ART
)

names(datasets) <- c(
  "MACS_species", "MACS_GO",
  "Fulcher_Pre-HIV_species", "Fulcher_Pre-HIV_GO", "Fulcher_Post-HIV_species", "Fulcher_Post-HIV_GO",
  "Garcia_species", "Garcia_GO",
  "Rocafort_Boston_MSM", "Rocafort_Boston_heterosexual", "Rocafort_Bostwana_heterosexual", "Rocafort_Uganda_heterosexual",
  "Armstrong_no_ART", "Armstrong_with_ART"
)


# Significant
datasets_sig <- lapply(datasets, function(x) {
  x_sig <- x %>% dplyr::filter(p < 0.01 & num_corr > max(x[, "num_corr"]) * 0.2 & BH_p < 0.1)
})

# For each dataset extract the genus names
for (i in c(1:8)) {
  print(i)
  datasets_sig[[i]] <- datasets_sig[[i]] %>%
    mutate(genus = unlist(purrr::map(stringi::stri_split(species, fixed = "_"), 3))) %>%
    group_by(genus) %>%
    summarise(num_genus = n()) %>%
    mutate(dataset = names(datasets)[i])
}
for (i in 9:12) {
  datasets_sig[[i]] <- datasets_sig[[i]] %>%
    mutate(genus = unlist(purrr::map(stringi::stri_split(species, fixed = " "), 1))) %>%
    group_by(genus) %>%
    summarise(num_genus = n()) %>%
    mutate(dataset = names(datasets)[i])
}

for (i in 13:14) {
  datasets_sig[[i]] <- datasets_sig[[i]] %>%
    mutate(genus = gsub("s__|g__|s__uncultured_", "", species)) %>%
    mutate(genus = unlist(purrr::map(stringi::stri_split(genus, fixed = "_"), 1))) %>%
    mutate(genus = gsub("\\[|]", "", genus)) %>%
    group_by(genus) %>%
    summarise(num_genus = n()) %>%
    mutate(dataset = names(datasets)[i])
}

datasets_sig <- bind_rows(datasets_sig)
datasets_sig_long <- datasets_sig %>%
  pivot_wider(id_cols = "genus", names_from = "dataset", values_from = "num_genus") %>%
  column_to_rownames("genus")
datasets_sig_long[is.na(datasets_sig_long)] <- 0
sample_type <- "Gut"
library(circlize)
library(ComplexHeatmap)
col_fun <- colorRamp2(c(min(datasets_sig_long), max(datasets_sig_long)), c("white", "red"))

datasets_sig_long <- datasets_sig_long[names(sort(rowSums(datasets_sig_long), decreasing = T)), ]
hm <- Heatmap(as.matrix(datasets_sig_long),
  row_names_gp = gpar(fontsize = 40, fontface = "italic"),
  column_names_gp = gpar(fontsize = 40),
  row_names_side = "left",
  row_dend_side = "right",
  column_title = paste(sample_type, "Genera"),
  column_title_gp = gpar(fontsize = 60),
  heatmap_height = unit(110, "cm"),
  heatmap_width = unit(45, "cm"),
  column_labels = names(datasets),
  border = TRUE,
  row_names_max_width = max_text_width(row.names(datasets_sig_long), gp = gpar(fontsize = 20)),
  cluster_columns = F, cluster_rows = F, col = col_fun,
  heatmap_legend_param = list(
    legend_height = unit(7, "cm"),
    legend_width = unit(5, "cm"),
    title_gp = gpar(fontsize = 40, fontface = "bold"),
    labels_gp = gpar(fontsize = 40),
    direction = "horizontal",
    title = "Number of dysbiotic \nspecies/genus"
  )
)

h <- draw(hm, heatmap_legend_side = "right")
pdf(file.path("Outputs/External_Dataset", "Genus_dysco_heatmap.pdf"), 30, 60)
print(h)
dev.off()
