rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
sample_type <- args[1]
source("code/load.R")
# Directory including DYSCO scores without controlling for sexual activity
output_dir <- file.path("Outputs", "Dysbiosis", sample_type)
# Directory including DYSCO scores with controlling for sexual activity
output_dir_sa <- file.path("Outputs", "Dysbiosis_Controlled_for_SA", sample_type)
out_fig <- file.path(output_dir, "Figures")
# Upload scores
df_diff <- read.csv(file.path(output_dir, "Scores", "scores_wide_diff_p.csv")) %>% column_to_rownames("X")
df_diff_sa <- read.csv(file.path(output_dir_sa, "Scores", "scores_wide_diff_p.csv")) %>% column_to_rownames("X")

#  Upload p-values
pval_diff <- read_csv(file.path(output_dir, "Scores", "pvals_p_diff.csv")) %>% column_to_rownames("...1")
pval_diff_sa <- read_csv(file.path(output_dir_sa, "Scores", "pvals_p_diff.csv")) %>% column_to_rownames("...1")
all.equal(rownames(pval_diff), rownames(df_diff))
all.equal(rownames(pval_diff_sa), rownames(df_diff_sa))

################### Summary
# gut species
mylist <- list()
myspecies <- list() # Species that were not significant after controlling for sexual activity
datasets_num <- ifelse(sample_type == "Gut", 4, 5)
for (index in 1:datasets_num) {
  sig_species_sc <- pval_diff[pval_diff[, index] < 0.05, ]
  sig_species_sa <- pval_diff_sa[pval_diff_sa[, index] < 0.05, ]

  # Proportion of dysbiotic species that were remained significant after controlling for sexual activity
  print(sum(rownames(sig_species_sc) %in% rownames(sig_species_sa)) / nrow(sig_species_sc))

  # dysbiotic species that were NOT remained significant after controlling for sexual activity
  non_sig <- rownames(sig_species_sc)[!rownames(sig_species_sc) %in% rownames(sig_species_sa)]

  mylist[[index]] <- list(sig_species_sc[non_sig, ], pval_diff_sa[non_sig, ])

  myspecies[[index]] <- non_sig
}
vec <- sort(unique(unlist(myspecies)))
writeLines(vec, file.path(out_fig, "sexual_activity_species.txt"))

################ Scatterplots comparing scores with and without controlling for sexual activity

# Merge data
df_merge <- df_diff %>%
  rownames_to_column("Taxa") %>%
  full_join(df_diff_sa %>% rownames_to_column("Taxa"), by = "Taxa", suffix = c("_sc", "_sa"))

df_merge_p <- pval_diff %>%
  rownames_to_column("Taxa") %>%
  full_join(pval_diff_sa %>% rownames_to_column("Taxa"), by = "Taxa", suffix = c("_sc", "_sa"))

df_merge_p <- df_merge_p %>% arrange(order(match(Taxa, df_merge$Taxa)))
all.equal(df_merge$Taxa, df_merge_p$Taxa)

mycol <- c(
  "Not explained by the number of sexual partners" = "#1B9E77",
  "Potentially due to the number of sexual partners" = "#E7298A",
  "Possibly not dysbiotic" = "grey"
)

df_list <- list()
plot_list <- list()

if (sample_type == "Gut") {
  datasets <- c("species", "go", "mtb", "plasma_mtb")
} else {
  datasets <- c("species", "go", "mtb", "plasma_mtb", "species2")
}


for (d in datasets) {
  if (d == "species") {
    d1 <- paste(sample_type, "Species")
  } else if (d == "go") {
    d1 <- paste(sample_type, "GO")
  } else if (d == "mtb") {
    d1 <- paste(sample_type, "Metabolites")
  } else if (d == "plasma_mtb") {
    d1 <- paste("Plamsa", "Metabolites")
  } else {
    d1 <- paste("Gut", "Species")
  }

  col_name_sc <- paste0(d, "_SC_sc")
  col_name_sa <- paste0(d, "_SC_sa")

  df_merge1 <- df_merge %>%
    mutate(significant_sa = ifelse(df_merge_p[, col_name_sa] < 0.05 & df_merge_p[, col_name_sc] < 0.05, "Not explained by the number of sexual partners",
      ifelse(df_merge_p[, col_name_sa] < 0.05 & df_merge_p[, col_name_sc] > 0.05, "Not explained by the number of sexual partners",
        ifelse(df_merge_p[, col_name_sa] > 0.05 & df_merge_p[, col_name_sc] < 0.05, "Potentially due to the number of sexual partners", "Possibly not dysbiotic")
      )
    )) %>%
    mutate(significant_sa = factor(significant_sa, levels = names(mycol))) %>%
    mutate(
      Taxa = gsub("s__", "", Taxa),
      Taxa = gsub("_", " ", Taxa),
      Taxa_itl = paste0("italic('", Taxa, "')")
    )


  # The index for the highest dysbiotic scores. Those species will be annotated in the figures
  index <- order(apply(abs(df_merge1[c(col_name_sc, col_name_sa)]), 1, max), decreasing = T)
  df_merge_sub_sa <- df_merge1 %>%
    filter(Taxa %in% df_merge1[index, ]$Taxa[1:5]) %>%
    filter(significant_sa != "Possibly not dysbiotic")


  cor_res <- cor.test(df_merge1[, col_name_sa], df_merge1[, col_name_sc], method = "pearson")

  p_sa <- ggplot(na.omit(df_merge1), aes(x = .data[[col_name_sa]], y = .data[[col_name_sc]])) +
    geom_jitter(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.7, aes(col = significant_sa)) +
    geom_text_repel(
      data = df_merge_sub_sa, aes(x = .data[[col_name_sa]], y = .data[[col_name_sc]], label = Taxa_itl),
      min.segment.length = 0.1,
      size = 3,
      parse = TRUE,
      segment.color = "grey50"
    ) +
    theme_bw() +
    theme(
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    ) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = d1, y = "DYSCO without adjustement ",
      x = "DYSCO with adjustment",
      color = "Dysbiosis",
      subtitle = paste0("R = ", round(cor_res$estimate, 3), "; ", "p ", rstatix::p_format(cor_res$p.value, digits = 4))
    ) +
    scale_color_manual(values = mycol)

  plot_list[[d]] <- p_sa

  # Save the table
  df_merge2 <- df_merge1 %>%
    filter(significant_sa != "not significant") %>%
    dplyr::select(Taxa, significant_sa)
  colnames(df_merge2)[2] <- d

  df_list[[d]] <- df_merge2
}

df_all <- Reduce(df_list, f = full_join)
write.csv(df_all, file.path(out_fig, paste0(sample_type, "_dysbiosis_type.csv")), row.names = F, quote = F)

lg <- cowplot::get_legend(plot_list[[1]])
plot_list <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
plot_list[[6]] <- lg

if (sample_type == "Gut") {
  lbs <- c("A", "B", "C", "D")
} else {
  lbs <- c("A", "B", "C", "D", "E")
}


pdf(file.path(out_fig, paste0(sample_type, "_scatterplots.pdf")), width = 11, height = 13)
print(cowplot::plot_grid(plotlist = plot_list, nrow = 3, ncol = 2, labels = lbs))
dev.off()
