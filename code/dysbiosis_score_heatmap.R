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

# Trim the species names
rownames(df_diff) <- gsub("s__", "", rownames(df_diff))
rownames(df_diff) <- gsub("_", " ", rownames(df_diff))
rownames(pval_diff) <- gsub("s__", "", rownames(pval_diff))
rownames(pval_diff) <- gsub("_", " ", rownames(pval_diff))
rownames(df_diff_sa) <- gsub("s__", "", rownames(df_diff_sa))
rownames(df_diff_sa) <- gsub("_", " ", rownames(df_diff_sa))
rownames(pval_diff_sa) <- gsub("s__", "", rownames(pval_diff_sa))
rownames(pval_diff_sa) <- gsub("_", " ", rownames(pval_diff_sa))


# Specify species that after adjustment for sexual activity are not significant anymore
species_sa <- read_lines(file.path(out_fig, "sexual_activity_species.txt"))
species_sa <- gsub("s__", "", species_sa)
species_sa <- gsub("_", " ", species_sa)

# Find species that are significant in at least one data modality
sig_species <- rownames(pval_diff)[rowSums(pval_diff < 0.05) > 0]


############# Heatmap before adjustment for SA
cor_mat <- df_diff[sig_species, ]
df_diff_sub <- df_diff[rownames(cor_mat), ]
pval_diff_sub <- pval_diff[rownames(cor_mat), ]

if (sample_type == "Gut") {
  col_names <- c(
    paste0(sample_type, "\nspecies"), "GO",
    paste0(sample_type, "\nmetabolites"), "Plasma\nmetabolites"
  )
} else {
  col_names <- c(
    paste0(sample_type, "\nspecies"), "GO",
    paste0(sample_type, "\nmetabolites"), "Plasma\nmetabolites",
    "Gut\nspecies"
  )
}

# Colors of the rownames
row_name_colors <- rep("#1B9E77", nrow(cor_mat))
row_name_colors[rownames(cor_mat) %in% species_sa] <- "#E7298A"
names(row_name_colors) <- rownames(cor_mat)

col_fun <- colorRamp2(c(min(df_diff), 0, max(df_diff)), c("darkgreen", "white", "red"))
# Create the first heatmap with row clustering
hm1 <- Heatmap(as.matrix(cor_mat),
  row_names_gp = gpar(fontsize = 30, fontface = "italic", col = row_name_colors),
  column_names_gp = gpar(fontsize = 30),
  row_names_side = "left",
  # column_title = ,
  column_title_gp = gpar(fontsize = 50),
  heatmap_height = unit(75, "cm"),
  heatmap_width = unit(50, "cm"),
  column_labels = col_names,
  border = TRUE,
  row_names_max_width = max_text_width(row.names(cor_mat), gp = gpar(fontsize = 30)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (pval_diff_sub[i, j] < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 35))
    } else if (pval_diff_sub[i, j] < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 35))
    } else if (pval_diff_sub[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 35))
    }
  },
  cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun,
  heatmap_legend_param = list(
    legend_height = unit(5, "cm"),
    legend_width = unit(5, "cm"),
    title_gp = gpar(fontsize = 30, fontface = "bold"),
    labels_gp = gpar(fontsize = 30),
    direction = "horizontal",
    title = "Dysbiosis Score"
  )
)


h <- draw(hm1, heatmap_legend_side = "bottom")
pdf(file.path(out_fig, paste0(sample_type, "_heatmap.pdf")), 20, 35)
print(h)
dev.off()


# Get the order of species in the heatmap
row_order_hm1 <- row_order(draw(hm1))
# Second heatmap will show the scores after controlling for sexual activity. The order of species
# between the two heatmaps should be the same and should be based on the ordering of the first heatmap.
cor_mat1 <- df_diff_sa[rownames(cor_mat), ]
pval_diff_sa <- pval_diff_sa[rownames(cor_mat), ]
col_fun <- colorRamp2(c(min(df_diff, df_diff_sa), 0, max(df_diff, df_diff_sa)), c("darkgreen", "white", "red"))

hm2 <- Heatmap(as.matrix(cor_mat1),
  row_names_gp = gpar(fontsize = 30, fontface = "italic"),
  show_heatmap_legend = FALSE,
  column_names_gp = gpar(fontsize = 30),
  row_names_side = "left",
  # column_title = "With adjustement for sexual acitivity",
  column_title_gp = gpar(fontsize = 50),
  heatmap_height = unit(75, "cm"),
  heatmap_width = unit(30, "cm"),
  column_labels = col_names,
  border = TRUE,
  row_names_max_width = max_text_width(row.names(cor_mat1), gp = gpar(fontsize = 30)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (pval_diff_sa[i, j] < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 35))
    } else if (pval_diff_sa[i, j] < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 35))
    } else if (pval_diff_sa[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 35))
    }
  },
  cluster_columns = FALSE, cluster_rows = FALSE, row_order = row_order_hm1, col = col_fun
)


all.equal(rownames(cor_mat)[row_order(draw(hm1))], rownames(cor_mat1)[row_order(draw(hm2))])

ht_list <- hm1 + hm2
h <- draw(ht_list, merge_legends = TRUE, heatmap_legend_side = "bottom", ht_gap = unit(5, "mm"), legend_gap = unit(10, "cm"))
pdf(file.path(out_fig, paste0(sample_type, "_heatmap_controlled_for_sa.pdf")), 40, 35)
print(h)
dev.off()
