# Load libraries
packages <- c(
  "TreeSummarizedExperiment", "phyloseq", "ggrepel", "ggpubr", "gridExtra",
  "ComplexHeatmap", "circlize", "viridis", "tibble", "tidyverse", "readr", "vegan", "rstatix",
  "NetCoMi", "RColorBrewer", "readxl", "kableExtra","microViz"
)
invisible(lapply(packages, library, character.only = TRUE))

output_dir <- "Outputs"
code_dir <- "code"

source(file.path(code_dir, "functions.R"))

# Set variables
status_cols <- c("#0072b5ff", "#bc3c29ff")
status <- c("nc", "sc")
group_names <- c("Non-HIV", "Pre-HIV")
names(status_cols) <- status
gut_phylum_levels <- c(
  "p__Firmicutes", "p__Bacteroidota", "p__Actinobacteria", "p__Proteobacteria",
  "p__Verrucomicrobia", "p__Spirochaetes", "p__Euryarchaeota", "p__Uroviricota", "p__Apicomplexa", "p__Unclassified", "p__Missing"
)
oral_phylum_levels <- c(
  "p__Firmicutes", "p__Bacteroidota", "p__Actinobacteria", "p__Proteobacteria", "p__Fusobacteria",
  "p__Uroviricota", "p__Candidatus_Saccharibacteria", "p__Unclassified", "p__Missing"
)
mycols <- c("#D95F02", "#7570B3", "#E7298A", "#00FF00", "#E31A1C", "#E6AB02", "#00FFFF", "#FFFF99", "#1B9E77", "#6A3D9A", "#666666", "#A6761D", "blue", "black", "brown")
