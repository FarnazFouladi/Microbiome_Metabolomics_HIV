# Create network for a subset of features
rm(list = ls())
sample_type <- "Gut"
source("code/load.R")
source("code/network_function.R")
threshold <- 0.3
pval_t <- 0.01


############################################# Gut Species
name = "Species"
sub_dir <- file.path(output_dir,"Dysbiosis/Gut","Gut_Species","sub-network")
dir.create(sub_dir)
res_g1 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Gut_Species", "species_g1.rds"))
res_g2 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Gut_Species", "species_g2.rds"))
cor_diff <- read.csv(file.path(output_dir,"Dysbiosis/Gut","Gut_Species", "species_differential_correlations.csv"))
sig_diff <- cor_diff[cor_diff$diff_cor > threshold & cor_diff$p < pval_t, ]
sig_diff_hold <- sig_diff %>% filter(grepl("Holdem",Taxa)|grepl("Holdem",name) )
species <- unique(c(sig_diff_hold$Taxa,sig_diff_hold$name))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% species,colnames(g1_mat) %in% species]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% species,colnames(g2_mat) %in% species]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% species,colnames(g1_mat_cooccur) %in% species]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% species,colnames(g2_mat_cooccur) %in% species]
phylum_levels <- c(
  "p__Firmicutes", "p__Bacteroidota", "p__Actinobacteria", "p__Proteobacteria",
  "p__Verrucomicrobia", "p__Spirochaetes", "p__Euryarchaeota", "p__Uroviricota", "p__Apicomplexa", "p__Unclassified", "p__Missing"
)
names(phylum_levels) <- mycols[1:length(phylum_levels)]
taxonomies <- read_rds("data/Processed_Data/Gut/taxonomies.rds")

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_hold,
  #taxonomies = taxonomies,
  datasets = c("Gut Species","Gut Species"),
  t = "Species", phylum_levels = phylum_levels,
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.8,
  cexNodes = 1.5, cexNodes2 = 1,
  cexHubs = 2, cexHubs2 = 1.5,
  cexLabels = 4, cexLabels2 = 4,
  cexHubLabels = 5, cexHubLabels2 = 5.5,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "Holdemanella", showLabel = T
)


############################################# Gut Metabolites
name = "Gut_Metabolites"
sub_dir <- file.path(output_dir,"Dysbiosis/Gut","Gut_MTB_POS","sub-network")
res_g1 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Gut_MTB_POS", "mtb_pos_g1.rds"))
res_g2 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Gut_MTB_POS", "mtb_pos_g2.rds"))
sig_diff <- read.table(file.path(output_dir,"Dysbiosis/Gut","Gut_MTB_POS", "mtb_pos_differential_correlations_sig.txt"),
                       sep = "\t",quote = "",check.names = FALSE,comment.char = "",header = T)


# Bases
sig_diff_sub <- sig_diff %>% filter(grepl("Thymine|Xanthine|Hypoxanthine|Histamine|Sphingosine",Metabolites)) %>%
  mutate(Taxa = paste0("data1 - ",Taxa),
         metabolite_ID = paste0("data2 - ",metabolite_ID),
         Metabolites = paste0("data2 - ",Metabolites))
features <- unique(c(sig_diff_sub$Taxa,sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features,colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features,colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features,colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features,colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>% left_join(sig_diff_sub[,c("metabolite_ID","Metabolites")] %>% distinct(metabolite_ID,.keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites),metabolite_ID,Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub),metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub),names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub),names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub),names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub),names_mat$metabolite_ID))

colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites

sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Gut MTB-POS"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.95,
  cexNodes = 1, cexNodes2 = 0.9,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 4,
  cexHubLabels = 2, cexHubLabels2 = 3,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "Bases_Histamine",showLabel = T
)

############################################# Plasma Metabolites
name = "Plasma_Metabolites"
sub_dir <- file.path(output_dir,"Dysbiosis/Gut","Plasma_MTB_NEG","sub-network")
res_g1 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Plasma_MTB_NEG", "plasma_mtb_neg_g1.rds"))
res_g2 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Plasma_MTB_NEG", "plasma_mtb_neg_g2.rds"))
sig_diff <- read.table(file.path(output_dir,"Dysbiosis/Gut","Plasma_MTB_NEG", "plasma_mtb_neg_differential_correlations_sig.txt"),
                       sep = "\t",quote = "",check.names = FALSE,comment.char = "",header = T)

bile_acids <- unique(sig_diff$Metabolites[grepl("Bile",sig_diff$sub_class)])
# Bile acids
sig_diff_sub <- sig_diff %>% filter(grepl("Bile",sub_class)) %>%
  mutate(Taxa = paste0("data1 - ",Taxa),
         metabolite_ID = paste0("data2 - ",metabolite_ID),
         Metabolites = paste0("data2 - ",Metabolites))
features <- unique(c(sig_diff_sub$Taxa,sig_diff_sub$metabolite_ID))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features,colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features,colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features,colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features,colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(metabolite_ID = colnames(g1_mat_sub)) %>% left_join(sig_diff_sub[,c("metabolite_ID","Metabolites")] %>% distinct(metabolite_ID,.keep_all = T)) %>%
  mutate(Metabolites = ifelse(is.na(Metabolites),metabolite_ID,Metabolites)) %>%
  arrange(order(match(colnames(g1_mat_sub),metabolite_ID)))
stopifnot(all.equal(colnames(g1_mat_sub),names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_sub),names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub),names_mat$metabolite_ID))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub),names_mat$metabolite_ID))


colnames(g1_mat_sub) <- names_mat$Metabolites
rownames(g1_mat_sub) <- names_mat$Metabolites
colnames(g2_mat_sub) <- names_mat$Metabolites
rownames(g2_mat_sub) <- names_mat$Metabolites
colnames(g1_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g1_mat_cooccur_sub) <- names_mat$Metabolites
colnames(g2_mat_cooccur_sub) <- names_mat$Metabolites
rownames(g2_mat_cooccur_sub) <- names_mat$Metabolites
sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Metabolites)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Plasma MTB-NEG"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.7,
  cexNodes = 1, cexNodes2 = 1,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 4,
  cexHubLabels = 2, cexHubLabels2 = 3,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "Bile_derivatives_neg",showLabel = T
)

############################################# Gut GO
name = "GO"
sub_dir <- file.path(output_dir,"Dysbiosis/Gut","Gut_GO","sub-network")
res_g1 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Gut_GO", "go_g1.rds"))
res_g2 <- readRDS(file.path(output_dir,"Dysbiosis/Gut","Gut_GO", "go_g2.rds"))
sig_diff <- read.table(file.path(output_dir,"Dysbiosis/Gut","Gut_GO", "go_differential_correlations_sig.txt"),
                       sep = "\t",quote = "",check.names = FALSE,comment.char = "",header = T)

sig_diff <- sig_diff %>% group_by(Description) %>%
  mutate(num = n()) 

degree <- sig_diff %>% distinct(Description,.keep_all = T) %>% arrange(desc(num))

# carbamoyl-phosphate synthase (glutamine-hydrolyzing) activity
# glutamine biosynthesis process
# glutamine-tRNA-ligase activity
sig_diff_sub <- sig_diff %>% filter(grepl("0004088|0006542|0004819",GO)) %>%
  mutate(Taxa = paste0("data1 - ",Taxa),
         GO = paste0("data2 - ",GO),
         Description = paste0("data2 - ",Description))
features <- unique(c(sig_diff_sub$Taxa,sig_diff_sub$GO))
g1_mat <- res_g1$corr
g1_mat_cooccur <- res_g1$mat_cooccur
g2_mat <- res_g2$corr
g2_mat_cooccur <- res_g2$mat_cooccur

g1_mat_sub <- g1_mat[rownames(g1_mat) %in% features,colnames(g1_mat) %in% features]
g2_mat_sub <- g2_mat[rownames(g2_mat) %in% features,colnames(g2_mat) %in% features]
g1_mat_cooccur_sub <- g1_mat_cooccur[rownames(g1_mat_cooccur) %in% features,colnames(g1_mat_cooccur) %in% features]
g2_mat_cooccur_sub <- g2_mat_cooccur[rownames(g2_mat_cooccur) %in% features,colnames(g2_mat_cooccur) %in% features]

names_mat <- data.frame(GO = colnames(g1_mat_sub)) %>% left_join(sig_diff_sub[,c("GO","Description")] %>% distinct(GO,.keep_all = T)) %>%
  mutate(Description = ifelse(is.na(Description),GO,Description)) %>%
  arrange(order(match(colnames(g1_mat_sub),GO)))
stopifnot(all.equal(colnames(g1_mat_sub),names_mat$GO))
stopifnot(all.equal(colnames(g2_mat_sub),names_mat$GO))
stopifnot(all.equal(colnames(g1_mat_cooccur_sub),names_mat$GO))
stopifnot(all.equal(colnames(g2_mat_cooccur_sub),names_mat$GO))


colnames(g1_mat_sub) <- names_mat$Description
rownames(g1_mat_sub) <- names_mat$Description
colnames(g2_mat_sub) <- names_mat$Description
rownames(g2_mat_sub) <- names_mat$Description
colnames(g1_mat_cooccur_sub) <- names_mat$Description
rownames(g1_mat_cooccur_sub) <- names_mat$Description
colnames(g2_mat_cooccur_sub) <- names_mat$Description
rownames(g2_mat_cooccur_sub) <- names_mat$Description
sig_diff_sub1 <- data.frame(Taxa = sig_diff_sub$Taxa, name = sig_diff_sub$Description)

create_networks(
  g1_mat = g1_mat_sub, g1_mat_cooccur = g1_mat_cooccur_sub,
  g2_mat = g2_mat_sub, g2_mat_cooccur = g2_mat_cooccur_sub,
  sig_diff = sig_diff_sub1,
  taxonomies = NULL,
  t = NULL, datasets = c("Gut Species", "Gut GO"),
  threshold = threshold, pval_th = pval_t,
  layoutGroup = "union", layoutGroup2 = "2",
  repulsion = 0.85, repulsion2 = 0.95,
  cexNodes = 1, cexNodes2 = 0.8,
  cexHubs = 1.5, cexHubs2 = 1.5,
  cexLabels = 1, cexLabels2 = 3.5,
  cexHubLabels = 2, cexHubLabels2 = 6,
  cexTitle = 3.8, cexTitle2 = 5,
  legend_cex = 3, legend_cex2 = 5,
  association_cex = 3, association_cex2 = 5,
  group_names = c("", ""), outputdir = sub_dir, prefix = "Glutamine",showLabel = T
)






