# Create network for a subset of features
rm(list = ls())
sample_type <- "Oral"
source("code/load.R")
source("code/create_network.R")
input_dir <- file.path("Outputs", "Differential_Correlations","Tables")
sub_dir <- file.path("Outputs","Differential_Correlations",  "Networks",sample_type)
dir.create(sub_dir,recursive = T)

############################################# Oral GO, Haemophilus_parainfluenzae
delta_r <- read.csv(file.path(input_dir,"oral_go_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"oral_go_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"oral_go_bh.csv"),check.names = F,row.names = 1) 
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0

sig_diff <- read.csv(file.path(input_dir,"oral_go_cor_diff_sig.csv"),check.names = F) 
sig_diff_f <- sig_diff %>% filter(grepl("Haemophilus_parainfluenzae",t1)) %>% arrange(t2)
annotation_go <- sig_diff_f %>% select(c("t2","Description")) %>% distinct()
features <- annotation_go$t2
net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]
colnames(net.r)[1:length(features)] <- annotation_go$Description
rownames(net.r)[1:length(features)] <- annotation_go$Description
colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))

create_network(
  mat = net.r,
  datasets = c("Oral Species","Oral GO"),
  repulsion = 0.85, 
  cexNodes = 1, 
  cexHubs = 15, 
  cexLabels = 1, 
  cexHubLabels = 15, 
  cexTitle = 3.8, 
  legend_cex = 3, 
  association_cex = 3, 
  outputdir = sub_dir, 
  prefix = "oral_go", 
)

############################################# Oral metabolites, LysoPC(16:0/0:0)
delta_r <- read.csv(file.path(input_dir,"oral_meta.n_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"oral_meta.n_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"oral_meta.n_bh.csv"),check.names = F,row.names = 1) 
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0

sig_diff <- read.csv(file.path(input_dir,"oral_meta.n_cor_diff_sig.csv"),check.names = F) 
features <- c("LysoPC(16:0/0:0)") 
sig_diff_f <- sig_diff %>% filter(Metabolites %in% features) %>% arrange(t2)
annotation <- sig_diff_f %>% select(c("t2","Metabolites")) %>% distinct()
net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]
colnames(net.r)[1:length(features)] <- annotation$Metabolites
rownames(net.r)[1:length(features)] <- annotation$Metabolites

colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))

create_network(
  mat = net.r,
  datasets = c("Oral Species","Oral MTB-N"),
  repulsion = 0.85, 
  cexNodes = 4, 
  cexHubs = 5, 
  cexLabels = 6, 
  cexHubLabels = 7, 
  cexTitle = 3.8, 
  legend_cex = 3, 
  association_cex = 3, 
  outputdir = sub_dir, 
  prefix = "oral_mtb.n", 
)

############################################# Plasma Metabolites
delta_r <- read.csv(file.path(input_dir,"oral_plasma.meta.n_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"oral_plasma.meta.n_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"oral_plasma.meta.n_bh.csv"),check.names = F,row.names = 1) 
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0

features <- c("HMDB0062551")
sig_diff <- read.csv(file.path(input_dir,"oral_plasma.meta.n_cor_diff_sig.csv"),check.names = F) 
sig_diff_f <- sig_diff %>% filter(t2 %in% features) %>% arrange(t2)
annotation <- sig_diff_f %>% select(c("t2","Metabolites")) %>% distinct()
net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]
colnames(net.r)[1:length(features)] <- annotation$Metabolites
rownames(net.r)[1:length(features)] <- annotation$Metabolites

colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))

create_network(
  mat = net.r,
  datasets = c("Oral Species","Plasma MTB-N"),
  repulsion = 0.9, 
  cexNodes = 4, 
  cexHubs = 5, 
  cexLabels = 7, 
  cexHubLabels = 8, 
  cexTitle = 3.8, 
  legend_cex = 3, 
  association_cex = 3, 
  outputdir = sub_dir, 
  prefix = "oral_plasma_mtb.n", 
)

############################################# Oral Species and Gut Species
delta_r <- read.csv(file.path(input_dir,"oral_gut_species_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"oral_gut_species_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"oral_gut_species_bh.csv"),check.names = F,row.names = 1) 
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("^s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0
sig_diff <- read.csv(file.path(input_dir,"oral_gut_species_cor_diff_sig.csv"),check.names = F) 
sig_diff_f <- sig_diff %>% filter(grepl("gut_s__Roseburia|^gut_s__Eubacterium|^gut_s__Bifidobacterium",t2) |
                                    grepl("s__Streptococcus|^s__Porphyromonas|^s__Prevotella",t1) )

net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]


colnames(net.r)[grep("^s__",colnames(net.r))] <- gsub("s__","taxa ",colnames(net.r)[grep("^s__",colnames(net.r))])
rownames(net.r)[grep("^s__",rownames(net.r))] <- gsub("s__","taxa ",rownames(net.r)[grep("^s__",rownames(net.r))])
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))

create_network(
  mat = net.r,
  datasets = c("Oral Species","Gut Species"),
  repulsion = 0.88, 
  cexNodes = 2.5, 
  cexHubs = 3.5, 
  cexLabels = 7, 
  cexHubLabels = 8, 
  cexTitle = 3.8, 
  legend_cex = 3, 
  association_cex = 3, 
  outputdir = sub_dir, 
  prefix = "oral_gut_species", 
)

