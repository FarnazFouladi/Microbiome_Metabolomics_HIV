# Create network for a subset of features
rm(list = ls())
sample_type <- "Gut"
source("code/load.R")
source("code/create_network.R")
input_dir <- file.path("Outputs", "Differential_Correlations","Tables")
sub_dir <- file.path("Outputs","Differential_Correlations", "Networks",sample_type)
dir.create(sub_dir,recursive = T)
############################################# Gut Species - Holdemanella spp.
delta_r <- read.csv(file.path(input_dir,"gut_species_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"gut_species_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"gut_species_bh.csv"),check.names = F,row.names = 1) 
colnames(adj.pval) <- colnames(pval)
rownames(adj.pval) <- rownames(pval)
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
sig_diff <- read.csv(file.path(input_dir,"gut_species_cor_diff_sig.csv"),check.names = F) 
sig_diff_f <- sig_diff %>% filter(grepl("Holdem",t1)|grepl("Holdem",t2) )
# Save source data
write.table(sig_diff_f %>% mutate(cor_diff = r2 - r1) %>%
  select(c("t1","t2","cor_diff","p","BH_p")),
  file = file.path(sub_dir,"Holdemanella_source_data.txt"),sep = "\t", row.names = F, quote = F)

species <- unique(c(sig_diff_f$t1,sig_diff_f$t2))
net.r <- delta_r[species,species]
adj.pval <- adj.pval[species,species]
adj.pval[is.na(adj.pval)] <- 1

colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))

create_network(
  mat = net.r,pval_mat = adj.pval,
  datasets = c("Gut Species","Gut Species"),
  repulsion = 0.83, 
  cexNodes = 2, 
  cexHubs = 3, 
  cexLabels = 0.5, 
  cexHubLabels = 0.8, 
  outputdir = sub_dir, 
  prefix = "Holdemanella", 
)

############################################# Gut GO, enzymes involved in purine and pyrimidine biosynthesis pathways 
delta_r <- read.csv(file.path(input_dir,"gut_go_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"gut_go_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"gut_go_bh.csv"),check.names = F,row.names = 1) 
colnames(adj.pval) <- colnames(pval)
rownames(adj.pval) <- rownames(pval)
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0

sig_diff <- read.csv(file.path(input_dir,"gut_go_cor_diff_sig.csv"),check.names = F) 
features <- c("GO:0004641","GO:0004088")
sig_diff_f <- sig_diff %>% filter(t2 %in% features) %>% arrange(t2)
# Save source data
write.table(sig_diff_f %>% mutate(cor_diff = r2 - r1) %>%
              select(c("t2","Description","t1","cor_diff","p","BH_p")),
            file = file.path(sub_dir,"go_source_data.txt"),sep = "\t", row.names = F, quote = F)

annotation_go <- sig_diff_f %>% select(c("t2","Description")) %>% distinct()
net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]
adj.pval <- adj.pval[rownames(net.r),colnames(net.r)]
#colnames(net.r)[1:length(features)] <- annotation_go$Description
#rownames(net.r)[1:length(features)] <- annotation_go$Description
colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))
adj.pval[is.na(adj.pval)] <- 1

create_network(
  mat = as.matrix(net.r),pval_mat = adj.pval,
  datasets = c("Gut Species","Gut GO"),
  repulsion = 0.85, 
  cexNodes = 2, 
  cexHubs = 3, 
  cexLabels = 0.5, 
  cexHubLabels = 0.8, 
  outputdir = sub_dir, 
  prefix = "go", 
)

############################################# Gut Metabolites, purine and pyrimidine derivatives 
delta_r <- read.csv(file.path(input_dir,"gut_meta.p_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"gut_meta.p_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"gut_meta.p_bh.csv"),check.names = F,row.names = 1) 
colnames(adj.pval) <- colnames(pval)
rownames(adj.pval) <- rownames(pval)
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0
sig_diff <- read.csv(file.path(input_dir,"gut_meta.p_cor_diff_sig.csv"),check.names = F) 
sig_diff_f <- sig_diff %>% filter(sub_class =="Pyrimidines and pyrimidine derivatives" | sub_class == 
                                    "Purines and purine derivatives")
# Save source data
write.table(sig_diff_f %>% mutate(cor_diff = r2 - r1) %>%
              select(c("Metabolites","t1","cor_diff","p","BH_p")),
            file = file.path(sub_dir,"gut_mtb.p_source_data.txt"),sep = "\t", row.names = F, quote = F)

annotation <- sig_diff_f %>% select(c("t2","Metabolites")) %>% distinct()
features <- unique(c(sig_diff_f$t2))
net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]
adj.pval <- adj.pval[rownames(net.r),colnames(net.r)]
colnames(net.r)[1:length(features)] <- annotation$Metabolites
rownames(net.r)[1:length(features)] <- annotation$Metabolites
colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))
adj.pval[is.na(adj.pval)] <- 1

create_network(
  mat = net.r,pval_mat = adj.pval,
  datasets = c("Gut Species","Gut MTB-P"),
  repulsion = 0.85, 
  cexNodes = 2, 
  cexHubs = 3, 
  cexLabels = 0.5, 
  cexHubLabels = 0.8, 
  outputdir = sub_dir, 
  prefix = "gut_mtb.p", 
)
############################################# Gut Metabolites2
delta_r <- read.csv(file.path(input_dir,"gut_meta.p_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"gut_meta.p_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"gut_meta.p_bh.csv"),check.names = F,row.names = 1) 
colnames(adj.pval) <- colnames(pval)
rownames(adj.pval) <- rownames(pval)
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0


sig_diff <- read.csv(file.path(input_dir,"gut_meta.p_cor_diff_sig.csv"),check.names = F) 
sig_diff_f <- sig_diff %>% filter(Metabolites %in% c("Histamine","Sphingosine(1+)"))
annotation <- sig_diff_f %>% select(c("t2","Metabolites")) %>% distinct()
features <- unique(c(sig_diff_f$t2))
net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]
adj.pval <- adj.pval[rownames(net.r),colnames(net.r)]
colnames(net.r)[1:length(features)] <- annotation$Metabolites
rownames(net.r)[1:length(features)] <- annotation$Metabolites


colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))
adj.pval[is.na(adj.pval)] <- 1

create_network(
  mat = net.r,pval_mat = adj.pval,
  datasets = c("Gut Species","Gut MTB-P"),
  repulsion = 0.85, 
  cexNodes = 2, 
  cexHubs = 3, 
  cexLabels = 0.5, 
  cexHubLabels = 0.8, 
  outputdir = sub_dir, 
  prefix = "gut_mtb.p_histamine", 
)

############################################# Plasma Metabolites, bile acids
delta_r <- read.csv(file.path(input_dir,"gut_plasma.meta.n_delta_cor.csv"),check.names = F,row.names = 1) 
pval <- read.csv(file.path(input_dir,"gut_plasma.meta.n_p.csv"),check.names = F,row.names = 1) 
adj.pval <- read.csv(file.path(input_dir,"gut_plasma.meta.n_bh.csv"),check.names = F,row.names = 1) 
colnames(adj.pval) <- colnames(pval)
rownames(adj.pval) <- rownames(pval)
adj.pval[lower.tri(adj.pval)] = t(adj.pval)[lower.tri(adj.pval)]
delta_r[is.na(pval) | !(pval<0.01 & adj.pval < 0.1 & abs(delta_r)  > 0.3)] <- 0
p <- grep("s__",rownames(delta_r))
q<- (length(p)+1):nrow(delta_r)
delta_r[p,p] <-0
delta_r[q,q] <-0

sig_diff <- read.csv(file.path(input_dir,"gut_plasma.meta.n_cor_diff_sig.csv"),check.names = F) 
sig_diff_f <- sig_diff %>% filter(grepl("Bile|bile",sub_class )) %>% arrange(t2)
# Save source data
write.table(sig_diff_f %>% mutate(cor_diff = r2 - r1) %>%
              select(c("Metabolites","t1","cor_diff","p","BH_p")),
            file = file.path(sub_dir,"plasma_mtb.n_source_data.txt"),sep = "\t", row.names = F, quote = F)

features <- unique(sig_diff_f$Metabolites)
annotation <- sig_diff_f %>% select(c("t2","Metabolites")) %>% distinct()
net.r <- delta_r[unique(c(sig_diff_f$t2,sig_diff_f$t1)),unique(c(sig_diff_f$t2,sig_diff_f$t1))]
adj.pval <- adj.pval[rownames(net.r),colnames(net.r)]
colnames(net.r)[1:length(features)] <- annotation$Metabolites
rownames(net.r)[1:length(features)] <- annotation$Metabolites

colnames(net.r) <- gsub("s__","taxa ",colnames(net.r))
rownames(net.r) <- gsub("s__","taxa ",rownames(net.r))
colnames(net.r) <- gsub("_"," ",colnames(net.r))
rownames(net.r) <- gsub("_"," ",rownames(net.r))
adj.pval[is.na(adj.pval)] <- 1

create_network(
  mat = net.r,pval_mat = adj.pval,
  datasets = c("Gut Species","Plasma MTB-N"),
  repulsion = 0.75, 
  cexNodes = 2, 
  cexHubs = 3, 
  cexLabels = 0.5, 
  cexHubLabels = 0.8, 
  outputdir = sub_dir, 
  prefix = "plasma_mtb.n", 
)

