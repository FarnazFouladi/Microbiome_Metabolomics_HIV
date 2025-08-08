# Create heatmaps for a subset of features
rm(list = ls())
source("code/load.R")
source("code/trend_heatmap.R")
input_dir <- file.path("Outputs", "Correlations_trend_analysis")
sample_type <- "Gut"
sub_dir <- file.path("Outputs","Differential_Correlations", "Networks",sample_type)
dir.create(sub_dir,recursive = T)
############################################# Gut Species
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("gut_species","_cor_diff_sig_with_trend_test.csv")),check.names = F)
sig_diff_f <- trend_res %>% filter(grepl("Holdem",t1)|grepl("Holdem",t2) )
sig_diff_ordered <- sig_diff_f
sig_diff_ordered$t1[grepl("Holdem",sig_diff_f$t2)] <- sig_diff_f$t2[grepl("Holdem",sig_diff_f$t2)]
sig_diff_ordered$t2[grepl("Holdem",sig_diff_f$t2)] <- sig_diff_f$t1[grepl("Holdem",sig_diff_f$t2)]

sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(t1,":",t2)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" ) 
trend_heatmap(mat1,annot,sub_dir,"Holdemanella")

############################################# Gut GO
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("gut_go","_cor_diff_sig_with_trend_test.csv")),check.names = F)
features <- c("GO:0004641","GO:0004088")
sig_diff_f <- trend_res %>% filter(t2 %in% features)
sig_diff_ordered <- sig_diff_f

sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(Description,":",t1)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"go")

############################################# Gut metabolites
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("gut_meta.p","_cor_diff_sig_with_trend_test.csv")),check.names = F)
sig_diff_f <- trend_res %>% filter(sub_class =="Pyrimidines and pyrimidine derivatives" | sub_class == 
                                    "Purines and purine derivatives")
sig_diff_ordered <- sig_diff_f
sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(Metabolites,":",t1)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"gut_mtb.p")

############################################# Gut metabolites2
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("gut_meta.p","_cor_diff_sig_with_trend_test.csv")),check.names = F)
sig_diff_f <- trend_res %>% filter(Metabolites %in% c("Histamine","Sphingosine(1+)"))
sig_diff_ordered <- sig_diff_f
sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(Metabolites,":",t1)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"gut_mtb.p_histamine")

############################################# Plasma metabolites
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("gut_plasma.meta.n","_cor_diff_sig_with_trend_test.csv")),check.names = F)
sig_diff_f <- trend_res %>% filter(grepl("Bile|bile",sub_class )) 
sig_diff_ordered <- sig_diff_f
sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(Metabolites,":",t1)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"gut_plasma.meta.n")

############################################ Oral GO
sample_type <- "Oral"
sub_dir <- file.path("Outputs","Differential_Correlations", "Networks",sample_type)
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("oral_go","_cor_diff_sig_with_trend_test.csv")),check.names = F)
sig_diff_f <- trend_res %>% filter(grepl("Haemophilus_parainfluenzae",t1))
sig_diff_ordered <- sig_diff_f

sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(t1,":",Description)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"oral_go")

############################################ Oral metabolites
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("oral_meta.n","_cor_diff_sig_with_trend_test.csv")),check.names = F)
features <- c("LysoPC(16:0/0:0)") 
sig_diff_f <- trend_res %>% filter(Metabolites %in% features) 
sig_diff_ordered <- sig_diff_f

sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(Metabolites,":",t1)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"oral_meta.n")

############################################ Plasma metabolites
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("oral_plasma.meta.n","_cor_diff_sig_with_trend_test.csv")),check.names = F)
features <- c("HMDB0062551")
sig_diff_f <- trend_res %>% filter(t2 %in% features) 
sig_diff_ordered <- sig_diff_f

sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(Metabolites,":",t1)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"oral_plasma.meta.n")

############################################ Oral gut metabolites
trend_res <- read.csv(file.path("Outputs/Correlations_trend_analysis",paste0("oral_gut_species","_cor_diff_sig_with_trend_test.csv")),check.names = F)
sig_diff_f <- trend_res %>% filter(grepl("gut_s__Roseburia|^gut_s__Eubacterium_|^gut_s__Bifidobacterium",t2))
sig_diff_ordered <- sig_diff_f

sig_diff_ordered <- sig_diff_ordered %>% arrange(`adj.p.trend`) %>% mutate(pairs = paste0(t1,":",t2)) %>%
  select(c("pairs","r1.g1" ,"r2.g2","r3.g3","r4.g4", "adj.p.trend"   )) %>% column_to_rownames("pairs")
mat1 <- sig_diff_ordered %>% select(c("r1.g1" ,"r2.g2","r3.g3","r4.g4"))
annot <- sig_diff_ordered %>% select("adj.p.trend" )

trend_heatmap(mat1,annot,sub_dir,"oral_gut_species")

