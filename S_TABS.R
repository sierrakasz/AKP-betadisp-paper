## Anna Karenina Kaszubinski et al.

#Supplemental Tables

#packages
library(BaylorEdPsych)
library(car)
library(foreign)
library(GGally)
library(ggforce)
library(ggfortify)
library(ggplot2)
library(grid)
library(ggpubr)
library(ggthemes)
library(lme4)
library(phyloseq)
library(plyr)
library(PMCMR)
library(microbiome)
library(mlogit)
library(nnet)
library(randomForest)
library(reshape2)
library(rsample)
library(tidyverse)
library(vegan)
#for reproducibilty
set.seed(1234)


# Load data ---------------------------------------------------------------
#otu table, taxonomy table, and tree file
otu <- read.csv("table.csv")
tax <- read.csv("taxonomy.csv")
tree <- read_tree('tree.nwk')

#load in metadata
metadata=(read.csv("HPMMMetadata_pt2_cod.csv",header=TRUE))
metadata$BMI <- as.numeric(metadata$BMI)
#format it into phyloseq format
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

#put otu table in phyloseq format
rownames(otu) <- otu$OTUID
otu <- otu[,-1]
OTU=otu_table(otu, taxa_are_rows=TRUE)
#merge tree and otu tables
physeq_otu.tree=phyloseq(OTU,tree, sampdat)

#put taxonomy table in phyloseq format
rownames(tax) <- tax$OTUID
tax_reduc <- merge(tax, otu, by = "row.names")
rownames(tax_reduc) <- tax_reduc$Row.names
tax_reduc <- tax_reduc[,-1]
tax_f <- tax_reduc[,1:7]
tax_f <- as.matrix(tax_f)
TAX=tax_table(tax_f)
taxa_names(TAX)=row.names(OTU)

#merge it all together into one phyloseq object
physeq <- merge_phyloseq(physeq_otu.tree, TAX)
#30906 taxa, 878 samples

#triming out taxa that are not representative of .01% of mean sequence number
physeq_trim <- prune_taxa(taxa_sums(physeq) > sum(otu) *.001 / 878, physeq)
#8692 taxa

# Table S1 ----------------------------------------------------------------

#table was created in excel


# Table S2 ----------------------------------------------------------------

#table created in excel


# Table S3 ----------------------------------------------------------------
#rarefy
#minimum library sizes: 3,000, 5,000, 7,000
physeq_3000 <- rarefy_even_depth(physeq_trim, sample.size = 3000)
#removed 16 samples, 90 OTUs
#8621 taxa, 862 samples
physeq_5000 <- rarefy_even_depth(physeq_trim, sample.size = 5000)
#removed 25 samples, 33 OTUs
#8664 taxa, 853 samples
physeq_7000 <- rarefy_even_depth(physeq_trim, sample.size = 7000)
#remove 45 samples, 83 OTUs
#8611 taxa, 833 samples

#normalize
# removing taxa not present in at least a certain percentage of samples
#cut offs: 1%, 3%, 10%
normalize_wout_rarefying <- function(physeq, level) {
  physeq_norm <- genefilter_sample(physeq, filterfun_sample(function(x) x >= 1), 
                                   A = level*nsamples(physeq_trim))
  ps_filtered <- prune_taxa(physeq_norm, physeq)
  
  return(ps_filtered)
}

physeq_1per <- normalize_wout_rarefying(physeq, .01)
#1500 taxa, 878 samples
physeq_3per <- normalize_wout_rarefying(physeq, .03)
#643 taxa, 878 samples
physeq_10per <- normalize_wout_rarefying(physeq, .1)
#216 taxa, 878 samples

#phy list - phyloseq objects which have rarefaction/normalization levels
phy_list <- list(physeq_3000, physeq_5000, physeq_7000,
                 physeq_1per, physeq_3per, physeq_10per)

names(phy_list) <- c('physeq_3000', 'physeq_5000', 'physeq_7000', 'physeq_1per',
                     'physeq_3per', 'physeq_10per')

#add this column to keep track of objects as they go into lists
sample_data(phy_list[[1]])$Method <- 'Rare_3000'
sample_data(phy_list[[2]])$Method <- 'Rare_5000'
sample_data(phy_list[[3]])$Method <- 'Rare_7000'
sample_data(phy_list[[4]])$Method <- 'Norm_1per'
sample_data(phy_list[[5]])$Method <- 'Norm_3per'
sample_data(phy_list[[6]])$Method <- 'Norm_10per'

#samp_area_list - all sample areas as unique names
samp_area_list <- c('Nose', 'Rectum', 'Ears', 'Mouth', 'Eyes')

#compare subset - list with all rarefaction/normalization level phyloseq objects
# subsetted by sample area
#total of 30 phyloseq objects
compare_subset <- vector('list')
for(i in 1:length(phy_list)) {
  for(j in 1:length(samp_area_list)) {
    a <- (phyloseq::subset_samples(phy_list[[i]], Sample_Area == samp_area_list[j]))
    compare_subset <- append(compare_subset, a)
  }
}

#changes names for easier access
names(compare_subset) <- c('physeq_3000_nos', 'physeq_3000_rec', 'physeq_3000_ear',
                           'physeq_3000_mou', 'physeq_3000_eye',
                           'physeq_5000_nos', 'physeq_5000_rec', 'physeq_5000_ear',
                           'physeq_5000_mou', 'physeq_5000_eye',
                           'physeq_7000_nos', 'physeq_7000_rec', 'physeq_7000_ear',
                           'physeq_7000_mou', 'physeq_7000_eye',
                           'physeq_1per_nos', 'physeq_1per_rec', 'physeq_1per_ear',
                           'physeq_1per_mou', 'physeq_1per_eye',
                           'physeq_3per_nos', 'physeq_3per_rec', 'physeq_3per_ear',
                           'physeq_3per_mou', 'physeq_3per_eye',
                           'physeq_10per_nos', 'physeq_10per_rec', 'physeq_10per_ear',
                           'physeq_10per_mou', 'physeq_10per_eye')


#function for taking subsetted phyloseq objects and finding corresponding 
#unifrac or weighted unifrac distances 
compare_dist_u <- vector('list')
compare_dist_w <- vector('list')
for(j in 1:length(compare_subset)) {
  a <- phyloseq::distance(compare_subset[[j]], "unifrac")
  compare_dist_u[[j]] <- a
  b <- phyloseq::distance(compare_subset[[j]], "wunifrac")
  compare_dist_w[[j]] <- b
}

#changes names for easier access
names(compare_dist_u) <- c('GPdist_3000_nos', 'GPdist_3000_rec', 'GPdist_3000_ear',
                           'GPdist_3000_mou', 'GPdist_3000_eye',
                           'GPdist_5000_nos', 'GPdist_5000_rec', 'GPdist_5000_ear',
                           'GPdist_5000_mou', 'GPdist_5000_eye',
                           'GPdist_7000_nos', 'GPdist_7000_rec', 'GPdist_7000_ear',
                           'GPdist_7000_mou', 'GPdist_7000_eye',
                           'GPdist_1per_nos', 'GPdist_1per_rec', 'GPdist_1per_ear',
                           'GPdist_1per_mou', 'GPdist_1per_eye',
                           'GPdist_3per_nos', 'GPdist_3per_rec', 'GPdist_3per_ear',
                           'GPdist_3per_mou', 'GPdist_3per_eye',
                           'GPdist_10per_nos', 'GPdist_10per_rec', 'GPdist_10per_ear',
                           'GPdist_10per_mou', 'GPdist_10per_eye')

names(compare_dist_w) <- c('GPdist_w_3000_nos', 'GPdist_w_3000_rec', 'GPdist_w_3000_ear',
                           'GPdist_w_3000_mou', 'GPdist_w_3000_eye',
                           'GPdist_w_5000_nos', 'GPdist_w_5000_rec', 'GPdist_w_5000_ear',
                           'GPdist_w_5000_mou', 'GPdist_w_5000_eye',
                           'GPdist_w_7000_nos', 'GPdist_w_7000_rec', 'GPdist_w_7000_ear',
                           'GPdist_w_7000_mou', 'GPdist_w_7000_eye',
                           'GPdist_w_1per_nos', 'GPdist_w_1per_rec', 'GPdist_w_1per_ear',
                           'GPdist_w_1per_mou', 'GPdist_w_1per_eye',
                           'GPdist_w_3per_nos', 'GPdist_w_3per_rec', 'GPdist_w_3per_ear',
                           'GPdist_w_3per_mou', 'GPdist_w_3per_eye',
                           'GPdist_w_10per_nos', 'GPdist_w_10per_rec', 'GPdist_w_10per_ear',
                           'GPdist_w_10per_mou', 'GPdist_w_10per_eye')

#make sample dataframes used for beta-dispersion function 
sample_df_list <- vector('list')
for(i in 1:length(compare_subset)) {
  df <- data.frame(sample_data(compare_subset[[i]]))
  sample_df_list[[i]] <- df 
}

names(sample_df_list) <- c('sampledf_3000_nos', 'sampledf_3000_rec', 'sampledf_3000_ear',
                           'sampledf_3000_mou', 'sampledf_3000_eye',
                           'sampledf_5000_nos', 'sampledf_5000_rec', 'sampledf_5000_ear',
                           'sampledf_5000_mou', 'sampledf_5000_eye',
                           'sampledf_7000_nos', 'sampledf_7000_rec', 'sampledf_7000_ear',
                           'sampledf_7000_mou', 'sampledf_7000_eye',
                           'sampledf_1per_nos', 'sampledf_1per_rec', 'sampledf_1per_ear',
                           'sampledf_1per_mou', 'sampledf_1per_eye',
                           'sampledf_3per_nos', 'sampledf_3per_rec', 'sampledf_3per_ear',
                           'sampledf_3per_mou', 'sampledf_3per_eye',
                           'sampledf_10per_nos', 'sampledf_10per_rec', 'sampledf_10per_ear',
                           'sampledf_10per_mou', 'sampledf_10per_eye')


#calculate beta-dispersion for each sample area, distance metric, MoD/CoD,
# and rarefaction/normalization level. 
beta_disp_list_u_mod <- vector('list')
beta_disp_list_w_mod <- vector('list')
beta_disp_list_u_cod <- vector('list')
beta_disp_list_w_cod <- vector('list')

for(k in 1:30) {
  a <- betadisper(compare_dist_u[[k]], sample_df_list[[k]]$MoD)
  beta_disp_list_u_mod[[k]] <- a
  b <- betadisper(compare_dist_w[[k]], sample_df_list[[k]]$MoD)
  beta_disp_list_w_mod[[k]] <- b
  c <- betadisper(compare_dist_u[[k]], sample_df_list[[k]]$CoD_Simple2)
  beta_disp_list_u_cod[[k]] <- c
  d <- betadisper(compare_dist_w[[k]], sample_df_list[[k]]$CoD_Simple2)
  beta_disp_list_w_cod[[k]] <- d
}

names(beta_disp_list_u_mod) <- c('betadisp_mod_3000_nos', 'betadisp_mod_3000_rec', 'betadisp_mod_3000_ear',
                                 'betadisp_mod_3000_mou', 'betadisp_mod_3000_eye',
                                 'betadisp_mod_5000_nos', 'betadisp_mod_5000_rec', 'betadisp_mod_5000_ear',
                                 'betadisp_mod_5000_mou', 'betadisp_mod_5000_eye',
                                 'betadisp_mod_7000_nos', 'betadisp_mod_7000_rec', 'betadisp_mod_7000_ear',
                                 'betadisp_mod_7000_mou', 'betadisp_mod_7000_eye',
                                 'betadisp_mod_1per_nos', 'betadisp_mod_1per_rec', 'betadisp_mod_1per_ear',
                                 'betadisp_mod_1per_mou', 'betadisp_mod_1per_eye',
                                 'betadisp_mod_3per_nos', 'betadisp_mod_3per_rec', 'betadisp_mod_3per_ear',
                                 'betadisp_mod_3per_mou', 'betadisp_mod_3per_eye',
                                 'betadisp_mod_10per_nos', 'betadisp_mod_10per_rec', 'betadisp_mod_10per_ear',
                                 'betadisp_mod_10per_mou', 'betadisp_mod_10per_eye')

names(beta_disp_list_w_mod) <- c('betadisp_w_mod_3000_nos', 'betadisp_w_mod_3000_rec', 'betadisp_w_mod_3000_ear',
                                 'betadisp_w_mod_3000_mou', 'betadisp_w_mod_3000_eye',
                                 'betadisp_w_mod_5000_nos', 'betadisp_w_mod_5000_rec', 'betadisp_w_mod_5000_ear',
                                 'betadisp_w_mod_5000_mou', 'betadisp_w_mod_5000_eye',
                                 'betadisp_w_mod_7000_nos', 'betadisp_w_mod_7000_rec', 'betadisp_w_mod_7000_ear',
                                 'betadisp_w_mod_7000_mou', 'betadisp_w_mod_7000_eye',
                                 'betadisp_w_mod_1per_nos', 'betadisp_w_mod_1per_rec', 'betadisp_w_mod_1per_ear',
                                 'betadisp_w_mod_1per_mou', 'betadisp_w_mod_1per_eye',
                                 'betadisp_w_mod_3per_nos', 'betadisp_w_mod_3per_rec', 'betadisp_w_mod_3per_ear',
                                 'betadisp_w_mod_3per_mou', 'betadisp_w_mod_3per_eye',
                                 'betadisp_w_mod_10per_nos', 'betadisp_w_mod_10per_rec', 'betadisp_w_mod_10per_ear',
                                 'betadisp_w_mod_10per_mou', 'betadisp_w_mod_10per_eye')

names(beta_disp_list_u_cod) <- c('betadisp_cod_3000_nos', 'betadisp_cod_3000_rec', 'betadisp_cod_3000_ear',
                                 'betadisp_cod_3000_mou', 'betadisp_cod_3000_eye',
                                 'betadisp_cod_5000_nos', 'betadisp_cod_5000_rec', 'betadisp_cod_5000_ear',
                                 'betadisp_cod_5000_mou', 'betadisp_cod_5000_eye',
                                 'betadisp_cod_7000_nos', 'betadisp_cod_7000_rec', 'betadisp_cod_7000_ear',
                                 'betadisp_cod_7000_mou', 'betadisp_cod_7000_eye',
                                 'betadisp_cod_1per_nos', 'betadisp_cod_1per_rec', 'betadisp_cod_1per_ear',
                                 'betadisp_cod_1per_mou', 'betadisp_cod_1per_eye',
                                 'betadisp_cod_3per_nos', 'betadisp_cod_3per_rec', 'betadisp_cod_3per_ear',
                                 'betadisp_cod_3per_mou', 'betadisp_cod_3per_eye',
                                 'betadisp_cod_10per_nos', 'betadisp_cod_10per_rec', 'betadisp_cod_10per_ear',
                                 'betadisp_cod_10per_mou', 'betadisp_cod_10per_eye')

names(beta_disp_list_w_cod) <- c('betadisp_w_cod_3000_nos', 'betadisp_w_cod_3000_rec', 'betadisp_w_cod_3000_ear',
                                 'betadisp_w_cod_3000_mou', 'betadisp_w_cod_3000_eye',
                                 'betadisp_w_cod_5000_nos', 'betadisp_w_cod_5000_rec', 'betadisp_w_cod_5000_ear',
                                 'betadisp_w_cod_5000_mou', 'betadisp_w_cod_5000_eye',
                                 'betadisp_w_cod_7000_nos', 'betadisp_w_cod_7000_rec', 'betadisp_w_cod_7000_ear',
                                 'betadisp_w_cod_7000_mou', 'betadisp_w_cod_7000_eye',
                                 'betadisp_w_cod_1per_nos', 'betadisp_w_cod_1per_rec', 'betadisp_w_cod_1per_ear',
                                 'betadisp_w_cod_1per_mou', 'betadisp_w_cod_1per_eye',
                                 'betadisp_w_cod_3per_nos', 'betadisp_w_cod_3per_rec', 'betadisp_w_cod_3per_ear',
                                 'betadisp_w_cod_3per_mou', 'betadisp_w_cod_3per_eye',
                                 'betadisp_w_cod_10per_nos', 'betadisp_w_cod_10per_rec', 'betadisp_w_cod_10per_ear',
                                 'betadisp_w_cod_10per_mou', 'betadisp_w_cod_10per_eye')

#correctly format beta-dispersion values into dataframe to do statistical tests with 
get_beta_formated_u_mod <- vector('list')
get_beta_formated_w_mod <- vector('list')
get_beta_formated_u_cod <- vector('list')
get_beta_formated_w_cod <- vector('list')

for(i in 1:length(beta_disp_list_u_mod)) {
  beta_obj <- as.data.frame(beta_disp_list_u_mod[[i]]$distances)
  beta_obj$SampleID <- rownames(beta_obj)
  colnames(beta_obj) <- c('distances', 'SampleID')
  beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 
  get_beta_formated_u_mod[[i]] <- beta_sa
  
  beta_obj <- as.data.frame(beta_disp_list_w_mod[[i]]$distances)
  beta_obj$SampleID <- rownames(beta_obj)
  colnames(beta_obj) <- c('distances', 'SampleID')
  beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 
  get_beta_formated_w_mod[[i]] <- beta_sa
  
  beta_obj <- as.data.frame(beta_disp_list_u_cod[[i]]$distances)
  beta_obj$SampleID <- rownames(beta_obj)
  colnames(beta_obj) <- c('distances', 'SampleID')
  beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 
  get_beta_formated_u_cod[[i]] <- beta_sa
  
  beta_obj <- as.data.frame(beta_disp_list_w_cod[[i]]$distances)
  beta_obj$SampleID <- rownames(beta_obj)
  colnames(beta_obj) <- c('distances', 'SampleID')
  beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 
  get_beta_formated_w_cod[[i]] <- beta_sa
}

names(get_beta_formated_u_mod) <- c('betavalues_mod_3000_nos', 'betavalues_mod_3000_rec', 'betavalues_mod_3000_ear',
                                    'betavalues_mod_3000_mou', 'betavalues_mod_3000_eye',
                                    'betavalues_mod_5000_nos', 'betavalues_mod_5000_rec', 'betavalues_mod_5000_ear',
                                    'betavalues_mod_5000_mou', 'betavalues_mod_5000_eye',
                                    'betavalues_mod_7000_nos', 'betavalues_mod_7000_rec', 'betavalues_mod_7000_ear',
                                    'betavalues_mod_7000_mou', 'betavalues_mod_7000_eye',
                                    'betavalues_mod_1per_nos', 'betavalues_mod_1per_rec', 'betavalues_mod_1per_ear',
                                    'betavalues_mod_1per_mou', 'betavalues_mod_1per_eye',
                                    'betavalues_mod_3per_nos', 'betavalues_mod_3per_rec', 'betavalues_mod_3per_ear',
                                    'betavalues_mod_3per_mou', 'betavalues_mod_3per_eye',
                                    'betavalues_mod_10per_nos', 'betavalues_mod_10per_rec', 'betavalues_mod_10per_ear',
                                    'betavalues_mod_10per_mou', 'betavalues_mod_10per_eye')

names(get_beta_formated_w_mod) <- c('betavalues_w_mod_3000_nos', 'betavalues_w_mod_3000_rec', 'betavalues_w_mod_3000_ear',
                                    'betavalues_w_mod_3000_mou', 'betavalues_w_mod_3000_eye',
                                    'betavalues_w_mod_5000_nos', 'betavalues_w_mod_5000_rec', 'betavalues_w_mod_5000_ear',
                                    'betavalues_w_mod_5000_mou', 'betavalues_w_mod_5000_eye',
                                    'betavalues_w_mod_7000_nos', 'betavalues_w_mod_7000_rec', 'betavalues_w_mod_7000_ear',
                                    'betavalues_w_mod_7000_mou', 'betavalues_w_mod_7000_eye',
                                    'betavalues_w_mod_1per_nos', 'betavalues_w_mod_1per_rec', 'betavalues_w_mod_1per_ear',
                                    'betavalues_w_mod_1per_mou', 'betavalues_w_mod_1per_eye',
                                    'betavalues_w_mod_3per_nos', 'betavalues_w_mod_3per_rec', 'betavalues_w_mod_3per_ear',
                                    'betavalues_w_mod_3per_mou', 'betavalues_w_mod_3per_eye',
                                    'betavalues_w_mod_10per_nos', 'betavalues_w_mod_10per_rec', 'betavalues_w_mod_10per_ear',
                                    'betavalues_w_mod_10per_mou', 'betavalues_w_mod_10per_eye')

names(get_beta_formated_u_cod) <- c('betavalues_cod_3000_nos', 'betavalues_cod_3000_rec', 'betavalues_cod_3000_ear',
                                    'betavalues_cod_3000_mou', 'betavalues_cod_3000_eye',
                                    'betavalues_cod_5000_nos', 'betavalues_cod_5000_rec', 'betavalues_cod_5000_ear',
                                    'betavalues_cod_5000_mou', 'betavalues_cod_5000_eye',
                                    'betavalues_cod_7000_nos', 'betavalues_cod_7000_rec', 'betavalues_cod_7000_ear',
                                    'betavalues_cod_7000_mou', 'betavalues_cod_7000_eye',
                                    'betavalues_cod_1per_nos', 'betavalues_cod_1per_rec', 'betavalues_cod_1per_ear',
                                    'betavalues_cod_1per_mou', 'betavalues_cod_1per_eye',
                                    'betavalues_cod_3per_nos', 'betavalues_cod_3per_rec', 'betavalues_cod_3per_ear',
                                    'betavalues_cod_3per_mou', 'betavalues_cod_3per_eye',
                                    'betavalues_cod_10per_nos', 'betavalues_cod_10per_rec', 'betavalues_cod_10per_ear',
                                    'betavalues_cod_10per_mou', 'betavalues_cod_10per_eye')

names(get_beta_formated_w_cod) <- c('betavalues_w_cod_3000_nos', 'betavalues_w_cod_3000_rec', 'betavalues_w_cod_3000_ear',
                                    'betavalues_w_cod_3000_mou', 'betavalues_w_cod_3000_eye',
                                    'betavalues_w_cod_5000_nos', 'betavalues_w_cod_5000_rec', 'betavalues_w_cod_5000_ear',
                                    'betavalues_w_cod_5000_mou', 'betavalues_w_cod_5000_eye',
                                    'betavalues_w_cod_7000_nos', 'betavalues_w_cod_7000_rec', 'betavalues_w_cod_7000_ear',
                                    'betavalues_w_cod_7000_mou', 'betavalues_w_cod_7000_eye',
                                    'betavalues_w_cod_1per_nos', 'betavalues_w_cod_1per_rec', 'betavalues_w_cod_1per_ear',
                                    'betavalues_w_cod_1per_mou', 'betavalues_w_cod_1per_eye',
                                    'betavalues_w_cod_3per_nos', 'betavalues_w_cod_3per_rec', 'betavalues_w_cod_3per_ear',
                                    'betavalues_w_cod_3per_mou', 'betavalues_w_cod_3per_eye',
                                    'betavalues_w_cod_10per_nos', 'betavalues_w_cod_10per_rec', 'betavalues_w_cod_10per_ear',
                                    'betavalues_w_cod_10per_mou', 'betavalues_w_cod_10per_eye')

#run KW test for every single comparison
#non-parametric means comparison test


for(i in 1:length(get_beta_formated_u_mod)) {
  print(kruskal.test(distances ~ MoD, data = get_beta_formated_u_mod[[i]]))
  print(kruskal.test(distances ~ MoD, data = get_beta_formated_w_mod[[i]]))
  print(kruskal.test(distances ~ CoD_Simple2, data = get_beta_formated_u_cod[[i]]))
  print(kruskal.test(distances ~ CoD_Simple2, data = get_beta_formated_w_cod[[i]]))
}

#run FK test for every single comparison
#non-parametric variance comparison test


for(i in 1:length(get_beta_formated_u_mod)) {
  print(fligner.test(distances ~ MoD, data = get_beta_formated_u_mod[[i]]))
  print(fligner.test(distances ~ MoD, data = get_beta_formated_w_mod[[i]]))
  print(fligner.test(distances ~ CoD_Simple2, data = get_beta_formated_u_cod[[i]]))
  print(fligner.test(distances ~ CoD_Simple2, data = get_beta_formated_w_cod[[i]]))
}

#final table organized in excel


# Table S4 ----------------------------------------------------------------

#pull out the OTU tables for each normalization level
physeq_1per_otu <- data.frame(otu_table(physeq_1per))
physeq_3per_otu <- data.frame(otu_table(physeq_3per))
physeq_10per_otu <- data.frame(otu_table(physeq_10per))

# 1 percent cut off
#transpose the data, summarize the library size for each sample
physeq_1per_totu <- data.frame(t(physeq_1per_otu))
physeq_1per_totu$Sums <- rowSums(physeq_1per_totu)
physeq_1per_totu$SampleID <- rownames(physeq_1per_totu)
p1per_totu <- merge(metadata, physeq_1per_totu, by = 'SampleID')

#test for normality
shapiro.test(p1per_totu$Sums)

#separate out by body site
p1per_totu_mou <- p1per_totu %>% filter(Sample_Area == 'Mouth')
p1per_totu_nos <- p1per_totu %>% filter(Sample_Area == 'Nose')

#test for differences among library sizes and MOD/COD
kruskal.test(Sums ~ MoD, data = p1per_totu_mou)
posthoc.kruskal.nemenyi.test(x = p1per_totu_mou$Sums, 
                             g = p1per_totu_mou$MoD, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ MoD, data = p1per_totu_nos)
posthoc.kruskal.nemenyi.test(x = p1per_totu_nos$Sums, 
                             g = p1per_totu_nos$MoD, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ CoD_Simple2, data = p1per_totu_mou)
posthoc.kruskal.nemenyi.test(x = p1per_totu_mou$Sums, 
                             g = p1per_totu_mou$CoD_Simple2, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ CoD_Simple2, data = p1per_totu_nos)
posthoc.kruskal.nemenyi.test(x = p1per_totu_nos$Sums, 
                             g = p1per_totu_nos$CoD_Simple2, p.adjust.method = 'bonf',
                             dist='Tukey')

#repeat this process for 3%
physeq_3per_totu <- data.frame(t(physeq_3per_otu))
physeq_3per_totu$Sums <- rowSums(physeq_3per_totu)
physeq_3per_totu$SampleID <- rownames(physeq_3per_totu)
p3per_totu <- merge(metadata, physeq_3per_totu, by = 'SampleID')

p3per_totu_mou <- p3per_totu %>% filter(Sample_Area == 'Mouth')
p3per_totu_nos <- p3per_totu %>% filter(Sample_Area == 'Nose')

kruskal.test(Sums ~ MoD, data = p3per_totu_mou)
posthoc.kruskal.nemenyi.test(x = p3per_totu_mou$Sums, 
                             g = p3per_totu_mou$MoD, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ MoD, data = p3per_totu_nos)
posthoc.kruskal.nemenyi.test(x = p3per_totu_nos$Sums, 
                             g = p3per_totu_nos$MoD, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ CoD_Simple2, data = p3per_totu_mou)
posthoc.kruskal.nemenyi.test(x = p3per_totu_mou$Sums, 
                             g = p3per_totu_mou$CoD_Simple2, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ CoD_Simple2, data = p3per_totu_nos)
posthoc.kruskal.nemenyi.test(x = p3per_totu_nos$Sums, 
                             g = p3per_totu_nos$CoD_Simple2, p.adjust.method = 'bonf',
                             dist='Tukey')

#and 10%
physeq_10per_totu <- data.frame(t(physeq_10per_otu))
physeq_10per_totu$Sums <- rowSums(physeq_10per_totu)
physeq_10per_totu$SampleID <- rownames(physeq_10per_totu)
p10per_totu <- merge(metadata, physeq_10per_totu, by = 'SampleID')

p10per_totu_mou <- p10per_totu %>% filter(Sample_Area == 'Mouth')
p10per_totu_nos <- p10per_totu %>% filter(Sample_Area == 'Nose')

kruskal.test(Sums ~ MoD, data = p10per_totu_mou)
posthoc.kruskal.nemenyi.test(x = p10per_totu_mou$Sums, 
                             g = p10per_totu_mou$MoD, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ MoD, data = p10per_totu_nos)
posthoc.kruskal.nemenyi.test(x = p10per_totu_nos$Sums, 
                             g = p10per_totu_nos$MoD, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ CoD_Simple2, data = p10per_totu_mou)
posthoc.kruskal.nemenyi.test(x = p10per_totu_mou$Sums, 
                             g = p10per_totu_mou$CoD_Simple2, p.adjust.method = 'bonf',
                             dist='Tukey')

kruskal.test(Sums ~ CoD_Simple2, data = p10per_totu_nos)
posthoc.kruskal.nemenyi.test(x = p10per_totu_nos$Sums, 
                             g = p10per_totu_nos$CoD_Simple2, p.adjust.method = 'bonf',
                             dist='Tukey')



# Table S5 ----------------------------------------------------------------

#check if diversity is being lost among normalization levels
alph_richnes <- vector('list')
for(i in 1:length(phy_list)) {
  erich <- estimate_richness(phy_list[[i]], measures = c('Chao1', "Shannon"))
  erich <- add_rownames(erich, "SampleID")
  erich_sums <- merge(erich, metadata)
  alph_richnes[[i]] <- erich_sums
}
alph_richnes[[1]]$Method <- 'Rare_3000'
alph_richnes[[2]]$Method <- 'Rare_5000'
alph_richnes[[3]]$Method <- 'Rare_7000'
alph_richnes[[4]]$Method <- 'Norm_1'
alph_richnes[[5]]$Method <- 'Norm_3'
alph_richnes[[6]]$Method <- 'Norm_10'

df_Alph <- ldply (alph_richnes, data.frame)

#summary statistics
df_Alph %>% group_by(Method) %>% summarise_at(c('Chao1', "Shannon"), funs(mean, sd))
df_Alph$Method <- as.factor(df_Alph$Method)

#statistical test 

print(kruskal.test(Chao1 ~ Method, data = df_Alph))
out <- posthoc.kruskal.nemenyi.test(x=df_Alph$Chao1, g=df_Alph$Method, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

print(kruskal.test(Shannon ~ Method, data = df_Alph))
out <- posthoc.kruskal.nemenyi.test(x=df_Alph$Shannon, g=df_Alph$Method, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

# Table S6 ----------------------------------------------------------------
physeq_5000 <- rarefy_even_depth(physeq_trim, sample.size = 5000)
#removed 25 samples, 33 OTUs
#8664 taxa, 853 samples

#samp_area_list - has all sample areas as unique names
#should be 5
samp_area_list <- c('Nose', 'Rectum', 'Ears', 'Mouth', 'Eyes')

# subsetted phyloseq object by sample area
physeq_by_bodysites <- vector('list')

for(j in 1:length(samp_area_list)) {
  a <- (phyloseq::subset_samples(physeq_5000, Sample_Area == samp_area_list[j]))
  physeq_by_bodysites <- append(physeq_by_bodysites, a)
}

names(physeq_by_bodysites) <- samp_area_list

#calculate unifrac distances 
dist_u <- vector('list')
for(j in 1:length(physeq_by_bodysites)) {
  a <- phyloseq::distance(physeq_by_bodysites[[j]], "unifrac")
  dist_u[[j]] <- a
}

names(dist_u) <- samp_area_list

#make sample dataframes from statistical comparisons
sample_df_list <- vector('list')
for(i in 1:length(physeq_by_bodysites)) {
  df <- data.frame(sample_data(physeq_by_bodysites[[i]]))
  sample_df_list[[i]] <- df 
}

#caluculate beta-dispersion 
#separate out MOD and COD
beta_disp_list_u_mod <- vector('list')
beta_disp_list_u_cod <- vector('list')

for(k in 1:5) {
  a <- betadisper(dist_u[[k]], sample_df_list[[k]]$MoD)
  beta_disp_list_u_mod[[k]] <- a
  c <- betadisper(dist_u[[k]], sample_df_list[[k]]$CoD_Simple2)
  beta_disp_list_u_cod[[k]] <- c
  
}

#correctly format beta-dispersion values into dataframe to do statistical tests with 
get_beta_formated_u_mod <- vector('list')
get_beta_formated_u_cod <- vector('list')

for(i in 1:length(beta_disp_list_u_mod)) {
  beta_obj <- as.data.frame(beta_disp_list_u_mod[[i]]$distances)
  beta_obj$SampleID <- rownames(beta_obj)
  colnames(beta_obj) <- c('distances', 'SampleID')
  beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 
  get_beta_formated_u_mod[[i]] <- beta_sa
  
  beta_obj <- as.data.frame(beta_disp_list_u_cod[[i]]$distances)
  beta_obj$SampleID <- rownames(beta_obj)
  colnames(beta_obj) <- c('distances', 'SampleID')
  beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 
  get_beta_formated_u_cod[[i]] <- beta_sa
  
}

#combine beta-dispersion values into one dataframe 
betadisp_mod_df <- rbind(get_beta_formated_u_mod[[1]], 
                         get_beta_formated_u_mod[[2]],
                         get_beta_formated_u_mod[[3]],
                         get_beta_formated_u_mod[[4]],
                         get_beta_formated_u_mod[[5]])

betadisp_cod_df <- rbind(get_beta_formated_u_cod[[1]], 
                         get_beta_formated_u_cod[[2]],
                         get_beta_formated_u_cod[[3]],
                         get_beta_formated_u_cod[[4]],
                         get_beta_formated_u_cod[[5]])

# run summary statistics
betadisp_mod_df %>% group_by(Sample_Area) %>% 
  summarize_at(c('distances'), funs(mean,sd))

betadisp_mod_df %>% group_by(MoD) %>% 
  summarize_at(c('distances'), funs(mean,sd))

betadisp_cod_df %>% group_by(Sample_Area) %>% 
  summarize_at(c('distances'), funs(mean,sd))

betadisp_cod_df %>% group_by(CoD_Simple2) %>% 
  summarize_at(c('distances'), funs(mean,sd))

# statistical tests
kruskal.test(distances ~ Sample_Area, data = betadisp_mod_df)
posthoc.kruskal.nemenyi.test(x = betadisp_mod_df$distances, 
                             g = betadisp_mod_df$Sample_Area, 
                             p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(distances ~ Sample_Area, data = betadisp_cod_df)
posthoc.kruskal.nemenyi.test(x = betadisp_cod_df$distances, 
                             g = betadisp_cod_df$Sample_Area, 
                             p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(distances ~ MoD, data = betadisp_mod_df)
posthoc.kruskal.nemenyi.test(x = betadisp_mod_df$distances, 
                             g = betadisp_mod_df$MoD, 
                             p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(distances ~ CoD_Simple2, data = betadisp_cod_df)
posthoc.kruskal.nemenyi.test(x = betadisp_cod_df$distances, 
                             g = betadisp_cod_df$CoD_Simple2, 
                             p.adjust.method = 'bonf', dist='Tukey')


# Table S7 ----------------------------------------------------------------
#collapse metadata into case rather than by sample area for case metrics
metadata_coll <- subset(metadata, select=-c(SampleID, X, Sample_Area, Description))
metadata_coll <- unique(metadata_coll)
metadata_coll <- metadata_coll[complete.cases(metadata_coll), ]

#test for normality
shapiro.test(metadata_coll$Age)
shapiro.test(metadata_coll$BMI)

#test to see if certain ages/ bmis are present in different MOD/CODs
kruskal.test(Age ~ MoD, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Age, 
                             g = metadata_coll$MoD, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Age ~ CoD_Simple2, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Age, 
                             g = metadata_coll$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(BMI ~ MoD, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$BMI, 
                             g = metadata_coll$MoD, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(BMI ~ CoD_Simple2, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$BMI, 
                             g = metadata_coll$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Sex ~ MoD, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Sex, 
                             g = metadata_coll$MoD, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Sex ~ CoD_Simple2, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Sex, 
                             g = metadata_coll$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Race ~ MoD, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Race, 
                             g = metadata_coll$MoD, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Race ~ CoD_Simple2, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Race, 
                             g = metadata_coll$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(BroadPMI ~ MoD, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$BroadPMI, 
                             g = metadata_coll$MoD, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(BroadPMI ~ CoD_Simple2, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$BroadPMI, 
                             g = metadata_coll$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Event_Location ~ MoD, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Event_Location, 
                             g = metadata_coll$MoD, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Event_Location ~ CoD_Simple2, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Event_Location, 
                             g = metadata_coll$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Season ~ MoD, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Season, 
                             g = metadata_coll$MoD, p.adjust.method = 'bonf', dist='Tukey')

kruskal.test(Season ~ CoD_Simple2, data = metadata_coll)
posthoc.kruskal.nemenyi.test(x = metadata_coll$Season, 
                             g = metadata_coll$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')



# Table S8 ----------------------------------------------------------------

#for this section, you will need to do some manipulation of the code to get
#the full supplemental table. 

#you want to run the script for every body site, M/COD, and model type

#body sites: change get_beta_formated_u_mod[[#]]
#1 - nose 
#2 - rectum
#3 - ears
#4 - mouth
#5 - eyes

#m/cod: 
#get_beta_formated_u_mod or _cod
#MoD or CoD_Simple2 for choice = and table($)
#reflevel = 'Natural' or 'Cardio'


#model type
#null model
#MoD or CoD_Simple2 ~ 1
#just beta-dispersion
# ~ distances
#full model
# ~ Sex + Race + Age + BroadPMI + Season + Event_Location + BMI
#significant model
# ~ include the metadata that came out as significant in the full model 
# p value < 0.1

                                      #number for body site
long_df <- mlogit.data(get_beta_formated_u_mod[[1]], 
                       choice = 'MoD', shape = 'wide')
                              #choice for m/cod
              #m/cod                                  m/cod
model1 <- mlogit(MoD ~ 1 | Sex + Race + Age + BroadPMI + Season + Event_Location + BMI, data = long_df, reflevel = 'Natural')
summary(model1)
correct = model1$probabilities
binarycorrect = colnames(correct)[apply(correct,1,which.max)]
                            #body site      m/cod
table(get_beta_formated_u_mod[[1]]$MoD, binarycorrect)


# Table S9 ----------------------------------------------------------------

## for this section, you will need to manipulate the code to get the full
## data table
## first, find non-core community beta-dispersion

#look inside function for specific changes that need to be made
preparing_data_for_core <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  #change MoD/ CoD_Simple2 when applicable
  met <- metadata[,c('SampleID', 'MoD')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  #change factors 'Natural', 'Suicide', 'Accident', 'Homicide'
  # or 'Cardio', 'Drug', 'Gunshot', 'Asphyx', 'BFT', 'Other', 'Unknown'
  #change MoD/ CoD_Simple2 when applicable
  mtotus$MoD <- factor(mtotus$MoD, levels = c('Natural', 'Suicide', 'Accident', 'Homicide'))
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  #change based on MOD/COD
  new_df <- mtotus %>% group_by(MoD) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(Upset)
}

#start with MoD
all_otus_mod_nose <- preparing_data_for_core(physeq_by_bodysites[['Nose']])
all_otus_mod_mou <- preparing_data_for_core(physeq_by_bodysites[['Mouth']])

#clean up the data
#mouth is shown as an example
#change variable name based on bodysite
all_otus_mod_mou <- all_otus_mod_mou %>% filter(total != 0)
all_otus_mod_mou$Natural <- as.numeric(as.character(all_otus_mod_mou$Natural))
all_otus_mod_mou$Suicide <- as.numeric(as.character(all_otus_mod_mou$Suicide))
all_otus_mod_mou$Accident <- as.numeric(as.character(all_otus_mod_mou$Accident))
all_otus_mod_mou$Homicide <- as.numeric(as.character(all_otus_mod_mou$Homicide))

#otus withins MODs only
#change variable names based on bodysite
nat_otus_mod_mou <- all_otus_mod_mou %>% filter(Natural != 0) %>% 
  filter(Suicide == 0) %>% filter(Accident == 0) %>% filter(Homicide == 0)
acc_otus_mod_mou <- all_otus_mod_mou %>% filter(Accident != 0) %>% 
  filter(Suicide == 0) %>% filter(Natural == 0) %>% filter(Homicide == 0)
hom_otus_mod_mou <- all_otus_mod_mou %>% filter(Natural == 0) %>% 
  filter(Suicide == 0) %>% filter(Accident == 0) %>% filter(Homicide != 0)
sui_otus_mod_mou <- all_otus_mod_mou %>% filter(Natural == 0) %>% 
  filter(Suicide != 0) %>% filter(Accident == 0) %>% filter(Homicide == 0)

#core otus merge it
#mouth is shown as an example. Repeat for additional body sites
core_merger_mod_mou <- rbind(nat_otus_mod_mou, acc_otus_mod_mou, 
                             hom_otus_mod_mou, sui_otus_mod_mou)


#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- core_merger_mod_mou$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
core_mod_mou_bd <- pop_taxa(physeq_by_bodysites[['Mouth']], impTaxa)
core_mod_mou_bd

#unifrac distances
betaphy_mod_mou_core <- phyloseq::distance(core_mod_mou_bd, "unifrac")
#beta-dispersion 
#change the number depending on the order listed in samp_area_list
betadisp_mod_mou_core <- betadisper(betaphy_mod_mou_core, sample_df_list[[4]]$MoD)

# get it correctly formatted
beta_obj <- as.data.frame(betadisp_mod_mou_core$distances)
beta_obj$SampleID <- rownames(beta_obj)
colnames(beta_obj) <- c('distances', 'SampleID')
betacore_df_mod_mou <- merge(beta_obj, metadata, by = 'SampleID') 

#COD
all_otus_cod_nose <- preparing_data_for_core(physeq_by_bodysites[['Nose']])
all_otus_cod_mou <- preparing_data_for_core(physeq_by_bodysites[['Mouth']])
all_otus_cod_ear <- preparing_data_for_core(physeq_by_bodysites[['Ears']])

#clean up the data
#mouth is shown as an example
#change variable name based on bodysite
all_otus_cod_mou <- all_otus_cod_mou %>% filter(total != 0)
all_otus_cod_mou$Cardio <- as.numeric(as.character(all_otus_cod_mou$Cardio))
all_otus_cod_mou$Drug <- as.numeric(as.character(all_otus_cod_mou$Drug))
all_otus_cod_mou$Gunshot <- as.numeric(as.character(all_otus_cod_mou$Gunshot))
all_otus_cod_mou$Asphyx <- as.numeric(as.character(all_otus_cod_mou$Asphyx))
all_otus_cod_mou$BFT <- as.numeric(as.character(all_otus_cod_mou$BFT))
all_otus_cod_mou$Other <- as.numeric(as.character(all_otus_cod_mou$Other))
#all_otus_cod_mou$Unknown <- as.numeric(as.character(all_otus_cod_mou$Unknown))

#otus withins cods only
#change variable names based on bodysite
#remove unknown for mouth, no unkonwn in that bodysite
car_otus_cod_mou <- all_otus_cod_mou %>% filter(Cardio != 0) %>% filter(Drug == 0) %>% 
  filter(Gunshot == 0) %>% filter(Asphyx == 0) %>% filter(BFT == 0) %>% filter(Other == 0) 
dru_otus_cod_mou <- all_otus_cod_mou %>% filter(Cardio == 0) %>% filter(Drug != 0) %>% 
  filter(Gunshot == 0) %>% filter(Asphyx == 0) %>% filter(BFT == 0) %>% filter(Other == 0) 
gun_otus_cod_mou <- all_otus_cod_mou %>% filter(Cardio == 0) %>% filter(Drug == 0) %>% 
  filter(Gunshot != 0) %>% filter(Asphyx == 0) %>% filter(BFT == 0) %>% filter(Other == 0) 
asp_otus_cod_mou <- all_otus_cod_mou %>% filter(Cardio == 0) %>% filter(Drug == 0) %>% 
  filter(Gunshot == 0) %>% filter(Asphyx != 0) %>% filter(BFT == 0) %>% filter(Other == 0) 
bft_otus_cod_mou <- all_otus_cod_mou %>% filter(Cardio == 0) %>% filter(Drug == 0) %>% 
  filter(Gunshot == 0) %>% filter(Asphyx == 0) %>% filter(BFT != 0) %>% filter(Other == 0)
oth_otus_cod_mou <- all_otus_cod_mou %>% filter(Cardio == 0) %>% filter(Drug == 0) %>% 
  filter(Gunshot == 0) %>% filter(Asphyx == 0) %>% filter(BFT == 0) %>% filter(Other != 0)
#unk_otus_cod_mou <- all_otus_cod_mou %>% filter(Cardio == 0) %>% filter(Drug == 0) %>% 
# filter(Gunshot == 0) %>% filter(Asphyx == 0) %>% filter(BFT == 0) %>% filter(Other == 0) %>% 
# filter(Unknown != 0)

#core otus merge it
core_merger_cod_mou <- rbind(car_otus_cod_mou, dru_otus_cod_mou, 
                             gun_otus_cod_mou, asp_otus_cod_mou, bft_otus_cod_mou,
                             oth_otus_cod_mou)


#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- core_merger_cod_mou$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
core_cod_mou_bd <- pop_taxa(physeq_by_bodysites[['Mouth']], impTaxa)
core_cod_mou_bd

#unifrac distances
betaphy_cod_mou_core <- phyloseq::distance(core_cod_mou_bd, "unifrac")
#beta-dispersion 
#change the number depending on the order listed in samp_area_list
betadisp_cod_mou_core <- betadisper(betaphy_cod_mou_core, sample_df_list[[4]]$CoD_Simple2)

# get it correctly formatted
beta_obj <- as.data.frame(betadisp_cod_mou_core$distances)
beta_obj$SampleID <- rownames(beta_obj)
colnames(beta_obj) <- c('distances', 'SampleID')
betacore_df_cod_mou <- merge(beta_obj, metadata, by = 'SampleID') 

## second, find the random forest indicator taxa

#function gets the data set up
#need to change depending on MoD/ CoD_Simple2
random_foresting_setup<- function(data_phy) {
  otu <- as.data.frame(t(otu_table(data_phy)))
  otu$SampleID <- rownames(otu)
  metadata <- data.frame(sample_data(data_phy))
  #change this MoD/CoD_Simple2
  meta_sa <- metadata %>% select(SampleID, MoD)
  meta_sa$SampleID <- as.character(meta_sa$SampleID)
  otu_f <- merge(meta_sa, otu, by = 'SampleID')
  otu_f <- otu_f[,-1]
  names(otu_f) <- make.names(names(otu_f))
  return(otu_f)
}

#example of MoD for mouth region.
#change depending on interest
rf_mod_nos <- random_foresting_setup(get_beta_formated_u_mod[[4]])

#find the indicator taxa
boruta.bank_train <- Boruta(MoD~., data = rf_mod_nos, doTrace = 2, pValue = .05)
print(boruta.bank_train)

# dataframe of important taxa
tax_of_importance <- data.frame(boruta.bank_train[["finalDecision"]])
tax_of_importance$Taxa <- rownames(tax_of_importance)
colnames(tax_of_importance) <- c('Decision', 'Taxa')
#keep confirmed and tenative
tax_of_importance <- tax_of_importance %>% filter(Decision != 'Rejected')

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- tax_of_importance$Taxa
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change following variable names based on what the data is (body site and MOD/COD)
rf_mod_nos_bd <- pop_taxa(physeq_by_bodysites[['Nose']], impTaxa)
rf_mod_nos_bd

#unifrac distances
betaphy_mod_nos_rf <- phyloseq::distance(rf_mod_nos_bd, "unifrac")
#beta-dispersion 
#change the number depending on the order listed in samp_area_list
betadisp_mod_nos_rf <- betadisper(betaphy_mod_nos_rf, sample_df_list[[1]]$MoD)
# get it correctly formatted
beta_obj <- as.data.frame(betadisp_mod_nos_rf$distances)
beta_obj$SampleID <- rownames(beta_obj)
colnames(beta_obj) <- c('distances', 'SampleID')
beta_rf_mod_nos <- merge(beta_obj, metadata, by = 'SampleID') 


#for this section, you will need to do some manipulation of the code to get
#the full supplemental table. 

#you want to run the script for every body site, M/COD, and model type
# but this time, for the non-core and random forest indicator taxa communities
# same as above 

#non-core: betacore_df_m/cod_body site dataframes
#random forest: beta_rf_m/cod_body site dataframes

#model type
#null model
#MoD or CoD_Simple2 ~ 1
#just beta-dispersion
# ~ distances
#full model
# ~ Sex + Race + Age + BroadPMI + Season + Event_Location + BMI
#significant model
# ~ include the metadata that came out as significant in the full model 
# p value < 0.1

#example for MOD
long_df <- mlogit.data(betacore_df_mod_mou, 
                       choice = 'MoD', shape = 'wide')

model1 <- mlogit(MoD ~ 1 | distances, data = long_df, reflevel = 'Natural')
summary(model1)
correct = model1$probabilities
binarycorrect = colnames(correct)[apply(correct,1,which.max)]
table(betacore_df_mod_mou$MoD, binarycorrect)

#example for COD
long_df <- mlogit.data(betacore_df_cod_mou, 
                       choice = 'CoD_Simple2', shape = 'wide')

model1 <- mlogit(CoD_Simple2 ~ 1 | distances, data = long_df, reflevel = 'Cardio')
summary(model1)
correct = model1$probabilities
binarycorrect = colnames(correct)[apply(correct,1,which.max)]
table(betacore_df_cod_mou$CoD_Simple2, binarycorrect)

# Table S10 ---------------------------------------------------------------

#this table was made in excel, which summarizes the data found in table
# S11. Please see the next section for the code


# Table S11 ---------------------------------------------------------------

#for this section, all comparisons have been coded out
#need to change modeling section to get full dataframe as above

#pairwise comparisons
#natural - accident
pc_mod_blr_nos <- get_beta_formated_u_mod[[1]]
pc_mod_blr_nos_N <- pc_mod_blr_nos %>% filter(MoD == 'Natural')
pc_mod_blr_nos_A <- pc_mod_blr_nos %>% filter(MoD == 'Accident')
pc_mod_blr_nos_NA <- rbind(pc_mod_blr_nos_N, pc_mod_blr_nos_A)

pc_mod_blr_mou <- get_beta_formated_u_mod[[4]]
pc_mod_blr_mou_N <- pc_mod_blr_mou %>% filter(MoD == 'Natural')
pc_mod_blr_mou_A <- pc_mod_blr_mou %>% filter(MoD == 'Accident')
pc_mod_blr_mou_NA <- rbind(pc_mod_blr_mou_N, pc_mod_blr_mou_A)

pc_mod_blr_rf_nos <- beta_rf_mod_nos
pc_mod_blr_rf_nos_N <- pc_mod_blr_rf_nos %>% filter(MoD == 'Natural')
pc_mod_blr_rf_nos_A <- pc_mod_blr_rf_nos %>% filter(MoD == 'Accident')
pc_mod_blr_rf_nos_NA <- rbind(pc_mod_blr_rf_nos_N, pc_mod_blr_rf_nos_A)

pc_mod_blr_rf_mou <- beta_rf_mod_mou
pc_mod_blr_rf_mou_N <- pc_mod_blr_rf_mou %>% filter(MoD == 'Natural')
pc_mod_blr_rf_mou_A <- pc_mod_blr_rf_mou %>% filter(MoD == 'Accident')
pc_mod_blr_rf_mou_NA <- rbind(pc_mod_blr_rf_mou_N, pc_mod_blr_rf_mou_A)

#natural - homicide
pc_mod_blr_nos <- get_beta_formated_u_mod[[1]]
pc_mod_blr_nos_N <- pc_mod_blr_nos %>% filter(MoD == 'Natural')
pc_mod_blr_nos_H <- pc_mod_blr_nos %>% filter(MoD == 'Homicide')
pc_mod_blr_nos_NA <- rbind(pc_mod_blr_nos_N, pc_mod_blr_nos_H)

pc_mod_blr_mou <- get_beta_formated_u_mod[[4]]
pc_mod_blr_mou_N <- pc_mod_blr_mou %>% filter(MoD == 'Natural')
pc_mod_blr_mou_H <- pc_mod_blr_mou %>% filter(MoD == 'Homicide')
pc_mod_blr_mou_NA <- rbind(pc_mod_blr_mou_N, pc_mod_blr_mou_H)

pc_mod_blr_rf_nos <- beta_rf_mod_nos
pc_mod_blr_rf_nos_N <- pc_mod_blr_rf_nos %>% filter(MoD == 'Natural')
pc_mod_blr_rf_nos_H <- pc_mod_blr_rf_nos %>% filter(MoD == 'Homicide')
pc_mod_blr_rf_nos_NA <- rbind(pc_mod_blr_rf_nos_N, pc_mod_blr_rf_nos_H)

pc_mod_blr_rf_mou <- beta_rf_mod_mou
pc_mod_blr_rf_mou_N <- pc_mod_blr_rf_mou %>% filter(MoD == 'Natural')
pc_mod_blr_rf_mou_H <- pc_mod_blr_rf_mou %>% filter(MoD == 'Homicide')
pc_mod_blr_rf_mou_NA <- rbind(pc_mod_blr_rf_mou_N, pc_mod_blr_rf_mou_H)

#suicide - everything else
pc_mod_blr_nos <- get_beta_formated_u_mod[[1]]
pc_mod_blr_nos[pc_mod_blr_nos$MoD == 'Suicide', "Suicide"] <- 'Suicide'
pc_mod_blr_nos[pc_mod_blr_nos$MoD != 'Suicide', "Suicide"] <- 'Not'
pc_mod_blr_nos$Suicide <- as.factor(pc_mod_blr_nos$Suicide)

pc_mod_blr_mou <- get_beta_formated_u_mod[[4]]
pc_mod_blr_mou[pc_mod_blr_mou$MoD == 'Suicide', "Suicide"] <- 'Suicide'
pc_mod_blr_mou[pc_mod_blr_mou$MoD != 'Suicide', "Suicide"] <- 'Not'
pc_mod_blr_mou$Suicide <- as.factor(pc_mod_blr_mou$Suicide)

pc_mod_blr_rf_nos <- beta_rf_mod_nos
pc_mod_blr_rf_nos[pc_mod_blr_rf_nos$MoD == 'Suicide', "Suicide"] <- 'Suicide'
pc_mod_blr_rf_nos[pc_mod_blr_rf_nos$MoD != 'Suicide', "Suicide"] <- 'Not'
pc_mod_blr_rf_nos$Suicide <- as.factor(pc_mod_blr_rf_nos$Suicide)

pc_mod_blr_rf_mou <- beta_rf_mod_nos
pc_mod_blr_rf_mou[pc_mod_blr_rf_mou$MoD == 'Suicide', "Suicide"] <- 'Suicide'
pc_mod_blr_rf_mou[pc_mod_blr_rf_mou$MoD != 'Suicide', "Suicide"] <- 'Not'
pc_mod_blr_rf_mou$Suicide <- as.factor(pc_mod_blr_rf_mou$Suicide)

#cardiovascular diease vs. drug use
pc_cod_blr_ear <- get_beta_formated_u_cod[[3]]
pc_cod_blr_ear_C <- pc_cod_blr_ear %>% filter(CoD_Simple2 == 'Cardio')
pc_cod_blr_ear_D <- pc_cod_blr_ear %>% filter(CoD_Simple2 == 'Drug')
pc_cod_blr_ear_CA <- rbind(pc_cod_blr_ear_C, pc_cod_blr_ear_D)

pc_cod_blr_mou <- get_beta_formated_u_cod[[4]]
pc_cod_blr_mou_C <- pc_cod_blr_mou %>% filter(CoD_Simple2 == 'Cardio')
pc_cod_blr_mou_D <- pc_cod_blr_mou %>% filter(CoD_Simple2 == 'Drug')
pc_cod_blr_mou_CA <- rbind(pc_cod_blr_mou_C, pc_cod_blr_mou_D)

pc_cod_blr_rf_ear <- beta_rf_cod_ear
pc_cod_blr_rf_ear_C <- pc_cod_blr_rf_ear %>% filter(CoD_Simple2 == 'Cardio')
pc_cod_blr_rf_ear_D <- pc_cod_blr_rf_ear %>% filter(CoD_Simple2 == 'Drug')
pc_cod_blr_rf_ear_CA <- rbind(pc_cod_blr_rf_ear_C, pc_cod_blr_rf_ear_D)

pc_cod_blr_rf_mou <- beta_rf_cod_mou
pc_cod_blr_rf_mou_C <- pc_cod_blr_rf_mou %>% filter(CoD_Simple2 == 'Cardio')
pc_cod_blr_rf_mou_D <- pc_cod_blr_rf_mou %>% filter(CoD_Simple2 == 'Drug')
pc_cod_blr_rf_mou_CA <- rbind(pc_cod_blr_rf_mou_C, pc_cod_blr_rf_mou_D)

#violent vs. non-violent
# already in metadata

#disease vs. non-disease state
pc_cod_blr_ear <- get_beta_formated_u_cod[[1]]
pc_cod_blr_ear[pc_cod_blr_ear$MoD == 'Natural', "Disease"] <- 'Diseased'
pc_cod_blr_ear[pc_cod_blr_ear$MoD != 'Natural', "Disease"] <- 'Not'
pc_cod_blr_ear$Disease <- as.factor(pc_cod_blr_ear$Disease)

pc_cod_blr_mou <- get_beta_formated_u_cod[[3]]
pc_cod_blr_mou[pc_cod_blr_mou$MoD == 'Natural', "Disease"] <- 'Diseased'
pc_cod_blr_mou[pc_cod_blr_mou$MoD != 'Natural', "Disease"] <- 'Not'
pc_cod_blr_mou$Disease <- as.factor(pc_cod_blr_mou$Disease)

pc_cod_blr_rf_ear <- beta_rf_cod_ear
pc_cod_blr_rf_ear[pc_cod_blr_rf_ear$MoD == 'Natural', "Disease"] <- 'Diseased'
pc_cod_blr_rf_ear[pc_cod_blr_rf_ear$MoD != 'Natural', "Disease"] <- 'Not'
pc_cod_blr_rf_ear$Disease <- as.factor(pc_cod_blr_rf_ear$Disease)

pc_cod_blr_rf_mou <- beta_rf_cod_mou
pc_cod_blr_rf_mou[pc_cod_blr_rf_mou$MoD == 'Natural', "Disease"] <- 'Diseased'
pc_cod_blr_rf_mou[pc_cod_blr_rf_mou$MoD != 'Natural', "Disease"] <- 'Not'
pc_cod_blr_rf_mou$Disease <- as.factor(pc_cod_blr_rf_mou$Disease)


#for this section, you will need to do some manipulation of the code to get
#the full supplemental table. 

#you want to run the script for every body site, M/COD, and model type
# but this time, for the non-core and random forest indicator taxa communities
# same as above 

#non-core: betacore_df_m/cod_body site dataframes
#random forest: beta_rf_m/cod_body site dataframes

#model type
#full model
# ~ Sex + Race + Age + BroadPMI + Season + Event_Location + BMI
#significant model
# ~ include the metadata that came out as significant in the full model 
# p value < 0.1

model = glm(MoD ~ distances, family = binomial(link = 'logit'),
            data = pc_mod_blr_nos_NA)
summary(model)

#check significance of reduction in error by adding predictors
chidiff = model$null.deviance - model$deviance
dfdiff = model$df.null - model$df.residual
pchisq(chidiff, dfdiff, lower.tail = F)

PseudoR2(model)

#correct predictions
correct = model$fitted.values
binarycorrect = ifelse(correct > 0.5, 1,0)
binarycorrectlab = factor(binarycorrect,
                          levels = c(0,1),
                          #change labels for what is being compared
                          labels = c('Accident', 'Natural'))
    #change data input and $MoD or CoD_Simple2
table(pc_mod_mlr_nos_NA$MoD, binarycorrectlab)


# Table S12 ---------------------------------------------------------------
#load in files for just case study metadata

#load in metadata
metadata=(read.csv("HPMMMeta_matched_design_nose.csv",header=TRUE))
metadata$BMI <- as.numeric(as.character((metadata$BMI)))
#format it into phyloseq format
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

#put otu table in phyloseq format
rownames(otu) <- otu$OTUID
otu <- otu[,-1]
OTU=otu_table(otu, taxa_are_rows=TRUE)
#merge tree and otu tables
physeq_otu.tree=phyloseq(OTU,tree, sampdat)

#put taxonomy table in phyloseq format
rownames(tax) <- tax$OTUID
tax_reduc <- merge(tax, otu, by = "row.names")
rownames(tax_reduc) <- tax_reduc$Row.names
tax_reduc <- tax_reduc[,-1]
tax_f <- tax_reduc[,1:7]
tax_f <- as.matrix(tax_f)
TAX=tax_table(tax_f)
taxa_names(TAX)=row.names(OTU)

#merge it all together into one phyloseq object
physeq <- merge_phyloseq(physeq_otu.tree, TAX)
physeq
#30906 taxa, 43 samples

#triming out taxa that are not representative of .01% of mean sequence number
physeq_trim <- prune_taxa(taxa_sums(physeq) > sum(otu) *.001 / 878, physeq)
physeq_trim
#952 taxa

#rarefy
physeq_5000 <- rarefy_even_depth(physeq_trim, sample.size = 5000)
physeq_5000
#10 OTUs removed. 951 taxa

###find some indicator taxa

#set up data format correctly from physeq object
random_foresting_setup<- function(data_phy) {
  otu <- as.data.frame(t(otu_table(data_phy)))
  otu$SampleID <- rownames(otu)
  metadata <- data.frame(sample_data(data_phy))
  meta_sa <- metadata %>% select(SampleID, Suicide)
  meta_sa$SampleID <- as.character(meta_sa$SampleID)
  otu_f <- merge(meta_sa, otu, by = 'SampleID')
  otu_f <- otu_f[,-1]
  names(otu_f) <- make.names(names(otu_f))
  return(otu_f)
}

cs_nose <- random_foresting_setup(physeq_5000)


boruta.bank_train <- Boruta(Suicide~., data = cs_nose, doTrace = 2, pValue = .05)
print(boruta.bank_train)

# dataframe of important taxa
tax_of_importance <- data.frame(boruta.bank_train[["finalDecision"]])
tax_of_importance$Taxa <- rownames(tax_of_importance)
colnames(tax_of_importance) <- c('Decision', 'Taxa')
#keep confirmed and tenative
tax_of_importance <- tax_of_importance %>% filter(Decision != 'Rejected')
tax_of_importance

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- tax_of_importance$Taxa
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change following variable names based on what the data is (body site and MOD/COD)
cs_nos_bd <- pop_taxa(physeq_5000, impTaxa)
tax_table(cs_nos_bd)


### test preliminarily if beta-disperion might be important
#test permutationally if different
#change beta variable based on comparison
beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Suicide)
  print(return(permutest(beta)))
}

beta_dispersion_calc(physeq_5000)

#calculate unifrac distances 
beta_df_cs <- phyloseq::distance(physeq_5000, "unifrac")
samp_df <- data.frame(sample_data(physeq_5000))

betadisp_df_cs <- betadisper(beta_df_cs, samp_df$Suicide)

beta_obj <- as.data.frame(betadisp_df_cs$distances)
beta_obj$SampleID <- rownames(beta_obj)
colnames(beta_obj) <- c('distances', 'SampleID')
beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 

#summary statistics
beta_sa %>% group_by(Suicide) %>% 
  summarize_at(c('distances'), funs(mean,sd))

# statistical tests
kruskal.test(distances ~ Suicide, data = beta_sa)

#logistic regression
model = glm(Suicide ~  distances, 
            family = binomial(link = 'logit'),
            data = beta_sa)
summary(model)

#check significance of reduction in error by adding predictors
chidiff = model$null.deviance - model$deviance
dfdiff = model$df.null - model$df.residual
pchisq(chidiff, dfdiff, lower.tail = F)

PseudoR2(model)

#correct predictions
correct = model$fitted.values
binarycorrect = ifelse(correct > 0.5, 1,0)
binarycorrect = factor(binarycorrect,
                       levels = c(0,1),
                       labels = c('Not', 'Suicide'))

table(beta_sa$Suicide, binarycorrect)


# Table S13 ---------------------------------------------------------------

##homicide vs. suicidal gunshots

#indicator taxa
random_foresting_setup<- function(data_phy) {
  otu <- as.data.frame(t(otu_table(data_phy)))
  otu$SampleID <- rownames(otu)
  metadata <- data.frame(sample_data(data_phy))
  meta_sa <- metadata %>% select(SampleID, MoD, CoD_Simple2)
  meta_sa$SampleID <- as.character(meta_sa$SampleID)
  otu_f <- merge(meta_sa, otu, by = 'SampleID')
  otu_f <- otu_f[,-1]
  names(otu_f) <- make.names(names(otu_f))
  return(otu_f)
}

rf_mod_nose <- random_foresting_setup(physeq_by_bodysites[[1]])

rf_mod_nos_gun <- rf_mod_nose %>%  filter(CoD_Simple2 == 'Gunshot')

boruta.bank_train <- Boruta(MoD~., data = rf_mod_nos_gun, doTrace = 2, pValue = .05)
print(boruta.bank_train)

# dataframe of important taxa
tax_of_importance <- data.frame(boruta.bank_train[["finalDecision"]])
tax_of_importance$Taxa <- rownames(tax_of_importance)
colnames(tax_of_importance) <- c('Decision', 'Taxa')
#keep confirmed and tenative
tax_of_importance <- tax_of_importance %>% filter(Decision != 'Rejected')

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- tax_of_importance$Taxa
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change following variable names based on what the data is (body site and MOD/COD)
rf_mod_nos_ind <- pop_taxa(physeq_by_bodysites[['Nose']], impTaxa)
rf_mod_nos_ind_tax <- data.frame(tax_table(rf_mod_nos_ind))

#test permutationally if different
#change beta variable based on comparison
beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$MoD)
  print(return(permutest(beta)))
}

phy_nos_gun <- subset_samples(physeq_by_bodysites[['Nose']], CoD_Simple2 == 'Gunshot')

beta_dispersion_calc(phy_nos_gun)

#calculate unifrac distances 
beta_df_cs <- phyloseq::distance(phy_nos_gun, "unifrac")
samp_df <- data.frame(sample_data(phy_nos_gun))

betadisp_df_cs <- betadisper(beta_df_cs, samp_df$MoD)

beta_obj <- as.data.frame(betadisp_df_cs$distances)
beta_obj$SampleID <- rownames(beta_obj)
colnames(beta_obj) <- c('distances', 'SampleID')
beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 

#summary statistics
beta_sa %>% group_by(MoD) %>% 
  summarize_at(c('distances'), funs(mean,sd))

# statistical tests
kruskal.test(distances ~ MoD, data = beta_sa)

#logistic regression
# Sex + BMI + Race + Season + Event_Location + BroadPMI + Age
model = glm(MoD ~  distances + Age, 
            family = binomial(link = 'logit'),
            data = beta_sa)
summary(model)

#check significance of reduction in error by adding predictors
chidiff = model$null.deviance - model$deviance
dfdiff = model$df.null - model$df.residual
pchisq(chidiff, dfdiff, lower.tail = F)

PseudoR2(model)

#correct predictions
correct = model$fitted.values
binarycorrect = ifelse(correct > 0.5, 1,0)
binarycorrect = factor(binarycorrect,
                       levels = c(0,1),
                       labels = c('Homicide', 'Suicide'))

table(beta_sa$MoD, binarycorrect)


## accidental vs. suicidal drug use
random_foresting_setup<- function(data_phy) {
  otu <- as.data.frame(t(otu_table(data_phy)))
  otu$SampleID <- rownames(otu)
  metadata <- data.frame(sample_data(data_phy))
  meta_sa <- metadata %>% select(SampleID, MoD, CoD_Simple2)
  meta_sa$SampleID <- as.character(meta_sa$SampleID)
  otu_f <- merge(meta_sa, otu, by = 'SampleID')
  otu_f <- otu_f[,-1]
  names(otu_f) <- make.names(names(otu_f))
  return(otu_f)
}

rf_mod_nose <- random_foresting_setup(physeq_by_bodysites[[1]])
rf_mod_nos_Drug <- rf_mod_nose %>%  filter(CoD_Simple2 == 'Drug')

boruta.bank_train <- Boruta(MoD~., data = rf_mod_nos_Drug, doTrace = 2, pValue = .05)
print(boruta.bank_train)

# dataframe of important taxa
tax_of_importance <- data.frame(boruta.bank_train[["finalDecision"]])
tax_of_importance$Taxa <- rownames(tax_of_importance)
colnames(tax_of_importance) <- c('Decision', 'Taxa')
#keep confirmed and tenative
tax_of_importance <- tax_of_importance %>% filter(Decision != 'Rejected')

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- tax_of_importance$Taxa
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change following variable names based on what the data is (body site and MOD/COD)
rf_mod_nos_ind <- pop_taxa(physeq_by_bodysites[['Nose']], impTaxa)
rf_mod_nos_ind_tax <- data.frame(tax_table(rf_mod_nos_ind))

#test permutationally if different
#change beta variable based on comparison
beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$MoD)
  print(return(permutest(beta)))
}


phy_nos_Drug <- subset_samples(physeq_by_bodysites[['Nose']], CoD_Simple2 == 'Drug')

beta_dispersion_calc(phy_nos_Drug)

#calculate unifrac distances 
beta_df_cs <- phyloseq::distance(phy_nos_Drug, "unifrac")
samp_df <- data.frame(sample_data(phy_nos_Drug))

betadisp_df_cs <- betadisper(beta_df_cs, samp_df$MoD)

beta_obj <- as.data.frame(betadisp_df_cs$distances)
beta_obj$SampleID <- rownames(beta_obj)
colnames(beta_obj) <- c('distances', 'SampleID')
beta_sa <- merge(beta_obj, metadata, by = 'SampleID') 
beta_sa <- beta_sa %>% filter(MoD != 'Natural')

#summary statistics
beta_sa %>% group_by(MoD) %>% 
  summarize_at(c('distances'), funs(mean,sd))

# statistical tests
kruskal.test(distances ~ MoD, data = beta_sa)

# Sex + BMI + Race + Season + Event_Location + BroadPMI + Age
model = glm(MoD ~  distances +Season, 
            family = binomial(link = 'logit'),
            data = beta_sa)
summary(model)

#check significance of reduction in error by adding predictors
chidiff = model$null.deviance - model$deviance
dfdiff = model$df.null - model$df.residual
pchisq(chidiff, dfdiff, lower.tail = F)

PseudoR2(model)

#correct predictions
correct = model$fitted.values
binarycorrect = ifelse(correct > 0.5, 1,0)
binarycorrect = factor(binarycorrect,
                       levels = c(0,1),
                       labels = c('Accident', 'Suicide'))

table(beta_sa$MoD, binarycorrect)


# Table S14 ---------------------------------------------------------------
random_foresting_setup<- function(data_phy) {
  otu <- as.data.frame(t(otu_table(data_phy)))
  otu$SampleID <- rownames(otu)
  metadata <- data.frame(sample_data(data_phy))
  #change MoD or CoD_Simple2 when needed
  meta_sa <- metadata %>% select(SampleID, MoD)
  meta_sa$SampleID <- as.character(meta_sa$SampleID)
  otu_f <- merge(meta_sa, otu, by = 'SampleID')
  otu_f <- otu_f[,-1]
  names(otu_f) <- make.names(names(otu_f))
  return(otu_f)
}


rf_mod_mou <- random_foresting_setup(physeq_bodysite[[4]])
rf_mod_nos <- random_foresting_setup(physeq_bodysite[[1]])
rf_cod_mou <- random_foresting_setup(physeq_bodysite[[4]])
rf_cod_nos <- random_foresting_setup(physeq_bodysite[[1]])
rf_cod_ear <- random_foresting_setup(physeq_bodysite[[3]])

#OOB error model
#change based on input data
m1 <- randomForest(
  #MoD/ CoD_Simple2
  formula = MoD ~ .,
  #data frame
  data    = rf_mod_acc_mou,
  ntree= 2000
)
