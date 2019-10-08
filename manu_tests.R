rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(car)
library(GGally)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lme4)
library(phyloseq)
library(plyr)
library(PMCMR)
library(microbiome)
library(randomForest)
library(rsample)
library(tidyverse)
library(vegan)
#for reproducibilty
set.seed(1234)

#load in files
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

#get some summary stats about metadata
#interested in how age and BMI is distributed throughout Manners of Death (MoD) and Causes of Death (CoD)
metadata %>% group_by(Sample_Area) %>% summarize_all(funs(length))
metadata %>% group_by(MoD) %>% summarize_at(c('Age'), funs(mean,median,sd))
metadata %>% group_by(CoD_Simple) %>% summarize_at(c('Age'), funs(mean,median,sd))
metadata %>% group_by(CoD_Simple2) %>% summarize_at(c('Age'), funs(mean,median,sd))
metadata %>% group_by(MoD) %>% summarize_at(c('BMI'), funs(mean,median,sd))
metadata %>% group_by(CoD_Simple) %>% summarize_at(c('BMI'), funs(mean,median,sd))
metadata %>% group_by(CoD_Simple2) %>% summarize_at(c('BMI'), funs(mean,median,sd))

#regress BMI and Age to see if there are differences
m1 <- lm(formula = as.numeric(BMI) ~ MoD, data = metadata)
summary(m1)
m1

m2 <- lm(formula = as.numeric(BMI) ~ CoD_Simple2, data = metadata)
summary(m2)
m2


m3 <- lm(formula = as.numeric(Age) ~ MoD, data = metadata)
summary(m3)
m3

m4 <- lm(formula = as.numeric(Age) ~ CoD_Simple2, data = metadata)
summary(m4)
m4


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


# How does methods/normalization affect beta dispersion?

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

#find average library sizes for normalized physeq objs. 
physeq_1per_otu <- as.data.frame(otu_table(physeq_1per))
physeq_3per_otu <- as.data.frame(otu_table(physeq_3per))
physeq_10per_otu <- as.data.frame(otu_table(physeq_10per))

physeq_1per_sums <- c()
for(i in 1:length(colnames(physeq_1per_otu))) {
 a <- sum(physeq_1per_otu[,i])
 physeq_1per_sums <- c(a, physeq_1per_sums)
}

mean(physeq_1per_sums)
sd(physeq_1per_sums)
median(physeq_1per_sums)
min(physeq_1per_sums)
max(physeq_1per_sums)

physeq_3per_sums <- c()
for(i in 1:length(colnames(physeq_3per_otu))) {
  a <- sum(physeq_3per_otu[,i])
  physeq_3per_sums <- c(a, physeq_3per_sums)
}

mean(physeq_3per_sums)
sd(physeq_3per_sums)
median(physeq_3per_sums)
min(physeq_3per_sums)
max(physeq_3per_sums)

physeq_10per_sums <- c()
for(i in 1:length(colnames(physeq_10per_otu))) {
  a <- sum(physeq_10per_otu[,i])
  physeq_10per_sums <- c(a, physeq_10per_sums)
}

mean(physeq_10per_sums)
sd(physeq_10per_sums)
median(physeq_10per_sums)
min(physeq_10per_sums)
max(physeq_10per_sums)

# normalization method comparison
#   rarefaction
#     three levels
#   'Core' normalization
#     three levels
# dissimilarity matrix comparison
#   unifrac and weighted unifrac for relative abundance/ binary comparisons 

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

#check if diversity is being lost among samples
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

df_Alph %>% group_by(Method) %>% summarise_at(c('Chao1', "Shannon"), funs(mean, sd))
df_Alph$Method <- as.factor(df_Alph$Method)

print(kruskal.test(Chao1 ~ Method, data = df_Alph))
out <- posthoc.kruskal.nemenyi.test(x=df_Alph$Chao1, g=df_Alph$Method, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

print(kruskal.test(Shannon ~ Method, data = df_Alph))
out <- posthoc.kruskal.nemenyi.test(x=df_Alph$Shannon, g=df_Alph$Method, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

df_Alph$Method <- factor(df_Alph$Method, levels = c('Rare_3000', 'Rare_5000', 'Rare_7000', 
                                                       'Norm_1', 'Norm_3', "Norm_10"))
p <- ggplot(df_Alph, aes(x=Method, y=Chao1, fill = Method)) + 
  geom_boxplot()
p

p2 <- ggplot(df_Alph, aes(x=Method, y=Shannon, fill = Method)) + 
  geom_boxplot()
p2

theme_set(theme_classic(base_size = 18))
tiff("alph_div.TIF", width = 3000, height = 3500, res=300)
ggarrange(p,p2,
          nrow = 2, ncol = 1)
dev.off()

#samp_area_list - has all sample areas as unique names
#should be 5
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


#function for taking subsetted phyloseq objects and finding corresponding unifrac or weighted 
#unifrac distances 
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


#calculate beta-dispersion for each sample area, distance metric, MoD/CoD, and rarefaction/
#normalization level. 
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


# Matrix ------------------------------------------------------------------
#matrix was cleaned up in excel with the results of every comparison
compare_matrix_data <- read.csv('question_1_cod2.csv')
#removes any NAs
compare_matrix_data <- compare_matrix_data[complete.cases(compare_matrix_data), ]

#specify order of rarefaction/normalization level for figure
compare_matrix_data$Level <- factor(compare_matrix_data$Level, levels = c('10%', '3%', '1%',
                                                                          '3,000', '5,000', '7,000'))
#convert p values from the tests into categorical variable
compare_matrix_data[compare_matrix_data$KW_pvalue > 0.1, "KW_pvalue_new"] <- '1'
compare_matrix_data[compare_matrix_data$KW_pvalue <= 0.1, "KW_pvalue_new"] <- '0.1'
compare_matrix_data[compare_matrix_data$KW_pvalue <= 0.05, "KW_pvalue_new"] <- '0.05'

compare_matrix_data[compare_matrix_data$FK_pvalue > 0.1, "FK_pvalue_new"] <- '1'
compare_matrix_data[compare_matrix_data$FK_pvalue <= 0.1, "FK_pvalue_new"] <- '0.1'
compare_matrix_data[compare_matrix_data$FK_pvalue <= 0.05, "FK_pvalue_new"] <- '0.05'

#filter out dataframe to make separate panels within the figure
comp_data_mod <- filter(compare_matrix_data, Comparison == 'Manner of Death')
comp_data_cod <- filter(compare_matrix_data, Comparison == 'Cause of Death')

comp_data_mod_w <- filter(comp_data_mod, Matrix == 'W_Unifrac')
comp_data_cod_W <- filter(comp_data_cod, Matrix == 'W_Unifrac')
comp_data_mod_un <- filter(comp_data_mod, Matrix == 'Unifrac')
comp_data_cod_un <- filter(comp_data_cod, Matrix == 'Unifrac')

#KW test results 
#each of the panels: a,b,c,d
a <- ggplot(comp_data_mod_un, aes(comp_data_mod_un$Body_Site, comp_data_mod_un$Level, 
                                fill= comp_data_mod_un$KW_pvalue_new)) + 
  geom_tile() + ylab('Method (Normalization Cutoff/Minimum Library Size)') + xlab('') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                      values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Unifrac', subtitle = 'Manner of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "none")

b <- ggplot(comp_data_mod_w, aes(comp_data_mod_w$Body_Site, comp_data_mod_w$Level, 
                                  fill= comp_data_mod_w$KW_pvalue_new)) + 
  geom_tile() + ylab('') + xlab('') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                    values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Weighted Unifrac', subtitle = 'Manner of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "none")

c <- ggplot(comp_data_cod_un, aes(comp_data_cod_un$Body_Site, comp_data_cod_un$Level, 
                                  fill= comp_data_cod_un$KW_pvalue_new)) + 
  geom_tile() + ylab('Method (Normalization Cutoff/Minimum Library Size)') + xlab('Body Site') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                    values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Unifrac', subtitle = 'Cause of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

d <- ggplot(comp_data_cod_W, aes(comp_data_cod_W$Body_Site, comp_data_cod_W$Level, 
                                  fill= comp_data_cod_W$KW_pvalue_new)) + 
  geom_tile() + ylab('') + xlab('Body Site') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                    values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Weighted Unifrac', subtitle = 'Cause of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

#final figures
theme_set(theme_classic(base_size = 18))
tiff("pvalues_KW_cod2.TIF", width = 4500, height = 4500, res=300)
ggarrange(a,b,c,d,
          labels = c('A', 'B', 'C', 'D'),
          nrow = 2, ncol = 2)
dev.off()

#FK test results 
#each of the panels: a,b,c,d
a <- ggplot(comp_data_mod_un, aes(comp_data_mod_un$Body_Site, comp_data_mod_un$Level, 
                                  fill= comp_data_mod_un$FK_pvalue_new)) + 
  geom_tile() + ylab('Method (Normalization Cutoff/Minimum Library Size)') + xlab('') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                    values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Unifrac', subtitle = 'Manner of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "none")

b <- ggplot(comp_data_mod_w, aes(comp_data_mod_w$Body_Site, comp_data_mod_w$Level, 
                                 fill= comp_data_mod_w$FK_pvalue_new)) + 
  geom_tile() + ylab('') + xlab('') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                    values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Weighted Unifrac', subtitle = 'Manner of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "none")

c <- ggplot(comp_data_cod_un, aes(comp_data_cod_un$Body_Site, comp_data_cod_un$Level, 
                                  fill= comp_data_cod_un$FK_pvalue_new)) + 
  geom_tile() + ylab('Method (Normalization Cutoff/Minimum Library Size)') + xlab('Body Site') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                    values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Unifrac', subtitle = 'Cause of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

d <- ggplot(comp_data_cod_W, aes(comp_data_cod_W$Body_Site, comp_data_cod_W$Level, 
                                 fill= comp_data_cod_W$FK_pvalue_new)) + 
  geom_tile() + ylab('') + xlab('Body Site') +
  scale_fill_manual(name = "P value", labels = c('Less than 0.05', 'Less than 0.1', 'Non-significant'),
                    values = c('#F1AC37', '#35CECE', '#E6EFEF')) + ggtitle(label= 'Weighted Unifrac', subtitle = 'Cause of Death') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + theme(legend.position = "bottom")

#final figure
theme_set(theme_classic(base_size = 18))
tiff("pvalues_FK_cod2.TIF", width = 4500, height = 4500, res=300)
ggarrange(a,b,c,d,
          labels = c('A', 'B', 'C', 'D'),
          nrow = 2, ncol = 2)
dev.off()


# pairwise comparisons ----------------------------------------------------

#these ended up being significant comparisons, look at pairwise comparisons 

# unifrac
#mod
# 10% eyes, mouth, nose
# 3% eyes, mouth
# 1% eyes, mouth
# 3,000 mouth, nose
# 5,000 mouth, nose
# 7,000 mouth, nose

#cod
#10% eyes, rectum
#3% eyes, rectum
#1% eyes, rectum
#3,000 eyes, nose, rectum
#5,000 eyes, rectum
#7,000 eyes, nose, rectum 

posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_3000_mou']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_3000_mou']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_3000_nos']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_3000_nos']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_5000_mou']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_5000_mou']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_5000_nos']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_5000_nos']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_7000_mou']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_7000_mou']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_7000_nos']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_7000_nos']]$MoD, p.adjust.method = 'bonf', dist='Tukey')

posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_cod[['betavalues_cod_1per_eye']]$distances, g = get_beta_formated_u_cod[['betavalues_cod_1per_eye']]$CoD_Simple2, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_1per_mou']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_1per_mou']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_3per_eye']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_3per_eye']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_3per_mou']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_3per_mou']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_10per_eye']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_10per_eye']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_10per_mou']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_10per_mou']]$MoD, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = get_beta_formated_u_mod[['betavalues_mod_10per_nos']]$distances, g = get_beta_formated_u_mod[['betavalues_mod_10per_nos']]$MoD, p.adjust.method = 'bonf', dist='Tukey')

posthoc.kruskal.nemenyi.test(x = betavalues_cod_3000_rec$distances, g = betavalues_cod_3000_rec$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_3000_eye$distances, g = betavalues_cod_3000_eye$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_3000_nos$distances, g = betavalues_cod_3000_nos$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_5000_rec$distances, g = betavalues_cod_5000_rec$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_5000_eye$distances, g = betavalues_cod_5000_eye$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_7000_rec$distances, g = betavalues_cod_7000_rec$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_7000_eye$distances, g = betavalues_cod_7000_eye$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_7000_nos$distances, g = betavalues_cod_7000_nos$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')

posthoc.kruskal.nemenyi.test(x = betavalues_cod_1per_eye$distances, g = betavalues_cod_1per_eye$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_1per_rec$distances, g = betavalues_cod_1per_rec$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_3per_eye$distances, g = betavalues_cod_3per_eye$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_3per_rec$distances, g = betavalues_cod_3per_rec$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_10per_eye$distances, g = betavalues_cod_10per_eye$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')
posthoc.kruskal.nemenyi.test(x = betavalues_cod_10per_rec$distances, g = betavalues_cod_10per_rec$CoD_Simple, p.adjust.method = 'bonf', dist='Tukey')



