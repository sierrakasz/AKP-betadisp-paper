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
set.seed(1234)

otu <- read.csv("table.csv")
tax <- read.csv("taxonomy.csv")
tree <- read_tree('tree.nwk')

#make a table in excel for rare data
metadata=(read.csv("HPMMMetadata_pt2_cod.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

rownames(otu) <- otu$OTUID
otu <- otu[,-1]
OTU=otu_table(otu, taxa_are_rows=TRUE)
physeq_otu.tree=phyloseq(OTU,tree, sampdat)

rownames(tax) <- tax$OTUID
tax_reduc <- merge(tax, otu, by = "row.names")
rownames(tax_reduc) <- tax_reduc$Row.names
tax_reduc <- tax_reduc[,-1]
tax_reduc <- tax_reduc[,-1]

tax_f <- tax_reduc[,1:7]
tax_f <- as.matrix(tax_f)
TAX=tax_table(tax_f)
taxa_names(TAX)=row.names(OTU)

physeq_otu.tax=phyloseq(OTU,TAX, sampdat)

#normalize
normalize_wout_rarefying <- function(physeq) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  sample_nums <- rowSums(otu == 0)
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  len <- length(colnames(otu)) -1
  core_perc <- sample_nums / len
  new_otu <- cbind(otu, core_perc)
  new_core_taxa <- new_otu$OTUID[new_otu$core_perc <= 0.99]
  new_core_otu <- filter(new_otu, new_otu$OTUID %in% new_core_taxa)
  new_core_otu <- subset(new_core_otu, select=-c(OTUID, core_perc))
  new_core_tax <- filter(tax, tax$OTUID %in% new_core_taxa)
  dis_core_tax <- filter(tax, !tax$OTUID %in% new_core_tax$OTUID)
  newphyseq <- prune_taxa(new_core_tax$OTUID, physeq)
  return(newphyseq)
}

physeq_norm <- normalize_wout_rarefying(physeq_otu.tax)
physeq_norm_tree <- merge_phyloseq(physeq_norm, tree)

#rarefy
physeq_beta <- rarefy_even_depth(physeq_otu.tree, sample.size = 5000)
physeq <- rarefy_even_depth(physeq_otu.tax, sample.size = 5000)


#check age effect
ord = ordinate(physeq_beta, method="PCoA", distance="wunifrac")
ordplot=plot_ordination(physeq_beta, ord, color="Age_bin_10")
ordplot
ordplot + geom_point(size = 4) + 
  stat_ellipse(alpha = 0.4, geom = "polygon", aes(fill = Age_bin_10))



  GPdist=phyloseq::distance(physeq_beta, "unifrac")
  sampledf <- data.frame(sample_data(physeq_beta))
  
  sampledf_new <- sampledf %>% filter(Age < 77)
  
  beta <- betadisper(GPdist, sampledf$Age_bin_10)
  permutest(beta)

  beta <- betadisper(GPdist, sampledf$Age_bin_25)
  permutest(beta)
  
  
  beta <- betadisper(GPdist, sampledf$Age_bin_35)
  permutest(beta)
  
  
  beta <- betadisper(GPdist, sampledf$MoD)
  permutest(beta)

  beta_df <- as.data.frame(beta$distances)
  beta_df$SampleID <- rownames(beta_df)
  colnames(beta_df) <- c('distances', 'SampleID')
  beta_sa <- merge(beta_df, metadata, by = 'SampleID')
  
  plot(y = beta_sa$distances, x = beta_sa$Age_bin_10)  
  