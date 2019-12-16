packages
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

#load in files
#otu table, taxonomy table, and tree file
otu <- read.csv("table.csv")
tax <- read.csv("taxonomy.csv")
tree <- read_tree('tree.nwk')

#load in metadata
metadata=(read.csv("HPMMMetadata_pt2_cod.csv",header=TRUE))
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
#30906 taxa, 878 samples

#triming out taxa that are not representative of .01% of mean sequence number
physeq_trim <- prune_taxa(taxa_sums(physeq) > sum(otu) *.001 / 878, physeq)
#8692 taxa

#rarefy

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

names(beta_disp_list_u_mod) <- c('Nose_MOD', 'Rectum_MOD', 'Ears_MOD',
                                 'Mouth_MOD', 'Eyes_MOD')

names(beta_disp_list_u_cod) <- c('Nose_COD', 'Rectum_COD', 'Ears_COD',
                                 'Mouth_COD', 'Eyes_COD')


beta_disp_formatted_combo_cod <- rbind(get_beta_formated_u_cod[[1]], 
                                       get_beta_formated_u_cod[[2]],
                                       get_beta_formated_u_cod[[3]],
                                       get_beta_formated_u_cod[[4]],
                                       get_beta_formated_u_cod[[5]])

beta_disp_formatted_combo_mod <- rbind(get_beta_formated_u_mod[[1]], 
                                      get_beta_formated_u_mod[[2]],
                                      get_beta_formated_u_mod[[3]],
                                      get_beta_formated_u_mod[[4]],
                                      get_beta_formated_u_mod[[5]])

write.csv(beta_disp_formatted_combo_mod, 'beta_disp_mod.csv')
write.csv(beta_disp_formatted_combo_cod, 'beta_disp_cod.csv')

#check for correlation
summary(beta_disp_formatted_combo_mod)
myvars <- names(beta_disp_formatted_combo_mod) %in% c('distances', 'Sample_Area', 
                                                      'Sex', 'Race', 'Age', 
                                                      'FinePMI', 'BMI', 'MoD')

#convert to numeric because corr function requires it
beta_disp_formatted_combo_mod_cor <- beta_disp_formatted_combo_mod[myvars] 
beta_disp_formatted_combo_mod_cor <- beta_disp_formatted_combo_mod_cor[complete.cases(beta_disp_formatted_combo_mod_cor),]
beta_disp_formatted_combo_mod_cor$Sample_Area <- as.numeric(beta_disp_formatted_combo_mod_cor$Sample_Area)
beta_disp_formatted_combo_mod_cor$Sex <- as.numeric(beta_disp_formatted_combo_mod_cor$Sex)
beta_disp_formatted_combo_mod_cor$Race <- as.numeric(beta_disp_formatted_combo_mod_cor$Race)
beta_disp_formatted_combo_mod_cor$FinePMI <- as.numeric(beta_disp_formatted_combo_mod_cor$FinePMI)
beta_disp_formatted_combo_mod_cor$MoD <- as.numeric(beta_disp_formatted_combo_mod_cor$MoD)
a <- cor(beta_disp_formatted_combo_mod_cor, use='pairwise.complete.obs')
#not correlated! 
symnum(a)

bdfcmod_NA <- beta_disp_formatted_combo_mod %>% filter(MoD != 'Suicide')

long_df <- mlogit.data(bdfcmod_NA, choice = 'MoD', shape = 'wide')

model1 <- mlogit(MoD ~ 1 | distances + Sample_Area + Age + Race + Sex + FinePMI, data = long_df, reflevel = 'Natural')
summary(model1)

correct = model1$probabilities
binarycorrect = colnames(correct)[apply(correct,1,which.max)]
table(bdfcmod_NA$MoD, binarycorrect)

#total, manner of death 
  # p value is significant ( < 0.05), but McFadden is really low : 0.009. 
  # log likelihood: -1098.1, chi-squared in 20.261
  # accident and homicide are different from natural, as beta dispersion goes down, more likely to
  # have died by accident, than homicide. 
#most of the model is incorrectly classified. Accidents: 81.92%, natural: 20%. H and S: 0
# overall 37.87% better than random chance, but not great. 

272/(272 + 60) * 100
51/(204 + 51) * 100
(272+51) / nrow(beta_disp_formatted_combo_mod) * 100

myvars <- names(beta_disp_formatted_combo_mod) %in% c('distances', 'Sample_Area', 
                                                      'Sex', 'Race', 'Age', 
                                                      'FinePMI', 'BMI', 'CoD_Simple2')


beta_disp_formatted_combo_cod_cor <- beta_disp_formatted_combo_cod[myvars] 
beta_disp_formatted_combo_cod_cor <- beta_disp_formatted_combo_cod_cor[complete.cases(beta_disp_formatted_combo_cod_cor),]
beta_disp_formatted_combo_cod_cor$Sample_Area <- as.numeric(beta_disp_formatted_combo_cod_cor$Sample_Area)
beta_disp_formatted_combo_cod_cor$Sex <- as.numeric(beta_disp_formatted_combo_cod_cor$Sex)
beta_disp_formatted_combo_cod_cor$Race <- as.numeric(beta_disp_formatted_combo_cod_cor$Race)
beta_disp_formatted_combo_cod_cor$FinePMI <- as.numeric(beta_disp_formatted_combo_cod_cor$FinePMI)
beta_disp_formatted_combo_cod_cor$CoD_Simple2 <- as.numeric(beta_disp_formatted_combo_cod_cor$CoD_Simple2)
a <- cor(beta_disp_formatted_combo_cod_cor, use='pairwise.complete.obs')
#not correlated! 
symnum(a)

long_df <- mlogit.data(beta_disp_formatted_combo_cod, choice = 'CoD_Simple2', shape = 'wide')

model1 <- mlogit(CoD_Simple2 ~ 1 | distances, data = long_df, reflevel = 'Cardio')
summary(model1)

correct = model1$probabilities
binarycorrect = colnames(correct)[apply(correct,1,which.max)]
table(beta_disp_formatted_combo_mod$CoD_Simple2, binarycorrect)

#total, cause of death 
# p value is significant ( < 0.05), but McFadden is really low : 0.024. 
# log likelihood: -1325.4, chi-squared in 65.012
# asphyxiation, bft, drug related, and gunshot deaths go down as cardio goes up 
#most of the model is incorrectly classified. and classified as drug-related deaths.  

#separate out body sites
#change model variables as needed
multi_log_reg_mod <- function(df) {
  long_df <- mlogit.data(df, choice = 'MoD', shape = 'wide')
  model1 <- mlogit(MoD ~ 1 | distances + Age + Race + Sex + FinePMI, data = long_df, reflevel = 'Natural')
  print(summary(model1))
  correct = model1$probabilities
  binarycorrect = colnames(correct)[apply(correct,1,which.max)]
  return(table(df$MoD, binarycorrect))
}

multi_log_reg_mod(get_beta_formated_u_mod[[1]])
multi_log_reg_mod(get_beta_formated_u_mod[[2]])
multi_log_reg_mod(get_beta_formated_u_mod[[3]])
multi_log_reg_mod(get_beta_formated_u_mod[[4]])
multi_log_reg_mod(get_beta_formated_u_mod[[5]])

multi_log_reg_cod <- function(df) {
  long_df <- mlogit.data(df, choice = 'CoD_Simple2', shape = 'wide')
  model1 <- mlogit(CoD_Simple2 ~ 1 | distances + Age + Race + Sex + FinePMI, data = long_df, reflevel = 'Cardio')
  print(summary(model1))
  correct = model1$probabilities
  binarycorrect = colnames(correct)[apply(correct,1,which.max)]
  return(table(df$CoD_Simple2, binarycorrect))
}

multi_log_reg_cod(get_beta_formated_u_cod[[1]])
multi_log_reg_cod(get_beta_formated_u_cod[[2]])
multi_log_reg_cod(get_beta_formated_u_cod[[3]])
multi_log_reg_cod(get_beta_formated_u_cod[[4]])
multi_log_reg_cod(get_beta_formated_u_cod[[5]])

#try adding PCoA value for MOD
test <- get_beta_formated_u_mod[[4]]
test_merge <- merge(test, PCoA_values_mod_acc_mou, by = 'SampleID')
test_merge1 <- merge(test, PCoA_values_mod_hom_mou, by = 'SampleID')
test_merge2 <- merge(test, PCoA_values_mod_sui_mou, by = 'SampleID')
test_merge3 <- merge(test, PCoA_values_mod_nat_mou, by = 'SampleID')
merger <- rbind(test_merge, test_merge1, test_merge2, test_merge3)

long_df <- mlogit.data(merger, choice = 'MoD', shape = 'wide')
model1 <- mlogit(MoD ~ 1 | distances + Sums+ Age + Race + Sex + FinePMI, data = long_df, reflevel = 'Natural')
print(summary(model1))
correct = model1$probabilities
binarycorrect = colnames(correct)[apply(correct,1,which.max)]
table(merger$MoD, binarycorrect)


#now we are going to try specific pairwise compairsons with body area that performed the best


###plots
### natural vs accident
pc_mod_mlr_nos <- get_beta_formated_u_mod[[1]]
pc_mod_mlr_nos_N <- pc_mod_mlr_nos %>% filter(MoD == 'Natural')
pc_mod_mlr_nos_A <- pc_mod_mlr_nos %>% filter(MoD == 'Accident')
pc_mod_mlr_nos_NA <- rbind(pc_mod_mlr_nos_N, pc_mod_mlr_nos_A)



model = glm(MoD ~ distances +  Race + Event_Location + Age, family = binomial(link = 'logit'),
            data = pc_mod_mlr_nos_NA)
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
                       labels = c('Accident', 'Natural'))

table(pc_mod_mlr_nos_NA$MoD, binarycorrectlab)


pc_mod_mlr_nos_NA$binarycorrect <- binarycorrectlab
pc_mod_mlr_nos_NA$correctnum <- binarycorrect
pc_mod_mlr_nos_NA$correct <- 1

for(i in 1:length(rownames(pc_mod_mlr_nos_NA))) {
  if(pc_mod_mlr_nos_NA$MoD[i] == pc_mod_mlr_nos_NA$binarycorrect[i]) {
    pc_mod_mlr_nos_NA$correct[i] <- 'Correct' } else {
      pc_mod_mlr_nos_NA$correct[i] <- 'Incorrect' } 
} 

grob <- grobTree(textGrob("Accident", x=0.1,  y=0.1, hjust=0,
                          gp=gpar(col="#FF5533", fontsize=13)))
grob2 <- grobTree(textGrob("Natural", x=0.7,  y=0.9, hjust=0,
                          gp=gpar(col="#FFBB33", fontsize=13)))

a <- ggplot(pc_mod_mlr_nos_NA, aes(x=distances, y=correctnum)) +
  geom_point(size = 3 , alpha = 1/5) + ylab('Predictions') + xlab('Beta-dispersion') +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = T) + annotation_custom(grob) + annotation_custom(grob2)
a

sampdat=sample_data(pc_mod_mlr_nos_NA)
sample_names(sampdat)=pc_mod_mlr_nos_NA$SampleID
sampdat$SampleID <- as.character(sampdat$SampleID)
physeq_N_nos <- subset_samples(physeq_by_bodysites[['Nose']], MoD == 'Natural')
physeq_A_nos <- subset_samples(physeq_by_bodysites[['Nose']], MoD == 'Accident')
physeq_NA_nos <- merge_phyloseq(physeq_N_nos, physeq_A_nos, sampdat)

#plot
ord = ordinate(physeq_NA_nos, method="PCoA", distance="unifrac")
ordplot=plot_ordination(physeq_NA_nos, ord, color="MoD")
d <- ordplot + geom_point(size = 3, aes(color = sample_data(physeq_NA_nos)$MoD, shape = sample_data(physeq_NA_nos)$correct)) + 
  scale_color_manual('Manner of Death' ,values = c("#FFBB33", "#FF5533")) + 
  scale_shape_manual('Prediction', values = c(15,19))
d


### cardio vs. drug-related
pc_cod_mlr_nos <- get_beta_formated_u_cod[[1]]
pc_cod_mlr_nos_C <- pc_cod_mlr_nos %>% filter(CoD_Simple2 == 'Cardio')  
pc_cod_mlr_nos_D <- pc_cod_mlr_nos %>% filter(CoD_Simple2 == 'Drug')  
pc_cod_mlr_nos_CD <- rbind(pc_cod_mlr_nos_C, pc_cod_mlr_nos_D)


model1 = glm(CoD_Simple2 ~ distances + BMI + Race + Event_Location + BroadPMI + Age, 
            family = binomial(link = 'logit'),
            data = pc_cod_mlr_nos_CD)
summary(model1)


#correct predictions
correct = model1$fitted.values
binarycorrect = ifelse(correct > 0.5, 1,0)
binarycorrectlab = factor(binarycorrect,
                          levels = c(0,1),
                          labels = c('Cardio', 'Drug'))

table(pc_cod_mlr_nos_CD$CoD_Simple2, binarycorrectlab)


pc_cod_mlr_nos_CD$binarycorrect <- binarycorrectlab
pc_cod_mlr_nos_CD$correctnum <- binarycorrect
pc_cod_mlr_nos_CD$correct <- 1

for(i in 1:length(rownames(pc_cod_mlr_nos_CD))) {
  if(pc_cod_mlr_nos_CD$CoD_Simple2[i] == pc_cod_mlr_nos_CD$binarycorrect[i]) {
    pc_cod_mlr_nos_CD$correct[i] <- 'Correct' } else {
      pc_cod_mlr_nos_CD$correct[i] <- 'Incorrect' } 
} 

grob3 <- grobTree(textGrob("Cardiovascular Disease", x=0.1,  y=0.1, hjust=0,
                          gp=gpar(col="#5DD06D", fontsize=13)))
grob4 <- grobTree(textGrob("Drug-related Deaths", x=0.6,  y=0.9, hjust=0,
                           gp=gpar(col="#5D8BD0", fontsize=13)))

b <- ggplot(pc_cod_mlr_nos_CD, aes(x=distances, y=correctnum)) +
  geom_point(size = 3 , alpha = 1/5) + ylab('Predictions') + xlab('Beta-dispersion') +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = T) + annotation_custom(grob3) + annotation_custom(grob4)
b

sampdat=sample_data(pc_cod_mlr_nos_CD)
sample_names(sampdat)=pc_cod_mlr_nos_CD$SampleID
sampdat$SampleID <- as.character(sampdat$SampleID)
physeq_C_nos <- subset_samples(physeq_by_bodysites[['Nose']], CoD_Simple2 == 'Cardio')
physeq_D_nos <- subset_samples(physeq_by_bodysites[['Nose']], CoD_Simple2 == 'Drug')
physeq_CD_nos <- merge_phyloseq(physeq_C_nos, physeq_D_nos, sampdat)

#plot
ord = ordinate(physeq_CD_nos, method="PCoA", distance="unifrac")
ordplot=plot_ordination(physeq_CD_nos, ord, color="CoD_Simple2")
e <- ordplot + geom_point(size = 3, aes(color = sample_data(physeq_CD_nos)$CoD_Simple2, 
                                   shape = sample_data(physeq_CD_nos)$correct)) + 
  scale_color_manual('Cause of Death' ,values = c("#5DD06D", "#5D8BD0")) + 
  scale_shape_manual('Prediction', values = c(15,19))
e

### disease
pc_cod_mlr_nos <- get_beta_formated_u_cod[[1]]
pc_cod_mlr_nos[pc_cod_mlr_nos$MoD == 'Natural', "Disease"] <- 'Diseased'
pc_cod_mlr_nos[pc_cod_mlr_nos$MoD != 'Natural', "Disease"] <- 'Not'
pc_cod_mlr_nos$Disease <- as.factor(pc_cod_mlr_nos$Disease)

model2 = glm(Disease ~ distances + BMI + Race + Event_Location + Age, 
             family = binomial(link = 'logit'),
             data = pc_cod_mlr_nos)
summary(model2)


#correct predictions
correct = model2$fitted.values
binarycorrect = ifelse(correct > 0.5, 1,0)
binarycorrectlab = factor(binarycorrect,
                          levels = c(0,1),
                          labels = c('Diseased', 'Not'))

table(pc_cod_mlr_nos$Disease, binarycorrectlab)


pc_cod_mlr_nos$binarycorrect <- binarycorrectlab
pc_cod_mlr_nos$correctnum <- binarycorrect
pc_cod_mlr_nos$correct <- 1

for(i in 1:length(rownames(pc_cod_mlr_nos))) {
  if(pc_cod_mlr_nos$Disease[i] == pc_cod_mlr_nos$binarycorrect[i]) {
    pc_cod_mlr_nos$correct[i] <- 'Correct' } else {
      pc_cod_mlr_nos$correct[i] <- 'Incorrect' } 
} 

grob5 <- grobTree(textGrob("Diseased", x=0.1,  y=0.1, hjust=0,
                           gp=gpar(col="#935DD0", fontsize=13)))
grob6 <- grobTree(textGrob("Not Diseased", x=0.7,  y=0.9, hjust=0,
                           gp=gpar(col="#D0A95D", fontsize=13)))

c <- ggplot(pc_cod_mlr_nos, aes(x=distances, y=correctnum)) +
  geom_point(size = 3 , alpha = 1/5) + ylab('Predictions') + xlab('Beta-dispersion') +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = T) + annotation_custom(grob5) + annotation_custom(grob6)
c

sampdat=sample_data(pc_cod_mlr_nos)
sample_names(sampdat)=pc_cod_mlr_nos$SampleID
sampdat$SampleID <- as.character(sampdat$SampleID)
physeq_nos <- merge_phyloseq(physeq_by_bodysites[['Nose']], sampdat)

#plot
ord = ordinate(physeq_nos, method="PCoA", distance="unifrac")
ordplot=plot_ordination(physeq_nos, ord, color="Disease")
f <- ordplot + geom_point(size = 3, aes(color = sample_data(physeq_nos)$Disease, 
                                        shape = sample_data(physeq_nos)$correct)) + 
  scale_color_manual('Disease Status' ,values = c("#935DD0", "#D0A95D")) + 
  scale_shape_manual('Prediction', values = c(15,19))
f

theme_set(theme_classic(base_size = 20))
tiff("pc_mlr_betadisp.TIF", width = 6000, height = 3500, res=300)
ggarrange(a,b,c,d,e,f,
          labels = c('A', 'B', 'C', 'D', 'E', 'F'),
          nrow = 2, ncol = 3)
dev.off()
