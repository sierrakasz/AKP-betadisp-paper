## Anna Karenina Kaszubinski et al.

#Figures


#packages
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
library(pBrackets)
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


# Load files --------------------------------------------------------------


#otu table, taxonomy table, and tree file
otu <- read.csv("table.csv")
tax <- read.csv("taxonomy.csv")
tree <- read_tree('tree.nwk')

#load in metadata
metadata=(read.csv("HPMMMetadata_pt2_cod.csv",header=TRUE))
metadata$BMI <- as.numeric(as.character((metadata$BMI)))
metadata$SampleID <- as.character((metadata$SampleID))
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


# Fig. 1 ------------------------------------------------------------------


#combine beta-dispersion values into a dataframe 
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

betadisp_mod_df$Sample_Area <- factor(betadisp_mod_df$Sample_Area, levels = c('Eyes', 'Ears', 'Nose',
                                                                              'Rectum', 'Mouth'))
a <- ggplot(betadisp_mod_df, aes(x=Sample_Area, y=distances)) + 
  geom_boxplot(lwd = 1, aes(fill = Sample_Area)) + ylab('Beta-dispersion among Manners of Death') + 
  scale_fill_manual(values= c('#C27C1D', '#0D6EBB', '#F2BE78','#18B3FE', '#95918B')) + 
  labs(color = 'Body Site') + theme(legend.position="none") + theme(axis.title.x=element_blank())

a

betadisp_mod_df$MoD <- factor(betadisp_mod_df$MoD, levels = c('Natural', 'Suicide', 
                                                              'Accident', 'Homicide'))
c <- ggplot(betadisp_mod_df, aes(x=MoD, y=distances)) + 
  geom_boxplot(lwd = 1, aes(fill = MoD)) + ylab('Beta-dispersion among Manners of Death') + 
  scale_fill_manual(values= c('#29B26A','#4CB47D', '#70B38F', '#8CB29D'), name = 'Manners of Death') +
  theme(legend.position="none") + theme(axis.title.x=element_blank())

c

betadisp_cod_df$Sample_Area <- factor(betadisp_cod_df$Sample_Area, levels = c('Eyes', 'Ears', 'Nose',
                                                                              'Rectum', 'Mouth'))
betadisp_cod_df <- betadisp_cod_df %>% filter(Pack_ID != '2014-S39')

b <- ggplot(betadisp_cod_df, aes(x=Sample_Area, y=distances)) + 
  geom_boxplot(lwd = 1, aes(fill = Sample_Area)) + ylab('Beta-dispersion among Causes of Death') + 
  scale_fill_manual(values= c('#C27C1D', '#0D6EBB', '#F2BE78','#18B3FE', '#95918B')) + 
  labs(color = 'Body Site') + theme(axis.title.x=element_blank()) + 
  theme(legend.position="none")

b


betadisp_cod_df %>% group_by(CoD_Simple2) %>% summarize_at(c('distances'), funs(mean,sd))

betadisp_cod_df$CoD_Simple2 <- factor(betadisp_cod_df$CoD_Simple2, levels = c('Cardio', 'Other', 
                                                                              'Drug', 'BFT', 'Gunshot',
                                                                              'Asphyx', 'Unknown'))
d <- ggplot(betadisp_cod_df, aes(x=CoD_Simple2, y=distances)) + 
  geom_boxplot(lwd=1, aes(fill = CoD_Simple2)) + ylab('Beta-dispersion among Causes of Death') +
  scale_fill_manual(values = c('#5B2FB7', '#6641B3', '#7556B4',
                               '#866FB4', '#998BB6', '#A39BB4', '#B0ADB7')) + 
  theme(legend.position="none") + theme(axis.title.x=element_blank())
d

#final figure
theme_set(theme_classic(base_size = 18))
tiff("FIG1.TIFF", width = 3500, height = 3500, res=300)
ggarrange(a,e,d,g,
          labels = c('A', 'B', 'C', 'D'),
          nrow = 2, ncol = 2)
dev.off()


# Fig 2 -------------------------------------------------------------------

#pulled data from Table S9
#data found in Table S10
#code reproduces figure, not full analyses. See supplemental tables code for 
#analyses 

sam_numb <- as.numeric(1:10)
beta_disp <- c(rep('Full communities', 5), rep('Random forest indicators', 5))
body_site <- c(rep(c('Nose', 'Mouth', 'Nose', 'Ears', 'Mouth'), 2))
of_Deaths <- c(rep(c('MOD', 'MOD', 'COD', 'COD', 'COD'), 2))
mcfadden <- c(0.277, 0.250, 0.361, 0.304, 0.298,
              0.310, 0.255, 0.377, 0.359, 0.291)
percent_correct <- c(61.0, 55.3, 62.8, 62.9, 60.6,
                     61.0, 54.1, 59.9, 41.3, 54.7)
n <- c(172, 170, 172, 167, 170,
       172, 170, 172, 167, 170)

df_for_graph <- as.data.frame(cbind(sam_numb, beta_disp, body_site, of_Deaths, mcfadden,
                                    percent_correct, n))
df_for_graph$sam_numb <- factor(df_for_graph$sam_numb, levels = c(1,2,3,4,5
                                                                  ,6,7,8,9,10))
df_for_graph$percent_correct <- as.numeric(as.character(df_for_graph$percent_correct))
df_for_graph$mcfadden <- as.numeric(as.character(df_for_graph$mcfadden))
df_for_graph$of_Deaths <- factor(df_for_graph$of_Deaths, levels = c('MOD', 'COD'))


grob <- grobTree(textGrob("Full Communities", x=0.07,  y=0.02, 1, hjust=0,
                          gp=gpar(col="black", fontsize=14)))
grob2 <- grobTree(textGrob("RF Indicators", x=0.62,  y=0.02, hjust=0,
                           gp=gpar(col="black", fontsize=14)))

theme_set(theme_classic(base_size = 20))
newlabs <- c("n=172", "n=170", "n=172", "n=167", "n=170",
             "n=172", "n=170", "n=172", "n=167", "n=170")

p <- ggplot(df_for_graph, aes(sam_numb, percent_correct, fill = body_site)) + geom_col() +
  scale_fill_manual("Body Sites", values = c("Nose" = "#F2BE78", "Mouth" = "#95918B", "Ears" = "#0D6EBB")) +
  facet_wrap(~of_Deaths, scales = 'free_x', strip.position = "bottom") + 
  ylab('Percent Correct') + xlab("") +
  annotation_custom(grob) + annotation_custom(grob2) + 
  scale_x_discrete(labels= newlabs) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

df_for_graph$sam_numb <- factor(df_for_graph$sam_numb, levels = c(1,2,6,7,3,4,5,8,9,10))
p2 <- ggplot(df_for_graph, aes(sam_numb, mcfadden)) + geom_line(group = 10) +
  geom_point(aes(color = body_site), size = 4) +
  scale_color_manual("Body Sites", 
                     values = c("Nose" = "#F2BE78", "Mouth" = "#95918B", "Ears" = '#0D6EBB')) +
  ylab('McFadden R-squared') + theme(axis.line.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                     axis.title.x=element_blank()
  )

p2

tiff("FIG2.TIF", width = 3500, height = 3500, res=300)
ggarrange(p2, p, nrow = 2)
dev.off()

# Fig. 3 ------------------------------------------------------------------
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
                          gp=gpar(col="#FF5533", fontsize=20)))
grob2 <- grobTree(textGrob("Natural", x=0.7,  y=0.9, hjust=0,
                           gp=gpar(col="#FFBB33", fontsize=20)))

a <- ggplot(pc_mod_mlr_nos_NA, aes(x=distances, y=correctnum)) +
  geom_point(size = 3 , alpha = 1/5) + ylab('Classification') + xlab('Beta-dispersion') +
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
  scale_shape_manual('Classification', values = c(15,19)) + theme(legend.position = "bottom",
                                                                  legend.box = "vertical")
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
                           gp=gpar(col="#5DD06D", fontsize=20)))
grob4 <- grobTree(textGrob("Drug-related", x=0.6,  y=0.9, hjust=0,
                           gp=gpar(col="#5D8BD0", fontsize=20)))

b <- ggplot(pc_cod_mlr_nos_CD, aes(x=distances, y=correctnum)) +
  geom_point(size = 3 , alpha = 1/5) + ylab('Classification') + xlab('Beta-dispersion') +
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
  scale_shape_manual('Classification', values = c(15,19)) +  theme(legend.position = "bottom",
                                                                   legend.box = "vertical")
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
                           gp=gpar(col="#935DD0", fontsize=20)))
grob6 <- grobTree(textGrob("Not Diseased", x=0.7,  y=0.9, hjust=0,
                           gp=gpar(col="#D0A95D", fontsize=20)))

c <- ggplot(pc_cod_mlr_nos, aes(x=distances, y=correctnum)) +
  geom_point(size = 3 , alpha = 1/5) + ylab('Classification') + xlab('Beta-dispersion') +
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
  scale_shape_manual('Classification', values = c(15,19)) +  theme(legend.position = "bottom",
                                                                   legend.box = "vertical")
f

theme_set(theme_classic(base_size = 20))
tiff("FIG3.TIF", width = 6000, height = 4000, res=300)
ggarrange(a,b,c,d,e,f,
          labels = c('A', 'B', 'C', 'D', 'E', 'F'),
          nrow = 2, ncol = 3)
dev.off()

# Fig. 4 ------------------------------------------------------------------

#final figure was made in powerpoint
