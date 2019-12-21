
load("C:/Users/sierr/Documents/ThesisProject_AKP/R Scripts/AKP_session.RData")
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

#plot beta-dispersion values
plot(beta_disp_list_u_mod[[1]])
boxplot(beta_disp_list_u_mod[[1]])

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

# variables correlated among themselves -----------------------------------

#linear regression of Age and BMI 
m1 <- lm(formula = as.numeric(BMI) ~ as.numeric(Age), data = metadata)
summary(m1)
cor.test(metadata$Age, metadata$BMI, method=c("spearman"))

theme_set(theme_classic(base_size = 18))
tiff("age-bmi-reg.TIF", width = 3000, height = 3000, res=300)
p <- ggplot(metadata, aes(Age, BMI)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  ggtitle(label = 'BMI ~ -0.0640 x Age + 31.6',
       subtitle = 'R-squared: 0.0127   P-value: < 0.001 rho: -0.0972') + theme(
         plot.title = element_text(hjust = 0.5),
         plot.subtitle = element_text(hjust = 0.5)
       ) + ylab('BMI (kg/m^2)')

p
dev.off()

# very weak, but significant negative correlation with Age and BMI 
#control for age and control for BMI


# variables correlated with MOD/COD ---------------------------------------
#collapse metadata into case rather than by sample area for case metrics
metadata_coll <- subset(metadata, select=-c(SampleID, X, Sample_Area, Description))
metadata_coll <- unique(metadata_coll)
metadata_coll <- metadata_coll[complete.cases(metadata_coll), ]

#interested in how age, race, sex, and BMI is distributed throughout Manners of Death (MoD) and Causes of Death (CoD)
metadata_coll %>% group_by(MoD) %>% summarize_at(c('Age'), funs(mean,median,sd))
metadata_coll %>% group_by(CoD_Simple2) %>% summarize_at(c('Age'), funs(mean,median,sd))
metadata_coll %>% group_by(MoD) %>% summarize_at(c('BMI'), funs(mean,median,sd))
metadata_coll %>% group_by(CoD_Simple2) %>% summarize_at(c('BMI'), funs(mean,median,sd))
metadata_coll %>% gather(observation, Val, Sex) %>%group_by(MoD,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, Sex) %>%group_by(CoD_Simple2,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, Race) %>%group_by(MoD,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, Race) %>%group_by(CoD_Simple2,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, BroadPMI) %>%group_by(MoD,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, BroadPMI) %>%group_by(CoD_Simple2,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, Event_Location) %>%group_by(MoD,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, Event_Location) %>%group_by(CoD_Simple2,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, Season) %>%group_by(MoD,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)
metadata_coll %>% gather(observation, Val, Season) %>%group_by(CoD_Simple2,observation, Val) %>% 
  summarise(n= n()) %>% ungroup() %>% spread(Val, n, fill=0)

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

# variables correlated with beta-disp -------------------------------------
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

#building linear models to get an idea of what variables are important for describing beta-dispersion
#includes variables above
l0 = lm(distances ~ 1, data=betadisp_mod_df)
l1 = lm(distances ~ Sample_Area, data=betadisp_mod_df)
l2 = lm(distances ~ Sample_Area + Age, data=betadisp_mod_df)
l3 = lm(distances ~ Sample_Area + BMI, data=betadisp_mod_df)
l4 = lm(distances ~ Sample_Area + Sex, data=betadisp_mod_df)
l5 = lm(distances ~ Sample_Area + Race, data=betadisp_mod_df)
l6 = lm(distances ~ Sample_Area + FinePMI, data=betadisp_mod_df)
l7 = lm(distances ~ Sample_Area + MoD, data=betadisp_mod_df)
l8 = lm(distances ~ Sample_Area + Age + BMI, data=betadisp_mod_df)
l9 = lm(distances ~ Sample_Area + Age + FinePMI, data=betadisp_mod_df)
l10 = lm(distances ~ Sample_Area + Age + FinePMI + MoD, data=betadisp_mod_df)

anova(l0,l1,l2,l4,l5,l6,l7,l9,l10)
anova(l3)
AIC(l0,l1,l2,l4,l5,l6,l7,l9,l10)
AIC(l3, l8)
#best model for MOD: l10. Includes body site, age, fine PMI, and MOD

p0 = lm(distances ~ 1, data=betadisp_cod_df)
p1 = lm(distances ~ Sample_Area, data=betadisp_cod_df)
p2 = lm(distances ~ Sample_Area + Age, data=betadisp_cod_df)
p3 = lm(distances ~ Sample_Area + BMI, data=betadisp_cod_df)
p4 = lm(distances ~ Sample_Area + Sex, data=betadisp_cod_df)
p5 = lm(distances ~ Sample_Area + Race, data=betadisp_cod_df)
p6 = lm(distances ~ Sample_Area + FinePMI, data=betadisp_cod_df)
p7 = lm(distances ~ Sample_Area + CoD_Simple2, data=betadisp_cod_df)
p8 = lm(distances ~ Sample_Area + FinePMI + CoD_Simple2, data=betadisp_cod_df)
p9 = lm(distances ~ Sample_Area + FinePMI + CoD_Simple2 + BMI, data=betadisp_cod_df)
p10 = lm(distances ~ Sample_Area + CoD_Simple2 + BMI, data=betadisp_cod_df)

anova(p0,p1,p2,p4,p5,p6,p7,p8)
anova(p3,p9,p10)
AIC(p0,p1,p2,p4,p5,p6,p7,p8)
AIC(p3,p9,p10)

#best mode for COD: includes body site, fine PMI, and COD 

#figures
betadisp_mod_df %>% group_by(Sample_Area) %>% summarize_at(c('distances'), funs(mean,sd))

betadisp_mod_df$Sample_Area <- factor(betadisp_mod_df$Sample_Area, levels = c('Eyes', 'Ears', 'Nose',
                                                                              'Rectum', 'Mouth'))
a <- ggplot(betadisp_mod_df, aes(x=Sample_Area, y=distances, color = Sample_Area)) + 
  geom_boxplot() + geom_point(position = 'jitter') + ylab('Beta-dispersion among Manners of Death') + 
  scale_color_manual(values= c('#37D6F1','#F1B337', '#37F196', '#8C37F1', '#F1375C')) + 
  labs(color = 'Body Site') + theme(legend.position="right") + theme(axis.title.x=element_blank(),
axis.text.x=element_blank())
  
a

b <- ggplot(betadisp_mod_df, aes(x=Age, y=distances)) + 
  geom_point() + xlab('Age') + ylab('Beta-dispersion among Manners of Death') + geom_smooth(method = "lm")

b

betadisp_mod_df$FinePMI <- factor(betadisp_mod_df$FinePMI, levels = c('Less24', '25-48', '49-72', 'Great73'))

c <- ggplot(betadisp_mod_df, aes(x=FinePMI, y=distances, fill = FinePMI)) + 
  geom_boxplot() + ylab('Beta-dispersion among Manners of Death') + 
  scale_fill_manual(values= c('#98F137','#FED718', '#FE7918', '#F1374B'),
                    name = "PMI Estimate", 
                    labels = c("< 24 hrs", "25-48 hrs", "49-72 hrs", '> 73 hrs')) + 
  theme(legend.position="bottom") + theme(axis.title.x=element_blank(),
                                          axis.text.x=element_blank())

c
betadisp_mod_df %>% group_by(MoD) %>% summarize_at(c('distances'), funs(mean,sd))

betadisp_mod_df$MoD <- factor(betadisp_mod_df$MoD, levels = c('Natural', 'Suicide', 
                                                              'Accident', 'Homicide'))
d <- ggplot(betadisp_mod_df, aes(x=MoD, y=distances, color = MoD)) + 
  geom_boxplot() + geom_point(position = 'jitter') + ylab('Beta-dispersion among Manners of Death') + 
  scale_color_manual(values= c('#B9FE18','#18B3FE', '#A118FE', '#FE4C18'), name = 'Manners of Death') +
  theme(legend.position="right") + theme(axis.title.x=element_blank(),
                                          axis.text.x=element_blank())

d

betadisp_cod_df %>% group_by(Sample_Area) %>% summarize_at(c('distances'), funs(mean,sd))

betadisp_cod_df$Sample_Area <- factor(betadisp_cod_df$Sample_Area, levels = c('Eyes', 'Ears', 'Nose',
                                                                              'Rectum', 'Mouth'))
betadisp_cod_df <- betadisp_cod_df %>% filter(Pack_ID != '2014-S39')

e <- ggplot(betadisp_cod_df, aes(x=Sample_Area, y=distances, color = Sample_Area)) + 
  geom_boxplot() + geom_point(position = 'jitter') + ylab('Beta-dispersion among Causes of Death') + 
  scale_color_manual(values= c('#37D6F1','#F1B337', '#37F196', '#8C37F1', '#F1375C')) + 
  labs(color = 'Body Site') + theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + 
  theme(legend.position="right")

e

betadisp_cod_df$FinePMI <- factor(betadisp_cod_df$FinePMI, levels = c('Less24', '25-48', '49-72', 'Great73'))

f <- ggplot(betadisp_cod_df, aes(x=FinePMI, y=distances, fill = FinePMI)) + 
  geom_boxplot() + ylab('Beta-dispersion among Causes of Death') + 
  scale_fill_manual(values= c('#98F137','#FED718', '#FE7918', '#F1374B'),
                    name = "PMI Estimate", 
                    labels = c("< 24 hrs", "25-48 hrs", "49-72 hrs", '> 73 hrs')) + 
  theme(legend.position="bottom") + theme(axis.title.x=element_blank(),
                                          axis.text.x=element_blank())

f
betadisp_cod_df %>% group_by(CoD_Simple2) %>% summarize_at(c('distances'), funs(mean,sd))

betadisp_cod_df$CoD_Simple2 <- factor(betadisp_cod_df$CoD_Simple2, levels = c('Cardio', 'Other', 
                                                              'Drug', 'BFT', 'Gunshot',
                                                              'Asphyx', 'Unknown'))
g <- ggplot(betadisp_cod_df, aes(x=CoD_Simple2, y=distances, color = CoD_Simple2)) + 
  geom_boxplot() + geom_point(position = 'jitter') + ylab('Beta-dispersion among Causes of Death') +
  scale_color_manual(values = c('#BE8111','#11BEB2','#117600','#BE9611',
                               '#BE11A1','#BE1E11', '#1C0117'),
    name = "Causes of Death", 
                    labels = c("Cardiovascular Disease", "Other", "Drug-related Deaths",
                               'Blunt Force Trauma', 'Gunshot Wound', 'Asphyxiation',
                               'Unknown')) + 
  theme(legend.position="right") + theme(axis.title.x=element_blank(),
                                          axis.text.x=element_blank())
g

#final figure
theme_set(theme_classic(base_size = 20))
tiff("betadisp_bs_modcod.TIF", width = 5000, height = 3500, res=300)
ggarrange(a,e,d,g,
          labels = c('A', 'B', 'C', 'D'),
          nrow = 2, ncol = 2)
dev.off()

# microbial taxa associated with beta-disp extremes -----------------------
#for each body site and comparison (mod/cod), going to label beta-dispersion distances categorically
#based on quantiles of data


#split out beta-dispersion values by MODs within body site 
mods_accidents <- vector('list')
mods_homicides <- vector('list')
mods_natural <- vector('list')
mods_suicides <- vector('list')

for(i in 1:5) {
  a <- get_beta_formated_u_mod[[i]] %>% filter(MoD == 'Accident')
  mods_accidents[[i]] <- a
  b <- get_beta_formated_u_mod[[i]] %>% filter(MoD == 'Homicide')
  mods_homicides[[i]] <- b
  c <- get_beta_formated_u_mod[[i]] %>% filter(MoD == 'Natural')
  mods_natural[[i]] <- c
  d <- get_beta_formated_u_mod[[i]] %>% filter(MoD == 'Suicide')
  mods_suicides[[i]] <- d
}

#split out beta-dispersion values by CODs within body site
cods_asp <- vector('list')
cods_bft <- vector('list')
cods_drug <- vector('list')
cods_gun <- vector('list')
cods_cardio <- vector('list')
cods_other <- vector('list')

for(i in 1:5) {
  a <- get_beta_formated_u_cod[[i]] %>% filter(CoD_Simple2 == 'Asphyx')
 cods_asp[[i]] <- a
 b <- get_beta_formated_u_cod[[i]] %>% filter(CoD_Simple2 == 'BFT')
 cods_bft[[i]] <- b
 c <- get_beta_formated_u_cod[[i]] %>% filter(CoD_Simple2 == 'Drug')
 cods_drug[[i]] <- c
 d <- get_beta_formated_u_cod[[i]] %>% filter(CoD_Simple2 == 'Gunshot')
 cods_gun[[i]] <- d
 e <- get_beta_formated_u_cod[[i]] %>% filter(CoD_Simple2 == 'Cardio')
 cods_cardio[[i]] <- e
 f <- get_beta_formated_u_cod[[i]] %>% filter(CoD_Simple2 == 'Other')
 cods_other[[i]] <- f
}

#going to add new metadata to correpsonding phyloseq object
physeq_bodysite_acc <- vector('list')
physeq_bodysite_hom <- vector('list')
physeq_bodysite_nat <- vector('list')
physeq_bodysite_sui <- vector('list')

for(i in 1:5) {
  sampdat=sample_data(mods_accidents[[i]])
  sample_names(sampdat)=mods_accidents[[i]]$SampleID
  phy_acc <- subset_samples(physeq_by_bodysites[[i]], MoD == 'Accident')
  physeq_bodysite_acc[[i]] <- merge_phyloseq(phy_acc, sampdat)
  
  sampdat=sample_data(mods_homicides[[i]])
  sample_names(sampdat)=mods_homicides[[i]]$SampleID
  phy_hom <- subset_samples(physeq_by_bodysites[[i]], MoD == 'Homicide')
  physeq_bodysite_hom[[i]] <- merge_phyloseq(phy_hom, sampdat)
  
  sampdat=sample_data(mods_natural[[i]])
  sample_names(sampdat)=mods_natural[[i]]$SampleID
  phy_nat <- subset_samples(physeq_by_bodysites[[i]], MoD == 'Natural')
  physeq_bodysite_nat[[i]] <- merge_phyloseq(phy_nat, sampdat)
  
  sampdat=sample_data(mods_suicides[[i]])
  sample_names(sampdat)=mods_suicides[[i]]$SampleID
  phy_sui <- subset_samples(physeq_by_bodysites[[i]], MoD == 'Suicide')
  physeq_bodysite_sui[[i]] <- merge_phyloseq(phy_sui, sampdat)
}

physeq_bodysite_asp <- vector('list')
physeq_bodysite_bft <- vector('list')
physeq_bodysite_gun <- vector('list')
physeq_bodysite_drug <- vector('list')
physeq_bodysite_cardio <- vector('list')
physeq_bodysite_other <- vector('list')

for(i in 1:5) {
  sampdat=sample_data(cods_asp[[i]])
  sample_names(sampdat)=cods_asp[[i]]$SampleID
  phy_asp <- subset_samples(physeq_by_bodysites[[i]], CoD_Simple2 == 'Asphyx')
  physeq_bodysite_asp[[i]] <- merge_phyloseq(phy_asp, sampdat)
  
  sampdat=sample_data(cods_bft[[i]])
  sample_names(sampdat)=cods_bft[[i]]$SampleID
  phy_bft <- subset_samples(physeq_by_bodysites[[i]], CoD_Simple2 == 'BFT')
  physeq_bodysite_bft[[i]] <- merge_phyloseq(phy_bft, sampdat)
  
  sampdat=sample_data(cods_gun[[i]])
  sample_names(sampdat)=cods_gun[[i]]$SampleID
  phy_gun <- subset_samples(physeq_by_bodysites[[i]], CoD_Simple2 == 'Gunshot')
  physeq_bodysite_gun[[i]] <- merge_phyloseq(phy_gun, sampdat)
  
  sampdat=sample_data(cods_cardio[[i]])
  sample_names(sampdat)=cods_cardio[[i]]$SampleID
  phy_cardio <- subset_samples(physeq_by_bodysites[[i]], CoD_Simple2 == 'Cardio')
  physeq_bodysite_cardio[[i]] <- merge_phyloseq(phy_cardio, sampdat)
  
  sampdat=sample_data(cods_drug[[i]])
  sample_names(sampdat)=cods_drug[[i]]$SampleID
  phy_drug <- subset_samples(physeq_by_bodysites[[i]], CoD_Simple2 == 'Drug')
  physeq_bodysite_drug[[i]] <- merge_phyloseq(phy_drug, sampdat)
  
  sampdat=sample_data(cods_other[[i]])
  sample_names(sampdat)=cods_other[[i]]$SampleID
  phy_oth <- subset_samples(physeq_by_bodysites[[i]], CoD_Simple2 == 'Other')
  physeq_bodysite_other[[i]] <- merge_phyloseq(phy_oth, sampdat)
  
}



#function for running random forest classification among 
#beta-dispersion values
#output is the classification error matrix and overall classification error
random_foresting_setup<- function(data_phy) {
  otu <- as.data.frame(t(otu_table(data_phy)))
  otu$SampleID <- rownames(otu)
  metadata <- data.frame(sample_data(data_phy))
  meta_sa <- metadata %>% select(SampleID, distances)
  meta_sa$SampleID <- as.character(meta_sa$SampleID)
  otu_f <- merge(meta_sa, otu, by = 'SampleID')
  otu_f <- otu_f[,-1]
  names(otu_f) <- make.names(names(otu_f))
  return(otu_f)
}

rf_mod_acc_mou <- random_foresting_setup(physeq_bodysite_acc[[4]])
rf_mod_nat_mou <- random_foresting_setup(physeq_bodysite_nat[[4]])
rf_mod_sui_mou <- random_foresting_setup(physeq_bodysite_sui[[4]])
rf_mod_hom_mou <- random_foresting_setup(physeq_bodysite_nat[[4]])

#OOB error model
m1 <- randomForest(
  formula = distances ~ .,
  data    = rf_mod_acc_mou,
  ntree= 500
)

#make test/training set
# 80% train,20% test
valid_split <- initial_split(rf_mod_acc_mou, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "distances")]
y_test <- otu_valid$distances

#test/training set model
rf_oob_comp <- randomForest(
  formula = distances~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree = 500
)

oob <- sqrt(m1$mse)
validation <- sqrt(rf_oob_comp$mse)

# compare error rates among models
tibble::tibble(
  `Out of Bag Error` = oob,
  `Test error` = validation,
  ntrees = 1:rf_oob_comp$ntree
) %>%
  gather(Metric, Error, -ntrees) %>%
  ggplot(aes(ntrees, Error, color = Metric)) +
  geom_line() +
  xlab("Number of trees")

#OOB error rate out performed test error for beta-disp with 2,000 trees, reduce to 500 trees.

m1
rf_oob_comp

#extract taxa names to ID indicator taxa
tax <- data.frame(tax_table(physeq_5000))
tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
tax$predictors <- rownames(tax)

#make RF models for each body site and for MOD/COD
#change formula and data
m1 <- randomForest(
  formula = distances ~ .,
  data    = rf_mod_hom_mou,
  ntree= 500
)

m1

#find most important indicator taxa, change based on SD of mean decrease in gini
imp <- importance(m1)
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(IncNodePurity))
imp.sort <- imp.sort %>% filter(IncNodePurity > sd(imp.sort$IncNodePurity))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

#change the name of this to save different taxa
imp_mod_hom_mou <- merge(imp.sort, tax)

#make our new phyloseq objects with the indicator taxa so that we can calculate PCoA components

#put taxonomy table in phyloseq format
#change name depending on body site/ MOD and COD
tax <- subset(imp_mod_hom_mou, select=-c(IncNodePurity))
rownames(tax) <- tax$predictors
tax_reduc <- merge(tax, otu, by = "row.names")
rownames(tax_reduc) <- tax_reduc$Row.names
tax_reduc <- tax_reduc[,-1]
tax_f <- tax_reduc[,1:7]
tax_f <- as.matrix(tax_f)
TAX=tax_table(tax_f)

#change this name depending on which indicator taxa are being used
physeq_notax <- physeq_5000
tax_table(physeq_notax) <- NULL
physeq_imp_mod_hom_mou <- merge_phyloseq(physeq_notax, TAX)
physeq_imp_mod_hom_mou <- subset_samples(physeq_imp_mod_hom_mou, Sample_Area == 'Mouth')
physeq_imp_mod_hom_mou <- subset_samples(physeq_imp_mod_hom_mou, MoD == 'Homicide')
df <- data.frame(otu_table(physeq_imp_mod_hom_mou))
a <- colSums(df)
#look for zeros
View(a)
#prune those samples
control <- prune_samples(sample_names(physeq_imp_mod_hom_mou) != c('WCME_176', 'WCME_69'),
                         physeq_imp_mod_hom_mou)
#or, no samples need to be pruned
control <- physeq_imp_mod_hom_mou

#what about alpha diversity of indicator taxa??
#Shannon: How difficult it is to predict the identity of a randomly chosen individual.
  erich <- estimate_richness(physeq_5000, measures = c("Chao1","Shannon"))
  erich <- add_rownames(erich, "SampleID")
  erich_sums <- merge(erich, metadata)
  alph_richnes[[i]] <- erich_sums
  
#eigen values from this - redo this code 
ord = ordinate(control, method="PCoA", distance="unifrac")
PCoA_values_mod_hom_mou <- data.frame(ord$vectors)
PCoA_values_mod_hom_mou$Sums <- rowSums(PCoA_values_mod_hom_mou)
PCoA_values_mod_hom_mou$SampleID <- rownames(PCoA_values_mod_hom_mou) 
PCoA_values_mod_hom_mou <- PCoA_values_mod_hom_mou[,c('SampleID', 'Sums')]



random_foresting_betadisp <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$distance_cat)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data)
  print(Forest)
  return(Forest)
}


forest_acc <- vector('list')
forest_hom <- vector('list')
forest_nat <- vector('list')
forest_sui <- vector('list')
forest_asp <- vector('list')
forest_bft <- vector('list')
forest_gun <- vector('list')
forest_car <- vector('list')
forest_dru <- vector('list')
forest_oth <- vector('list')

#look at error rate with number of observations

for(i in 1:5) {
  a <-random_foresting_betadisp(physeq_bodysite_acc[[i]])
  forest_acc[[i]] <- a
  b <- random_foresting_betadisp(physeq_bodysite_hom[[i]])
  forest_hom[[i]] <- b
  c <- random_foresting_betadisp(physeq_bodysite_nat[[i]])
  forest_nat[[i]] <- c
  d <- random_foresting_betadisp(physeq_bodysite_sui[[i]])
  forest_sui[[i]] <- d
  e <- random_foresting_betadisp(physeq_bodysite_asp[[i]])
  forest_asp[[i]] <- e
  f <- random_foresting_betadisp(physeq_bodysite_bft[[i]])
  forest_bft[[i]] <- f
  g <- random_foresting_betadisp(physeq_bodysite_gun[[i]])
  forest_gun[[i]] <- g
  h <- random_foresting_betadisp(physeq_bodysite_cardio[[i]])
  forest_car[[i]] <- h
  ii <- random_foresting_betadisp(physeq_bodysite_other[[i]])
  forest_oth[[i]] <- ii
  j <- random_foresting_betadisp(physeq_bodysite_drug[[i]])
  forest_dru[[i]] <- j
}  
  

#function for determing indicator taxa for random forest classification among categorical
#beta-dispersion values
#output is top 10 predictors based on mean decrease in Gini
#function writes predictors to a csv file

for_pred_acc <- vector('list')
for_pred_hom <- vector('list')
for_pred_nat <- vector('list')
for_pred_sui <- vector('list')
for_pred_asp <- vector('list')
for_pred_bft <- vector('list')
for_pred_gun <- vector('list')
for_pred_car <- vector('list')
for_pred_dru <- vector('list')
for_pred_oth <- vector('list')

forest_predictors_betadisp <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.10 <- imp.sort[1:10, ]
  tax <- data.frame(tax_table(physeq_5000))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.10 <- merge(imp.10, tax)
  write.table(imp.10, sep=",", "random_forest_predictors_betadisp.csv", append = TRUE)
  return(imp.10)
}

#resulting output went into table
for(i in 1:5) {
  a <-forest_predictors_betadisp(forest_acc[[i]])
  for_pred_acc[[i]] <- a
  b <- forest_predictors_betadisp(forest_hom[[i]])
  for_pred_hom[[i]] <- b
  c <- forest_predictors_betadisp(forest_nat[[i]])
  for_pred_nat[[i]] <- c
  d <- forest_predictors_betadisp(forest_sui[[i]])
  for_pred_sui[[i]] <- d
  e <- forest_predictors_betadisp(forest_asp[[i]])
  for_pred_asp[[i]] <- e
  f <- forest_predictors_betadisp(forest_bft[[i]])
  for_pred_bft[[i]] <- f
  g <- forest_predictors_betadisp(forest_gun[[i]])
  for_pred_gun[[i]] <- g
  h <- forest_predictors_betadisp(forest_car[[i]])
  for_pred_car[[i]] <- h
  i <- forest_predictors_betadisp(forest_oth[[i]])
  for_pred_oth[[i]] <- i
  j <- forest_predictors_betadisp(forest_dru[[i]])
  for_pred_dru[[i]] <- j
}  

save.image("C:/Users/sierr/Documents/ThesisProject_AKP/R Scripts/AKP_session.RData")


