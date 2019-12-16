####RF: continuous response####

#Here we will perform Random Forest using a continuous phenotype (e.g. RF using regression trees)
library(randomForest)
#for each body site and comparison (mod/cod), going to label beta-dispersion distances categorically
#based on quantiles of data

#change meta_sa to MoD/ CoD_Simple2 depending on what you are looking for
random_foresting_setup<- function(data_phy) {
  otu <- as.data.frame(t(otu_table(data_phy)))
  otu$SampleID <- rownames(otu)
  metadata <- data.frame(sample_data(data_phy))
  meta_sa <- metadata %>% select(SampleID, CoD_Simple2)
  meta_sa$SampleID <- as.character(meta_sa$SampleID)
  otu_f <- merge(meta_sa, otu, by = 'SampleID')
  otu_f <- otu_f[,-1]
  names(otu_f) <- make.names(names(otu_f))
  return(otu_f)
}

rf_mod_nose <- random_foresting_setup(physeq_by_bodysites[[1]])
rf_mod_rec <- random_foresting_setup(physeq_by_bodysites[[2]])
rf_mod_ear <- random_foresting_setup(physeq_by_bodysites[[3]])
rf_mod_mou <- random_foresting_setup(physeq_by_bodysites[[4]])
rf_mod_eye <- random_foresting_setup(physeq_by_bodysites[[5]])

rf_cod_nose <- random_foresting_setup(physeq_by_bodysites[[1]])
rf_cod_rec <- random_foresting_setup(physeq_by_bodysites[[2]])
rf_cod_ear <- random_foresting_setup(physeq_by_bodysites[[3]])
rf_cod_mou <- random_foresting_setup(physeq_by_bodysites[[4]])
rf_cod_eye <- random_foresting_setup(physeq_by_bodysites[[5]])

#choose nose to get parameters
#play with OOB or test/train set error and # of trees

#for MOD: 

#OOB error model
m1 <- randomForest(
  formula = CoD_Simple2 ~ .,
  data    = rf_cod_nose,
  ntree= 2000
)

#make test/training set
# 80% train,20% test
valid_split <- initial_split(rf_cod_nose, .8)
otu_train <- analysis(valid_split)
otu_valid <- assessment(valid_split)
x_test <- otu_valid[setdiff(names(otu_valid), "CoD_Simple2")]
y_test <- otu_valid$CoD_Simple2

#test/training set model
rf_oob_comp <- randomForest(
  formula = CoD_Simple2~ .,
  data    = otu_train,
  xtest   = x_test,
  ytest   = y_test,
  ntree = 2000
)

oob <- sqrt(m1$err.rate)
validation <- sqrt(rf_oob_comp$err.rate)

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

#OOB error rate out performed test error for MOD with 2,000 trees
#OOB error rate out performed test error for COD with 2,000 trees
m1
rf_oob_comp

#extract taxa names to ID indicator taxa
tax <- data.frame(tax_table(physeq_5000))
tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
tax$predictors <- rownames(tax)

#change order of categories to reduce bias
rf_mod_nose$MoD <- factor(rf_mod_nose$MoD, levels = c('Homicide', 'Accident', 'Suicide', 'Natural'))
#didn't change the outcome.

#make RF models for each body site and for MOD/COD
#change formula and data
m1 <- randomForest(
  formula = MoD ~ .,
  data    = rf_mod_nose,
  ntree= 2000
)

m1

#find most important indicator taxa, change based on SD of mean decrease in gini
  imp <- importance(m1)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort <- imp.sort %>% filter(MeanDecreaseGini > sd(imp.sort$MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  
#change the name of this to save different taxa
  imp_mod_ear <- merge(imp.sort, tax)

#make our new phyloseq objects with the indicator taxa so that we can calculate PCoA components

  #put taxonomy table in phyloseq format
  #change name depending on body site/ MOD and COD
  tax <- subset(imp_mod_mou, select=-c(MeanDecreaseGini))
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
  physeq_imp_mod_mou <- merge_phyloseq(physeq_notax, TAX)
  physeq_imp_mod_mou <- subset_samples(physeq_imp_mod_mou, Sample_Area == 'Mouth')
  df <- data.frame(otu_table(physeq_imp_mod_mou))
  a <- colSums(df)
  #look for zeros
  View(a)
  #prune those samples
  control <- prune_samples(sample_names(physeq_imp_mod_eye) != c('WCME_512', 'WCME_513'),
                           physeq_imp_mod_eye)
  #or, no samples need to be pruned
  control <- physeq_imp_mod_mou
  ord = ordinate(control, method="PCoA", distance="unifrac")
  PCoA_values_mod_mou <- data.frame(ord$vectors)
  PCoA_values_mod_mou$Sums <- rowSums(PCoA_values_mod_mou)
  PCoA_values_mod_mou$SampleID <- rownames(PCoA_values_mod_mou) 
  PCoA_values_mod_mou <- PCoA_values_mod_mou[,c('SampleID', 'Sums')]
