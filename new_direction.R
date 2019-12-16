#testing with low sample sizes

#test three body areas: mouth, nose, ears

# Gunshot deaths ----------------------------------------------------------
rf_cod_ear$Gun[rf_cod_ear$CoD_Simple2 == 'Gunshot'] <- 'Gun'
rf_cod_ear$Gun[rf_cod_ear$CoD_Simple2 != 'Gunshot'] <- 'Not'
rf_cod_ear$Gun <- as.factor(rf_cod_ear$Gun)

rf_cod_mou$Gun[rf_cod_mou$CoD_Simple2 == 'Gunshot'] <- 'Gun'
rf_cod_mou$Gun[rf_cod_mou$CoD_Simple2 != 'Gunshot'] <- 'Not'
rf_cod_mou$Gun <- as.factor(rf_cod_mou$Gun)

rf_cod_nose$Gun[rf_cod_nose$CoD_Simple2 == 'Gunshot'] <- 'Gun'
rf_cod_nose$Gun[rf_cod_nose$CoD_Simple2 != 'Gunshot'] <- 'Not'
rf_cod_nose$Gun <- as.factor(rf_cod_nose$Gun)

boruta.bank_train <- Boruta(Gun~., data = rf_cod_nose, doTrace = 2, pValue = .05)
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
rf_cod_nos_ind <- pop_taxa(physeq_by_bodysites[['Nose']], impTaxa)
rf_cod_nos_ind_tax <- data.frame(tax_table(rf_cod_nos_ind))

#test permutationally if different
#change beta variable based on comparison
beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Gun)
  print(return(permutest(beta)))
}

sample_data(physeq_by_bodysites[['Ears']])$Gun[sample_data(physeq_by_bodysites[['Ears']])$CoD_Simple2 == 'Gunshot'] <- 'Gun'
sample_data(physeq_by_bodysites[['Ears']])$Gun[sample_data(physeq_by_bodysites[['Ears']])$CoD_Simple2 != 'Gunshot'] <- 'Not'
sample_data(physeq_by_bodysites[['Ears']])$Gun <- as.factor(sample_data(physeq_by_bodysites[['Ears']])$Gun)

sample_data(physeq_by_bodysites[['Mouth']])$Gun[sample_data(physeq_by_bodysites[['Mouth']])$CoD_Simple2 == 'Gunshot'] <- 'Gun'
sample_data(physeq_by_bodysites[['Mouth']])$Gun[sample_data(physeq_by_bodysites[['Mouth']])$CoD_Simple2 != 'Gunshot'] <- 'Not'
sample_data(physeq_by_bodysites[['Mouth']])$Gun <- as.factor(sample_data(physeq_by_bodysites[['Mouth']])$Gun)

sample_data(physeq_by_bodysites[['Nose']])$Gun[sample_data(physeq_by_bodysites[['Nose']])$CoD_Simple2 == 'Gunshot'] <- 'Gun'
sample_data(physeq_by_bodysites[['Nose']])$Gun[sample_data(physeq_by_bodysites[['Nose']])$CoD_Simple2 != 'Gunshot'] <- 'Not'
sample_data(physeq_by_bodysites[['Nose']])$Gun <- as.factor(sample_data(physeq_by_bodysites[['Nose']])$Gun)


beta_dispersion_calc(physeq_by_bodysites[['Ears']])
beta_dispersion_calc(physeq_by_bodysites[['Mouth']])
beta_dispersion_calc(physeq_by_bodysites[['Nose']])


# Hom vs. Sui gunshot deaths ----------------------------------------------------------
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
rf_mod_ear <- random_foresting_setup(physeq_by_bodysites[[3]])
rf_mod_mou <- random_foresting_setup(physeq_by_bodysites[[4]])

rf_mod_ear_gun <- rf_mod_ear %>%  filter(CoD_Simple2 == 'Gunshot')
rf_mod_mou_gun <- rf_mod_mou %>%  filter(CoD_Simple2 == 'Gunshot')
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

phy_ear_gun <- subset_samples(physeq_by_bodysites[['Ears']], CoD_Simple2 == 'Gunshot')
phy_mou_gun <- subset_samples(physeq_by_bodysites[['Mouth']], CoD_Simple2 == 'Gunshot')
phy_nos_gun <- subset_samples(physeq_by_bodysites[['Nose']], CoD_Simple2 == 'Gunshot')

beta_dispersion_calc(phy_ear_gun)
beta_dispersion_calc(phy_mou_gun)
beta_dispersion_calc(phy_nos_gun)

# Blunt force trauma ----------------------------------------------------------
rf_cod_ear$BFT[rf_cod_ear$CoD_Simple2 == 'BFT'] <- 'BFT'
rf_cod_ear$BFT[rf_cod_ear$CoD_Simple2 != 'BFT'] <- 'Not'
rf_cod_ear$BFT <- as.factor(rf_cod_ear$BFT)

rf_cod_mou$BFT[rf_cod_mou$CoD_Simple2 == 'BFT'] <- 'BFT'
rf_cod_mou$BFT[rf_cod_mou$CoD_Simple2 != 'BFT'] <- 'Not'
rf_cod_mou$BFT <- as.factor(rf_cod_mou$BFT)

rf_cod_nose$BFT[rf_cod_nose$CoD_Simple2 == 'BFT'] <- 'BFT'
rf_cod_nose$BFT[rf_cod_nose$CoD_Simple2 != 'BFT'] <- 'Not'
rf_cod_nose$BFT <- as.factor(rf_cod_nose$BFT)

boruta.bank_train <- Boruta(BFT~., data = rf_cod_nose, doTrace = 2, pValue = .05)
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
rf_cod_nos_ind <- pop_taxa(physeq_by_bodysites[['Nose']], impTaxa)
rf_cod_nos_ind_tax <- data.frame(tax_table(rf_cod_nos_ind))

#test permutationally if different
#change beta variable based on comparison
beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$BFT)
  print(return(permutest(beta)))
}

sample_data(physeq_by_bodysites[['Ears']])$BFT[sample_data(physeq_by_bodysites[['Ears']])$CoD_Simple2 == 'BFT'] <- 'BFT'
sample_data(physeq_by_bodysites[['Ears']])$BFT[sample_data(physeq_by_bodysites[['Ears']])$CoD_Simple2 != 'BFT'] <- 'Not'
sample_data(physeq_by_bodysites[['Ears']])$BFT <- as.factor(sample_data(physeq_by_bodysites[['Ears']])$BFT)

sample_data(physeq_by_bodysites[['Mouth']])$BFT[sample_data(physeq_by_bodysites[['Mouth']])$CoD_Simple2 == 'BFT'] <- 'BFT'
sample_data(physeq_by_bodysites[['Mouth']])$BFT[sample_data(physeq_by_bodysites[['Mouth']])$CoD_Simple2 != 'BFT'] <- 'Not'
sample_data(physeq_by_bodysites[['Mouth']])$BFT <- as.factor(sample_data(physeq_by_bodysites[['Mouth']])$BFT)

sample_data(physeq_by_bodysites[['Nose']])$BFT[sample_data(physeq_by_bodysites[['Nose']])$CoD_Simple2 == 'BFT'] <- 'BFT'
sample_data(physeq_by_bodysites[['Nose']])$BFT[sample_data(physeq_by_bodysites[['Nose']])$CoD_Simple2 != 'BFT'] <- 'Not'
sample_data(physeq_by_bodysites[['Nose']])$BFT <- as.factor(sample_data(physeq_by_bodysites[['Nose']])$BFT)


beta_dispersion_calc(physeq_by_bodysites[['Ears']])
beta_dispersion_calc(physeq_by_bodysites[['Mouth']])
beta_dispersion_calc(physeq_by_bodysites[['Nose']])

# Acc vs. Hom BFT
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
rf_mod_ear <- random_foresting_setup(physeq_by_bodysites[[3]])
rf_mod_mou <- random_foresting_setup(physeq_by_bodysites[[4]])

rf_mod_ear_BFT <- rf_mod_ear %>%  filter(CoD_Simple2 == 'BFT')
rf_mod_mou_BFT <- rf_mod_mou %>%  filter(CoD_Simple2 == 'BFT')
rf_mod_nos_BFT <- rf_mod_nose %>%  filter(CoD_Simple2 == 'BFT')

boruta.bank_train <- Boruta(MoD~., data = rf_mod_nos_BFT, doTrace = 2, pValue = .05)
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

phy_ear_BFT <- subset_samples(physeq_by_bodysites[['Ears']], CoD_Simple2 == 'BFT')
phy_mou_BFT <- subset_samples(physeq_by_bodysites[['Mouth']], CoD_Simple2 == 'BFT')
phy_nos_BFT <- subset_samples(physeq_by_bodysites[['Nose']], CoD_Simple2 == 'BFT')

beta_dispersion_calc(phy_ear_BFT)
beta_dispersion_calc(phy_mou_BFT)
beta_dispersion_calc(phy_nos_BFT)

# Drug use ----------------------------------------------------------
rf_cod_ear$Drug[rf_cod_ear$CoD_Simple2 == 'Drug'] <- 'Drug'
rf_cod_ear$Drug[rf_cod_ear$CoD_Simple2 != 'Drug'] <- 'Not'
rf_cod_ear$Drug <- as.factor(rf_cod_ear$Drug)

rf_cod_mou$Drug[rf_cod_mou$CoD_Simple2 == 'Drug'] <- 'Drug'
rf_cod_mou$Drug[rf_cod_mou$CoD_Simple2 != 'Drug'] <- 'Not'
rf_cod_mou$Drug <- as.factor(rf_cod_mou$Drug)

rf_cod_nose$Drug[rf_cod_nose$CoD_Simple2 == 'Drug'] <- 'Drug'
rf_cod_nose$Drug[rf_cod_nose$CoD_Simple2 != 'Drug'] <- 'Not'
rf_cod_nose$Drug <- as.factor(rf_cod_nose$Drug)

boruta.bank_train <- Boruta(Drug~., data = rf_cod_nose, doTrace = 2, pValue = .05)
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
rf_cod_nos_ind <- pop_taxa(physeq_by_bodysites[['Nose']], impTaxa)
rf_cod_nos_ind_tax <- data.frame(tax_table(rf_cod_nos_ind))

#test permutationally if different
#change beta variable based on comparison
beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "unifrac")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Drug)
  print(return(permutest(beta)))
}

sample_data(physeq_by_bodysites[['Ears']])$Drug[sample_data(physeq_by_bodysites[['Ears']])$CoD_Simple2 == 'Drug'] <- 'Drug'
sample_data(physeq_by_bodysites[['Ears']])$Drug[sample_data(physeq_by_bodysites[['Ears']])$CoD_Simple2 != 'Drug'] <- 'Not'
sample_data(physeq_by_bodysites[['Ears']])$Drug <- as.factor(sample_data(physeq_by_bodysites[['Ears']])$Drug)

sample_data(physeq_by_bodysites[['Mouth']])$Drug[sample_data(physeq_by_bodysites[['Mouth']])$CoD_Simple2 == 'Drug'] <- 'Drug'
sample_data(physeq_by_bodysites[['Mouth']])$Drug[sample_data(physeq_by_bodysites[['Mouth']])$CoD_Simple2 != 'Drug'] <- 'Not'
sample_data(physeq_by_bodysites[['Mouth']])$Drug <- as.factor(sample_data(physeq_by_bodysites[['Mouth']])$Drug)

sample_data(physeq_by_bodysites[['Nose']])$Drug[sample_data(physeq_by_bodysites[['Nose']])$CoD_Simple2 == 'Drug'] <- 'Drug'
sample_data(physeq_by_bodysites[['Nose']])$Drug[sample_data(physeq_by_bodysites[['Nose']])$CoD_Simple2 != 'Drug'] <- 'Not'
sample_data(physeq_by_bodysites[['Nose']])$Drug <- as.factor(sample_data(physeq_by_bodysites[['Nose']])$Drug)


beta_dispersion_calc(physeq_by_bodysites[['Ears']])
beta_dispersion_calc(physeq_by_bodysites[['Mouth']])
beta_dispersion_calc(physeq_by_bodysites[['Nose']])

# Acc vs. Hom BFT
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
rf_mod_ear <- random_foresting_setup(physeq_by_bodysites[[3]])
rf_mod_mou <- random_foresting_setup(physeq_by_bodysites[[4]])

rf_mod_ear_Drug <- rf_mod_ear %>%  filter(CoD_Simple2 == 'Drug')
rf_mod_mou_Drug <- rf_mod_mou %>%  filter(CoD_Simple2 == 'Drug')
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

phy_ear_Drug <- subset_samples(physeq_by_bodysites[['Ears']], CoD_Simple2 == 'Drug')
phy_mou_Drug <- subset_samples(physeq_by_bodysites[['Mouth']], CoD_Simple2 == 'Drug')
phy_nos_Drug <- subset_samples(physeq_by_bodysites[['Nose']], CoD_Simple2 == 'Drug')

beta_dispersion_calc(phy_ear_Drug)
beta_dispersion_calc(phy_mou_Drug)
beta_dispersion_calc(phy_nos_Drug)

# Suicide deaths ----------------------------------------------------------
rf_mod_ear$Suicide[rf_mod_ear$MoD == 'Suicide'] <- 'Suicide'
rf_mod_ear$Suicide[rf_mod_ear$MoD != 'Suicide'] <- 'Not'
rf_mod_ear$Suicide <- as.factor(rf_mod_ear$Suicide)

rf_mod_mou$Suicide[rf_mod_mou$MoD == 'Suicide'] <- 'Suicide'
rf_mod_mou$Suicide[rf_mod_mou$MoD != 'Suicide'] <- 'Not'
rf_mod_mou$Suicide <- as.factor(rf_mod_mou$Suicide)

rf_mod_nose$Suicide[rf_mod_nose$MoD == 'Suicide'] <- 'Suicide'
rf_mod_nose$Suicide[rf_mod_nose$MoD != 'Suicide'] <- 'Not'
rf_mod_nose$Suicide <- as.factor(rf_mod_nose$Suicide)

boruta.bank_train <- Boruta(Suicide~., data = rf_mod_nose, doTrace = 2, pValue = .05)
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
  beta <- betadisper(GPdist, sampledf$Suicide)
  print(return(permutest(beta)))
}

sample_data(physeq_by_bodysites[['Ears']])$Suicide[sample_data(physeq_by_bodysites[['Ears']])$MoD == 'Suicide'] <- 'Suicide'
sample_data(physeq_by_bodysites[['Ears']])$Suicide[sample_data(physeq_by_bodysites[['Ears']])$MoD != 'Suicide'] <- 'Not'
sample_data(physeq_by_bodysites[['Ears']])$Suicide <- as.character(sample_data(physeq_by_bodysites[['Ears']])$Suicide)
sample_data(physeq_by_bodysites[['Ears']])$Suicide[is.na(sample_data(physeq_by_bodysites[['Ears']])$Suicide)] <- 'Suicide'

sample_data(physeq_by_bodysites[['Mouth']])$Suicide[sample_data(physeq_by_bodysites[['Mouth']])$MoD == 'Suicide'] <- 'Suicide'
sample_data(physeq_by_bodysites[['Mouth']])$Suicide[sample_data(physeq_by_bodysites[['Mouth']])$MoD != 'Suicide'] <- 'Not'
sample_data(physeq_by_bodysites[['Mouth']])$Suicide <- as.character(sample_data(physeq_by_bodysites[['Mouth']])$Suicide)
sample_data(physeq_by_bodysites[['Mouth']])$Suicide[is.na(sample_data(physeq_by_bodysites[['Mouth']])$Suicide)] <- 'Suicide'

sample_data(physeq_by_bodysites[['Nose']])$Suicide[sample_data(physeq_by_bodysites[['Nose']])$MoD == 'Suicide'] <- 'Suicide'
sample_data(physeq_by_bodysites[['Nose']])$Suicide[sample_data(physeq_by_bodysites[['Nose']])$MoD != 'Suicide'] <- 'Not'
sample_data(physeq_by_bodysites[['Nose']])$Suicide <- as.character(sample_data(physeq_by_bodysites[['Nose']])$Suicide)
sample_data(physeq_by_bodysites[['Nose']])$Suicide[is.na(sample_data(physeq_by_bodysites[['Nose']])$Suicide)] <- 'Suicide'

beta_dispersion_calc(physeq_by_bodysites[['Ears']])
beta_dispersion_calc(physeq_by_bodysites[['Mouth']])
beta_dispersion_calc(physeq_by_bodysites[['Nose']])

