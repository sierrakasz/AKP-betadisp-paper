## core

#change this depending on the MOD/COD
preparing_data_for_core <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  #change MoD/ CoD_Simple2 when applicable
  met <- metadata[,c('SampleID', 'CoD_Simple2')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  #change factors 'Natural', 'Suicide', 'Accident', 'Homicide'
  # or 'Cardio', 'Drug', 'Gunshot', 'Asphyx', 'BFT', 'Other', 'Unknown'
  mtotus$CoD_Simple2 <- factor(mtotus$CoD_Simple2, levels = c('Cardio', 'Drug', 'Gunshot', 'Asphyx', 'BFT', 'Other', 'Unknown'))
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  #change based on MOD/COD
  new_df <- mtotus %>% group_by(CoD_Simple2) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(Upset)
}

all_otus_mod_nose <- preparing_data_for_core(physeq_by_bodysites[['Nose']])
all_otus_mod_mou <- preparing_data_for_core(physeq_by_bodysites[['Mouth']])

#clean up the data
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