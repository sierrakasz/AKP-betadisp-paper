library(Boruta)


#tested .1, .05, and .01 p value. Not really different in the number of atrributes changing
# ranging from 11-13 important attributes. 
boruta.bank_train <- Boruta(MoD~., data = rf_mod_nose, doTrace = 2, pValue = .05)
print(boruta.bank_train)

boruta.bank_train <- Boruta(MoD~., data = rf_mod_mou, doTrace = 2, pValue = .05)
print(boruta.bank_train)

boruta.bank_train <- Boruta(CoD_Simple2~., data = rf_cod_nose, doTrace = 2, pValue = .05)
print(boruta.bank_train)

boruta.bank_train <- Boruta(CoD_Simple2~., data = rf_cod_mou, doTrace = 2, pValue = .05)
print(boruta.bank_train)

boruta.bank_train <- Boruta(CoD_Simple2~., data = rf_cod_ear, doTrace = 2, pValue = .05)
print(boruta.bank_train)

#next, grab out indicator taxa from RF model
#beta-dispersion from those taxa. 

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
rf_cod_ear_bd <- pop_taxa(physeq_by_bodysites[['Ears']], impTaxa)
rf_cod_ear_bd

#unifrac distances
betaphy_cod_ear_rf <- phyloseq::distance(rf_cod_ear_bd, "unifrac")
#beta-dispersion 
#change the number depending on the order listed in samp_area_list
betadisp_cod_ear_rf <- betadisper(betaphy_cod_ear_rf, sample_df_list[[3]]$CoD_Simple2)
# get it correctly formatted
beta_obj <- as.data.frame(betadisp_cod_ear_rf$distances)
beta_obj$SampleID <- rownames(beta_obj)
colnames(beta_obj) <- c('distances', 'SampleID')
beta_df_cod_ear <- merge(beta_obj, metadata, by = 'SampleID') 
