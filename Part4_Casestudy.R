### case study


#load in files
#otu table, taxonomy table, and tree file
otu <- read.csv("table.csv")
tax <- read.csv("taxonomy.csv")
tree <- read_tree('tree.nwk')

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



