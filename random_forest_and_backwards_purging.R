#Overall workflow:
#		1. prepare input data - convert to matrix format, add any covariates, linear regression (Zhao et al. 2012) correction of genotypes and/or phenotypes
#		2. optimize parameters for mtry and ntree
#		3. perform rf on subsets of loci to id set for backwards purging
#		4. backwards purging using set of conservative subset of loci





####RF: categorical response####

#Here we will perform Random Forest using a categorical phenotype (e.g. RF using classification trees)

#Set up environment
library(randomForest)
library(vcfR)

#randomForest encodes genotypes in the following way:
#   0 = homozygote
#   1 = heterozygote
#   2 = homozygote 2

#the input file format is a data matrix with one row per individual, and one column per locus. genotypes are therefore
#represented by single digits corresponding with the above values.
#genlight objects (adegenet) are encoded in the following manner (single digit genotypes reflecting the number of 
#copies of the second allele):
#   0 = homozygote major allele
#   1 = heterozygote
#   2 = homozygote minor allele
#   NA = missing data

#so this means we can simply transform our genlight object into a data matrix, and use this for random forest. 

#however, RF does not tolerate missing data, so before importing data, impute missing values using a program such as Beagle (or FastPhase, LinkImpute...)
#if use beagle, first sort markers with bedtools, then add vcf header lines back in (bedtools removes them). 

#import imputed dataset (in vcf file format) and convert to matrix format
rf.vcf = vcfR::read.vcfR("input.impute.vcf", convertNA = TRUE, verbose = TRUE)
rf.gl = vcfR::vcfR2genlight(rf.vcf)
rf.mat = as.matrix(rf.gl)

#save genotypes matrix
write.csv(rf.mat, "input.impute.matrix.csv")
#now manipulate to add columns for any covariates we want to correct for (col 2:n), and our phenotype of interest (col n+1)

#my input data typically look like this:
#	sample  Pop	 LifeHist  loc1  loc2  loc3
#	sample1	pop1 coaster   0     1	   1
#	sample2	pop1 coaster   2     0	   1

#where Pop is a covariate that represents population structure for the dataset (typically either a categorical population assignment or a continuous value from principal component axis (or axes))
#and LifeHist is (in this case) a categorical respose variable for which we are identifying predictor loci

#import data
rf.df.1 = read.csv("input.impute.matrix.pheno.csv", header = TRUE)
rf.df.1[1:10,1:5]

#make sample names (col 1) as row names
rf.df = rf.df.1[,-1]
rf.df[1:10,1:5]
rownames(rf.df) = rf.df.1[,1]
rf.df[1:10,1:5]

# If implementing correction for population stratification (and/or other covariates): 
# Here we will correct the genotypes using the approach of Zhao et al.(2012). 
# We will not correct the phenotype itself, as it is binary and we want to maintain its categorical distribution.
rf.df.corrected = rf.df #store as new object so still have original on hand
ncol(rf.df.corrected) #id last col of data frame
rf.df.corrected[1:10,1:5]
rf.df.corrected[,3:9884] <- NA  # Keep columns 1-2 with population ID and phenotype, but then replace genotypes with NAs, we will then overwrite this with residuals for covariate correction (Zhao et al. 2012)
rf.df.corrected[1:10,1:5]

# Now correct the genotypes using the regression/residual method (Zhao et al. 2012). We're using a standard linear regression because Zhao et al. 2012 found that the correction procedure is robust to selection of the link function
for (i in 3:ncol(rf.df)){ #start with col 3 because that is the first col of genotypes
  #LM_SNP_i <- lm(rf.df[,i] ~ factor(rf.df$Pop)) # apply linear model to all loci and the response; categorical covariate
  LM_SNP_i <- lm(rf.df[,i] ~ rf.df$Pop) # apply linear model to all loci and the response; continuous covariate
  rf.df.corrected[,i] <- LM_SNP_i$residuals
  colnames(rf.df.corrected)[i]<-colnames(rf.df)[i] 
  if(i%%50==0) print(i)
}


# Verify that the residuals have been written to the data frame properly, using the last column as an example
rf.df.corrected[,9884]-LM_SNP_i$residuals  #Should all be zero if correct (open file and check entire thing anyway)
rf.df.corrected[1:10,1:5]

# Export a copy of the corrected data for future reference
write.csv(rf.df.corrected,file="infile.impute.matrix.pheno.corrected.csv")

# Left this code in in case you need to worry about unbalanced sample sizes across phenotypes:
# Before running Random Forest, let's also check for an imbalance in the response variable 
# because any imbalances can bias the results.
length(which(rf.df$LifeHist==1)) 
length(which(rf.df$LifeHist==2))

# the phenotypes are pretty evenly distributed, so we shouldn't need to worry about a correction
# BUT if we needed to, we would over-sample the underrepresented class
# and under-sample the overrepresented class. This is performed by setting the sampsize parameter to 2/3
# of the class with the lower sample size and using those sample sizes in conjunction with the strata option
# 2/3 * 59 = 39 -> sub-sample each phenotype category to 35 observations
x = (2/3) * 30; x
sample_size <- c(20,20)
# You would then  use these sample sizes in conjunction with the strata option for RF analyses below




####################################
#########Ready for Random Forest!!!!
####################################

# Now run Random Forest analysis. Since this is a binary trait, we need to conduct a classification RF

# First, we need to optimize mtry by running different values of mtry at different values of ntree. 

# We will run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p, where p is the number of loci
# We will initially run each of these mtry values at ntree=100 to 1000 (by increments of 100). 
# We are looking for a plateau where the out-of-bag error rate (OOB-ER) stops decreasing with larger values of ntree
# Once we reach the plateau, we will choose the mtry value that minimizes the OOB-ER.

sqrt(9882) #change to fit your number of loci
2*sqrt(9882)
0.1*(9882)
0.2*(9882)
(9882)/3

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(99,199,988,1976,3294,9882)){    #values of mtry based on total # loci
    rf_ij <- randomForest(x = rf.df.corrected[,3:9884], y = as.factor(rf.df.corrected$LifeHist), importance=TRUE, 
                          proximity=TRUE, ntree=i, mtry=j) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
head(results_optimization)
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")
head(results_optimization)

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 99],results_optimization$OOB_ER[results_optimization$mtry == 99], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1.0), lwd = 2)
lines(results_optimization$ntree[results_optimization$mtry == 199],results_optimization$OOB_ER[results_optimization$mtry == 199], col="blue", lwd = 2)
lines(results_optimization$ntree[results_optimization$mtry == 988],results_optimization$OOB_ER[results_optimization$mtry == 988], col="green", lwd = 2)
lines(results_optimization$ntree[results_optimization$mtry == 1976],results_optimization$OOB_ER[results_optimization$mtry == 1976], col="purple", lwd = 2)
lines(results_optimization$ntree[results_optimization$mtry == 3294],results_optimization$OOB_ER[results_optimization$mtry == 3294], col="orange", lwd = 2)
lines(results_optimization$ntree[results_optimization$mtry == 9882],results_optimization$OOB_ER[results_optimization$mtry == 9882], col="red", lwd = 2)
legend(200, 0.8, c("mtry99","mtry199","mtry988","mtry1976","mtry3294","mtry9882"),
       pt.bg=c("black","blue","green","purple","orange","red"),
       pch=c(21,21,21,21,21,21), bty='n', x.intersp=.5, y.intersp = .6, text.font=1, 
       pt.cex=2, cex = 1)


# remember, want the value for mtry associated with the LOWEST OOB-ER!!!!


# Now begin the full Random Forest analyses

# Recall that even though we optimized mtry, we must now run a larger number of trees in order to achieve convergence of importance values between forests.
# As a starting point, you could try growing 25,000 trees and then increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# since the OOB-ER reached a plateau for a given mtry value after about 400 trees.
rf_all_1 = randomForest(x = rf.df.corrected[,3:9884], y = as.factor(rf.df.corrected$LifeHist), importance=TRUE, proximity=TRUE, mtry=9882, ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size) #turn on strata options if using correction for unbalanced sample sizes
save(rf_all_1,file="rf_all_1.Rdata")

rf_all_2 = randomForest(x = rf.df.corrected[,3:9884], y = as.factor(rf.df.corrected$LifeHist), importance=TRUE, proximity=TRUE, mtry=9882, ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_all_2,file="rf_all_2.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")
head(importance_rf_all_1)
nrow(importance_rf_all_1)
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")
head(importance_rf_all_2)
nrow(importance_rf_all_2)

cor(importance_rf_all_1,importance_rf_all_2) #I typically look for a correlation of >= 0.95 between forests. Increase ntree to get near that target.


rf_all_3 = randomForest(x = rf.df.corrected[,3:9884], y = as.factor(rf.df.corrected$LifeHist), importance=TRUE, proximity=TRUE, mtry=9882, ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_all_3,file="rf_all_3.Rdata")

importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")
head(importance_rf_all_3)

# The predictive ability of classification trees is measured by the out-of-bag error rate.  An error rate is calculated for each tree within a forest. 
# We will use the error rate from the last tree in the forest, which takes all previous trees into account and thus represents the error rate after the model stabilizes/converges

rf_all_1_err.rate <- rf_all_1$err.rate[150000]; rf_all_1_err.rate 
rf_all_2_err.rate <- rf_all_2$err.rate[150000]; rf_all_2_err.rate
rf_all_3_err.rate <- rf_all_3$err.rate[150000]; rf_all_3_err.rate

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")
head(importance_rf_all)

# Export importance values for future reference
write.csv(importance_rf_all, file="rf_importance_values_mtry9882_ntree150000_popstratcorr.csv", row.names=FALSE)



# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of life history phenotype. 
# For each subset, we will use the mtry identified above since that is the optimal setting that we previously found.

##### Best 2% 
names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes 
genotypes_2perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_2perc_unique]

# Now conduct RF on this subset
rf_2perc_1 = randomForest(x = genotypes_2perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_2perc_1,file="rf_2perc_1.Rdata")

rf_2perc_2 = randomForest(x = genotypes_2perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_2perc_2,file="rf_2perc_2.Rdata")

rf_2perc_3 = randomForest(x = genotypes_2perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_2perc_3,file="rf_2perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_2perc_1_err.rate <- rf_2perc_1$err.rate[150000]; rf_2perc_1_err.rate 
rf_2perc_2_err.rate <- rf_2perc_2$err.rate[150000]; rf_2perc_2_err.rate 
rf_2perc_3_err.rate <- rf_2perc_3$err.rate[150000]; rf_2perc_3_err.rate 

rm(rf_2perc_1,rf_2perc_2,rf_2perc_3) # remove the objects to save memory in R 


##### Best 3% 

names_best_3perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.97))]
names_best_3perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.97))]
names_best_3perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.97))]
names_best_3perc_unique<-unique(c(names_best_3perc_1,names_best_3perc_2,names_best_3perc_3))

# Extract genotypes
genotypes_3perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_3perc_unique]

# Now conduct RF on this subset
rf_3perc_1 = randomForest(x = genotypes_3perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_3perc_1,file="rf_3perc_1.Rdata")

rf_3perc_2 = randomForest(x = genotypes_3perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_3perc_2,file="rf_3perc_2.Rdata")

rf_3perc_3 = randomForest(x = genotypes_3perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_3perc_3,file="rf_3perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_3perc_1_err.rate <- rf_3perc_1$err.rate[150000]; rf_3perc_1_err.rate
rf_3perc_2_err.rate <- rf_3perc_2$err.rate[150000]; rf_3perc_2_err.rate
rf_3perc_3_err.rate <- rf_3perc_3$err.rate[150000]; rf_3perc_3_err.rate

rm(rf_3perc_1,rf_3perc_2,rf_3perc_3) # remove the objects to save memory in R 

##### Best 4% 

names_best_4perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.96))]
names_best_4perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.96))]
names_best_4perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.96))]
names_best_4perc_unique<-unique(c(names_best_4perc_1,names_best_4perc_2,names_best_4perc_3))

# Extract genotypes
genotypes_4perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_4perc_unique]

# Now conduct RF on this subset
rf_4perc_1 = randomForest(x = genotypes_4perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_4perc_1,file="rf_4perc_1.Rdata")

rf_4perc_2 = randomForest(x = genotypes_4perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_4perc_2,file="rf_4perc_2.Rdata")

rf_4perc_3 = randomForest(x = genotypes_4perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_4perc_3,file="rf_4perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_4perc_1_err.rate <- rf_4perc_1$err.rate[150000]; rf_4perc_1_err.rate
rf_4perc_2_err.rate <- rf_4perc_2$err.rate[150000]; rf_4perc_2_err.rate
rf_4perc_3_err.rate <- rf_4perc_3$err.rate[150000]; rf_4perc_3_err.rate

rm(rf_4perc_1,rf_4perc_2,rf_4perc_3) # remove the objects to save memory in R 

##### Best 5% 

names_best_5perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.95))]
names_best_5perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.95))]
names_best_5perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.95))]
names_best_5perc_unique<-unique(c(names_best_5perc_1,names_best_5perc_2,names_best_5perc_3))

# Extract genotypes
genotypes_5perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_5perc_unique]

# Now conduct RF on this subset
rf_5perc_1 = randomForest(x = genotypes_5perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_5perc_1,file="rf_5perc_1.Rdata")

rf_5perc_2 = randomForest(x = genotypes_5perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_5perc_2,file="rf_5perc_2.Rdata")

rf_5perc_3 = randomForest(x = genotypes_5perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_5perc_3,file="rf_5perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_5perc_1_err.rate <- rf_5perc_1$err.rate[150000]
rf_5perc_2_err.rate <- rf_5perc_2$err.rate[150000]
rf_5perc_3_err.rate <- rf_5perc_3$err.rate[150000]

rm(rf_5perc_1,rf_5perc_2,rf_5perc_3) # remove the objects to save memory in R 

##### Best 10% 

names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))


# Extract genotypes
genotypes_10perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_10perc_unique]

# Now conduct RF on this subset
rf_10perc_1 = randomForest(x = genotypes_10perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_10perc_1,file="rf_10perc_1.Rdata")

rf_10perc_2 = randomForest(x = genotypes_10perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_10perc_2,file="rf_10perc_2.Rdata")

rf_10perc_3 = randomForest(x = genotypes_10perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_10perc_3,file="rf_10perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_10perc_1_err.rate <- rf_10perc_1$err.rate[150000]
rf_10perc_2_err.rate <- rf_10perc_2$err.rate[150000]
rf_10perc_3_err.rate <- rf_10perc_3$err.rate[150000]

rm(rf_10perc_1,rf_10perc_2,rf_10perc_3) # remove the objects to save memory in R 

##### Best 20% 

names_best_20perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.80))]
names_best_20perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.80))]
names_best_20perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.80))]
names_best_20perc_unique<-unique(c(names_best_20perc_1,names_best_20perc_2,names_best_20perc_3))

# Extract genotypes
genotypes_20perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_20perc_unique]

# Now conduct RF on this subset
rf_20perc_1 = randomForest(x = genotypes_20perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_20perc_1,file="rf_20perc_1.Rdata")

rf_20perc_2 = randomForest(x = genotypes_20perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_20perc_2,file="rf_20perc_2.Rdata")

rf_20perc_3 = randomForest(x = genotypes_20perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_20perc_3,file="rf_20perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_20perc_1_err.rate <- rf_20perc_1$err.rate[150000]
rf_20perc_2_err.rate <- rf_20perc_2$err.rate[150000]
rf_20perc_3_err.rate <- rf_20perc_3$err.rate[150000]

rm(rf_20perc_1,rf_20perc_2,rf_20perc_3) # remove the objects to save memory in R 

##### Best 30% 

names_best_30perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.70))]
names_best_30perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.70))]
names_best_30perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.70))]
names_best_30perc_unique<-unique(c(names_best_30perc_1,names_best_30perc_2,names_best_30perc_3))

# Extract genotypes
genotypes_30perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_30perc_unique]

# Now conduct RF on this subset
rf_30perc_1 = randomForest(x = genotypes_30perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_30perc_1,file="rf_30perc_1.Rdata")

rf_30perc_2 = randomForest(x = genotypes_30perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_30perc_2,file="rf_30perc_2.Rdata")

rf_30perc_3 = randomForest(x = genotypes_30perc, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_30perc_3,file="rf_30perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_30perc_1_err.rate <- rf_30perc_1$err.rate[150000]
rf_30perc_2_err.rate <- rf_30perc_2$err.rate[150000]
rf_30perc_3_err.rate <- rf_30perc_3$err.rate[150000]

rm(rf_30perc_1,rf_30perc_2,rf_30perc_3) # remove the objects to save memory in R 



# Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach

All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
                              cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
                              cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
                              cbind(rf_4perc_1_err.rate,rf_4perc_2_err.rate,rf_4perc_3_err.rate),
                              cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate),
                              cbind(rf_10perc_1_err.rate,rf_10perc_2_err.rate,rf_10perc_3_err.rate),
                              cbind(rf_20perc_1_err.rate,rf_20perc_2_err.rate,rf_20perc_3_err.rate),
                              cbind(rf_30perc_1_err.rate,rf_30perc_2_err.rate,rf_30perc_3_err.rate))

# Plot error rates for the various subsets
All_initial_err.rate<-data.frame(All_initial_err.rate)
All_initial_err.rate$Number_loci<-c(9882,length(names_best_2perc_unique),length(names_best_3perc_unique),length(names_best_4perc_unique),length(names_best_5perc_unique),length(names_best_10perc_unique),length(names_best_20perc_unique),length(names_best_30perc_unique))
rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_err.rate, file="All_initial_err_rate_mtryallloci_ntree150000_popstratcorr.csv")

# Plot error rates as well
pdf("All_initial_err_rate_mtryallloci_ntree150000_popstratcorr.pdf", width=10, height=8)
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", 
     ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1, ylim=c(0.08,0.28), yaxt="n")
axis(side=2, at=seq(0.08, 0.28, by=0.02))
dev.off()

#use this information to select a subset of loci for backwards purging
#save names of loci that will be used for backwards purging
write.csv(names_best_4perc_unique, "names_best_4perc_unique.csv")

#select set of loci with lowest OOB-ER (be conservative) to continue with backwards purging
#now that we know which subset of loci seems most relevant to our phenotype, let's perform 
#backwards purging to identify optimized set of explanatory loci


#################### Backward purging approach


#this part takes a really long time so will likely have to run on hpcc 

#Import matrix of genotypes corrected for population stratification
data1 = read.csv("input.impute.matrix.pheno.corrected.csv", header = TRUE)
rf.df.corrected = data1[,-1]
rownames(rf.df.corrected) = data1[,1]


#Import locus names for desired set of loci
data2 = read.csv("names_best_4perc_unique.csv", header = TRUE)
names_best_4perc_unique = as.vector(data2[['x']])


#Backward purging approach
names_purging <- names_best_4perc_unique
genotypes_purging<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_purging]


rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_purging_1,file="rf_purging_4pct_1.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_purging_2,file="rf_purging_4pct_2.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
save(rf_purging_3,file="rf_purging_4pct_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[length(names_purging)]]<-names_purging
error_rate_best<-data.frame(V1=1:length(names_purging),V2=1:length(names_purging),V3=1:length(names_purging))
rownames(error_rate_best)<-1:length(names_purging)
error_rate_best[length(names_purging),] <- c(rf_purging_1$err.rate[150000],rf_purging_2$err.rate[150000],rf_purging_3$err.rate[150000])


for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1<-data.frame(importance(rf_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging)
  colnames(imp_purging_1)<-"Mean_Decrease_Accuracy1"
  rownames(imp_purging_2)<-colnames(genotypes_purging)
  colnames(imp_purging_2)<-"Mean_Decrease_Accuracy2"
  rownames(imp_purging_3)<-colnames(genotypes_purging)
  colnames(imp_purging_3)<-"Mean_Decrease_Accuracy3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]]<-names_keep
  genotypes_purging<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_keep]
  rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
  rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
  rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(rf.df.corrected$LifeHist), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=150000) #, strata=as.factor(rf.df.corrected$LifeHist), sampsize=sample_size)
  error_rate_best[length(names_purging)-i,] <- c(rf_purging_1$err.rate[150000],rf_purging_2$err.rate[150000],rf_purging_3$err.rate[150000])
}


error_rate_best$Average<-apply(error_rate_best,1,mean)
write.csv(error_rate_best, file="Backward_purging_OOB-ER_ntree150000_popstratcorr_4pct.csv") # Save the error rates

# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus
plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5,pch=16)

# Which group of loci yields the lowest error rate?
lowest.ER = which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) 
lowest.ER
write.csv(lowest.ER, file="Lowest_ER_loci_ntree150000_popstratcorr_4pct.csv")

# Export the names of the predictor loci
write.csv(names_all_iterations[[lowest.ER]],file="Predictor_loci_ntree150000_popstratcorr_purgewall_4pct.csv")





####RF: continuous response####

#Here we will perform Random Forest using a continuous phenotype (e.g. RF using regression trees)
library(randomForest)
library(vcfR)

#see above notes under categorical RF for information on formatting input data

#input vcf file, convert to adgenet genlight object then matrix
rf.vcf = vcfR::read.vcfR("input.impute.vcf", convertNA = TRUE, verbose = TRUE)
rf.gl = vcfR::vcfR2genlight(rf.vcf)
rf.mat = as.matrix(rf.gl)
write.csv(rf.mat, "input.impute.matrix.csv")
#now manipulate to add columns for any covariates we want to correct for (col 2:n), and our phenotype of interest (col n+1)
#again see notes for categorical response variable, above

#import data with covariates included
rf.df.1 = read.csv("input.impute.matrix.del13C.csv", header = TRUE)

#make indiv names as row names
rf.df = rf.df.1[,-1]
rf.df[1:10,1:5]
rownames(rf.df) = rf.df.1[,1]
rf.df[1:10,1:5]

# If implementing correction for population stratification (and/or other covariates): 
# Here we will correct the genotypes AND phenotypes using the approach of Zhao et al.(2012). 
# Note that for our categorical phenotype (classification RF) we only corrected the genotypes since the phenotype was binary
ncol(rf.df)
rf.df.corrected <- rf.df # create another data frame for the corrected phenotypes and genotypes
rf.df.corrected[,2:9884] <- NA  # Keep the first column with population ID, but then replace with NA's over which you can write the residuals.


# Now correct the genotypes AND phenotypes using the regression/residual method (Zhao et al. 2012). We're using a standard linear regression because Zhao et al. 2012 found that the correction procedure is robust to selection of the link function
# here, performing simple linear regression starting with col2 (to include phenotype), whereas with categorical phenotype we start
# at col3 (first col of genotypes)
for (i in 2:ncol(rf.df)){
  #LM_SNP_i <- lm(rf.df[,i] ~ factor(rf.df$Pop)) # categorical population variable: apply linear model to all loci and the response
  LM_SNP_i <- lm(rf.df[,i] ~ rf.df$Pop) # continuous population variable: apply linear model to all loci and the response
  rf.df.corrected[,i] <- LM_SNP_i$residuals 
  colnames(rf.df.corrected)[i]<-colnames(rf.df)[i] 
  if(i%%50==0) print(i)
}


# Verify that the residuals have been written to the data frame properly, using the last column as an example
rf.df.corrected[,9884]-LM_SNP_i$residuals  #Should all be zero if correct

# Export a copy of the corrected genotypes and phenotypes for future reference
write.csv(rf.df.corrected,file="input.impute.matrix.pheno.corrected.csv")





####################################
#########Ready for Random Forest!!!!
####################################

# Now run Random Forest analysis. Since this is a continuous trait, we need to conduct a regression RF

# First, we need to optimize mtry by running different values of mtry at different values of ntree. 

# We will run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p, where p is the number of loci
# We will initially run each of these mtry values at ntree=100 to 1000 (by increments of 100). 
# We are looking for a plateau where the proportion variation explained (PVE) stops increasing with larger values of ntree
# Once we reach the plateau, we will choose the mtry value that maximizes PVE.

sqrt(9884) #change this value based on the number of loci in your dataset
2*sqrt(9884)
0.1*9884
0.2*9884
9884/3


results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)
for (i in seq(from = 100, to = 1000, by = 100)){  # values for ntree; note did not reach plateau with ntree set to 100-2500 in increments of 100.
  print(i)
  for (j in c(99,199,988,1977,3295,9884)){    #values of mtry based on total no. loci
    rf_ij <- randomForest(x = rf.df[,3:9884], y = rf.df$LifeHist, importance=TRUE ,proximity=TRUE, ntree=i, mtry=j)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$rsq,1)))
  }
}


# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","PVE")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 99],results_optimization$PVE[results_optimization$mtry == 99], type="l", col="black", xlab="ntree",ylab="PVE", ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 199],results_optimization$PVE[results_optimization$mtry == 199], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 988],results_optimization$PVE[results_optimization$mtry == 988], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 1977],results_optimization$PVE[results_optimization$mtry == 1977], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 3295],results_optimization$PVE[results_optimization$mtry == 3295], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 9884],results_optimization$PVE[results_optimization$mtry == 9884], col="red")
legend(200, 0.8, c("mtry99","mtry199","mtry988","mtry1977","mtry3295","mtry9884"),
       pt.bg=c("black","blue","green","purple","orange","red"),
       pch=c(21,21,21,21,21,21), bty='n', x.intersp=.5, y.intersp = .6, text.font=1, 
       pt.cex=2, cex = 1)

# remember, want the mtry that has the GREATEST PVE!


# Now begin the full Random Forest analyses

# Recall that even though we optimized mtry, we must now run a larger number of trees in order to achieve convergence of importance values between forests.
# Start with ntree = 25,000 then increase from there if correlation between forests is not high enough
rf_all_1 = randomForest(x = rf.df.corrected[,3:9884], y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=9884, ntree=200000)
save(rf_all_1,file="rf_all_1.Rdata")
load("rf_all_1.Rdata")

rf_all_2 = randomForest(x = rf.df.corrected[,3:9884], y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=9884, ntree=200000)
save(rf_all_2,file="rf_all_2.Rdata")
load("rf_all_2.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) # Type 1 is mean decrease in accuracy 
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2) #I typically look for a correlation of >= 0.95 between forests. Increase ntree to get near that target.

rf_all_3 = randomForest(x = rf.df.corrected[,3:9884], y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=9884, ntree=200000)
save(rf_all_3,file="rf_all_3.Rdata")
load("rf_all_3.Rdata")

importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")

# The predictive ability of regression trees is measured by the proportion of variation explained.  A value is calculated for each tree within a forest.
# For each forest, we will take the proportion variation explained from the last tree (after convergence)

rf_all_1_rsq <- tail(rf_all_1$rsq,1)
rf_all_2_rsq <- tail(rf_all_2$rsq,1)
rf_all_3_rsq <- tail(rf_all_3$rsq,1)

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")

# Export importance values for future reference
write.csv(importance_rf_all,file="rf_importance_values_all_loci.csv",row.names=FALSE)

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of return time. 
# For each subset, we will use mtry=p since that is the optimal setting that we previously found.

##### Best 2% 

names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes
genotypes_2perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_2perc_unique]

# Now conduct RF on this subset
rf_2perc_1 = randomForest(x=genotypes_2perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=200000)
save(rf_2perc_1,file="rf_2perc_1.Rdata")
rf_2perc_2 = randomForest(x=genotypes_2perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=200000)
save(rf_2perc_2,file="rf_2perc_2.Rdata")
rf_2perc_3 = randomForest(x=genotypes_2perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=200000)
save(rf_2perc_3,file="rf_2perc_3.Rdata")

# Extract and save the r-squared values for this subset of loci
rf_2perc_1_rsq <- tail(rf_2perc_1$rsq,1)
rf_2perc_2_rsq <- tail(rf_2perc_2$rsq,1)
rf_2perc_3_rsq <- tail(rf_2perc_3$rsq,1)

rm(rf_2perc_1,rf_2perc_2,rf_2perc_3) # remove the objects to save memory in R 

##### Best 3% 

names_best_3perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.97))]
names_best_3perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.97))]
names_best_3perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.97))]
names_best_3perc_unique<-unique(c(names_best_3perc_1,names_best_3perc_2,names_best_3perc_3))

# Extract genotypes
genotypes_3perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_3perc_unique]

# Now conduct RF on this subset
rf_3perc_1 = randomForest(x=genotypes_3perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=200000)
save(rf_3perc_1,file="rf_3perc_1.Rdata")
rf_3perc_2 = randomForest(x=genotypes_3perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=200000)
save(rf_3perc_2,file="rf_3perc_2.Rdata")
rf_3perc_3 = randomForest(x=genotypes_3perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=200000)
save(rf_3perc_3,file="rf_3perc_3.Rdata")

# Extract and save the r-squared values for this subset of loci
rf_3perc_1_rsq <- tail(rf_3perc_1$rsq,1)
rf_3perc_2_rsq <- tail(rf_3perc_2$rsq,1)
rf_3perc_3_rsq <- tail(rf_3perc_3$rsq,1)

rm(rf_3perc_1,rf_3perc_2,rf_3perc_3) # remove the objects to save memory in R 

##### Best 4% 

names_best_4perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.96))]
names_best_4perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.96))]
names_best_4perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.96))]
names_best_4perc_unique<-unique(c(names_best_4perc_1,names_best_4perc_2,names_best_4perc_3))
write.csv(names_best_4perc_unique, "names_best_4perc_unique.csv")

# Extract genotypes
genotypes_4perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_4perc_unique]

# Now conduct RF on this subset
rf_4perc_1 = randomForest(x=genotypes_4perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=200000)
save(rf_4perc_1,file="rf_4perc_1.Rdata")
rf_4perc_2 = randomForest(x=genotypes_4perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=200000)
save(rf_4perc_2,file="rf_4perc_2.Rdata")
rf_4perc_3 = randomForest(x=genotypes_4perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=200000)
save(rf_4perc_3,file="rf_4perc_3.Rdata")

# Extract and save the r-squared values for this subset of loci
rf_4perc_1_rsq <- tail(rf_4perc_1$rsq,1)
rf_4perc_2_rsq <- tail(rf_4perc_2$rsq,1)
rf_4perc_3_rsq <- tail(rf_4perc_3$rsq,1)

rm(rf_4perc_1,rf_4perc_2,rf_4perc_3) # remove the objects to save memory in R 


##### Best 5% 

names_best_5perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.95))]
names_best_5perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.95))]
names_best_5perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.95))]
names_best_5perc_unique<-unique(c(names_best_5perc_1,names_best_5perc_2,names_best_5perc_3))

# Extract genotypes
genotypes_5perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_5perc_unique]

# Now conduct RF on this subset
rf_5perc_1 = randomForest(x=genotypes_5perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=200000)
save(rf_5perc_1,file="rf_5perc_1.Rdata")
rf_5perc_2 = randomForest(x=genotypes_5perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=200000)
save(rf_5perc_2,file="rf_5perc_2.Rdata")
rf_5perc_3 = randomForest(x=genotypes_5perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=200000)
save(rf_5perc_3,file="rf_5perc_3.Rdata")

# Extract and save the r-squared values for this subset of loci
rf_5perc_1_rsq <- tail(rf_5perc_1$rsq,1)
rf_5perc_2_rsq <- tail(rf_5perc_2$rsq,1)
rf_5perc_3_rsq <- tail(rf_5perc_3$rsq,1)

rm(rf_5perc_1,rf_5perc_2,rf_5perc_3) # remove the objects to save memory in R 


##### Best 10% 

names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))

# Extract genotypes
genotypes_10perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_10perc_unique]

# Now conduct RF on this subset
rf_10perc_1 = randomForest(x=genotypes_10perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=200000)
save(rf_10perc_1,file="rf_10perc_1.Rdata")
rf_10perc_2 = randomForest(x=genotypes_10perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=200000)
save(rf_10perc_2,file="rf_10perc_2.Rdata")
rf_10perc_3 = randomForest(x=genotypes_10perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=200000)
save(rf_10perc_3,file="rf_10perc_3.Rdata")

# Extract and save the r-squared values for this subset of loci
rf_10perc_1_rsq <- tail(rf_10perc_1$rsq,1)
rf_10perc_2_rsq <- tail(rf_10perc_2$rsq,1)
rf_10perc_3_rsq <- tail(rf_10perc_3$rsq,1)

rm(rf_10perc_1,rf_10perc_2,rf_10perc_3) # remove the objects to save memory in R 

##### Best 20% 

names_best_20perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.80))]
names_best_20perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.80))]
names_best_20perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.80))]
names_best_20perc_unique<-unique(c(names_best_20perc_1,names_best_20perc_2,names_best_20perc_3))

# Extract genotypes
genotypes_20perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_20perc_unique]

# Now conduct RF on this subset
rf_20perc_1 = randomForest(x=genotypes_20perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=200000)
save(rf_20perc_1,file="rf_20perc_1.Rdata")
rf_20perc_2 = randomForest(x=genotypes_20perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=200000)
save(rf_20perc_2,file="rf_20perc_2.Rdata")
rf_20perc_3 = randomForest(x=genotypes_20perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=200000)
save(rf_20perc_3,file="rf_20perc_3.Rdata")

# Extract and save the r-squared values for this subset of loci
rf_20perc_1_rsq <- tail(rf_20perc_1$rsq,1)
rf_20perc_2_rsq <- tail(rf_20perc_2$rsq,1)
rf_20perc_3_rsq <- tail(rf_20perc_3$rsq,1)

rm(rf_20perc_1,rf_20perc_2,rf_20perc_3) # remove the objects to save memory in R 


##### Best 30% 

names_best_30perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.70))]
names_best_30perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.70))]
names_best_30perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.70))]
names_best_30perc_unique<-unique(c(names_best_30perc_1,names_best_30perc_2,names_best_30perc_3))

# Extract genotypes
genotypes_30perc<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_best_30perc_unique]

# Now conduct RF on this subset
rf_30perc_1 = randomForest(x=genotypes_30perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=200000)
save(rf_30perc_1,file="rf_30perc_1.Rdata")
rf_30perc_2 = randomForest(x=genotypes_30perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=200000)
save(rf_30perc_2,file="rf_30perc_2.Rdata")
rf_30perc_3 = randomForest(x=genotypes_30perc, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=200000)
save(rf_30perc_3,file="rf_30perc_3.Rdata")

# Extract and save the r-squared values for this subset of loci
rf_30perc_1_rsq <- tail(rf_30perc_1$rsq,1)
rf_30perc_2_rsq <- tail(rf_30perc_2$rsq,1)
rf_30perc_3_rsq <- tail(rf_30perc_3$rsq,1)

rm(rf_30perc_1,rf_30perc_2,rf_30perc_3) # remove the objects to save memory in R 


# Now combine all of the r-squared values from the subsets and for all loci to identify a group for the backward purging approach

All_initial_rsq <- rbind(cbind(rf_all_1_rsq,rf_all_2_rsq,rf_all_3_rsq),
                         cbind(rf_2perc_1_rsq,rf_2perc_2_rsq,rf_2perc_3_rsq),
                         cbind(rf_3perc_1_rsq,rf_3perc_2_rsq,rf_3perc_3_rsq),
                         cbind(rf_4perc_1_rsq,rf_4perc_2_rsq,rf_4perc_3_rsq),
                         cbind(rf_5perc_1_rsq,rf_5perc_2_rsq,rf_5perc_3_rsq),
                         cbind(rf_10perc_1_rsq,rf_10perc_2_rsq,rf_10perc_3_rsq),
                         cbind(rf_20perc_1_rsq,rf_20perc_2_rsq,rf_20perc_3_rsq),
                         cbind(rf_30perc_1_rsq,rf_30perc_2_rsq,rf_30perc_3_rsq))

# Plot r-squared values for the various subsets
All_initial_rsq<-data.frame(All_initial_rsq)
All_initial_rsq$Number_loci<-c(9882,length(names_best_2perc_unique),length(names_best_3perc_unique),length(names_best_4perc_unique),length(names_best_5perc_unique),length(names_best_10perc_unique),length(names_best_20perc_unique),length(names_best_30perc_unique))
rownames(All_initial_rsq)<-c("All","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_rsq$Average<-apply(All_initial_rsq[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_rsq,file="All_initial_rsq.csv")

# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_rsq$Number_loci,All_initial_rsq$Average,log="x", pch=19,xlab="Number of Loci", ylab="Proportion Variation Explained",cex.lab=1.5,cex.axis=1.5)

#select set of loci with highest PVE (be conservative) to continue with backwards purging



#################### Backward purging approach

#this part takes a long time, probably best to run on hpcc

names_purging <- names_best_4perc_unique
write.csv(names_best_4perc_unique, "names_best_4perc_unique.csv")

genotypes_purging<-rf.df.corrected[,colnames(rf.df.corrected) %in% names_purging]

rf_purging_1 = randomForest(x=genotypes_purging, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=200000)
save(rf_purging_1,file="rf_purging_1.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=200000)
save(rf_purging_2,file="rf_purging_2.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=200000)
save(rf_purging_3,file="rf_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[length(names_purging)]]<-names_purging
rsq_best<-data.frame(V1=1:length(names_purging),V2=1:length(names_purging),V3=1:length(names_purging))
rownames(rsq_best)<-1:length(names_purging)
rsq_best[length(names_purging),] <- c(tail(rf_purging_1$rsq,1),tail(rf_purging_2$rsq,1),tail(rf_purging_3$rsq,1))


for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1<-data.frame(importance(rf_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging)
  colnames(imp_purging_1)<-"X.IncMSE1"
  rownames(imp_purging_2)<-colnames(genotypes_purging)
  colnames(imp_purging_2)<-"X.IncMSE2"
  rownames(imp_purging_3)<-colnames(genotypes_purging)
  colnames(imp_purging_3)<-"X.IncMSE3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(x=dont_keep,n=1),]]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]]<-names_keep
  genotypes_purging <- rf.df.corrected[,colnames(rf.df.corrected) %in% names_keep]
  rf_purging_1 = randomForest(x=genotypes_purging, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=200000)
  rf_purging_2 = randomForest(x=genotypes_purging, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=200000)
  rf_purging_3 = randomForest(x=genotypes_purging, y = rf.df.corrected$LifeHist, importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=200000)
  rsq_best[length(names_purging)-i,] <- c(tail(rf_purging_1$rsq,1),tail(rf_purging_2$rsq,1),tail(rf_purging_3$rsq,1))
  
}


rsq_best$Average<-apply(rsq_best,1,mean)
write.csv(rsq_best, file="Backward_purging_variation_explained_regression.csv") # Save the r-squared values
save(rsq_best,file="rsq_best.Rdata")

# Which group of loci explains the most variation?
best = which(rsq_best$Average==max(rsq_best$Average[-c(1)]))
best

# Export the names of the predictor loci
write.csv(names_all_iterations[[best]],file="Predictor_loci_regression.csv")

# Now plot the backward purging results. Omit the row from one locus since RF cannot be conducted with just one locus
plot(seq(1,169,1), rsq_best$Average[-c(1)], xlab="Number of Loci", ylab="Proportion Variation Explained",cex.lab=1.5,cex.axis=1.5,pch=16)
#note that length of x axis will vary depending on number of loci retained for backwards purging






