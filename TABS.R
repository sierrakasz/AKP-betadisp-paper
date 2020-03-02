## Anna Karenina Kaszubinski et al.

#Tables


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


# Table 1 -----------------------------------------------------------------
#collapse metadata into case rather than by sample area for case metrics
metadata_coll <- subset(metadata, select=-c(SampleID, X, Sample_Area, Description))
metadata_coll <- unique(metadata_coll)
metadata_coll <- metadata_coll[complete.cases(metadata_coll), ]

#age, race, sex, and BMI distributed throughout M/COD
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

#final table was made in excel


# Table 2 -----------------------------------------------------------------

# final table was made in excel from the results in Table S9. Please see
# supplemental table code for full code. 
