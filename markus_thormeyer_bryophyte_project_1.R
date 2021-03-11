#### USEFUL FUNCTIONS BY THE COURSE ####
### Observed Richness ###
# A simple count of the number of species found.
# Usually denoted by R or S
richness = function(x){
  return(sum(x>0))
}

### Abundance ###
# Simply the number of sequences in the sample
# Because of sampling and sequencing biases this is less useful than one would hope.
abundance = function(x){
  return(sum(x,na.rm=TRUE))
}

### Shannon's Diversity Index (entropy) ###
# Usually denoted by H'
# Does not have a specific mathematical explanation, but increases as diversity increases.
# Assumes infinite, well-mixed population.
# Also a measure of both richness and evenness.
shannons = function(x){
  present = x[x>0]
  p = present/sum(present)
  -sum(p*log(p))
}

### Simpson's Diversity Index ###
# Usually denoted by lambda (or l)
# Equals the probability any two randomly chosen individuals in the population are the same.
# Therefore, actually a similarity index. To get diversity, D = 1/l
# Measure of both richness and evenness.
simpsons = function(x){
  p = x/sum(x,na.rm=TRUE)
  sum(p^2)
}

### Pielou's evennenss function ###
# Often denoted J
# Measures community evenness such that 1 is completely even and 0 is uneven.
# Completely even means that all taxa have exactly the same abundance. 
# An uneven community is one dominated by one or a few large taxa and many rare taxa.
# Often used in conjunction with shannons or simpsons.
evenness = function(x){
  present = x[x>0]
  p = present/sum(present)
  -sum(p*log(p))/log(sum(x>0))
}

### Chao1 ###
# Chao1 is a non-parametric estimator of richness. That means it makes no assumptions (infinite population,
# well mixed, specific model to fit and find an asymptote). It is observed richness, corrected
# with the number of singletons (species only seen once) and doubletons (species seen twice).
# There is actually some fancy statistics behind using singletons and doubletons as correctors.
# This may just return observed richness if your dataset was already cleaned to remove rare ASVs
chao1 = function(x){
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  S = sum(x > 0)
  return(S + f1*(f1-1)/(2*(f2+1)))
} 

# Takes an otu table and turns it into an abundance table for a different rank of taxonomy.
group_by_rank = function(otu_table,taxonomy,metadata,rank){
  n_sample = nrow(otu_table)
  otu_table = cbind(t(otu_table),(taxonomy))
  if (rank>1){
    otu_table$taxa = apply(otu_table[,(n_sample+1):(n_sample+rank)],1,paste,collapse="")
  }else{
    otu_table$taxa = otu_table[,n_sample+1]
  }
  otu_table = group_by(otu_table,taxa)
  grouped_table = summarise_at(otu_table,1:n_sample,sum)
  taxa_names = grouped_table$taxa
  grouped_table = grouped_table[,-1]
  grouped_table = t(grouped_table)
  colnames(grouped_table) = taxa_names
  taxa_cols = 1:ncol(grouped_table)
  grouped_table = cbind(grouped_table,metadata)
  taxa_rel_abundance = colSums(grouped_table[,taxa_cols]/rowSums(grouped_table[,taxa_cols],na.rm=T),na.rm=T)
  grouped_table = group_by(grouped_table,SampleID)
  grouped_table = pivot_longer(grouped_table,names_to="taxa",values_to="abundance",cols=all_of(taxa_cols))
  grouped_table$taxa = factor(grouped_table$taxa,levels=names(sort(taxa_rel_abundance,decreasing=TRUE)))
  return(grouped_table)
}

################################################################################

# load in every package ever in existence
library(tidyverse)
library(dplyr)
library(data.table)

# set working directory
setwd("~/Desktop/BIOL403/Project")

# load in data
metaasv = read.csv("FINAL_master_table_melissa_bryophyte_NONRAREFIED.csv",header=TRUE,sep=",")
taxonomy = read.csv("taxonomy_bryophyte_propagated.csv",header=TRUE,sep=",")

otu_cols = c(40:1147)

# Create a dataset of stuff Im interested in, lab
metaasv = subset(metaasv, metaasv$Enivironment == "Lab")

# Write csv for later I guess
write.csv(metaasv, "bryophyte_lab_data.csv")

# GENERAL DIVERSITY METRICS

metaasv$abundance = apply(metaasv[,(40:1147)],1,abundance)
metaasv$shannons = apply(metaasv[,(40:1147)],1,shannons)
metaasv$richness = apply(metaasv[,(40:1147)],1,richness)
metaasv$simpsons = apply(metaasv[,(40:1147)],1,simpsons)
metaasv$evenness = apply(metaasv[,(40:1147)],1,evenness)
metaasv$chao1 = apply(metaasv[,(40:1147)],1,chao1)

sum(metaasv$abundance)
# total reads = 545423

## Total metrics
diversity_total <- c(mean(metaasv$shannons),
                    mean(metaasv$richness),
                    mean(metaasv$simpsons),
                    mean(metaasv$evenness))

## Bryophyte specific
metaasv_b = subset(metaasv, metaasv$Description == "Prot")
diversity_b <- c(mean(metaasv_b$shannons),
                mean(metaasv_b$richness),
                mean(metaasv_b$simpsons),
                mean(metaasv_b$evenness))

## Water specific
metaasv_w = subset(metaasv, metaasv$Description == "Water")
diversity_w <- c(mean(metaasv_w$shannons),
                 mean(metaasv_w$richness),
                 mean(metaasv_w$simpsons),
                 mean(metaasv_w$evenness))

diversity <- data.frame(diversity_total, diversity_b, diversity_w)

