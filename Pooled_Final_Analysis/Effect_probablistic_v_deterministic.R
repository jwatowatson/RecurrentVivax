###############################################
# Script to explore the different impact of 
# deterministic versus probablistic phasing
# To-do: understand the what is wrong with VHX_214_3
###############################################

rm(list = ls())
load('../RData/RPackages_List.RData')
new.pkg <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(new.pkg) > 0){ # Install using .R script
  stop('Please run the script Install_and_load_required_packages.R before returning to this script. Thanks.')}else{
    sapply(pkgs, require, character.only = TRUE) # Load all packages
  }

source("../Genetic_Model/Data_functions.R")
source('../Genetic_Model/iGraph_functions.R')
source("../Genetic_Model/post_prob_CLI.R") 
source("../Genetic_Model/test_Rn_compatible.R") 
source("../Genetic_Model/Data_Inflation_Functions.R")
source("../Genetic_Model/hap_combinations.R")

load('../RData/Data_for_relatedness.RData') # For Fs_combined
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # Pooled MS data from BPD and VHX

# Run models
Results_MaxHap0 = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 0)
Results_MaxHap20 = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 20) 

all(rownames(Results_MaxHap0) == rownames(Results_MaxHap20)) # Check the same
inds_to_compare = which(Results_MaxHap0$Phased == 'P_P' & Results_MaxHap20$Phased == "D_D")

# Plot
par(mfrow = c(2,2))

plot(Results_MaxHap0$L[inds_to_compare], 
     Results_MaxHap20$L[inds_to_compare], 
     xlab = 'Probablistic', ylab = 'Deterministic')
abline(a = 0, b = 1)

plot(Results_MaxHap0$I[inds_to_compare], 
     Results_MaxHap20$I[inds_to_compare], 
     xlab = 'Probablistic', ylab = 'Deterministic')
abline(a = 0, b = 1)

plot(Results_MaxHap0$C[inds_to_compare], 
     Results_MaxHap20$C[inds_to_compare],
     xlab = 'Probablistic', ylab = 'Deterministic')
abline(a = 0, b = 1)



