###############################################
# Script to explore the different impact of 
# deterministic versus probablistic phasing
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
RUN = F

#=====================================================
# Star with simple cases with one or two recurrence
#=====================================================
if(RUN){ # Run models
  Results_MaxHap0 = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 0)
  Results_MaxHap20 = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 20)
  #save(Results_MaxHap0, Results_MaxHap20, file = 'Results_prob_v_det_simple.RData')
}

#============================================
# Now explore inflated
#============================================
ind_calculated = which(MS_pooled$Episode_Identifier %in% Results_MaxHap0)
IDs_calculated = unique(MS_pooled$ID[ind_calculated])
IDs_remaining = unique(MS_pooled$ID[!MS_pooled$ID %in% IDs_calculated])
writeLines(sprintf('individuals with more than two recurrences: %s',length(IDs_remaining)))
MS_inflate = reformat_MSdata(filter(MS_pooled, ID %in% IDs_remaining), MSs = names(Fs_Combined)) # Remove all but remaining
MS_inflated = Inflate_into_pairs(MS_data = MS_inflate) # Inflate

if(RUN){ # Run models
  Results_inflate_MaxHap0 = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 0)
  Results_inflate_MaxHap20 = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 20)
  #save(Results_inflate_MaxHap0, Results_inflate_MaxHap20, file = 'Results_inflate_prob_v_det_simple.RData')
}


#============================================
# Now plot results
#============================================
load('Results_prob_v_det_simple.RData')
load('Results_inflate_prob_v_det_simple.RData')

# Check all the same
all(rownames(Results_MaxHap0) == rownames(Results_MaxHap20)) 
all(rownames(Results_inflate_MaxHap0) == rownames(Results_inflate_MaxHap20)) 

# ++++++++++++++++++ check row names  ++++++++++++++++++++
# James I cannot find the misnamed samples? Please can you
# provide an example? 
Rownames = list(rownames(Results_MaxHap0),
                rownames(Results_MaxHap20),
                rownames(Results_inflate_MaxHap0), 
                rownames(Results_inflate_MaxHap20))
lapply(Rownames, function(X){
  not_string = !(grepl('VHX', X) | grepl('BPD', X))  
  X[not_string]
})
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Extract comparisons usign probablistic vs deterministic approach
inds_to_compare = which(Results_MaxHap0$Phased == 'P_P' & Results_MaxHap20$Phased == "D_D")
inds_to_compare_inflate = which(Results_inflate_MaxHap0$Phased == 'P_P' & Results_inflate_MaxHap20$Phased == "D_D")

# Plot
par(mfrow = c(3,2))
for(state in c("C","L","I")){
  plot(Results_MaxHap0[inds_to_compare, state], 
       Results_MaxHap20[inds_to_compare, state], 
       xlab = 'Probablistic', ylab = 'Deterministic', main = state, sub = 'simple')
  abline(a = 0, b = 1)
  
  plot(Results_inflate_MaxHap0[inds_to_compare_inflate, state], 
       Results_inflate_MaxHap20[inds_to_compare_inflate, state], 
       xlab = 'Probablistic', ylab = 'Deterministic', main = state, sub = 'inflate')
  abline(a = 0, b = 1)
}


