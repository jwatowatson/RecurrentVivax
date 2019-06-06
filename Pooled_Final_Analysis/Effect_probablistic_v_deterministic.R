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
source("../Genetic_Model/PlottingFunction.R")
load('../RData/Data_for_relatedness.RData') # For Fs_combined
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # Pooled MS data from BPD and VHX
RUN = T
MSs = names(Fs_Combined)

#=====================================================
# Star with simple cases with one or two recurrence
#=====================================================
if(RUN){ # Run models
  Results_prob = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 0)
  Results_det = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T)
  save(Results_prob, Results_det, file = '../RData/Results_prob_v_det_simple.RData')
} else {
  load('Results_prob_v_det_simple.RData')
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
  Results_inflate_prob = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 0)
  Results_inflate_det = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T)
  save(Results_inflate_prob, Results_inflate_det, file = '../RData/Results_inflate_prob_v_det_simple.RData')
} else {
  load('Results_inflate_prob_v_det_simple.RData')
}

# By changing the code we can return numbers of compatable haplotypes
# no_haps_simple = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 0)
# no_haps_inflate = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, Max_Hap_genotypes = 0)
# 
# tail(sort(no_haps_simple)) #  two biggest: 648 1296 
# table(tail(sort(no_haps_inflate), 25)) # if we Max_Hap_genotypes = 800 we could get all but 13
# 288  648 1296 1728 
# 2    1    9   13

#============================================
# Now plot results
#============================================

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
inds_to_compare = (Results_MaxHap0$Phased == 'P_P' & Results_MaxHap20$Phased == "D_D")
inds_to_compare_inflate = (Results_inflate_MaxHap0$Phased == 'P_P' & Results_inflate_MaxHap20$Phased == "D_D")

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

#===========================================================
# When Max_Hap_genotypes = 20, how many are probablistically
# phased, thus potentially impacted by recrudescene frailty
#===========================================================
# First let's look at those that were not deterministically phased
Results_MaxHap20[Results_MaxHap20$Phased != 'D_D', ] # 5 different IDs
Results_inflate_MaxHap20[Results_inflate_MaxHap20$Phased != 'D_D', ] 

# How many were probablistically phased for both
Results_MaxHap20[Results_MaxHap20$Phased == 'P_P', ] # 0 
Results_inflate_MaxHap20[Results_inflate_MaxHap20$Phased == 'P_P', ] #6 

# Extract rownames (!= episode name for inflated)
rownames_simple_prob = rownames(Results_MaxHap20)[Results_MaxHap20$Phased != 'D_D']
rownames_inflate_prob = rownames(Results_inflate_MaxHap20)[Results_inflate_MaxHap20$Phased != 'D_D']

# Extract episode names
Episodes_prob_simple = unique(rownames_simple_prob[!is.na(rownames_simple_prob)])
Episodes_prob_inflate = unique(sapply(strsplit(rownames_inflate_prob[!is.na(rownames_inflate_prob)], 
                                               split = "%"), function(x){
                                                 ID = gsub('TID', '', x[1])
                                                 episode = do.call(rbind, strsplit(x[5], split = "_"))[,2]
                                                 paste(ID, episode, sep = '_')
                                               }))

IDs_prob_simple = apply(do.call(rbind, strsplit(Episodes_prob_simple, '_'))[,1:2], 1, paste, collapse = "_")
IDs_prob_inflate = apply(do.call(rbind, strsplit(Episodes_prob_inflate, '_'))[,1:2], 1, paste, collapse = "_")

# Extract IDs
Episodes_prob = c(Episodes_prob_simple, Episodes_prob_inflate)
IDs_prob = c(IDs_prob_simple, IDs_prob_inflate)

# How many IDs might be impacted? 
length(IDs_prob) # 24 IDs
mean(unique(MS_pooled$ID) %in% IDs_prob)*100 # 8%
MS_pooled_prob = MS_pooled[MS_pooled$ID %in% IDs_prob,]

# COI patterns and number of markers typed of those potentially impacted:
# most have some potential for recrudescence (decreasing diversity)
# all very high number of markers, so very potentially impacted
ddply(MS_pooled_prob, .variables = 'ID', .fun = function(x){
  coi_pattern = paste(table(x$Episode_Identifier), collapse = '_')
  num_markers = sum(sapply(MSs, function(ms) any(!is.na(x[,ms]))))
  data.frame(coi_pattern = coi_pattern, num_markers = num_markers)
  })

#-----------------------------------------------------------------
# Lets visually inspect these data to see if we can find evidence of 
# clones in the data that are not captured in the results
#-----------------------------------------------------------------
# Inspect simple first 
ColorPlot_MSdata(MS_data = reformat_MSdata(MS_pooled[MS_pooled$ID %in% IDs_prob_simple,])) # Plot of all data 
X = Results_MaxHap20[Results_MaxHap20$Phased != 'D_D', ] # 5 different IDs (501 at top of plot)
X[nrow(X):1,] # Reverse order to match plot
# VHX_501 results make sense
# VHX_352 results make sense
# VHX_239 results make sense
# BPD_598 results make sense
# BPD_562 results make sense 

# We need to process the results to check for inflated (not yet done)
ColorPlot_MSdata(MS_data = reformat_MSdata(MS_pooled[MS_pooled$ID %in% IDs_prob_inflate,])) # Plot of all data 
