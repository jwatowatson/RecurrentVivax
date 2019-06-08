##############################################################################
# Script to explore the impact of deterministic versus probablistic phasing
# Two motivations: 
# 1) Assess if faster to probabilistically phase and with which Max_Hap_comb
# 2) Assess if probablistic phasing prevents recrudescence phasing in the real
# data (i.e. the needle in the haystack problem seen in the simulated 
# stranger-graph setting)
##############################################################################
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
source("../Genetic_Model/hap_combinations_functions.R")
source("../Genetic_Model/PlottingFunction.R")
load('../RData/Data_for_relatedness.RData') # For Fs_combined
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # Pooled MS data from BPD and VHX
RUN = T
MSs = names(Fs_Combined)


#=====================================================
# Start with simple cases with one or two recurrence
#=====================================================
if(RUN){ # Run models compare det vs prob with Max_Hap_genotypes and Max_Hap_comb 800 vs 300
  tic.clearlog()
  
  #==============================================================================
  tic(msg = 'prob_800_simple')
  Results_prob = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T,
                               Max_Hap_genotypes = 0,
                               Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE) 
  
  tic(msg = 'pdet_800_simple')
  Results_pdet = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T,
                              Max_Hap_genotypes = 70,
                              Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE) 
  
  tic(msg = 'det_800_simple')
  Results_det = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T,
                              Max_Hap_genotypes = 100,
                              Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE) 

  save(Results_prob, Results_pdet, Results_det, file = '../RData/Results_prob_v_det_simple_800.RData')
  rm(c('Results_prob', 'Results_pdet', 'Results_det'))
  
  #==============================================================================
  tic(msg = 'prob_300_simple')
  Results_prob = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, 
                               Max_Hap_genotypes = 0,
                               Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE) 
  
  tic(msg = 'pdet_300_simple')
  Results_pdet = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T,
                               Max_Hap_genotypes = 70,
                               Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE) 
  
  
  tic(msg = 'det_300_simple')
  Results_det = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, 
                              Max_Hap_genotypes = 100,
                              Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE)
  
  save(Results_prob, Results_pdet, Results_det, file = '../RData/Results_prob_v_det_simple_300.RData')
  rm(c('Results_prob', 'Results_pdet', 'Results_det'))
  
} 


#============================================
# Now explore inflated
#============================================
ind_calculated = which(MS_pooled$Episode_Identifier %in% Results_prob)
IDs_calculated = unique(MS_pooled$ID[ind_calculated])
IDs_remaining = unique(MS_pooled$ID[!MS_pooled$ID %in% IDs_calculated])
writeLines(sprintf('individuals with more than two recurrences: %s',length(IDs_remaining)))
MS_inflate = reformat_MSdata(filter(MS_pooled, ID %in% IDs_remaining), MSs = names(Fs_Combined)) # Remove all but remaining
MS_inflated = Inflate_into_pairs(MS_data = MS_inflate) # Inflate

if(RUN){ # Run models
  
  #==============================================================================
  tic(msg = 'prob_800_inflate')
  Results_inflate_prob = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T,
                                       Max_Hap_genotypes = 0,
                                       Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE)
  
  tic(msg = 'pdet_800_inflate')
  Results_inflate_pdet = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T,
                               Max_Hap_genotypes = 70, 
                               Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE) 

  tic(msg = 'det_800_inflate')
  Results_inflate_det = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T,
                                      Max_Hap_genotypes = 100, 
                                      Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE) 

  save(Results_inflate_prob, Results_inflate_pdet, Results_inflate_det, file = '../RData/Results_inflate_prob_v_det_simple_800.RData')
  rm('Results_inflate_prob', 'Results_inflate_pdet', 'Results_inflate_det')

  #==============================================================================
  tic(msg = 'prob_300_inflate')
  Results_inflate_prob = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, 
                                       Max_Hap_genotypes = 0,
                                       Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE) # with Max_Hap_genotypes = 100 and Max_Hap_comb = 800, 
  
  tic(msg = 'pdet_300_inflate')
  Results_inflate_pdet = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T,
                                       Max_Hap_genotypes = 70, 
                                       Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE) 
  
  
  tic(msg = 'det_300_inflate') 
  Results_inflate_det = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, 
                                      Max_Hap_genotypes = 100, 
                                      Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE) # with Max_Hap_genotypes = 100 and Max_Hap_comb = 800, 
  
  save(Results_inflate_prob, Results_inflate_pdet, Results_inflate_det, file = '../RData/Results_inflate_prob_v_det_simple_300.RData')
  rm('Results_inflate_prob', 'Results_inflate_pdet', 'Results_inflate_det')
} 







#============================================
# Now plot results: this bit needs finishing 
#============================================

# Check all the same
all(rownames(Results_prob) == rownames(Results_det)) 
all(rownames(Results_inflate_prob) == rownames(Results_inflate_det)) 

# ++++++++++++++++++ check row names  ++++++++++++++++++++
# James I cannot find the misnamed samples? Please can you
# provide an example? 
Rownames = list(rownames(Results_prob),
                rownames(Results_det), 
                rownames(Results_inflate_prob),
                rownames(Results_inflate_det))
lapply(Rownames, function(X){
  not_string = !(grepl('VHX', X) | grepl('BPD', X))  
  X[not_string]
})
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Extract comparisons usign probablistic vs deterministic approach
inds_to_compare = (Results_prob$Phased == 'P1_P1' & Results_det$Phased == "D_D")
inds_to_compare_inflate = (Results_inflate_prob$Phased == 'P1_P1' & Results_inflate_det$Phased == "D_D")

# Plot
par(mfrow = c(3,2))
for(state in c("C","L","I")){
  plot(Results_prob[inds_to_compare, state], 
       Results_det[inds_to_compare, state], 
       xlab = 'Probablistic', ylab = 'Deterministic', main = state, sub = 'simple')
  abline(a = 0, b = 1)
  
  plot(Results_inflate_prob[inds_to_compare_inflate, state],
       Results_inflate_det[inds_to_compare_inflate, state],
       xlab = 'Probablistic', ylab = 'Deterministic', main = state, sub = 'inflate')
  abline(a = 0, b = 1)
}

#===========================================================
# When Max_Hap_genotypes = 100, how many are probablistically
# phased and do not converge, thus potentially impacted by recrudescene frailty
#===========================================================
# First let's look at those that were not deterministically phased
Results_det[Results_det$Phased != 'D_D', ] # Zero
Results_inflate_det[Results_inflate_det$Phased != 'D_D', ] # Not so many, moreover many converged

# Extract rownames for inflated (!= episode name for inflated)
rownames_inflate_prob = rownames(Results_inflate_det)[Results_inflate_det$Phased != 'D_D']

# Extract episode names
Episodes_prob_inflate = unique(sapply(strsplit(rownames_inflate_prob[!is.na(rownames_inflate_prob)], 
                                               split = "%"), function(x){
                                                 ID = gsub('TID', '', x[1])
                                                 episode = do.call(rbind, strsplit(x[5], split = "_"))[,2]
                                                 paste(ID, episode, sep = '_')
                                               }))

IDs_prob_inflate = apply(do.call(rbind, strsplit(Episodes_prob_inflate, '_'))[,1:2], 1, paste, collapse = "_")

# How many IDs might be impacted? 
length(IDs_prob_inflate) # 5 IDs
mean(unique(MS_pooled$ID) %in% IDs_prob_inflate)*100 # 2% episodes
MS_pooled_prob = MS_pooled[MS_pooled$ID %in% IDs_prob_inflate,]

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
# We need to process the results to check for inflated (not yet done)
ColorPlot_MSdata(MS_data = reformat_MSdata(MS_pooled[MS_pooled$ID %in% IDs_prob_inflate,])) # Plot of all data 
