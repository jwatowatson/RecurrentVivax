##############################################################################
# Script in which UpperComplexityIDs are extracted and saved
#
# Original motivation: explore impact of deterministic vs probablistic phasing
# 1) Assess if faster to probabilistically phase and with which Max_Hap_comb
# 2) Assess if probablistic phasing impacts recrudescence inference
# 
# ------- Finding ----------------------------------------------------------
# Summary: IDs that do not pass UpperComplexity considerably slow computation. Since no 
# results are never obtained from these, preferable to remove ahead of computation
#
# Here, these IDs are defined as any for which one or more episode comparison 
# (either simple or inflated) returns NA when Max_Hap_genotypes = 100 and Max_Hap_comb = 800. 
# Note that, when Max_Hap_genotypes = 0 for VHX_16, and when Max_Hap_comb = 300 for VHX_239, 
# both are probabilistically phased s.t. the hap combination count is capped 
# and they pass UpperComplexity. However, ultimately we opt for Max_Hap_genotypes > 0
# (deterministic phasing is faster) so inclusion of VHX_16 is pointless. 
# Moreover, there is little point in setting Max_Hap_comb = 300 to enable analysis of VHX_239 
# since doing so is prohibitively slow and returns P0 results, while increasing Max_Hap_comb 
# with inclusion of VHX_239 returns NA results for VHX_239; we thus exclude VHX_239. 
#
# James, this strategy results in 9 "too complex" whereas before we had 7
# Do you know how these additional 2 may have evaded being considered "too complex" previously? 
# e.g. were NAs in adjacency matrices ignored? (see lines 241-242)
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
RUN = F
MSs = names(Fs_Combined)


#=======================================================================================
# Run Max_Hap_genotypes = 0 (prob) and 100 (det) with Max_Hap_comb 800 vs 300
#=======================================================================================
if(RUN){ 
  tic.clearlog()
  
  tic(msg = 'prob_800_simple')
  Results_prob = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T,
                               Max_Hap_genotypes = 0, # All probablistically phased
                               Max_Hap_comb = 800) # Applies to all episodes 
  toc(log = TRUE, quiet = TRUE) 
  
  tic(msg = 'det_800_simple')
  Results_det = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T,
                              Max_Hap_genotypes = 100, # Most deterministically phased
                              Max_Hap_comb = 800) # Applies to episodes with > Max_Hap_genotypes only
  toc(log = TRUE, quiet = TRUE) 
  
  save(Results_prob, Results_det, file = '../RData/Results_prob_v_det_simple_800.RData')
  rm(list = c('Results_prob', 'Results_det'))
  
  tic(msg = 'prob_300_simple')
  Results_prob = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, 
                               Max_Hap_genotypes = 0,
                               Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE)
  
  
  tic(msg = 'det_300_simple')
  Results_det = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, verbose = T, 
                              Max_Hap_genotypes = 100,
                              Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE)
  
  save(Results_prob, Results_det, file = '../RData/Results_prob_v_det_simple_300.RData')
  rm(list = c('Results_prob', 'Results_det'))
} 


#======================================================================================
# Run inflated Max_Hap_genotypes = 0 (prob) and 100 (det) with Max_Hap_comb 800 vs 300
#======================================================================================
load('../RData/Results_prob_v_det_simple_800.RData')
ind_calculated = which(MS_pooled$Episode_Identifier %in% rownames(Results_prob))
IDs_calculated = unique(MS_pooled$ID[ind_calculated])
IDs_remaining = unique(MS_pooled$ID[!MS_pooled$ID %in% IDs_calculated])
writeLines(sprintf('individuals with more than two recurrences: %s',length(IDs_remaining)))
MS_inflate = reformat_MSdata(filter(MS_pooled, ID %in% IDs_remaining), MSs = names(Fs_Combined)) # Remove all but remaining
MS_inflated = Inflate_into_pairs(MS_data = MS_inflate) # Inflate

if(RUN){
  
  tic(msg = 'prob_800_inflate')
  Results_inflate_prob = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T,
                                       Max_Hap_genotypes = 0,
                                       Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE)
  
  tic(msg = 'det_800_inflate')
  Results_inflate_det = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T,
                                      Max_Hap_genotypes = 100, 
                                      Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE) 
  
  save(Results_inflate_prob, Results_inflate_det, file = '../RData/Results_inflate_prob_v_det_simple_800.RData')
  rm(list = 'Results_inflate_prob', 'Results_inflate_det')
  
  tic(msg = 'prob_300_inflate')
  Results_inflate_prob = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, 
                                       Max_Hap_genotypes = 0,
                                       Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE) # with Max_Hap_genotypes = 100 and Max_Hap_comb = 800, 
  
  tic(msg = 'det_300_inflate') 
  Results_inflate_det = post_prob_CLI(MSdata = MS_inflated, Fs = Fs_Combined, verbose = T, 
                                      Max_Hap_genotypes = 100, 
                                      Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE) # with Max_Hap_genotypes = 100 and Max_Hap_comb = 800, 
  
  save(Results_inflate_prob, Results_inflate_det, file = '../RData/Results_inflate_prob_v_det_simple_300.RData')
  rm(list = 'Results_inflate_prob', 'Results_inflate_det')
} 

if(RUN){
  log.txt <- tic.log(format = TRUE)
  save(log.txt, file = '../RData/log.txt_prob_vs_det.RData')
} else {
  load(file = '../RData/log.txt_prob_vs_det.RData')
}

#=================================================================================
# Times
# In spite of simulated results (see Explore_soln_prob_phasing_C_frailty.R),
# For simple, Max_Hap_comb = 800 faster than 300 
# For inflate, Max_Hap_comb = 300 faster than 800
#=================================================================================
log.txt # Load times

# Load results
load('../RData/Results_prob_v_det_simple_800.RData')
Results_prob_800 = Results_prob; 
Results_det_800 = Results_det 

rm(list = c('Results_prob', 'Results_det'))
load('../RData/Results_prob_v_det_simple_300.RData')
Results_prob_300 = Results_prob; 
Results_det_300 = Results_det

rm(list = c('Results_prob', 'Results_det'))
load('../RData/Results_inflate_prob_v_det_simple_800.RData')
Results_inflate_prob_800 = Results_inflate_prob; 
Results_inflate_det_800 = Results_inflate_det 

rm(list = c('Results_inflate_prob', 'Results_inflate_det'))
load('../RData/Results_inflate_prob_v_det_simple_300.RData')
Results_inflate_prob_300 = Results_inflate_prob; 
Results_inflate_det_300 = Results_inflate_det

# Check all the same
all(rownames(Results_det_800) == rownames(Results_prob_800) & 
      rownames(Results_prob_300) == rownames(Results_det_300) & 
      rownames(Results_prob_800) == rownames(Results_prob_300))

# Extract episodes that do not pass UpperComplexity in simple case regardless of parameters
UpperComplexityEps = rownames(Results_det_800)[is.na(Results_det_800$C)]
UpperComplexityIDs = sapply(strsplit(UpperComplexityEps, split = '_'), 
                            function(x){paste(x[1:2], collapse = '_')})

# VHX_239 (almost certainly a relapse based on data) is computed for Max_Hap_comb = 300 
# but only because it is phased probabilistically thereby capping the no. 
# of haplotype combinations and allowing it to evade the UpperComplexity test
Results_det_800[grepl("VHX_239", rownames(Results_det_800)), ]
Results_prob_800[grepl("VHX_239", rownames(Results_prob_800)), ]
Results_det_300[grepl("VHX_239", rownames(Results_det_300)), ]
Results_prob_300[grepl("VHX_239", rownames(Results_prob_300)), ]

# These are IDs that are too complex for simple and do not feature in inflate
UpperComplexityIDs %in% IDs_calculated 
UpperComplexityIDs %in% IDs_remaining

# Check methods give the same results: Yes 
Results_prob_800[rowSums(Results_prob_800 == Results_prob_300) != 5,]
Results_prob_300[rowSums(Results_prob_800 == Results_prob_300) != 5,]

# In general, Time also scales with Max_Hap_comb using the real data 
max_hap_combs = c(100, 300, 500, 800) # When we allow more than 1000 with sim data, start to skip problems as too complex
time_store = array(dim = c(length(IDs_calculated), length(max_hap_combs)), 
                   dimnames = list(IDs_calculated, max_hap_combs))

if(RUN){
  for(ID in IDs_calculated){
    inds = MS_pooled$ID %in% ID
    
    tic.clearlog()
    Results = t(sapply(max_hap_combs, function(x){
      tic()
      TH = post_prob_CLI(MSdata = MS_pooled[inds,], Fs = Fs_Combined, 
                         Max_Hap_genotypes = 0,
                         Max_Hap_comb = x)
      toc(log = T)
      TH 
    }))
    Time_msg = unlist(tic.log(format = T))
    Times = as.numeric(do.call(rbind, strsplit(Time_msg, split = ' '))[,1])
    time_store[ID,] = Times
  }
  save(time_store, file = '../RData/time_per_ID_prob.RData')
} else {
  load('../RData/time_per_ID_prob.RData') 
}

# Inspect time_store manually: 
time_store 
matplot(t(time_store), type = 'l', xaxt = 'n')
axis(side = 1, at = 1:length(max_hap_combs), labels = max_hap_combs)

# Samples for which time does not scale with max_hap_combs
odd_samples = names(which(apply(time_store, 1, function(x)any(diff(x) < 0)))) 
time_store[odd_samples,] # Main culprit: VHX_239 also VHX_52 in UpperComplexityIDs
odd_samples %in% UpperComplexityIDs
axis(side = 1, at = 1:length(max_hap_combs), labels = max_hap_combs)

#================================================================================================
# Let's save the names of the samples for which results are NA 
# and exclude from future computations ahead of applying function
#================================================================================================
UpperComplexityEpsInflate = names(which(apply(is.na(Results_inflate_det_800), 1, all)))
Results_inflate_prob[UpperComplexityEpsInflate, ]
UpperComplexityIDsInflate = unique(gsub('TID', '', do.call(rbind, strsplit(UpperComplexityEpsInflate, split = '%'))[,1]))
UpperComplexityIDs %in% UpperComplexityIDsInflate
UpperComplexityIDs = unique(c(UpperComplexityIDs, UpperComplexityIDsInflate))
save(UpperComplexityIDs, file = '../RData/UpperComplexityIDs.RData')

# Based on notes, complex pre-revision: "VHX_239_2" "VHX_33_2"  "VHX_39_2"  "VHX_461_2" "VHX_52_2"  "VHX_583_2"
pre_revision = c("VHX_239", "VHX_33", "VHX_39", "VHX_461", "VHX_52", "VHX_583")

# How were these below previously calculated? Were NAs in adjacency matrices previously ignored? 
UpperComplexityIDs[!UpperComplexityIDs %in% pre_revision] # Sample deemed not too complex previously
UpperComplexityIDs[!UpperComplexityIDs %in% pre_revision] %in% IDs_remaining # All in inflated

writeLines(sprintf("There are %s IDs that are too complex when Max_Hap_genotypes = 100 and Max_Hap_comb = 800: %s", 
                   length(UpperComplexityIDs), paste(UpperComplexityIDs, collapse = " ")))


# ================================================================================================
# Let's check probablistic converged and deterministic give the same results
# Focus on 800 since more likely P1
# Conlcusion: confident that if converged, probablistic gives same result as deterministic
# ================================================================================================
# First, what proportion of probablistic converged can be evaluated (i.e. can be compared to D_D)
# Simple uninflated data:
denom = sum(Results_prob_800$Phased == 'P1_P1', na.rm = T)
numtr = sum(Results_det_800$Phased[Results_prob_800$Phased == 'P1_P1'] == 'D_D', na.rm = T)
writeLines(sprintf('Of %s probablistic converged (noninflated data), %s percent can be evaluated', numtr, numtr/denom*100))

# Inflated data:
denom = sum(Results_inflate_prob_800$Phased == 'P1_P1', na.rm = T)
numtr = sum(Results_inflate_det_800$Phased[Results_inflate_prob_800$Phased == 'P1_P1'] == 'D_D', na.rm = T)
writeLines(sprintf('Of %s probablistic converged (inflated data), %s percent can be evaluated', numtr, round(numtr/denom*100)))

# Comparison: perfect agreement visually
inds_to_compare = (Results_prob_800$Phased == 'P1_P1' & Results_det_800$Phased == "D_D")
inds_to_compare_inflate = (Results_inflate_prob_800$Phased == 'P1_P1' & Results_inflate_det_800$Phased == "D_D")

# Plot
par(mfrow = c(3,2))
for(state in c("C","L","I")){
  plot(Results_prob_800[inds_to_compare, state], 
       Results_det_800[inds_to_compare, state], 
       xlab = 'Probablistic', ylab = 'Deterministic', main = state, sub = 'simple')
  abline(a = 0, b = 1)
  
  plot(Results_inflate_prob_800[inds_to_compare_inflate, state],
       Results_inflate_det_800[inds_to_compare_inflate, state],
       xlab = 'Probablistic', ylab = 'Deterministic', main = state, sub = 'inflate')
  abline(a = 0, b = 1)
}

# Zero un-evaluated simple 800 (since 239 excluded)
Results_prob_800[Results_prob_800$Phased != 'P1_P1', ] 
Results_det_800[Results_det_800$Phased != 'D_D', ]

# Several un-evaluated inflate, those with P0 in UpperComplexityIDs
Results_inflate_prob_800[Results_inflate_prob_800$Phased != 'P1_P1', ]
Results_inflate_det_800[Results_inflate_det_800$Phased != 'D_D', ]

#====================================================================
# When using det Max_Hap_genotypes = 800 vs 300, how many are probablistically
# phased and do not converge, thus potentially impacted by recrudescene frailty
# When Max_Hap_genotypes = 100: "VHX_33"  "VHX_551" "VHX_646" 
# When Max_Hap_genotypes = 0: "VHX_16"  "VHX_225" "VHX_33"  "VHX_551" "VHX_583" "VHX_646" 
# Regardless of Max_Hap_genotypes, all in UpperComplexityIDs and so once these 
# are removed, zero included will be potentially impacted by recrudescene
#===========================================================
# First let's look at those that did not converge
ind_P0_800 = grepl('P0', Results_inflate_prob_800$Phased) 
ind_P0_300 = grepl('P0', Results_inflate_prob_300$Phased) 

# Look at the results 
Results_inflate_prob_800[ind_P0_800,]
Results_inflate_prob_300[ind_P0_300,]

# Extract rownames for inflated (!= episode name for inflated)
names_P0_800 = rownames(Results_inflate_prob_800)[ind_P0_800]
names_P0_300 = rownames(Results_inflate_prob_300)[ind_P0_300]

# Extract episode names
IDs_P0_800 = unique(gsub('TID', '', do.call(rbind, strsplit(names_P0_800, split = "%"))[,1]))
IDs_P0_300 = unique(gsub('TID', '', do.call(rbind, strsplit(names_P0_300, split = "%"))[,1]))
setdiff(IDs_P0_800, IDs_P0_300) # Makes no difference if 800 or 300: 

# Having removed complexity, and given 300 and 800 give same P0s
# Following suggest best Max_Hap_genotypes could be 300 
colSums(time_store[!rownames(time_store) %in% UpperComplexityIDs, ]) 

# Data summary for potentially impacted IDs  
# COI patterns and number of markers typed of those potentially impacted:
# most have some potential for recrudescence (decreasing diversity)
# all very high number of markers, so very potentially impacted
MS_pooled_prob = MS_pooled[MS_pooled$ID %in% IDs_P0_800,]
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
ColorPlot_MSdata(MS_data = reformat_MSdata(MS_pooled[MS_pooled$ID %in% "VHX_16",])) # Plot of some data 



#-----------------------------------------------------------------
# Could some have been probablistically phased with a dynamical approach
# That just caps Max_hap_comb? Maybe VHX_16, VHX_225, VHX_583
#-----------------------------------------------------------------
MS_pooled_prob = MS_pooled[MS_pooled$ID %in% UpperComplexityIDs,]
Summary = ddply(MS_pooled_prob, .variables = 'ID', .fun = function(x){
  coi_pattern = paste(table(x$Episode_Identifier), collapse = '_')
  num_markers = sum(sapply(MSs, function(ms) any(!is.na(x[,ms]))))
  data.frame(coi_pattern = coi_pattern, num_markers = num_markers)
})

Summary$hapcnt = sapply(Summary$ID, function(id){
indsimple = grepl(paste(id,'_',sep=''),  rownames(Results_prob_300))
  indinflate = grepl(paste(id,'%',sep=''),  rownames(Results_inflate_prob_300))
  if(any(indinflate)){
    z = max(as.numeric(unlist(strsplit(Results_inflate_prob_300$Hapcnt[indinflate], split = '_'))), na.rm = T)
  } 
  if(any(indsimple)){
    z = max(as.numeric(unlist(strsplit(Results_prob_300$Hapcnt[indsimple], split = '_'))), na.rm = T)
  }
  if(all(!c(indinflate, indsimple))){z = NA}
  z
})

print(Summary)