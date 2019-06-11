##############################################################################
# Script to explore the impact of deterministic versus probablistic phasing
# AFTER EXCLUSION OF UpperComplexityIDs
#
# Two motivations: 
# 1) Assess if faster to probabilistically phase: no
# 2) Assess if Max_Hap_comb impacts recrudescence inference: no
# 
# Conclusion: run models by first excluding UpperComplexityIDs with
# Max_Hap_genotypes = 50 and Max_Hap_comb = 500
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
load('../RData/UpperComplexityIDs.RData') # IDs to avoid
RUN = F
MSs = names(Fs_Combined)


#=====================================================
# Start with simple cases with one or two recurrence
#=====================================================
if(RUN){ 
  tic.clearlog()
  
  tic(msg = 'prob_300_simple') # Slowest simple inc. complexity
  Results_prob = post_prob_CLI(MSdata = MS_pooled[!MS_pooled$ID %in% UpperComplexityIDs,], 
                               Fs = Fs_Combined, verbose = T,
                               Max_Hap_genotypes = 0, # All probablistically phased
                               Max_Hap_comb = 300) # Applies to all episodes 
  toc(log = TRUE, quiet = TRUE) 

  tic(msg = 'pdet_500_simple') # Intermediate 
  Results_pdet = post_prob_CLI(MSdata = MS_pooled[!MS_pooled$ID %in% UpperComplexityIDs,], 
                               Fs = Fs_Combined, verbose = T,
                              Max_Hap_genotypes = 50,  
                              Max_Hap_comb = 500) # Applies to episodes with > Max_Hap_genotypes only
  toc(log = TRUE, quiet = TRUE) 
  
  tic(msg = 'det_800_simple') # Fastest simple inc. complexity
  Results_det = post_prob_CLI(MSdata = MS_pooled[!MS_pooled$ID %in% UpperComplexityIDs,], 
                              Fs = Fs_Combined, verbose = T,
                              Max_Hap_genotypes = 100,  
                              Max_Hap_comb = 800) # Applies to episodes with > Max_Hap_genotypes only
  toc(log = TRUE, quiet = TRUE) 

  save(Results_prob, Results_pdet, Results_det, file = '../RData/Results_prob_v_det_simple_exc.NA.IDs.RData')
  rm(list = c('Results_prob', 'Results_pdet', 'Results_det'))
} 


#============================================
# Now explore inflated
#============================================
load('../RData/Results_prob_v_det_simple_exc.NA.IDs.RData')
ind_calculated = which(MS_pooled$Episode_Identifier %in% rownames(Results_prob))
IDs_calculated = unique(MS_pooled$ID[ind_calculated])
IDs_remaining = unique(MS_pooled$ID[!MS_pooled$ID %in% IDs_calculated])
writeLines(sprintf('individuals with more than two recurrences: %s',length(IDs_remaining)))
IDs_remaining_excNAs = IDs_remaining[!IDs_remaining %in% UpperComplexityIDs]
MS_inflate = reformat_MSdata(filter(MS_pooled, ID %in% IDs_remaining_excNAs), MSs = names(Fs_Combined)) # Remove all but remaining
MS_inflated = Inflate_into_pairs(MS_data = MS_inflate) # Inflate

if(RUN){ # Run models

  tic(msg = 'prob_300_inflate') # Fastest inflate inc. complex
  Results_inflate_prob = post_prob_CLI(MSdata = MS_inflated, 
                                       Fs = Fs_Combined, verbose = T,
                                       Max_Hap_genotypes = 0,
                                       Max_Hap_comb = 300)
  toc(log = TRUE, quiet = TRUE)
  
  tic(msg = 'pdet_500_inflate') # Intermediate
  Results_inflate_pdet = post_prob_CLI(MSdata = MS_inflated, 
                                       Fs = Fs_Combined, verbose = T,
                               Max_Hap_genotypes = 50, 
                               Max_Hap_comb = 500)
  toc(log = TRUE, quiet = TRUE) 

  tic(msg = 'det_800_inflate') # Slowest inflate inc. complex
  Results_inflate_det = post_prob_CLI(MSdata = MS_inflated, 
                                      Fs = Fs_Combined, verbose = T,
                                      Max_Hap_genotypes = 100, 
                                      Max_Hap_comb = 800)
  toc(log = TRUE, quiet = TRUE) 

  save(Results_inflate_prob, Results_inflate_pdet, Results_inflate_det, file = '../RData/Results_prob_v_det_inflate_exc.NA.IDs.RData')
  rm(list = 'Results_inflate_prob', 'Results_inflate_pdet', 'Results_inflate_det')
} 

if(RUN){
  log.txt <- tic.log(format = TRUE)
  save(log.txt, file = '../RData/log.txt_prob_vs_det_exc.NA.IDs.RData')
} else {
  load(file = '../RData/log.txt_prob_vs_det_exc.NA.IDs.RData')
}




#============================================
# Now let's consider the timing results 
#============================================
log.txt 
# Both deterministic faster:
# Max_Hap_genotypes = 50 and Max_Hap_comb = 500 
# Max_Hap_genotypes = 100 and Max_Hap_comb = 800 
# Note, Max_Hap_comb = 500 and 800 indifferent since all converge: 

#============================================================================
# Now let's look if there are any P0s, thus potential impact of probablistic 
# phasing on recrudescence inference: no
#============================================================================
load('../RData/Results_prob_v_det_simple_exc.NA.IDs.RData')
load('../RData/Results_prob_v_det_inflate_exc.NA.IDs.RData')

# First check NAs: no
any(is.na(Results_det$Phased))
any(is.na(Results_pdet$Phased))
any(is.na(Results_prob$Phased))
any(is.na(Results_inflate_det$Phased))
any(is.na(Results_inflate_pdet$Phased))
any(is.na(Results_inflate_prob$Phased))

# No P0s, so no concern about recrudesence inference
table(Results_det$Phased) # All are D_D
table(Results_pdet$Phased) # All besides two are D_D 
table(Results_prob$Phased) # All are "P1_P1" 
table(Results_inflate_det$Phased) # All besides 17 D_D
table(Results_inflate_pdet$Phased) # All besides 25 are D_D 
table(Results_inflate_prob$Phased) # All are "P1_P1" 

