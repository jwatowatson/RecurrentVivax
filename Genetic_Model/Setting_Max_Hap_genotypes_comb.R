##############################################################################################
# Initially, this script was written to choose reasonable limits for both Max_Hap_genotypes 
# and Max_Hap_comb based on the BPD and VHX data. Important because prob_post_CLI is slow
# when Max_Hap_comb is high as in some cases it needs to sum over:
# (a lot of graphs TIMES a lot of haplotypes). 
# It is now used to better understand interplay between Max_Hap_genotype and Max_Hap_comb 
# in combination with results generated in the following scripts: 
# Generate_analyse_probablistic_v_deterministic_inc.NA.IDs.R 
# Generate_analyse_probablistic_v_deterministic_exc.NA.IDs.R
#
# ********************** !!!!IMPORTANT!!!! ******************************
# 1) Not always true that no. compatible haploid genotypes > no. of compatible combinations
# 2) Mapping between no. compatible haploid genotypes and combinations is NOT one-to-one 
# 3) No. combinations does not necessarily scale with no. haploid genotypes compatible, e.g. 
# BPD_402_1: 16 genotypes comp with 196 compatible combinations
# VHX_501_1: 64 genotypes comp with 32 compatible combinations
# ***********************************************************************
#
# Re 3) it thus makes little sense to set Max_Hap_comb based on Max_Hap_genotypes, e.g.
# 97.5th percentile of no. of compatible haploid genotypes per episode BPD and VHX simple: 16 
# max no. of hap. genotype combinations given 97.5th percentile: 196 
# The 99th percentile, 50.24, is between 32 and 64, 
# for which max no. of hap. genotype combinations are 32 and 64 (less than for 16)
#
# Nevertheless, initial approach: set Max_Hap_genotypes to highest reasonable, 
# then set Max_Hap_comb to match highest in the set that pass Max_Hap_genotypes. 
# For Max_Hap_genotypes = 72, highest count of combinations = 1296 (both VHX_461_1)
# However, setting Max_hap_comb > 1000 is prohibitively slow when applied to many episodes. 
# If we keep Max_Hap_genotypes = 100, we will deterministically phase all with <= Max_Hap_genotypes 
# inc. VHX_461_1 with 72 comp genotypes and 1296 possible combinations. 
# However, from Generate_analyse_probablistic_v_deterministic_inc.NA.IDs.R, we learn that VHX_461
# is inc. in UpperComplexityIDs due to VHX_461_2. A.s. no reason to set Max_Hap_genotypes > 70. 
#
# If probablistic faster than deterministic, preferable to set 
# Max_Hap_genotypes to low in the main code, s.t. some that could be phased deterministically
# are phased probabilistically but with Max_Hap_comb set to sufficiently high
# to capture all that could be otherwise phased deterministically
# Doing so requires at least one high Max_Hap_genotypes run
# to show that deterministically phased is the same as probablistically phased. 
# This is shown in Generate_analyse_probablistic_v_deterministic_inc.NA.IDs.R 
# However, from the above mentioned scripts we find probablistic phasing is slower.
##############################################################################################
rm(list = ls())
load('../RData/RPackages_List.RData')
new.pkg <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(new.pkg) > 0){ # Install using .R script
  stop('Please run the script Install_and_load_required_packages.R before returning to this script. Thanks.')}else{
    sapply(pkgs, require, character.only = TRUE) # Load all packages
  }

load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # For MS_pooled
load('../RData/Data_for_relatedness.RData') # For Fs_Combined
source('../Genetic_Model/iGraph_functions.R')
source('../Genetic_Model/Likelihood_function.R')
source('../Genetic_Model/Data_functions.R')
source('../Genetic_Model/test_Rn_compatible.R')
source('../Genetic_Model/post_prob_CLI.R')
source('../Genetic_Model/hap_combinations_functions.R')
source('../Simulation_Study/BuildSimData.R')

# Remove MS data for which there are no recurrent data
N_episodes_typed = table(MS_pooled$ID[!duplicated(MS_pooled$Episode_Identifier)])
MSdata = filter(MS_pooled, ID %in% names(N_episodes_typed[N_episodes_typed>1]))
Episodes = unique(MSdata$Episode_Identifier)
Fs = Fs_Combined
Max_Eps = 3 # Limit on number of episodes (due to test_Rn_compatible)
Max_Tot_Vtx = 6 # Limit on number of vertices = cumulative COI
MSdata = reformat_MSdata(MSdata = MSdata, MSs=names(Fs)) # Reformat the data s.t. there are no NA gaps 
MSdata = arrange(MSdata, ID, Episode, MOI_id) # Make sure that the data are sorted correctly
MSs = names(Fs)
M = length(MSs)
IDs_all = as.character(unique(MSdata$ID)) # Character vector
yns = dlply(MSdata, 'ID') # Transform into a list
cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
sum_cns = sapply(cns, sum)
Tns = lapply(cns, length) # Number of episodes (also used in check below)
Tns_chr = lapply(Tns, as.character) # Character version needed for indexes in do.par
Tns = unlist(Tns)
all(sapply(yns, function(x) any(unique(x$Episode) == 1)))

#==========================================================================
# Number of haploid genotypes per person: this includes all episodes
#==========================================================================
n_haps_per_episode = lapply(yns, function(yn){ # For a given individual
  ynts = dlply(yn, 'Episode_Identifier') # Break into episodes
  sapply(ynts, function(ynt){ # For a given episode
    ynmt = ynt[,MSs,drop = FALSE] # Extract microsatellite data 
    ynmt_size = apply(ynmt, 2, function(x){length(unique(x))})  
    num_haps = prod(ynmt_size) # Calculate the number of haploid genotypes
  })
})


# Remove complex and no recurrence
complex = (Tns > Max_Eps) | (sum_cns > Max_Tot_Vtx) # 54 indiv. whose data are too complex 
no_recurrence = unlist(Tns) < 2 # individuals that have Tn = 1 (no recurrence)
IDs = IDs_all[(!complex & !no_recurrence)] 
N = length(IDs)

# Number of genotypes before and after removing complex
Num_with_wout = cbind("WithComplex" = tail(sort(unlist(n_haps_per_episode)),10), 
      "WOutComplex" = tail(sort(unlist(n_haps_per_episode[IDs])),10)) 
rownames(Num_with_wout) = NULL 
Num_with_wout

Limit = 72 # Limit (3rd hightest after removal of complex)


#==========================================================================
# Set limit on haploid genotypes
#==========================================================================
Max_Hap_genotypes = Limit # Set max to 95th percentil Limits based on script
number_comp = number_geno = array(dim = length(Episodes), dimnames = list(Episodes)) # Includes all


for(i in 1:N){
  
  id = IDs[i] 
  
  #==========================================================================
  # Extract data and processed data for the nth individual
  #==========================================================================
  yn = yns[[id]] 
  cn = cns[[id]]
  Tn = Tns[[id]]
  Tn_chr = Tns_chr[[id]]
  infections = unique(yn$Episode_Identifier) # IDs of infections for the nth individual
  
  #==========================================================================          
  # Generate all possible vertex haploid genotype label mappings of data onto a graph 
  #==========================================================================
  for(inf in infections){
    
    # Extract data for the tth infection
    ynt = filter(yn, Episode_Identifier == inf)[,MSs,drop = FALSE] 
    
    # Summarise data for compatibility check and Hap_combinations_probabilistic() 
    # alply ensures it's always as a list, important for Hap_combinations_probabilistic() 
    Y = alply(ynt, 2, function(x){sort(unique(x[!is.na(x)]))}, .dims = T) 
    
    # All haploid genotypes compatible with ynt (`unique` collapses repeats due to row per MOI)
    Hnt = expand.grid((lapply(ynt, unique)))
    total_haps_count = nrow(Hnt)
    number_geno[inf] = total_haps_count
    
    # If very many compatible haploid genotypes, adopt probablistic phasing approach
    if(total_haps_count > Max_Hap_genotypes){
      next()
    } else {
      hap_comps = hap_combinations_deterministic(Hnt, cnt = cn[inf], ynt)
      number_comp[inf] = length(hap_comps)
    }
  }
}


# Remove skipped due to total_haps_count > Max_Hap_genotypes
number_comp = number_comp[!is.na(number_comp)]
number_geno = number_geno[!is.na(number_geno)]

# Extract percentile before removing ID names
PROB = 0.99
Percentile = quantile(unlist(n_haps_per_episode[IDs]), probs = PROB)
names(n_haps_per_episode) = NULL # Remove id names (for next step)

# Check no. haploid genotypes compatible calculation gives same result as enumeration
all(number_geno == unlist(n_haps_per_episode)[names(number_geno)]) 

# Inspect relationship between no. haploid genotypes compatible and 
# no. of combinations of haploid genotypes compatible
cbind('No_geno' = number_geno[names(sort(number_comp))], 
      'No_geno_comb' = sort(number_comp))
plot(y = number_comp, x = number_geno[names(number_comp)])

# Report precentile
writeLines(sprintf('%sth percentile of no. of haploid genotypes per episode: %s', PROB*100, Percentile))
writeLines(sprintf('Max number of hap. genotype combinations given limit of %s: %s', Limit, max(number_comp)))



