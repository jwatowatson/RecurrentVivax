##############################################################################################
# This script was used to choose reasonable limits for Max_Hap_genotypes and Max_Hap_comb 
# based on the BPD and VHX data. Note that, the code takes a long time when Max_Hap_comb is 
# high as it needs to sum over a lot of graphs
#
# If probablistic is faster than deterministic, may be preferable to set 
# Max_Hap_genotypes to low, s.t. some that could be phased deterministically
# are phased probabilistically but with Max_Hap_comb set to sufficiently high
# to capture all that could be otherwise phased deterministically
#
# Doing so requires at least one high Max_Hap_genotypes run
# to show that deterministically phased is the same as probablistically phased. 
#
#
# Based on this script: 
# 97.5th percentile of no. of haploid genotypes per episode: 64
# 99th percentile of no. of haploid genotypes per episode: 128
# Max number of hap. genotype combinations given 97.5th percentile: 216
# Max number of hap. genotype combinations given 99th percentile: 1296
#
# However, Setting Max_hap_comb > 1000 can lead to UpperComplexity > 10^6
# Combinations is excessively slow when nrow(Hnt) = 288 choose cnt = 4 (280720440 rows)
#
# In clonclusion
# Let's set Max_Hap_genotypes = 100 for deterministic
# Let's set Max_Hap_comb = 800 always
##############################################################################################

rm(list = ls())
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # For MS_pooled
load('../RData/Data_for_relatedness.RData') # For Fs_Combined
source('../Genetic_Model/iGraph_functions.R')
source('../Genetic_Model/Likelihood_function.R')
source('../Genetic_Model/Data_functions.R')
source('../Genetic_Model/test_Rn_compatible.R')
source('../Genetic_Model/post_prob_CLI.R')
source('../Genetic_Model/hap_combinations.R')
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
yns = dlply(MSdata, 'ID') # transform into a list
cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
sum_cns = sapply(cns, sum)
Tns = lapply(cns, length) # Number of episodes (also used in check below)
Tns_chr = lapply(Tns, as.character) # Character version needed for indexes in do.par

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
names(n_haps_per_episode) = NULL # Remove IDs since episode ID sufficient 
n_haps_per_episode = unlist(n_haps_per_episode) # Convert to vector and then string
PROB = 0.975
Percentile = quantile(n_haps_per_episode, probs = PROB)


#==========================================================================
# Set limit on haploid genotypes
#==========================================================================
Max_Hap_genotypes = Percentile # Set max to 95th percentil Limits based on script
number_comp = number_geno = array(dim = length(Episodes), dimnames = list(Episodes)) # Includes all

# Loop over not too complex
complex = (Tns > Max_Eps) | (sum_cns > Max_Tot_Vtx) # 54 indiv. whose data are too complex 
no_recurrence = Tns < 2 # individuals that have Tn = 1 (no recurrence)
IDs = IDs_all[(!complex & !no_recurrence)] 
N = length(IDs)

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
      hap_comps = hap_combinations_deterministic(Hnt, cnt = cn[inf], ynt, Y)
      number_comp[inf] = length(hap_comps)
    }
  }
}


# Condense
number_comp = number_comp[!is.na(number_comp)]
number_geno = number_geno[!is.na(number_geno)]
all(number_geno ==  n_haps_per_episode[names(number_geno)]) # Check the same answer: yes

writeLines(sprintf('%sth percentile of no. of haploid genotypes per episode: %s', PROB*100, Percentile))
writeLines(sprintf('Max number of hap. genotype combinations given %sth percentile: %s', PROB*100, max(number_comp)))



