############################################################################
# This script defines a function that returns probablistic phasing combinations 
# compatible with the observed data
# 
# When a polyclonal infection is compatible with more than Max_Hap_genotypes, 
# we avoid using the deterministic approach, i.e. combinations() followed by scores, 
# where score = T for a combination of haploid genotypes (phasing) that is 
# compatible with the observed data.
#
# This is because the number of combinations is very large when the data are
# compatible with many haploid genotypes, even if sum(scores = T) is small
#
# Instead we permute and bootstrap the observed data to create combinations that 
# are guaranteed compatible with the data. Since this is a probablistic approach
# some of the combinations are duplicates. Rather than generating a fixed number of 
# combinations, we repeat until either the number found stabalises or exceeds
# some arbitarily high cut off, Max_Hap_comb. Without the cut off, the while loop 
# is liable to go on forever for very complex infections. 
#
# Ideally, the arbitary cute off should exceed the number of combinations compatible
# with infections whose haploid genotype count < Max_Hap_genotypes. Based on two simulated 
# infections whose number of total_haps_count = 74 and 96, the number of combinations, 
# nrow(Vt_Hnt_inds_comp) = 1296 and 7776, respectively (obtained setting Max_Hap_genotypes = 100)
# 
# Specifically, we permute saturated het. markers (markers with num. alleles obs. = COI)
# and bootstrap markers with reduncancy (markers with num. alleles observed < COI)
#
# Using the probablistic approach, we are likely to some over all possible phasings 
# if the number of compatable haploid genotypes is small, but with no guarantee. 
# The deterministic approach has the advantage of ensuring we sum over all 
# possible combinations of haplotypes compatible with the data, 
# but is not efficient since many combinations are not compatible. 
#
# Future iterations of the model with likely adopt the probablistic approach.
# Could also move deterministic process here so all phasing in one place.
############################################################################

Hap_combinations_probabilistic = function(Max_Hap_comb, cn = cn[inf], ynt, Y){

  names(Y) = colnames(ynt) # S.t. return combinations have names
  diff_unique = c(1,1,1) # Set a trio of non-zero differences s.t. while loop starts
  num_unique = 0 # Number of initial combinations 
  nrep = 1000 # We don't know how many are needed to capture most unique combinations, 
  # so...
  # either repeat until the last three nrep increases resulted in no change, or if number
  # unique exceeds Max_Hap_comb
  while(!all(diff_unique[1:3] %in% 0) & num_unique < Max_Hap_comb){ 
    
    nrep_comp = lapply(1:nrep, function(i){sapply(Y, function(as){
      num_obs = length(as) # Number of alleles at m-th marker
      if(num_obs > 1 & num_obs < cn){ # if het but not saturated, bootstrap...
        extra = sample(as, cn - num_obs, replace = T) 
        return(c(as, extra))}
      if(num_obs == cn){ # if het and saturated, permute
        return(sample(as, cn, replace = F))}
      if(num_obs == 1){ # Otherwise, return only value observed
        return(rep(as,cn))}
    })})
    
    # Second, remove duplicates (can take a while for a long list)
    unique_comp = unique(nrep_comp)
    
    # Third updated diff_unique with most recent at start
    diff_unique = c(length(unique_comp) - num_unique, diff_unique)
    if(diff_unique[1] < 0){next()} # If the current attempt is worse than before, try again
    num_unique = length(unique_comp) # update number that are unique
    nrep = nrep + 500 # increase number of comp to try
  }
  return(unique_comp)
}


