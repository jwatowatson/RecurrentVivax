############################################################################
# This script defines a function that returns probabilistic phasing combinations 
# compatible with the observed data (hap_combinations_probabilistic); a function
# that returns all possible phasing combinations compatible with the observed data 
# (hap_combinations_deterministic) and a function that returns NAs if the data are 
# missing (hap_combinations_missing_data)
# 
# hap_combinations_probabilistic: 
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
# combinations, we repeat until either the number found stabilises or exceeds
# some arbitarily high cut-off, Max_Hap_comb. Without the cut-off, the while loop 
# could be close to infinite for highly complex infections. 
#
# Ideally, the arbitary cut-off, Max_Hap_comb, should exceed the number of combinations 
# compatible with infections whose haploid genotype count < Max_Hap_genotypes. We choose 
# Max_Hap_genotypes and Max_Hap_comb in Setting_Max_Hap_genotypes_comb
# 
# Specifically, we permute saturated het. markers (markers with num. alleles obs. = COI)
# and bootstrap markers with redundancy (markers with num. alleles observed < COI)
#
# Using the probablistic approach, we are likely to some over all possible phasings 
# if the number of compatable haploid genotypes is small, but with no guarantee. 
# The deterministic approach has the advantage of ensuring we sum over all 
# possible combinations of haplotypes compatible with the data, 
# but is not efficient since many combinations are not compatible. 
#
# Future iterations of the model with likely adopt the probablistic approach.
############################################################################


### Why does this function not take as argument MSs??? The others do (e.g. deterministic one)
hap_combinations_probabilistic = function(Max_Hap_comb, cnt, ynt, Y){
  
  names(Y) = colnames(ynt) # S.t. return combinations have names
  diff_unique = c(1,1,1) # Set a trio of non-zero differences s.t. while loop starts
  num_unique = 0 # Number of initial combinations 
  nrep = Max_Hap_comb # We don't know how many are needed to capture most unique combinations, 
  # set to Max_Hap_comb so that 
  # the unique_comp returned doesn't vastly exceed Max_Hap_comb
  # so...
  # either repeat until the last three nrep increases resulted in no change, 
  # or if number unique exceeds Max_Hap_comb
  while(!all(diff_unique[1:3] %in% 0) & num_unique < Max_Hap_comb){ 
    
    nrep_comp = lapply(1:nrep, function(i){sapply(Y, function(as){
      num_obs = length(as) # Number of alleles at m-th marker
      if(num_obs > 1 & num_obs < cnt){ # if het but not saturated, bootstrap...
        extra = sample(as, cnt - num_obs, replace = T) 
        return(c(as, extra))}
      if(num_obs == cnt){ # if het and saturated, permute
        return(sample(as, cnt, replace = F))}
      if(num_obs == 1){ # Otherwise, return only value observed
        return(rep(as,cnt))}
      if(num_obs == 0){ # Otherwise, return only value observed
        return(rep(NA,cnt))}
    })})
    
    # Second, remove duplicates (can take a while for a long list)
    unique_comp = unique(nrep_comp)
    
    # Third updated diff_unique with most recent at start
    diff_unique = c(length(unique_comp) - num_unique, diff_unique)
    if(diff_unique[1] < 0){next()} # If the current attempt is worse than before, try again
    num_unique = length(unique_comp) # update number that are unique
    nrep = nrep + 0.5*Max_Hap_comb # increase (0.5*Max_Hap_comb s.t. unique_comp doesn't vastly exceed Max_Hap_comb)
    print(c(num_unique, nrep))
  }
  return(unique_comp)
}


hap_combinations_deterministic = function(Hnt, cnt, ynt, Y, MSs, M){
  
  # Indices of all combinations of nrow(Hnt) choose cn[inf] haploid genotypes for the tth infection, inf 
  Vt_Hnt_inds = combinations(nrow(Hnt), r = cnt, v = 1:nrow(Hnt))  
  
  # Check each combination to see if compatible with ynt 
  # (this is a bit like an ABC step with epsilon = 0)
  scores = apply(Vt_Hnt_inds, 1, function(x, Y){ 
    X = alply(Hnt[x,,drop=FALSE], 2, function(x){as.character(sort(unique(x[!is.na(x)])))})
    score = identical(X,Y) # Test that they are the same
    return(score)
  }, Y)
  
  # Extract only those indices that are compatible 
  Vt_Hnt_inds_comp = Vt_Hnt_inds[scores,,drop = FALSE]
  
  # Convert to character since used to index 
  Hnt_chr = matrix(sapply(Hnt, as.character), ncol = M)
  colnames(Hnt_chr) = MSs
  
  # Return all possible compatible combinations of haploid genotypes as a list (apply returns array)
  all_comp = lapply(apply(Vt_Hnt_inds_comp,1,function(i){list(Hnt_chr[i,,drop = F])}), function(x) x[[1]])
  
  return(all_comp)
}


hap_combinations_missing_data = function(ynt, MSs, M){
  
  # Convert to character since used to index 
  Hnt_chr = matrix(sapply(ynt, as.character), ncol = M)
  colnames(Hnt_chr) = MSs
  
  # Return data as a list
  only_comp = list(Hnt_chr)
  
  return(only_comp)
}

