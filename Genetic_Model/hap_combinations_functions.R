############################################################################
# This script defines a function that returns probabilistic phasing combinations 
# compatible with the observed data (hap_combinations_probabilistic); a function
# that returns all possible phasing combinations compatible with the observed data 
# (hap_combinations_deterministic) and a function that returns NAs if the data are 
# missing (hap_combinations_missing_data)
# 
# hap_combinations_probabilistic: 
# When a polyclonal infection is compatible with more than Max_Hap_genotypes 
# (this cut off is not ideal - see Setting_Max_Hap_genotypes_combs.R), 
# we avoid using the deterministic approach, i.e. combinations() and scores, 
# where score = T for each combination (phasing) compatible with the observed data.
#
# This is because, in general, the number of combinations is very large when the data are
# compatible with many haploid genotypes (or COI is large), even if sum(scores = T) is small.
# However, no. of combinations also scales with cnt, s.t. Max_Hap_genotypes is NOT a very
# robust way to avoid deterministic phasing
# e.g., no. combinations when nrow(Hnt) = 648 choose cnt = 3 (45139896 rows)
# e.g., no. combinations when nrow(Hnt) = 288 choose cnt = 4 (280720440 rows)
# 
# Instead we permute and bootstrap the observed data to create combinations that 
# are guaranteed compatible with the data. 
# Specifically, we permute saturated het. markers (markers with num. alleles obs. = COI)
# and bootstrap markers with redundancy (markers with num. alleles observed < COI)
# Since this is a probablistic approach some of the combinations are duplicates. 
# Rather than generating a fixed number of combinations, we repeat until either the 
# number found stabilises or exceeds a cut-off, Max_Hap_comb. Without the cut-off, the while loop 
# is liable to loop forever for highly complex infections. 
#
# Using the probablistic approach, we are likely to recover all possible phasings 
# if the number of compatable haploid genotypes is small and the output "converge" is true, 
# but without guarantee. 
# The deterministic approach has the advantage of ensuring we sum over all 
# possible combinations of haplotypes compatible with the data, 
# but is not efficient since many combinations are not compatible. 
#
# Future more-efficient version of deterministic phasing would use dynamic programming
# Future iterations of the model with likely adopt the probablistic approach.
# Future iterations could also recursively phase conditional on inferred state
# (doing so would upweight inrelated parsites if L and thus reduce state space)
############################################################################
hap_combinations_probabilistic = function(Max_Hap_comb, cnt, ynt){
  
  # Summarise data for compatibility check and Hap_combinations_probabilistic()
  # alply ensures it's always as a list, important for Hap_combinations_probabilistic()
  Y = alply(ynt, 2, function(x){sort(unique(x[!is.na(x)]))}, .dims = T)
  Y = lapply(Y, as.character) # Needs to be character s.t. sample('9') not interpreted as sample(1:9)
  MSs = colnames(ynt)
  M = length(MSs)
  diff_unique = c(1,1,1) # Set a trio of non-zero differences s.t. while loop starts
  num_unique = 0 # Number of initial combinations 
  nrep = Max_Hap_comb # We don't know how many are needed to capture most unique combinations, 
  # set to Max_Hap_comb so that 
  # the unique_comp returned doesn't vastly exceed Max_Hap_comb
  # so...
  # either repeat until the last three nrep increases resulted in no change, 
  # or if number unique exceeds Max_Hap_comb
  while(!all(diff_unique[1:3] %in% 0) & num_unique < Max_Hap_comb){ 
    
    nrep_comp = lapply(1:nrep, function(i){
      
      z = sapply(Y, function(as){
        num_obs = length(as) # Number of alleles at m-th marker
        if(num_obs == cnt){ # if het and saturated
          return(sample(as, cnt, replace = F))} # permute order
        if(num_obs > 1 & num_obs < cnt){ # if het but not saturated
          extra = sample(as, cnt - num_obs, replace = T) # bootstrap
          inds = sample(1:cnt, cnt, replace = F) # permute order
          return(c(as,extra)[inds])}
        if(num_obs == 1){ # Otherwise, return only value observed
          return(rep(as,cnt))}
        if(num_obs == 0){ # Otherwise, return only value observed
          return(rep(NA,cnt))}
      })
      
      # Sort s.t. matrices do not differ by row order (important for unique step below)
      z = matrix(z, ncol = M, dimnames = list(NULL,MSs))
      inds_sorted = sort.int(apply(z, 1, paste, collapse= ""), index.return = T)$ix
      comp = z[inds_sorted, , drop = F]
      return(comp)
    })
    
    # Second, remove exact duplicates (can take a while for a long list)
    unique_comp = unique(nrep_comp)
    
    # Third updated diff_unique with most recent at start
    diff_unique = c(length(unique_comp) - num_unique, diff_unique)
    if(diff_unique[1] < 0){next()} # If the current attempt is worse than before, try again
    num_unique = length(unique_comp) # update number that are unique
    nrep = nrep + 0.5*Max_Hap_comb # increase (0.5*Max_Hap_comb s.t. unique_comp doesn't vastly exceed Max_Hap_comb)
    print(c(num_unique, nrep))
  }
  
  to_return = list(hap_combs = unique_comp, 
                   coverge = all(diff_unique[1:3] %in% 0)) # Return logical to know if diff_unique converged or not
  return(to_return)
}



hap_combinations_deterministic = function(Hnt, cnt, ynt){
  
  MSs = colnames(ynt)
  M = length(MSs)
  
  # Indices of all combinations of nrow(Hnt) choose cn[inf] haploid genotypes for the tth infection, inf 
  Vt_Hnt_inds = combinations(nrow(Hnt), r = cnt, v = 1:nrow(Hnt))  
  
  # Check each combination to see if compatible with ynt 
  # (this is a bit like an ABC step with epsilon = 0)
  # This is essentially a huge for loop and so could be sped up with Rcpp
  all_comp = alply(Vt_Hnt_inds, 1, function(x){ 
    comb = Hnt[x,,drop=FALSE]
    for(MS in MSs){
      if(setequal(comb[,MS], ynt[,MS])){
        score = T
      } else {
        score = F; break()
      }
    }
    if(score){
      Hnt_chr = matrix(sapply(comb, as.character), ncol = M)
      colnames(Hnt_chr) = MSs
      return(Hnt_chr)
    } else {
      return(NULL)
    }})
  
  to_return = all_comp[!sapply(all_comp, is.null)] 
  
  return(to_return)
}





hap_combinations_missing_data = function(ynt){
  
  MSs = colnames(ynt)
  M = length(MSs)
  
  # Convert to character since used to index 
  Hnt_chr = matrix(sapply(ynt, as.character), ncol = M)
  colnames(Hnt_chr) = MSs
  
  # Return data as a list
  only_comp = list(Hnt_chr)
  
  return(only_comp)
}

