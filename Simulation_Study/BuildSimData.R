#' Make a matrix of simulated microsatellite marker data
#' 
#' @param Tn The total number of episodes (infections) for nth person. Tn has to be greater than 1
#' @paramCOIs The complexities of each infection C_1, ... , C_Tn (COI)
#' @param M The number of microsatellite markers
#' @param N The number of individuals
#' @param N_alleles The number of alleles per microsatellite 
#' (this is the same across all markers for simplicity)
#' @param K_poly_markers The number of markers which are polyallelic 
#' with number of different alleles equal to theCOI_j
#' in the simulated data, this controls the complexity of the haplotype reconstruction 
#' problem
#' @param relatedness The type of relatedness between episodes in the simulated data, 
#' this is either: 'clone': clonal (exact replicates); 'sibling': IBD of a half, 
#' e.g. probability of a half of being the same per marker; 
#' 'stranger': new draw at random from background frequencies.
#' 
#' @return A matrix of dimensions N*C_total x M+4
#' The first M columns are for the M markers
#' The data for each individual is on C_total rows, withCOI[i] for the ith episode


BuildSimData = function(Tn,COIs, M, N, N_alleles=10, K_poly_markers,
                        relatedness = c('Clone','Sibling','Stranger')){
  # Check the inputs
  if(length(COIs) != Tn & Tn > 1){
    stop('wrong inputs:COIs need to be of length Tn and Tn needs to be greater than 1')
  }
  if(!relatedness %in% c('Clone','Sibling','Stranger')){ 
    stop('Please specify relatedness from the list of: Clone, Sibling, Stranger')
  }
  # Make the vector of markers
  MS_markers = sapply(1:M, function(x) paste0('MS',x))
  FS = lapply(MS_markers, function(x) table(1:N_alleles)/N_alleles)
  names(FS) = MS_markers
  # the number of polyallelelic markers: make into vector of length N if not already
  if(length(K_poly_markers) == 1) K_poly_markers = rep(K_poly_markers,N)
  
  # Construct a dataframe for the simulated data of correct size
  C_total = sum(COIs)
  MS_data_sim = as.data.frame(array(dim = c(N * C_total, M+4)))
  colnames(MS_data_sim) = c('ID','MOI_id', MS_markers, 
                            'Episode', 'Episode_Identifier')
  
  for(n in 1:N){
    # Number of polyallelic markers: this determines the complexity of problem
    K_poly_allele = K_poly_markers[n]
    
    # indices of infection 1
    inf_1_ind = ((n-1)*C_total + 1) : ((n-1)*C_total +COIs[1])
    
    # Generate marker data for infection 1
    for(ms in 1:K_poly_allele){ # First we do the polyallelic markers
      # sample without replacement
      MS_data_sim[inf_1_ind,MS_markers[ms]] = sample(1:N_alleles, size =COIs[1], replace = F)
    }
    if(K_poly_allele < M){
      for(ms in (K_poly_allele+1):M){ # and then the non-polyallelic markers
        MS_data_sim[inf_1_ind, MS_markers[ms]] = sample(1:N_alleles, size = 1)
      }
    }
    
    # Add extra information needed for the Likelihood functions
    MS_data_sim$ID[inf_1_ind] = paste('SIM',n,sep = '_')
    MS_data_sim$Episode[inf_1_ind] = 1
    MS_data_sim$MOI_id[inf_1_ind] = 1:COIs[1]
    
    # Generate marker data for recurrent infections
    for(inf in 2:Tn){
      # indices of recurrent infections
      inf_Rec_ind = ((n-1)*C_total + sum(COIs[1:(inf-1)]) + 1) : ((n-1)*C_total + sum(COIs[1:(inf)]))
      # Add extra information needed for the Likelihood functions
      MS_data_sim$ID[inf_Rec_ind] = paste('SIM',n,sep = '_')
      MS_data_sim$Episode[inf_Rec_ind] = inf
      MS_data_sim$MOI_id[inf_Rec_ind] = 1:COIs[inf]
      ###############################################
      # We copy one line from infection 1
      if(relatedness == 'Clone'){
        # For the second infection we copy the first row of the first infection
        MS_data_sim[inf_Rec_ind[1], MS_markers] = MS_data_sim[inf_1_ind[1], MS_markers]
        if(COIs[inf]>1){
          # Generate alleles at random for infection 2
          for(ms in 1:K_poly_allele){
            # sample without replacement: remove x from the sampling
            x = MS_data_sim[inf_Rec_ind[1],MS_markers[ms]]
            MS_data_sim[inf_Rec_ind[2:COIs[inf]],MS_markers[ms]] = sample((1:N_alleles)[-x], size = COIs[inf]-1, replace = F)
          }
          if(K_poly_allele<M){
            for(ms in (K_poly_allele+1):M){
              # just copy from previous
              x = MS_data_sim[inf_Rec_ind[1],MS_markers[ms]]
              MS_data_sim[inf_Rec_ind, MS_markers[ms]] = x
            }
          }
        }
      }
      ###############################################
      # We copy with probability 0.5 from infection 1
      if(relatedness == 'Sibling'){
        IBD = sample(c(T,F), size = M, replace = T)
        MS_data_sim[inf_Rec_ind[1], MS_markers[IBD]] = MS_data_sim[inf_1_ind[1], MS_markers[IBD]]
        MS_data_sim[inf_Rec_ind[1], MS_markers[!IBD]] = sample(1:N_alleles, size =  sum(!IBD), replace = T)
        if(COIs[inf]>1){
          # Generate alleles at random for infection 2
          for(ms in 1:K_poly_allele){
            # sample without replacement: remove x from the sampling
            x = MS_data_sim[inf_Rec_ind[1],MS_markers[ms]]
            MS_data_sim[inf_Rec_ind[2:COIs[inf]],MS_markers[ms]] = sample((1:N_alleles)[-x], size = COIs[inf]-1, replace = F)
          }
          if(K_poly_allele<M){
            for(ms in (K_poly_allele+1):M){
              # just copy from previous
              x = MS_data_sim[inf_Rec_ind[1],MS_markers[ms]]
              MS_data_sim[inf_Rec_ind, MS_markers[ms]] = x
            }
          }
        }
      }
      ###############################################
      # Random data: same as for infection 1
      if(relatedness == 'Stranger'){
        for(ms in 1:K_poly_allele){ # First we do the polyallelic markers
          # sample without replacement
          MS_data_sim[inf_Rec_ind,MS_markers[ms]] = sample(1:N_alleles, size = COIs[inf], replace = F)
        }
        if(K_poly_allele<M){
          for(ms in (K_poly_allele+1):M){ # and then the non-polyallelic markers
            MS_data_sim[inf_Rec_ind, MS_markers[ms]] = sample(1:N_alleles, size = 1)
          }
        }
      }
    }
  }
  MS_data_sim$Episode_Identifier = apply(MS_data_sim, 1, function(x) paste(x['ID'],x['Episode'],sep='_'))
  
  return(MS_data_sim)
}

# Some examples of how to use this function
# xsclone=BuildSimData(Tn = 3,COIs = c(3,3,1), M = 3, N = 10, 
#                 relatedness = 'Clone', 
#                 N_alleles = 10, K_poly_markers = 1)
# xs_sib=BuildSimData(Tn = 3,COIs = c(3,3,1), M = 3, N = 10, 
#                 relatedness = 'Sibling', 
#                 N_alleles = 10, K_poly_markers = 1)
# xs=BuildSimData(Tn = 3,COIs = c(3,3,1), M = 3, N = 10,
#                 relatedness = 'Stranger', N_alleles = 10, K_poly_markers = 1)
