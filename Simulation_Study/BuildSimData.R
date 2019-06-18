#' Make a matrix of simulated microsatellite marker data
#' 
#' @param Tn The total number of episodes (infections) for nth person. Tn has to be greater than 1
#' @paramCOIs The complexities of each infection C_1, ... , C_Tn (COI)
#' @param M The number of microsatellite markers
#' @param N The number of individuals
#' @param N_alleles The number of alleles per microsatellite 
#' (this is the same across all markers for simplicity)
#' @param relationship The type of relationship between episodes in the simulated data, 
#' this is either: 'clone': clonal (exact replicates); 'sibling': IBD of a half, 
#' e.g. probability of a half of being the copied with probability one per marker; 
#' 'stranger': new draw at random from background frequencies.
#' 
#' @return A matrix of dimensions N*C_total x M+4
#' The first M columns are for the M markers
#' The data for each individual is on C_total rows, withCOI[i] for the ith episode


BuildSimData = function(Tn = 2, # Currently we consider only 2 
                        COIs, M, N, N_alleles=10, 
                        relationship = c('Clone','Sibling','Stranger')){
  # Check the inputs
  if(length(COIs) != Tn & Tn > 1){
    stop('wrong inputs:COIs need to be of length Tn and Tn needs to be greater than 1')
  }
  if(!relationship %in% c('Clone','Sibling','Stranger')){ 
    stop('Please specify relationship from the list of: Clone, Sibling, Stranger')
  }
  
  # Make the vector of markers and their allele frequencies
  MS_markers = sapply(1:M, function(x) paste0('MS',x))
  FS = lapply(MS_markers, function(x) table(1:N_alleles)/N_alleles)
  names(FS) = MS_markers
  
  # Construct a dataframe for the simulated data of correct size
  C_total = sum(COIs)
  MS_data_sim = as.data.frame(array(dim = c(N * C_total, M+4)))
  colnames(MS_data_sim) = c('ID','MOI_id', MS_markers, 'Episode', 'Episode_Identifier')
  
  for(n in 1:N){
    
    # indices of infection 1
    inf_1_ind = ((n-1) * C_total + 1) : ((n-1) * C_total + COIs[1])
    
    # Generate marker data for infection 1
    for(ms in 1:M){MS_data_sim[inf_1_ind,MS_markers[ms]]=sample(1:N_alleles, size=COIs[1], replace = T)}
    
    # Add extra information needed for the Likelihood functions
    MS_data_sim$ID[inf_1_ind] = paste('SIM',n,sep = '_') # Name simulated individual
    MS_data_sim$Episode[inf_1_ind] = 1 # Name episode
    MS_data_sim$MOI_id[inf_1_ind] = 1:COIs[1] # Name COI row
    
    # Generate marker data for recurrent infections
    for(inf in 2:Tn){
      
      # indices of recurrent infections
      inf_Rec_ind = ((n-1)*C_total + sum(COIs[1:(inf-1)]) + 1) : ((n-1)*C_total + sum(COIs[1:(inf)]))
      
      # Add extra information needed for the Likelihood functions
      MS_data_sim$ID[inf_Rec_ind] = paste('SIM',n,sep = '_') # Name simulated individual
      MS_data_sim$Episode[inf_Rec_ind] = inf # Name episode
      MS_data_sim$MOI_id[inf_Rec_ind] = 1:COIs[inf] # Name COI row
      
      ###############################################
      # We copy one line from infection 1
      if(relationship == 'Clone'){
        
        # For the second infection we copy the first row of the first infection
        MS_data_sim[inf_Rec_ind[1], MS_markers] = MS_data_sim[inf_1_ind[1], MS_markers]
        
        if(COIs[inf]>1){
          # Generate alleles at random for infection 2
          for(ms in 1:M){
            # sample without replacement: remove x from the sampling (prevent collapse to clone)
            MS_data_sim[inf_Rec_ind[2:COIs[inf]],MS_markers[ms]] = sample(1:N_alleles, size = COIs[inf]-1, replace = T)
          }
        }
      }
      
      ###############################################
      # We copy with probability 0.5 from infection 1
      if(relationship == 'Sibling'){
        IBD = sample(c(T,F), size = M, replace = T)
        MS_data_sim[inf_Rec_ind[1], MS_markers[IBD]] = MS_data_sim[inf_1_ind[1], MS_markers[IBD]]
        MS_data_sim[inf_Rec_ind[1], MS_markers[!IBD]] = sample(1:N_alleles, size =  sum(!IBD), replace = T)
        
        if(COIs[inf]>1){
          # Generate alleles at random for infection 2
          for(ms in 1:M){
            MS_data_sim[inf_Rec_ind[2:COIs[inf]],MS_markers[ms]] = sample((1:N_alleles), size = COIs[inf]-1, replace = T)
          }
        }
      }
      
      ###############################################
      # Random data: same as for infection 1
      if(relationship == 'Stranger'){
        for(ms in 1:M){ # First we do the polyallelic markers
          MS_data_sim[inf_Rec_ind,MS_markers[ms]] = sample(1:N_alleles, size = COIs[inf], replace = T)
        }
      }
    }
  }
  MS_data_sim$Episode_Identifier = apply(MS_data_sim, 1, function(x) paste(x['ID'],x['Episode'],sep='_'))
  
  return(list(MS_data_sim=MS_data_sim, FS=FS))
}


# # Some examples of how to use this function
xsclone=BuildSimData(Tn = 3,COIs = c(3,3,1),M = 3,N = 10,
                     relationship = 'Clone',
                     N_alleles = 10)

xs_sib=BuildSimData(Tn = 3,COIs = c(3,3,1), M = 3, N = 10,
                    relationship = 'Sibling',
                    N_alleles = 10)

xs_str=BuildSimData(Tn = 3,COIs = c(3,3,1), M = 3, N = 10,
                    relationship = 'Stranger',
                    N_alleles = 10)
