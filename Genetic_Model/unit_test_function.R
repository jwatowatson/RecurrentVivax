#================================================================================================
# Independently verify the code in compute_theta_reLapse() via example calculation
# Verification is not exhaustive: restricted to examples that are easily independently calculable 
# Ideally, we would make this into a unit test that is run each time we make a change or 
# install R package.
# Could be adapted for use with sim data also
#================================================================================================

rm(list = ls())
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # For MS_pooled
load('../RData/Data_for_relatedness.RData') # For Fs_Combined
source('../Genetic_Model/iGraph_functions.R')
source('../Genetic_Model/Likelihood_function.R')
source('../Genetic_Model/Data_functions.R')
source('../Genetic_Model/test_Rn_compatible.R')
source('../Genetic_Model/post_prob_CLI.R')
source('../Genetic_Model/hap_combinations.R')
source('BuildSimData.R')


#================================================================================================
# Function to calculate cases where two monoclonal episodes
#================================================================================================
by_hand = function(MS_data_reformated, Fs, alpha = 0, Max_Eps = 3){ 
  
  #===========================================================================
  # Define functions for internal use
  #===========================================================================
  if(alpha == 0){
    p_ym_hm_str = function(fs_him, fs_hjm, alpha, I_hij){fs_him * fs_hjm}
  } else {
    p_ym_hm_str = function(fs_him, fs_hjm, alpha, I_hij){
      ifelse(I_hij[m], alpha*fs_him+(1-alpha)*fs_him*fs_hjm, (1-alpha)*fs_him*fs_hjm)}
  }
  p_ym_hm_sib = function(fs_him, fs_hjm, alpha, I_hij){
    ifelse(I_hij[m],(.5+alpha)*fs_him + (.5-alpha)*fs_him*fs_hjm, (.5-alpha)*fs_him*fs_hjm)}
  p_ym_hm_clo = function(fs_him, fs_hjm, alpha, I_hij){
    ifelse(I_hij[m],1*fs_him,0)}
  
  #===========================================================================
  # Summarise data 
  #===========================================================================
  writeLines('Please be aware, function assumes all microsatellite names are preceded by \"PV.\"')
  MSs = colnames(MS_data_reformated)[grepl('PV.', colnames(MS_data_reformated))] # Extract microsatellites
  M = length(MSs) # Number of microsatellites
  IDs_all = as.character(unique(MS_data_reformated$ID)) # Character vector for indexing
  yns = lapply(IDs_all, function(x){yn = filter(MS_data_reformated, ID == x)})
  names(yns) = IDs_all
  cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
  IDs = names(which(sapply(cns, length) <= Max_Eps))
  
  sum_cns = sapply(cns[IDs], sum) # Size of Gn = sum(cnt) for check below
  Tns = lapply(cns, length) # Also used in check below
  Tns_chr = lapply(Tns, as.character) # Character version needed for indexes in do.par
  vtx_count_strs = sapply(cns[IDs], function(cn,  Max_Eps){
    vtx_count = c(cn, rep(0, Max_Eps-(length(cn)))) # Add 0 if Tn < Max_Eps
    vtx_count_str = paste(vtx_count, collapse = "_") # Create string for graph lookup
  },  Max_Eps = 3) 
  
  #===========================================================================
  # Check example ids conform to 1_2_0 expectation
  #===========================================================================
  writeLines('Please be aware, function assumes input ID (if any) are for individuals \n
             who experience a monoclonal primary infection followed by a single recurrence \n
             with complexity of infection (COI) equal to two and two sets of viable vertex haplotype labels.')
  
  #===========================================================================
  # Test 1: all 1_1_0 cases 
  #===========================================================================
  IDs_1_1_0 = names(which(vtx_count_strs == '1_1_0')) # Extract only on 1_1_0 cases
  MS_data_1_1_0 = MS_data_reformated[MS_data_reformated$ID %in% IDs_1_1_0, ] 
  Thetas_1_1_0 = array(dim = c(length(IDs_1_1_0), 2), dimnames = list(IDs_1_1_0, c('by_hand', 'by_model'))) # Create store
  G_1_1_0 = array(dim = c(length(IDs_1_1_0), 3), dimnames = list(IDs_1_1_0, c('G_str', 'G_sib', 'G_clo')))
  
  # Calculate thetas without using the model nor log-domain 
  for(id in IDs_1_1_0){
    
    yn = yns[[id]] # Extract data
    h1 = yn[yn$Episode == 1, MSs] # Extract haplotype on first vertex
    h2 = yn[yn$Episode == 2, MSs] # Extract haplotype on second vertex   
    I_h12 = h1 == h2 # Indicator 
    
    # Create store for probability of ms data give edge
    P_yn_w12m_str = rep(NA, M)
    P_yn_w12m_sib = rep(NA, M)
    P_yn_w12m_clo = rep(NA, M)
    
    for(m in 1:M){
      # Extract frequencies
      fs_h1m = Fs[[MSs[m]]][as.character(h1[MSs[m]])]
      fs_h2m = Fs[[MSs[m]]][as.character(h2[MSs[m]])]
      
      # Calculate probability of ms data give edge
      P_yn_w12m_str[m] = p_ym_hm_str(fs_h1m, fs_h2m, alpha, I_h12)
      P_yn_w12m_sib[m] = p_ym_hm_sib(fs_h1m, fs_h2m, alpha, I_h12)
      P_yn_w12m_clo[m] = p_ym_hm_clo(fs_h1m, fs_h2m, alpha, I_h12)
    }
    
    # Three graphs (product over microsatellites)
    G_str = prod(P_yn_w12m_str, na.rm = TRUE)
    G_sib = prod(P_yn_w12m_sib, na.rm = TRUE)
    G_clo = prod(P_yn_w12m_clo, na.rm = TRUE)
    
    # Save for verbose return
    G_1_1_0[id, c('G_str', 'G_sib', 'G_clo')] = c(G_str, G_sib, G_clo)
    
    p_yn_I = G_str # Probability of data given reinfection
    p_yn_C = G_clo # Probability of data given recrudesence 
    p_yn_L = (1/3) * (G_clo + G_str + G_sib)
    p_yn =  (1/3) * (p_yn_I + p_yn_C + p_yn_L)
    P_L_yn = (p_yn_L/3)/p_yn
    
    Thetas_1_1_0[id, 'by_hand'] = P_L_yn
    return(Thetas_1_1_0)
  }
}


load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # For MS_pooled
load('../RData/Data_for_relatedness.RData') # For Fs_Combined
N_episodes_typed = table(MS_pooled$ID[!duplicated(MS_pooled$Episode_Identifier)])
MSdata = filter(MS_pooled, ID %in% names(N_episodes_typed[N_episodes_typed>1]))
by_hand(MS_data_reformated = MSdata, Fs = Fs_Combined, alpha = 0)


# Calculate thetas under the model and unpackage
Thetas_by_model_1_1_0 = post_prob_CLI(MS_data_reformated = MS_data_1_1_0, Fs = Fs, alpha = alpha)







# Remove MS data for which there are no recurrent data
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

which(sum_cns) == 2 # These are the ones to keep
