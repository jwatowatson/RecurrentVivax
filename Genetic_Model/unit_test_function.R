#================================================================================================
# Independently verify the code in compute_theta_reLapse() via example calculation
# Verification is not exhaustive: restricted to examples that are easily independently calculable 
# Ideally, we would make this into a unit test that is run each time we make a change or 
# install R package.
# Could be adapted for use with sim data also
#================================================================================================

rm(list = ls())
load('../RData/RPackages_List.RData')
# new.pkg <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
# if(length(new.pkg) > 0){ # Install using .R script
#   stop('Please run the script Install_and_load_required_packages.R before returning to this script. Thanks.')}else{
#     sapply(pkgs, require, character.only = TRUE) # Load all packages
#   }
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # For MS_pooled
load('../RData/Data_for_relatedness.RData') # For Fs_Combined
source('./iGraph_functions.R')
source('./Likelihood_function.R')
source('./Data_functions.R')
source('./test_Rn_compatible.R')
source('./post_prob_CLI.R')
source('./hap_combinations.R')


#================================================================================================
# Function to calculate cases where two monoclonal episodes
#================================================================================================
by_hand = function(MSdata, Fs, alpha = 0, Max_Eps = 3){ 
  

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
  MSdata = reformat_MSdata(MSdata = MSdata, MSs=names(Fs)) # Reformat the data s.t. there are no NA gaps 
  MSdata = arrange(MSdata, ID, Episode, MOI_id) # Make sure that the data are sorted correctly
  MSs = names(Fs)
  M = length(MSs) # Number of microsatellites
  IDs = as.character(unique(MSdata$ID)) # Character vector for indexing
  yns = lapply(IDs, function(x){yn = filter(MSdata, ID == x)})
  names(yns) = IDs
  cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})

  vtx_count_strs = sapply(cns[IDs], function(cn,  Max_Eps){
    vtx_count = c(cn, rep(0, Max_Eps-(length(cn)))) # Add 0 if Tn < Max_Eps
    vtx_count_str = paste(vtx_count, collapse = "_") # Create string for graph lookup
  },  Max_Eps = 3) 
  
  #===========================================================================
  # Test 1: all 1_1_0 cases 
  #===========================================================================
  Thetas_1_1_0 = array(dim = length(IDs), dimnames = list(IDs)) # Create store
  G_1_1_0 = array(dim = c(length(IDs), 3), dimnames = list(IDs, c('G_str', 'G_sib', 'G_clo')))
  
  # Calculate thetas without using the model nor log-domain 
  for(id in IDs){
    
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
    p_yn = (1/3) * (p_yn_I + p_yn_C + p_yn_L)
    P_L_yn = (p_yn_L/3)/p_yn
    
    Thetas_1_1_0[id] = P_L_yn
  }
  return(Thetas_1_1_0)
}

# Recover allele frequencies
load('../RData/Data_for_relatedness.RData') # For Fs_Combined
Fs = Fs_Combined



load('../RData/GeneticModel/MS_data_PooledAnalysis.RData') # For MS_pooled
yns = dlply(MS_pooled, 'ID') # transform into a list
cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
cns_chr = sapply(cns, paste, collapse = "_")
MSdata_1_1_0 = filter(MS_pooled, ID %in% names(which(cns_chr == "1_1"))) # Keep individuals with two monoclonal episodes only

by_hand(MSdata = MSdata_1_1_0, Fs = Fs_Combined, alpha = 0)

post_prob_CLI(MSdata = MSdata_1_1_0, Fs = Fs, alpha = 0)[,'L']











