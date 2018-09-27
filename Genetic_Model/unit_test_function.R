#================================================================================================
# Function to independently verify the code in compute_theta_reLapse() via example calculation
# Verification is not exhaustive: restricted to examples that are easily independently calculable 
# Ideally we need this to be a stand alone script using a test data set that is run each time
# we make a change or install R package.
#================================================================================================
pass_unit_test = function(MS_data_reformated, Fs, 
                          id_input = c('54', '103'), 
                          precision = .Machine$double.eps*100, 
                          verbose_unit_test = FALSE, 
                          p = 0.2, alpha = 0){ 
  
  source('../Probabilistic_Model_Functions/Likelihood_function.R')
  source('../Probabilistic_Model_Functions/iGraph_functions.R')
  
  #===========================================================================
  # Define functions for internal use
  #===========================================================================
  if(alpha == 0){
    p_yn_wm_str = function(fs_him, fs_hjm, alpha, I_hij){fs_him * fs_hjm}
  } else {
    p_yn_wm_str = function(fs_him, fs_hjm, alpha, I_hij){
      ifelse(I_hij[m], alpha*fs_him+(1-alpha)*fs_him*fs_hjm, (1-alpha)*fs_him*fs_hjm)}
  }
  p_yn_wm_sib = function(fs_him, fs_hjm, alpha, I_hij){
    ifelse(I_hij[m],(.5+alpha)*fs_him + (.5-alpha)*fs_him*fs_hjm, (.5-alpha)*fs_him*fs_hjm)}
  p_yn_wm_clo = function(fs_him, fs_hjm, alpha, I_hij){
    ifelse(I_hij[m],1*fs_him,0)}
  
  #===========================================================================
  # Summarise data 
  #===========================================================================
  writeLines('Please be aware, function assumes all microsatellite names are preceded by \"PV.\"')
  MSs = colnames(MS_data_reformated)[grepl('PV.', colnames(MS_data_reformated))]
  M = length(MSs) # Number of microsatellites
  IDs_all = as.character(unique(MS_data_reformated$ID)) # Character vector for indexing
  yns = lapply(IDs_all, function(x){yn = filter(MS_data_reformated, ID == x)})
  names(yns) = IDs_all
  cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
  sum_cns = sapply(cns, sum) # Size of Gn = sum(cnt) for check below
  Tns = lapply(cns, length) # Also used in check below
  Tns_chr = lapply(Tns, as.character) # Character version needed for indexes in do.par
  vtx_count_strs = sapply(cns[IDs_all], function(cn,  Max_Eps){
    vtx_count = c(cn, rep(0, Max_Eps-(length(cn)))) # Add 0 if Tn < Max_Eps
    vtx_count_str = paste(vtx_count, collapse = "_") # Create string for graph lookup
  },  Max_Eps = 3)
  
  #===========================================================================
  # Check example ids conform to 1_2_0 expectation
  #===========================================================================
  writeLines('Please be aware, function assumes input ID (if any) are for individuals who experience a monoclonal primary infection followed by a single recurrence with complexity of infection (COI) equal to two and two sets of viable vertex haplotype labels.')
  if(!all(sapply(id_input, function(id)all(cns[[id]] == c(1,2))))){stop('Assumption of input ID having primary infection with COI = 1 followed by a single recurrence with COI = 2 not met.')}
  
  #===========================================================================
  # Test 1: all 1_1_0 cases 
  #===========================================================================
  IDs_1_1_0 = names(which(vtx_count_strs == '1_1_0')) # Extract only on 1_1_0 cases
  MS_data_1_1_0 = MS_data_reformated[MS_data_reformated$ID %in% IDs_1_1_0, ] 
  Thetas_1_1_0 = array(dim = c(length(IDs_1_1_0), 2), dimnames = list(IDs_1_1_0, c('by_hand', 'by_model'))) # Create store
  G_1_1_0 = array(dim = c(length(IDs_1_1_0), 3), dimnames = list(IDs_1_1_0, c('G_str', 'G_sib', 'G_clo')))
    
  # Calculate thetas without using the model and thus without using the log-domain 
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
      P_yn_w12m_str[m] = p_yn_wm_str(fs_h1m, fs_h2m, alpha, I_h12)
      P_yn_w12m_sib[m] = p_yn_wm_sib(fs_h1m, fs_h2m, alpha, I_h12)
      P_yn_w12m_clo[m] = p_yn_wm_clo(fs_h1m, fs_h2m, alpha, I_h12)
    }
    
    # Three graphs (product over microsatellites)
    G_str = prod(P_yn_w12m_str, na.rm = TRUE)
    G_sib = prod(P_yn_w12m_sib, na.rm = TRUE)
    G_clo = prod(P_yn_w12m_clo, na.rm = TRUE)
    
    # Save for verbose return
    G_1_1_0[id, c('G_str', 'G_sib', 'G_clo')] = c(G_str, G_sib, G_clo)
    
    p_yn_0 = G_str # Probability of data given no relatedness between
    p_yn_1 = 0.5*(G_sib + G_clo) # Probability of data given some relatedness between
    p_yn = p_yn_0 * (1-p) + p_yn_1 * p # Probability of data
    p_1_yn = (p_yn_1 * p)/p_yn # Probability of some relatedness between given data 
    
    Thetas_1_1_0[id, 'by_hand'] = p_1_yn
  }
  
  # Calculate thetas under the model and unpackage
  Thetas_by_model_1_1_0 = compute_theta_reLapse(MS_data_reformated = MS_data_1_1_0, Fs = Fs, alpha = alpha, p = p, 
                                                cores = 4, Max_Eps = 3, Max_Tot_Vtx = 6,
                                                UpperComplexity = 50000, #2*10^6, # Assuming 10ms per operation -> 55 hours
                                                verbose = FALSE)
  X_1_1_0 = Thetas_by_model_1_1_0[!is.na(Thetas_by_model_1_1_0)]
  Thetas_1_1_0[do.call(rbind, strsplit(names(X_1_1_0), split = '_'))[,1], 'by_model'] = X_1_1_0 
  
  # Compare (don't expect them to be exactly the same due to log sum exp approximation)
  Diff_1_1_0 = Thetas_1_1_0[,'by_hand'] - Thetas_1_1_0[,'by_model']
  max_Diff_1_1_0 = max(abs(Diff_1_1_0)) 
  PASS_1_1_0 = max_Diff_1_1_0 < precision 
  
  #===========================================================================
  # Test 1_2_0 cases (Default IDs are 54 and 103)
  #===========================================================================
  IDs_1_2_0 = id_input
  MS_data_1_2_0 = MS_data_reformated[MS_data_reformated$ID %in% IDs_1_2_0, ] 
  Thetas_1_2_0 = array(dim = c(length(IDs_1_2_0), 2), dimnames = list(IDs_1_2_0, c('by_hand', 'by_model'))) # Create store
  G_1_2_0 = vector('list', length = length(IDs_1_2_0))
  names(G_1_2_0) = IDs_1_2_0
    
  for(id in IDs_1_2_0){
    
    yn = yns[[id]] 
    t0_ind = yn$Episode == 1
    t1_ind = yn$Episode == 2
    
    # Each row is a haplotype
    sorted_yn1 = sapply(yn[t1_ind, MSs], function(y)(unique(sort(y))))
    H_matrix = unique(expand.grid(yn[t1_ind, MSs]))
    Poss_comb = combinations(n = 4, r = 2)
    Viable_comb_ind = apply(Poss_comb, 1, function(x){
      sorted_h = sapply(H_matrix[x, ], function(y)(unique(sort(y))))
      z = identical(sorted_h, sorted_yn1)
      return(z)
    })
    Viable_comb = Poss_comb[Viable_comb_ind, ]
    
    if(sum(Viable_comb_ind)!=2){stop('No. of viable sets of vertex haplotype labels not equal to 2')}
    
    # Vertex haplotype labels
    h_w1 = yn[t0_ind, MSs]
    h_w2 = H_matrix[Viable_comb[,1],]
    h_w3 = H_matrix[Viable_comb[,2],]
    
    # Create store for pr_yn_Glw
    G = array(dim = c(2,9))
    
    # check there are no clones 
    if(any(all(h_w1 == h_w2[1,]), 
           all(h_w1 == h_w3[1,]), 
           all(h_w1 == h_w2[2,]), 
           all(h_w1 == h_w3[2,]))){stop('Examples 1,2,0 include clones across')}
    else { # Go on
      
      # Since no clones, pr_yn_Glw = 0 for w = 1,2,5,7, l = 1,2
      G[,c(1,2,5,7)] = 0
      
      #----------------------------------------------------------------
      # Calculate pr_yn_Glw for remaining 5 graphs
      #----------------------------------------------------------------
      # Create store for probability of ms data give edge 
      pr_y_wm = array(dim = c(3,2,2,M), dimnames = list(c('12','13','23'), c('str', 'sib'), c('l1', 'l2'), NULL))
      
      for(l in 1:2){
        
        I_h12 = h_w1[1,] == h_w2[l,]
        I_h13 = h_w1[1,] == h_w3[l,]
        I_h23 = h_w2[l,] == h_w3[l,]
        
        for(m in 1:M){
          
          # Extract frequencies
          fs_h1m = Fs[[MSs[m]]][as.character(h_w1[1, m])]
          fs_h2m = Fs[[MSs[m]]][as.character(h_w2[l, m])]
          fs_h3m = Fs[[MSs[m]]][as.character(h_w3[l, m])]
          
          # Calculate probability of ms data give edge (don't need clone as already ascertained above that no compatible)
          pr_y_wm['12','str',l,m] = p_yn_wm_str(fs_h1m, fs_h2m, alpha, I_h12)
          pr_y_wm['13','str',l,m] = p_yn_wm_str(fs_h1m, fs_h3m, alpha, I_h13)
          pr_y_wm['23','str',l,m] = p_yn_wm_str(fs_h2m, fs_h3m, alpha, I_h23)
          pr_y_wm['12','sib',l,m] = p_yn_wm_sib(fs_h1m, fs_h2m, alpha, I_h12)
          pr_y_wm['13','sib',l,m] = p_yn_wm_sib(fs_h1m, fs_h3m, alpha, I_h13)
          pr_y_wm['23','sib',l,m] = p_yn_wm_sib(fs_h2m, fs_h3m, alpha, I_h23)
        }
      }
      
      # Drop the marker dimension
      pr_y_w <- apply(pr_y_wm, c(1,2,3), prod, na.rm = T)
      
      # Create adjacency matrices for remaining 5 graphs
      A = rbind(c('12' = 'sib', '13' = 'sib', '23' = 'sib'), # A3
                c('12' = 'str', '13' = 'str', '23' = 'sib'), # A4
                c('12' = 'sib', '13' = 'str', '23' = 'str'), # A6 
                c('12' = 'str', '13' = 'sib', '23' = 'str'), # A8
                c('12' = 'str', '13' = 'str', '23' = 'str')) # A9
      
      G[1,c(3,4,6,8,9)] = apply(A, 1, function(x){prod(diag(pr_y_w[,x,'l1']), na.rm = T)})
      G[2,c(3,4,6,8,9)] = apply(A, 1, function(x){prod(diag(pr_y_w[,x,'l2']), na.rm = T)})
      
      G_1_2_0[[id]] = G # Store for verbose return
      
      Thetas_1_2_0[id, 'by_hand'] <- (1/7 * sum(G[,c(1:3,5:8)]) * p) / ((1/7 * sum(G[,c(1:3,5:8)]) * p) + (.5 * sum(G[,c(4,9)]) * (1-p)))
    }
  }
  
  # Calculate under the model and unpackage
  Thetas_by_model_1_2_0 = compute_theta_reLapse(MS_data_reformated = MS_data_1_2_0, Fs = Fs, alpha = alpha, p = p)
  X_1_2_0 = Thetas_by_model_1_2_0[!is.na(Thetas_by_model_1_2_0)]
  Thetas_1_2_0[do.call(rbind, strsplit(names(X_1_2_0), split = '_'))[,1], 'by_model'] = X_1_2_0 
  
  # Compare (don't expect them to be exactly the same due to log sum exp approximation)
  Diff_1_2_0 = Thetas_1_2_0[,'by_hand'] - Thetas_1_2_0[,'by_model']
  max_Diff_1_2_0 = max(abs(Diff_1_2_0)) 
  PASS_1_2_0 = max_Diff_1_2_0 < precision 
  
  #===========================================================================
  # End of tests
  #===========================================================================
  if(verbose_unit_test){
    to_return = list(PASS_1_1_0 = PASS_1_1_0,
                     max_Diff_1_1_0 == max_Diff_1_1_0, 
                     G_1_1_0 = G_1_1_0, 
                     Thetas_1_1_0 = Thetas_1_1_0, 
                     PASS_1_2_0 = PASS_1_2_0,
                     max_Diff_1_2_0 == max_Diff_1_2_0, 
                     G_1_2_0 = G_1_2_0, 
                     Thetas_1_2_0 = Thetas_1_2_0)
    return(to_return)
  } else {
    return(PASS_1_1_0 & PASS_1_2_0)
  }
}
