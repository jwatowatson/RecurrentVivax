################################################################################
# To-do list:
# Remove MSs, don't loop over NAs if all NAs, 
# Bug: returning massively hight thetas always
# Match all notation before release
# need to work out what an OK upper limit should be and consolidate duplicate complexity checks
################################################################################
library(doParallel)

compute_theta_reLapse = function(MS_data_reformated, Fs, 
                                 p = 0.2, alpha = 0, 
                                 cores = 4, Max_Eps = 3, 
                                 Max_Tot_Vtx = 6,
                                 UpperComplexity = 10^6,# Assuming 10ms per operation -> 55 hours
                                 verbose = FALSE){
  
  if(verbose) writeLines('Setting up parameters to do computation....')
  
  #==========================================================================
  # Create store Thetas with length = num. of all infections inc. Tn = 1 etc.
  #==========================================================================
  all_infections = unique(MS_data_reformated$Episode_Identifier)
  Thetas = array(dim = length(all_infections), dimnames = list(all_infections))
  
  #==========================================================================
  # Calculate the prior probability or reinfection, log probs. and log frequencies
  #==========================================================================
  MSs = names(Fs)  
  M = length(MSs) 
  q = 1 - p 
  log_p = log(p) 
  log_q = log(q) 
  log_Fs = lapply(Fs, log)
  
  #==========================================================================
  # Calculate alpha terms required for logSumExp in nested functions
  #==========================================================================
  alpha_terms = c(alpha = alpha, 
                  log_alpha = log(alpha), 
                  log_1_alpha = log(1-alpha), 
                  log_0.5alpha = log(0.5 + alpha),
                  log_0.5_alpha = log(0.5 - alpha))
  
  #==========================================================================
  # Extract yn for all individuals for checks ahead of do.par and Gn calculation
  # It is important that IDs_all is a character vector since id used to index throughout
  #==========================================================================
  IDs_all = as.character(unique(MS_data_reformated$ID)) # Character vector
  yns = lapply(IDs_all, function(x){
    yn = filter(MS_data_reformated, ID == x)})
  names(yns) = IDs_all
  
  #==========================================================================
  # Extract COIs (cns) and num. of infections (Tns) for each individual
  #==========================================================================
  cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
  sum_cns = sapply(cns, sum) # Size of Gn = sum(cnt) for check below
  Tns = lapply(cns, length) # Also used in check below
  Tns_chr = lapply(Tns, as.character) # Character version needed for indexes in do.par
  
  #==========================================================================
  # Check for individuals whose data are too complex for theta calculation
  #==========================================================================
  complex = (Tns > Max_Eps) | (sum_cns >Max_Tot_Vtx)
  if(any(complex)){
    writeLines(sprintf('\nSkipping ID: %s, which are too complex either because they have more than %s episodes or more than %s vertices, sorry.', 
                       paste(IDs_all[complex], collapse = ', '), Max_Eps,Max_Tot_Vtx))}
  
  #==========================================================================
  # Check for individuals that have Tn = 1 
  # (inadequate data for theta calculation)
  #==========================================================================
  no_recurrence = Tns < 2
  if(any(no_recurrence)){
    writeLines(sprintf('\nSkipping ID: %s, since no recurrent data.', 
                       paste(IDs_all[no_recurrence], collapse = ', ')))
  }
  
  #==========================================================================
  # Extract IDs with adequate data that are not too complex
  # (we only calculate theta for these IDs)
  #==========================================================================
  IDs = IDs_all[(!complex & !no_recurrence)]
  N = length(IDs) # N is the number of individuals for whom theta is calculable
  if(N == 0) {
    stop('No infections with calculable theta')
  } else {
    writeLines(sprintf('\nTotal number of IDs with calculable theta: %s', N))
  } 
  
  #==========================================================================
  # Generate all poss. Lns given unique (Tns) for IDs with calculable theta
  #==========================================================================
  unique_Tns = unique(Tns[IDs])
  Lns = lapply(unique_Tns, function(x){
    Ln = permutations(n = 2, r = x-1, v = c(0,1), repeats.allowed = TRUE)
    rownames(Ln) = apply(Ln, 1, paste, collapse = '') 
    return(Ln)
  })
  names(Lns) = unique_Tns
  
  #==========================================================================
  # Calculate log Pr( Ln | p )
  #==========================================================================
  log_Pr_Lns = lapply(Lns, function(Ln){
    apply(Ln, 1, function(x){sum(x * log_p + (1-x) * log_q)})})
  
  #==========================================================================
  # Count num. of vertices (COIs) for each individual's tth infection 
  #==========================================================================
  vtx_count_strs = sapply(cns[IDs], function(cn,  Max_Eps){
    vtx_count = c(cn, rep(0, Max_Eps-(length(cn)))) # Add 0 if Tn < Max_Eps
    vtx_count_str = paste(vtx_count, collapse = "_") # Create string for graph lookup
  },  Max_Eps)
  unique_vtx_count_str = unique(vtx_count_strs)
  
  #==========================================================================
  # Check all the requisite graph look-ups exist before going on
  #==========================================================================
  graph_lookup_check = sapply(unique_vtx_count_str, function(x){
    if(!file.exists(sprintf('../../RData/graph_lookup/graph_lookup_%s.Rdata', x))){
      stop(sprintf('graph_lookup %s does not exist', x))} else {
        'graph lookup exists...'
      }
  })
  
  #==========================================================================
  # For each unique_vtx_count_str generate matrix log_Pr_G_Ln 
  # where each log_Pr_G_Ln has rows per Gnw; cols per Lnt
  #==========================================================================
  log_Pr_G_Lns = lapply(unique_vtx_count_str, function(x, log_p, log_q){
    load(sprintf('../../RData/graph_lookup/graph_lookup_%s.Rdata', x)) # loads all Gnw for Gn
    G_Ln_comp = sapply(graph_lookup, test_Rn_compatible) # For each Gnw, test compatibility with Rn
    cn = as.numeric(strsplit(x, split = '_')[[1]]) # Reconstruct cn
    Tn_chr = as.character(length(cn)) # Reconstruct Tn_chr
    Ln = Lns[[Tn_chr]] # Extract Ln 
    log_Pr_Ln = log_Pr_Lns[[Tn_chr]] # Extract log Pr( Rn | p)
    log_Pr_G_Ln = log(t(G_Ln_comp*(1/rowSums(G_Ln_comp)))) # log Pr( Gnw | Rn ) = 1 / W in matrix
    log_Pr_G_Ln[is.infinite(log_Pr_G_Ln)] = NA # Set -Inf due to log(0) to NA
    return(log_Pr_G_Ln)
  }, log_p, log_q)
  names(log_Pr_G_Lns) = unique_vtx_count_str
  
  if(verbose){
    writeLines('\nIn the following viable graphs include only those with sufficient numbers of vertices 
               and appropriate vertex labels to fully represent each individual\'s data.')}
  
  
  # If cores=1 then this does normal sequential computation
  if(cores>1) registerDoParallel(cores = cores)

  #***********************************************
  # Computation of likelihood values
  #***********************************************
  Thetas = foreach(i=1:N, .combine = c) %dopar% {
    # set id = '70' for vtx_counts_str = "1_2_2" if checking by hand
    # set id = '402' for vtx_counts_str = "4_1_0" if checking by hand
    id = IDs[i] 
    
    #==========================================================================
    # Extract data and processed data for the nth individual
    #==========================================================================
    yn = yns[[id]] 
    cn = cns[[id]]
    sum_cn = sum(cn)
    Tn = Tns[[id]]
    Tn_chr = Tns_chr[[id]]
    infections = unique(yn$Episode_Identifier) # IDs of infections for the nth individual
    Ln = Lns[[Tn_chr, drop = FALSE]] # Extract relevant Ln
    colnames(Ln) = infections[-1] # Extract relevant Ln and name
    vtx_count_str = vtx_count_strs[id] # Extract vertex sizes of Gn
    log_Pr_Ln = log_Pr_Lns[[Tn_chr]] # Extract log P(Ln)
    log_Pr_G_Ln = log_Pr_G_Lns[[vtx_count_str]] # Extract log P(Gnw | Ln)
    
    #==========================================================================          
    # Generate all possible vertex haplotype label mappings of data onto a graph 
    #==========================================================================
    # Store for haplotypes and indeces of compatible combinations thereof
    Haplotypes_and_combinations = vector('list', length = 0)
    
    for(inf in infections){
      # Extract data for the tth infection
      ynt = filter(yn, Episode_Identifier == inf)[,MSs,drop = FALSE] 
      # All haplotypes compatible with ynt (unique collapses repeats due to row per MOI)
      Hnt = expand.grid((lapply(ynt, unique)))
      # Indices of all combinations of nrow(Hnt) choose cn[inf] haplotypes for the tth infection, inf 
      Vt_Hnt_inds = combinations(nrow(Hnt), r = cn[inf], v = 1:nrow(Hnt))       
      # Summarise data for compatiblility check below 
      Y = apply(ynt, 2, function(x){sort(unique(x[!is.na(x)]))}) 
      # Check each combination to see if compatible with ynt 
      # (this is a bit like an ABC step with epsilon = 0)
      scores = apply(Vt_Hnt_inds, 1, function(x, Y){ 
        X = apply(Hnt[x,,drop=FALSE], 2, function(x){sort(unique(x[!is.na(x)]))}) 
        score = identical(X,Y) # Test that they are the same
        return(score)
      }, Y)
      # Extract only those indices that are compatible 
      Vt_Hnt_inds_comp = Vt_Hnt_inds[scores,,drop = FALSE]
      # Convert to character since used to index 
      Hnt_chr = matrix(sapply(Hnt, as.character), ncol = M)
      colnames(Hnt_chr) = MSs
      # Return all possible haplotypes and indices of compatible combinations of haplotypes 
      Haplotypes_and_combinations[[inf]] = list(Hnt = Hnt_chr, 
                                                Vt_Hnt_inds_comp = Vt_Hnt_inds_comp)
    }
    
    
    # Extract the number of compatible combinations for each infection
    num_comp_combs_Vt = lapply(Haplotypes_and_combinations, function(x){1:nrow(x$Vt_Hnt_inds_comp)})
    # Create a matrix of inds to create all possible Gnl for a given Gn
    labelled_G_ind = as.matrix(expand.grid(num_comp_combs_Vt))
    # Total number of possible mappings, L (same as prod(sapply(num_comp_combs_Vt, length)))
    L = nrow(labelled_G_ind)
    log_L = log(L) # Needed for log domain calculation
    # From labelled_G_ind create all Gnl 
    labelled_Gs = vector('list', L) # store in list
    for(label_ind in 1:nrow(labelled_G_ind)){ # labelled_G_ind inc. compatible labelled graphs only
      # For each labbeling in labelled_G_ind extract vertex data matrix 
      x <- labelled_G_ind[label_ind, ,drop=FALSE]
      vertex_data_matrix = array(dim = c(sum_cn, M))
      for(t in 1:Tn){ # Since each infection has its own Hnt and Vt_Hnt_inds_comp need to loop over infections
        Z = Haplotypes_and_combinations[[infections[t]]] 
        # Extract haplotypes from Hnt in compatible combinations whose indices are listed in Vt_Hnt_inds_comp
        labels = as.matrix(Z$Hnt[Z$Vt_Hnt_inds_comp[x[t],,drop=FALSE],,drop=FALSE]) 
        vertex_data_matrix[cumsum(c(1,cn))[t]:cumsum(cn)[t], ] = labels 
      }
      labelled_Gs[[label_ind]] = vertex_data_matrix 
    }
    
    
    #==========================================================================          
    # Load all Gnw and check complexity
    #==========================================================================
    load(sprintf('../../RData/graph_lookup/graph_lookup_%s.Rdata', vtx_count_str))
    length_graph_lookup = length(graph_lookup)
    # James: added as.numeric to avoid getting integer overflow leading to an error
    complexity_problem = as.numeric(L)*as.numeric(length_graph_lookup)
    
    if(verbose){
      writeLines(sprintf('\nComplexity of problem, ID: %s\nThe number of different viable vertex labeled graphs is: %s\nThe number of different viable edge labeled graphs is: %s\nThe number of graphs in viable graph space is therefore: %s', 
                         id, L, length_graph_lookup, complexity_problem))
    }
    
    # We do a complexity check - need to work out what an OK upper limit should be
    # Aimee: Rather than having two complexity checks, I wonder if we might do this check
    # outside do.par/forloop (I will give it some thought)
    if(complexity_problem > UpperComplexity){ # Return NAs and skip to next ID
      
      writeLines(sprintf('\nSkipping this problem (ID %s), too complex.', id))
      thetas = array(NA, dim = Tn, dimnames = list(infections)) # Return NAs
      #Thetas[names(thetas)] = thetas # Need to change this for parallel computation
      
      # return NA vector of thetas
      thetas
      
    } else { # Go on to calculate theta
      
      #==========================================================================          
      # Calculate probabilities of all Gnw summed over all Gnlw and record run time
      #==========================================================================
      tic() # Start clock to calculate time take per graph
      
      # Calculate probabilities of Gnw summed over all Gnlw 
      log_Pr_yn_Gnws_unnormalised = sapply(graph_lookup, Log_Pr_yn_Gnw_unnormalised, 
                                           labelled_Gs=labelled_Gs, cn=cn, Tn=Tn, 
                                           log_Fs=log_Fs, MSs=MSs, alpha_terms=alpha_terms) 
      
      z = toc(quiet = TRUE) # Stop clock
      ms_per_graph_space = 1000*(z$toc-z$tic) # Total time taken
      ms_per_graph = ms_per_graph_space/complexity_problem # Time taken per graph
      
      if(verbose){writeLines(sprintf('Run time (ms) over all graphs in graph space: %s\nRun time (ms) per graph in graph space: %s', round(ms_per_graph_space), round(ms_per_graph)))}
      #complexity_time[i,] = c(complexity_problem, ms_per_graph) # Store time 
      log_Pr_yn_Gnws = log_Pr_yn_Gnws_unnormalised - log_L # Normalise probabilities
      
      #==========================================================================          
      # Calculate P(yn | Ln) by summing over all Gnw 
      #==========================================================================
      log_Pr_yn_Ln = array(dim = dim(log_Pr_G_Ln), dimnames = dimnames(log_Pr_G_Ln))
      for(col_i in 1:ncol(log_Pr_G_Ln)){
        log_Pr_yn_Ln[,col_i] = log_Pr_yn_Gnws + log_Pr_G_Ln[,col_i]
      }
      log_Pr_yn_Ln = apply(log_Pr_yn_Ln, 2, logSumExp, na.rm = TRUE)
      
      #==========================================================================          
      # Calculate P(Ln | yn) 
      #==========================================================================
      # Pr_yn_Ln[] = 1 # Likelihood check: returns the prior
      log_Pr_yn_and_Ln = log_Pr_yn_Ln[names(log_Pr_Ln)] + log_Pr_Ln
      log_Pr_yn = logSumExp(log_Pr_yn_and_Ln)
      Pr_Ln_yn = exp(log_Pr_yn_and_Ln - log_Pr_yn) 
      
      #==========================================================================          
      # Calculate P(Lnt | yn) and return 
      #==========================================================================
      if(Tn == 2){
        thetas =  c(NA, Pr_Ln_yn['1'])
      }
      if(Tn == 3){
        thetas = c(NA, sum(Pr_Ln_yn[c('11','10')]), sum(Pr_Ln_yn[c('11','01')]))
      }
      names(thetas) = infections
      
      # return vector
      thetas
    }
  }
    
  
  # # Store to check relationship between problem complexity and time
  # complexity_time = array(NA, c(N,2), dimnames = list(IDs, c('Number of graphs in viable graph space', 
  #                                                            'Run time (ms) per graph')))
  
  ## Uncomment for doParallel
  # Thetas[names(thetas)] = thetas
  # Aimee I'm saving this so plot can be generated outside of function
  # James: need to think how to incorporate this into the parallel version
  # save(complexity_time, file = '../../RData/complexity_time.RData')
  return(Thetas)
  
}


