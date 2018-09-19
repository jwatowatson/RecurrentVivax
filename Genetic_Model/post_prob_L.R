# Notation matches manuscript aside from Tn, which counts from 1 here and 0 in the notation
# Variables names: conditioning bar is represented by an underscore while 'and' is explicit
# e.g. log_Pr_yn_Rn denotes log Pr(yn | Rn) while log_Pr_yn_and_Rn denotes log Pr(yn, Rn)

# # Can we collapse the following checks? 
# Max_Eps = 3, 
# Max_Tot_Vtx = 6,
# UpperComplexity = 10^6

post_prob_L = function(MS_data, # MS data (assumes no NA gaps in mixed infections)
                       Fs, # MS population frequencies 
                       p = c('C' = 0.01, 'L' = 0.2, 'I' = 0.79), # Population constant prior over C, L, I
                       alpha = 0, # Additative inbreeding constant
                       cores = 4, 
                       Max_Eps = 3, 
                       Max_Tot_Vtx = 6,
                       UpperComplexity = 10^6, # Assuming 10ms per operation -> 55 hours
                       verbose = FALSE){
  

  if(verbose) writeLines('Setting up parameters to do computation....')
  
  #==========================================================================
  # Retrieve population prior and recurrent eps identifiers
  #==========================================================================
  p_pop = sapply(c('C','L','I'), function(x)as.list(formals(post_prob_L)$p)[[x]])
  recurrent_eps_ind = as.numeric(do.call(rbind, strsplit(MS_data$Episode_Identifier, split = '_'))[,3]) > 1
  recurrent_eps = unique(MS_data$Episode_Identifier[recurrent_eps_ind])
  
  #==========================================================================
  # Check to see if using pop prior or prior from time-to-event and comment
  #==========================================================================
  if(is.null(dim(p)) & length(p) == 3){
    p_pop_ind = TRUE
    writeLines('Using population prior probabilities of recurrence states')
    if(sum(p) != 1) stop('Population prior probabilities do not sum to one')
  } else {
    p_pop_ind = FALSE
    writeLines('Using time-to-event prior probabilities of recurrence states')
    recurrent_eps_no_prior = recurrent_eps[which(!recurrent_eps %in% p$Episode_Identifier)]
    no_recurrent_eps_no_prior = length(recurrent_eps_no_prior)
    if(no_recurrent_eps_no_prior == 0){
      writeLines(sprintf('Using time-to-event prior probabilities for all recurrent infections', recurrent_eps_no_prior))
    } else {
      writeLines(sprintf('Using time-to-event prior probabilities for all but the following for which population prior probabilities are used: %s', 
                         paste(recurrent_eps_no_prior, collapse = ', ')))
    }
    for(x in recurrent_eps_no_prior){p = rbind(c(x, p_pop), p)} # Add to episodes based on the prior to the top of the list
  }
  
  #==========================================================================
  # Create Post_probs store with length = num. of all infections inc. those without recurrence
  #==========================================================================
  all_infections = unique(MS_data$Episode_Identifier)
  Post_probs = array(dim = length(all_infections), dimnames = list(all_infections))
  
  #==========================================================================
  # Extract microsatellites, related quantities and calculate log pop prior
  #==========================================================================
  MSs = names(Fs)  
  M = length(MSs) 
  log_Fs = lapply(Fs, log)
  if(p_pop_ind){
    log_p = log(p)
  } else {
    log_p = t(apply(p, 1, function(x)log(as.numeric(x[-1]))))
    rownames(log_p) = p[,1]
    colnames(log_p) = colnames(p)[-1]
  }
  
  #==========================================================================
  # Calculate alpha terms required for logSumExp in nested functions
  #==========================================================================
  alpha_terms = c(alpha = alpha, 
                  log_alpha = log(alpha), 
                  log_1_alpha = log(1-alpha), 
                  log_0.5alpha = log(0.5 + alpha),
                  log_0.5_alpha = log(0.5 - alpha))
  
  #==========================================================================
  # Extract yn for all individuals for checks ahead of do.par and Gn_ab calculation
  # It is important that IDs_all is a character vector since id used to index
  #==========================================================================
  IDs_all = as.character(unique(MS_data$ID)) # Character vector
  yns = lapply(IDs_all, function(x){yn = filter(MS_data, ID == x)})
  names(yns) = IDs_all
  
  #==========================================================================
  # Extract COIs (cns) and num. of infections (Tns) for each individual
  # Note that, in contrast to notation, here Tn counts from 1 
  #==========================================================================
  cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
  sum_cns = sapply(cns, sum) # Size of graph = sum(cnt) for check below
  Tns = lapply(cns, length) # Also used in check below
  Tns_chr = lapply(Tns, as.character) # Character version needed for indexes in do.par
  
  #==========================================================================
  # Check for individuals whose data are too complex for theta calculation
  #++++++++++++++++++++++++++++
  # May need to change this given added complexity of additional state
  #++++++++++++++++++++++++++++
  #==========================================================================
  complex = (Tns > Max_Eps) | (sum_cns >Max_Tot_Vtx)
  if(any(complex)){
    writeLines(sprintf('\nSkipping ID: %s, which are too complex either because they have more than %s episodes or more than %s vertices, sorry.', 
                       paste(IDs_all[complex], collapse = ', '), Max_Eps,Max_Tot_Vtx))}
  
  #==========================================================================
  # Check for individuals that have Tn = 1 (no recurrence)
  # (inadequate data for posterior probability calculation)
  #==========================================================================
  no_recurrence = Tns < 2
  if(any(no_recurrence)){
    writeLines(sprintf('\nSkipping ID: %s, since no recurrent data.', 
                       paste(IDs_all[no_recurrence], collapse = ', ')))
  }
  
  #==========================================================================
  # Extract IDs with adequate data that are not too complex
  # (we only calculate posterior probabilities for these IDs)
  #==========================================================================
  IDs = IDs_all[(!complex & !no_recurrence)]
  N = length(IDs) # N is the number of individuals for whom post. probs. are calculable
  if(N == 0) {stop('No infections with calculable posterior probabilities')
  } else {
    writeLines(sprintf('\nTotal number of IDs with calculable posterior probabilities: %s', N))
  } 
  
  #==========================================================================
  # Generate all poss. Rns given unique (Tns) for IDs with calculable posteriors
  #==========================================================================
  unique_Tns = unique(Tns[IDs])
  Rns = lapply(unique_Tns, function(x){
    Rn = permutations(n = 3, r = x-1, v = c('C','L','I'), repeats.allowed = TRUE)
    rownames(Rn) = apply(Rn, 1, paste, collapse = '') 
    return(Rn)
  })
  names(Rns) = unique_Tns
  
  #==========================================================================
  # Calculate log Pr(Rn | p)
  # If p_pop_ind = T: log_Pr_Rns is a list of length unique(Tn) with an entry by Rn^(t)
  # If p_pop_ind = F: log_Pr_Rns list of length N
  #==========================================================================
  if(p_pop_ind){ 
    log_Pr_Rns = lapply(Rns, function(Rn){apply(Rn, 1, function(x){sum(log_p[x])})})
  } else {
    log_p_IDs = lapply(IDs, function(x){log_p[grepl(paste(x,'_', sep = ''), rownames(log_p)),]})
    names(log_p_IDs) = IDs
    log_Pr_Rns = log_p_IDs # Set log_p equal to log Pr(Rn | p)
    for(x in names(which(unlist(Tns) > 2))){ # For those with more than recurrence sum over recurrences 
      log_Pr_Rns[[x]] = apply(Rns[[Tns_chr[[x]]]], 1, function(z) sum(diag(log_p_IDs[[x]][,z])))
    }
  }
  
  #==========================================================================
  # Count num. of vertices (COIs) for each individual's tth infection 
  # (save as a vector of strings as used to index)
  #==========================================================================
  vtx_count_strs = sapply(cns[IDs], function(cn,  Max_Eps){
    vtx_count = c(cn, rep(0, Max_Eps-(length(cn)))) # Add 0 if Tn < Max_Eps
    vtx_count_str = paste(vtx_count, collapse = "_") # Create string for graph lookup
  },  Max_Eps)
  unique_vtx_count_str = unique(vtx_count_strs)
  
  #==========================================================================
  # Check all the requisite graph look-ups exist before proceeding
  #==========================================================================
  graph_lookup_check = sapply(unique_vtx_count_str, function(x){
    if(!file.exists(sprintf('../../RData/graph_lookup/graph_lookup_%s.Rdata', x))){
      stop(sprintf('graph_lookup %s does not exist', x))} else {
        'graph lookup exists...'
      }
  })
  
  #==========================================================================
  # log_Pr_G_Rn: matrix with rows per Gnb and cols per Rnt
  #==========================================================================
  log_Pr_G_Rns = lapply(unique_vtx_count_str, function(x){
    load(sprintf('../../RData/graph_lookup/graph_lookup_%s.Rdata', x)) # loads all Gnb for given vtx_count_str
    G_Rn_comp = sapply(graph_lookup, test_Rn_compatible, Rns) # For each Gnb, test compatibility with Rn
    cn = as.numeric(strsplit(x, split = '_')[[1]]) # Reconstruct cn 
    Tn_chr = as.character(length(cn)) # Reconstruct Tn_chr
    Rn = Rns[[Tn_chr]] # Extract Rn 
    log_Pr_Rn = log_Pr_Rns[[Tn_chr]] # Extract log Pr( Rn | p)
    log_Pr_G_Rn = log(t(G_Rn_comp*(1/rowSums(G_Rn_comp)))) # log Pr( Gnb | Rn ) = 1 / B in matrix
    log_Pr_G_Rn[is.infinite(log_Pr_G_Rn) | is.nan(log_Pr_G_Rn)] = NA # Set -Inf due to log(0) and NAN due to division by 0 to NA
    # plot_Vivax_model(G) # uncomment when manually checking graphs make sense
    return(log_Pr_G_Rn)
  })
  names(log_Pr_G_Rns) = unique_vtx_count_str
  
  if(verbose){
    writeLines('\nIn the following viable graphs include only those with sufficient numbers of vertices 
               and appropriate vertex labels to fully represent each individual\'s data.')}
  
  # If cores=1 then this does normal sequential computation
  if(cores>1) registerDoParallel(cores = cores)
  
  #***********************************************
  # Computation of per-individual values
  #***********************************************
  Post_probs = foreach(i=1:N, .combine = c) %dopar% {
    # set id = 'BPD_70' for vtx_counts_str = "1_2_2" when checking by hand
    # set id = 'BPD_402' for vtx_counts_str = "4_1_0" when checking by hand
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
    Rn = Rns[[Tn_chr, drop = FALSE]] # Extract relevant Rn
    colnames(Rn) = infections[-1] # Extract relevant Rn and name
    vtx_count_str = vtx_count_strs[id] # Extract vertex sizes of Gn
    log_Pr_Rn = if(p_pop_ind){log_Pr_Rns[[Tn_chr]]}else{log_Pr_Rns[[id]]} # Extract log P(Rn)
    log_Pr_G_Rn = log_Pr_G_Rns[[vtx_count_str]] # Extract log P(Gnb | Rn)
    
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
      Haplotypes_and_combinations[[inf]] = list(Hnt = Hnt_chr, Vt_Hnt_inds_comp = Vt_Hnt_inds_comp)
    }
    
    
    # Extract the number of compatible combinations for each infection
    num_comp_combs_Vt = lapply(Haplotypes_and_combinations, function(x){1:nrow(x$Vt_Hnt_inds_comp)})
    # Create a matrix of inds to create all possible Gna for a given Gn
    labelled_G_ind = as.matrix(expand.grid(num_comp_combs_Vt))
    # Total number of possible mappings, A (same as prod(sapply(num_comp_combs_Vt, length)))
    A = nrow(labelled_G_ind)
    log_A = log(A) # Needed for log domain calculation
    # From labelled_G_ind create all Gna 
    Gabs = vector('list', A) # store in list
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
      Gabs[[label_ind]] = vertex_data_matrix 
    }
    
    
    #==========================================================================          
    # Load all Gnb and check complexity
    #==========================================================================
    load(sprintf('../../RData/graph_lookup/graph_lookup_%s.Rdata', vtx_count_str))
    length_graph_lookup = length(graph_lookup)
    # as.numeric to avoid getting integer overflow leading to an error
    complexity_problem = as.numeric(A)*as.numeric(length_graph_lookup)
    
    if(verbose){
      writeLines(sprintf('\nComplexity of problem, ID: %s\nThe number of different viable vertex labeled graphs is: %s\nThe number of different viable edge labeled graphs is: %s\nThe number of graphs in viable graph space is therefore: %s', 
                         id, A, length_graph_lookup, complexity_problem))
    }
    
    # We do a complexity check - need to work out what an OK upper limit should be
    # Aimee: Rather than having two complexity checks, I wonder if we might do this check
    # outside do.par/forloop (I will give it some thought)
    if(complexity_problem > UpperComplexity){ # Return NAs and skip to next ID
      
      writeLines(sprintf('\nSkipping this problem (ID %s), too complex.', id))
      Post_probs = array(NA, dim = Tn, dimnames = list(infections)) # Return NAs
      #Post_probs[names(Post_probs)] = Post_probs # Need to change this for parallel computation
      
      # return NA vector of Post_probs
      Post_probs
      
    } else { # Go on to calculate posterior probabilities
      
      #==========================================================================          
      # Calculate Pr(yn | Gnb) = sum from a = 1 to A over Pr(yn | Gnab) and record run time
      #==========================================================================
      tic() # Start clock to calculate time take per graph
      
      # Calculate probabilities of Gnb summed over all Gnab 
      log_Pr_yn_Gnbs_unnormalised = sapply(graph_lookup, Log_Pr_yn_Gnb_unnormalised, 
                                           Gabs=Gabs, cn=cn, Tn=Tn, 
                                           log_Fs=log_Fs, MSs=MSs, alpha_terms=alpha_terms) 
      
      z = toc(quiet = TRUE) # Stop clock
      ms_per_graph_space = 1000*(z$toc-z$tic) # Total time taken
      ms_per_graph = ms_per_graph_space/complexity_problem # Time taken per graph
      
      if(verbose){writeLines(sprintf('Run time (ms) over all graphs in graph space: %s\nRun time (ms) per graph in graph space: %s', round(ms_per_graph_space), round(ms_per_graph)))}
      #complexity_time[i,] = c(complexity_problem, ms_per_graph) # Store time 
      log_Pr_yn_Gnbs = log_Pr_yn_Gnbs_unnormalised - log_A # Normalise probabilities
      
      #==========================================================================          
      # Calculate P(yn | Rn) by summing over all Gnb
      #==========================================================================
      log_Pr_yn_Rn = array(dim = dim(log_Pr_G_Rn), dimnames = dimnames(log_Pr_G_Rn))
      for(col_i in 1:ncol(log_Pr_G_Rn)){
        log_Pr_yn_Rn[,col_i] = log_Pr_yn_Gnbs + log_Pr_G_Rn[,col_i]
      }
      log_Pr_yn_Rn = apply(log_Pr_yn_Rn, 2, logSumExp, na.rm = TRUE)
      
      #==========================================================================          
      # Calculate P(Rn | yn) 
      #==========================================================================
      # Pr_yn_Rn[] = 1 # Likelihood check: returns the prior
      log_Pr_yn_and_Rn = log_Pr_yn_Rn[names(log_Pr_Rn)] + log_Pr_Rn
      log_Pr_yn = logSumExp(log_Pr_yn_and_Rn)
      Pr_Rn_yn = exp(log_Pr_yn_and_Rn - log_Pr_yn) 
      
      #==========================================================================          
      # Calculate P(Rnt | yn) and return 
      #==========================================================================
      if(Tn == 2){
        Post_probs =  c(NA, Pr_Rn_yn['L'])
      }
      if(Tn == 3){ # Need to change 
        Post_probs = c(NA, sum(Pr_Rn_yn[c('LL','LI','LC')]), sum(Pr_Rn_yn[c('LL','CL','IL')]))
      }
      names(Post_probs) = infections
      
      Post_probs # return vector
    }
  }
  
  # # Store to check relationship between problem complexity and time
  # complexity_time = array(NA, c(N,2), dimnames = list(IDs, c('Number of graphs in viable graph space', 
  #                                                            'Run time (ms) per graph')))
  
  ## Uncomment for doParallel
  # Post_probs[names(Post_probs)] = Post_probs
  # Aimee I'm saving this so plot can be generated outside of function
  # James: need to think how to incorporate this into the parallel version
  # save(complexity_time, file = '../../RData/complexity_time.RData')
  return(Post_probs)
}