##############################################################################################
# Notation matches manuscript aside from Tn, which counts from 1 here and 0 in the notation
# Variables names: conditioning bar is represented by an underscore while 'and' is explicit
# e.g. log_Pr_yn_Rn denotes log Pr(yn | Rn) while log_Pr_yn_and_Rn denotes log Pr(yn, Rn)
#
# Currently, we have set various limits (Max_Eps, Max_Tot_Vtx, Max_Hap_genotypes, UpperComplexity)
# allowing us to analyses and interpret the data with which we are working. Ultimately, we aim 
# to generalise the model to accommodate different data types and release it for general use. 
# The general-purpose release will likely retain the same statistical framework, but with more
# computationally sophistication, i.e. packaged as software versus statistical model code, e.g
# in optimal code upperComplexity = 10^6 may be redundant 
##############################################################################################

post_prob_CLI = function(MSdata, # MS data 
                         Fs, # MS population frequencies 
                         p = c('C' = 1/3, 'L' = 1/3, 'I' = 1/3), # Uniform prior over C, L, I
                         alpha = 0, # Additive inbreeding constant
                         cores = 4, # Number of cores for parallel computation
                         Max_Eps = 3, # Limit on number of episodes (due to test_Rn_compatible) 
                         Max_Tot_Vtx = 6, # Limit on number of vertices = cumulative COI
                         Max_Hap_genotypes = 50, # This is a limit on deterministic phasing 
                         Max_Hap_comb = 5000, # Ideally exceeds that of an episide whose hap count < Max_Hap_genotypes
                         UpperComplexity = 10^6, # Assuming 1ms per operation -> 5 hours
                         verbose = FALSE){ # Set to true to return all messages
  
  
  source('../Genetic_Model/Hap_combinations_probabilistic.R')
  
  # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # # Simulate code for internal checks. Comment when not doing internal checks
  # # Ultimate goal: integrate into unit test
  # set.seed(1)
  # sim_output = BuildSimData(COIs = c(3,1), # COI from settings
  #                           M = 6, # Number of markers from settings
  #                           N = 2, # Number of individuals
  #                           N_alleles = 4,
  #                           relatedness = 'Clone')
  # MSdata = sim_output$MS_data_sim # MS data
  # Fs = sim_output$FS # MS population frequencies
  # p = c('C' = 1/3, 'L' = 1/3, 'I' = 1/3) # Uniform prior over C, L, I
  # alpha = 0 # Additive inbreeding constant
  # cores = 4 # Number of cores for parallel computation
  # Max_Eps = 3 # Limit on number of episodes (due to test_Rn_compatible) 
  # Max_Tot_Vtx = 6 # Limit on number of vertices = cumulative COI
  # Max_Hap_genotypes = 50 # Need to discuss this with Aimee: hack insert
  # Max_Hap_comb = 5000 
  # UpperComplexity = 10^6 # Assuming 1ms per operation -> 5 hours
  # verbose = FALSE
  # #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #==========================================================================
  # Check the MSdata and p have the correct structures
  #========================================================================== 
  if(any(!c('Episode','ID','Episode_Identifier','MOI_id') %in% colnames(MSdata))){
    stop('MSdata needs to include the following columns: ID,Episode,Episode_Identifier,MOI_id')
  }
  
  #==========================================================================
  # Reformat the data s.t. there are no NA gaps in mixed infections
  #==========================================================================
  MSdata = reformat_MSdata(MSdata = MSdata, MSs=names(Fs))
  
  # Make sure that the data are sorted correctly
  MSdata = arrange(MSdata, ID, Episode, MOI_id)
  
  if(verbose) writeLines('Setting up parameters to do computation....')
  
  #==========================================================================
  # Retrieve population prior and recurrent eps identifiers
  #==========================================================================
  p_pop = sapply(c('C','L','I'), function(x) as.list(formals(post_prob_CLI)$p)[[x]])
  recurrent_eps = unique(MSdata$Episode_Identifier[MSdata$Episode>1]) # Based on genetic data
  
  #==========================================================================
  # Check to see if using pop prior or prior from time-to-event and comment
  #==========================================================================
  if(is.null(dim(p)) & length(p) == 3){
    p_pop_ind = TRUE
    writeLines('Using population prior probabilities of recurrence states')
    if(sum(p) != 1) stop('Population prior probabilities do not sum to one')
    
  } else {
    p_pop_ind = FALSE
    p = arrange(p, Episode_Identifier)
    # make sure p contains the correct columns
    if(any(!c('Episode_Identifier','C','L','I') %in% names(p))){
      stop('p needs to include the following columns: Episode_Identifier,C,L,I')
    }
    
    writeLines('Using time-to-event prior probabilities of recurrence states')
    recurrent_eps_no_prior = recurrent_eps[which(!recurrent_eps %in% p$Episode_Identifier)]
    prior_with_genetic_ind = p$Episode_Identifier %in% recurrent_eps
    prior_no_genetic = p$Episode_Identifier[!prior_with_genetic_ind]
    if(length(prior_no_genetic) > 0){
      writeLines(sprintf('No genetic data for following time-to-event priors (thus removed): %s', 
                         paste(prior_no_genetic, collapse = ', ')))
      p = p[prior_with_genetic_ind, ]
    }
    num_recurrent_eps_no_prior = length(recurrent_eps_no_prior)
    if(num_recurrent_eps_no_prior == 0){
      writeLines('Using time-to-event prior probabilities for all recurrent infections with genetic data')
    } else {
      writeLines(sprintf('Using time-to-event prior probabilities but for the 
                         following for which population prior probabilities are used: %s', 
                         paste(recurrent_eps_no_prior, collapse = ', ')))
      # Add to episodes based on the prior to the top of the list
      for(x in recurrent_eps_no_prior){p = rbind(c(x, p_pop), p)} 
    }
    # check time-to-event and genetic now agree
    if(!all(p$Episode_Identifier %in% recurrent_eps) & all(recurrent_eps %in% p$Episode_Identifier)){
      stop('Disagreement between epsiodes with available time-to-event and genetic data')
    }
  }
  
  #==========================================================================
  # Create Post_probs store with length = num. of all infections inc. those without recurrence
  #==========================================================================
  all_infections = unique(MSdata$Episode_Identifier)
  
  #==========================================================================
  # Extract microsatellites, related quantities and calculate log pop prior
  #==========================================================================
  MSs = names(Fs)  
  M = length(MSs) 
  log_Fs = lapply(Fs, log)
  if(p_pop_ind){
    log_p = log(p)
  } else {
    log_p = t(apply(p, 1, function(x)log(as.numeric(x[ names(x) %in% c('C','L','I')]))))
    rownames(log_p) = p[,'Episode_Identifier']
    colnames(log_p) = colnames(p[,names(p) %in% c('C','L','I')])
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
  IDs_all = as.character(unique(MSdata$ID)) # Character vector
  yns = dlply(MSdata, 'ID') # transform into a list
  
  #==========================================================================
  # Extract COIs (cns) and num. of infections (Tns) for each individual
  # Note that, in contrast to notation, here Tn counts from 1 
  #==========================================================================
  cns = lapply(yns, function(x){cn = table(x$Episode_Identifier)})
  sum_cns = sapply(cns, sum) # Size of graph = sum(cnt) for check below
  Tns = lapply(cns, length) # Number of episodes (also used in check below)
  Tns_chr = lapply(Tns, as.character) # Character version needed for indexes in do.par
  
  #==========================================================================
  # Check for individuals whose data are too complex for theta calculation
  #==========================================================================
  complex = (Tns > Max_Eps) | (sum_cns > Max_Tot_Vtx)
  if(any(complex)){
    writeLines(sprintf('\nSkipping ID: %s, which are too complex either because they 
                       have more than %s episodes or more than %s vertices, sorry.', 
                       paste(IDs_all[complex], collapse = ', '), Max_Eps, Max_Tot_Vtx))}
  
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
    for(x in names(which(unlist(Tns[IDs]) > 2))){ # For those with more than recurrence sum over recurrences 
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
    if(!file.exists(sprintf('../RData/graph_lookup/graph_lookup_%s.Rdata', x))){
      stop(sprintf('graph_lookup %s does not exist', x))} else {
        'graph lookup exists...'
      }
  })
  
  #==========================================================================
  # log_Pr_G_Rn: matrix with rows per Gnb and cols per Rnt 
  # (haploide-genotyped unlabelled graphs)
  #==========================================================================
  log_Pr_G_Rns = lapply(unique_vtx_count_str, function(x){
    load(sprintf('../RData/graph_lookup/graph_lookup_%s.Rdata', x)) # loads all Gnb for given vtx_count_str
    G_Rn_comp = sapply(graph_lookup, test_Rn_compatible, Rns) # For each Gnb, test compatibility with Rn
    log_Pr_G_Rn = log(t(G_Rn_comp*(1/rowSums(G_Rn_comp)))) # log Pr( Gnb | Rn ) = 1 / B in matrix
    log_Pr_G_Rn[is.infinite(log_Pr_G_Rn) | is.nan(log_Pr_G_Rn)] = NA# Set -Inf due to log(0) and NAN due to division by 0 to NA  
    return(log_Pr_G_Rn)
  })
  names(log_Pr_G_Rns) = unique_vtx_count_str
  
  # #===================================================
  # # Aside: uncomment when manually checking log_Pr_G_Rns make sense
  # x = unique_vtx_count_str[10]
  # load(sprintf('../../RData/graph_lookup/graph_lookup_%s.Rdata', x)) # loads all Gnb for given vtx_count_str
  # G_Rn_comp = sapply(graph_lookup, test_Rn_compatible, Rns)
  # RR_comp_Gs = which(G_Rn_comp['C',])
  # par(mfrow = c(6,3), mar = c(1,1,1,1))
  # sapply(RR_comp_Gs, function(z) plot_Vivax_model(graph_lookup[[z]]))
  # #===================================================
  
  if(verbose){
    writeLines('\nIn the following viable graphs include only those with sufficient numbers of vertices 
               and appropriate vertex labels to fully represent each individual\'s data.')}
  
  #==========================================================================
  # Identify samples with too many compatable haploid genotypes 
  # Since warnings do not work inside %dopar% 
  # Added by Aimee May 24th
  #==========================================================================
  n_haps_per_episode = lapply(yns, function(yn){ # For a given individual
    ynts = dlply(yn, 'Episode_Identifier') # Break into episodes
    sapply(ynts, function(ynt){ # For a given episode
      ynmt = ynt[,MSs,drop = FALSE] # Extract microsatellite data 
      ynmt_size = apply(ynmt, 2, function(x){length(unique(x))})  
      num_haps = prod(ynmt_size) # Calculate the number of haploid genotypes
    })
  })
  names(n_haps_per_episode) = NULL # Remove IDs since episode ID sufficient 
  n_haps_per_episode = unlist(n_haps_per_episode) # Convert to vector and then string
  if(any(n_haps_per_episode > Max_Hap_genotypes)){ # If there are any too complex, print names to screen
    episodes_with_too_many_hap_to_phase = paste0(names(which(n_haps_per_episode > Max_Hap_genotypes)), collapse = ', ')
    writeLines(sprintf('Hack will be used on the following epsisodes with more than %s haploid genotypes: \n %s',
                       Max_Hap_genotypes,episodes_with_too_many_hap_to_phase))
  }
  
  
  #***********************************************
  # Computation of per-individual values
  #***********************************************
  # dopar in Windows needs the packages and functions explicitly passed to foreach command
  # otherwise it doesn't work for more than one core. This is still compatible with Mac
  if(cores>1) registerDoParallel(cores = cores) # If cores=1 then this does normal sequential computation
  Post_probs = foreach(i=1:N, .combine = rbind, 
                       .packages = c('dplyr','igraph','gtools',
                                     'matrixStats','Matrix','tictoc'), 
                       .export = c('Log_Pr_yn_Gab','Log_Pr_yn_Gnb_unnormalised',
                                   'test_Rn_compatible','test_cln_incompatible')
  ) %dopar%
  {
    
    # set id = 'BPD_91' for vtx_counts_str = "2_1_0" when checking by hand
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
    # Generate all possible vertex haploid genotype label mappings of data onto a graph 
    #==========================================================================
    Hap_combinations = vector('list', length = 0) # Store of haploid genotype combinations 
    
    for(inf in infections){
      
      # Extract data for the tth infection
      ynt = filter(yn, Episode_Identifier == inf)[,MSs,drop = FALSE] 
      
      # Summarise data for compatibility check and Hap_combinations_probabilistic() 
      # alply ensures it's always as a list, important for Hap_combinations_probabilistic() 
      Y = alply(ynt, 2, function(x){sort(unique(x[!is.na(x)]))}) 
      
      # All haploid genotypes compatible with ynt (`unique` collapses repeats due to row per MOI)
      Hnt = expand.grid((lapply(ynt, unique)))
      total_haps_count = nrow(Hnt)
      
      # This line allows the model to process entirely NA data 
      if(total_haps_count < cn[inf]){Hnt = ynt} 
      
      # If very many compatible haploid genotypes, adopt probablistic phasing approach
      if(total_haps_count > Max_Hap_genotypes){ 
        Hap_combinations[[inf]] = Hap_combinations_probabilistic(Max_Hap_comb, cn[inf], ynt, Y)}
      
      # If moderate haploid genotypes, adopt deterministic phasing approach
      if(total_haps_count >= cn[inf] & total_haps_count <= Max_Hap_genotypes){  
        
        # Indices of all combinations of nrow(Hnt) choose cn[inf] haploid genotypes for the tth infection, inf 
        Vt_Hnt_inds = combinations(nrow(Hnt), r = cn[inf], v = 1:nrow(Hnt))  
        
        # Check each combination to see if compatible with ynt 
        # (this is a bit like an ABC step with epsilon = 0)
        scores = apply(Vt_Hnt_inds, 1, function(x, Y){ 
          X = alply(Hnt[x,,drop=FALSE], 2, function(x){sort(unique(x[!is.na(x)]))}) 
          score = identical(X,Y) # Test that they are the same
          return(score)
        }, Y)
        
        # Extract only those indices that are compatible 
        Vt_Hnt_inds_comp = Vt_Hnt_inds[scores,,drop = FALSE]
        
        # Convert to character since used to index 
        Hnt_chr = matrix(sapply(Hnt, as.character), ncol = M)
        colnames(Hnt_chr) = MSs
        
        # Return all possible compatible combinations of haploid genotypes as a list (apply returns array)
        Hap_combinations[[inf]] = alply(Vt_Hnt_inds_comp,1,function(i){Hnt_chr[i,,drop = F]})
      }
    }
    
    
    # Extract vector of compatible combinations for each episode
    num_comp_combs_Vt = lapply(Hap_combinations, function(x){1:length(x)})
    
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
      for(t in 1:Tn){ # Since each infection has its own set of combinations, we need to loop over infections
        vertex_data_matrix[cumsum(c(1,cn))[t]:cumsum(cn)[t], ] = Hap_combinations[[infections[t]]][[x[t]]] 
      }
      Gabs[[label_ind]] = vertex_data_matrix # Characters since used to extract allele frequ. 
    }
    
    #==========================================================================          
    # Load all Gnb and check complexity
    #==========================================================================
    load(sprintf('../RData/graph_lookup/graph_lookup_%s.Rdata', vtx_count_str))
    length_graph_lookup = length(graph_lookup)
    # as.numeric to avoid getting integer overflow leading to an error
    complexity_problem = as.numeric(A)*as.numeric(length_graph_lookup)
    
    if(verbose){
      writeLines(sprintf('\nComplexity of problem, ID: 
                         %s\nThe number of different viable vertex labeled graphs is: 
                         %s\nThe number of different viable edge labeled graphs is: 
                         %s\nThe number of graphs in viable graph space is therefore: %s', 
                         id, A, length_graph_lookup, complexity_problem))
    }
    
    # Create recurrences names (avoid infections[-1] incase misordered)
    # I've added an `arrange` line at the start of the function so shouldn't be a problem
    inf_no = yn$Episode[!duplicated(yn$Episode_Identifier)] #do.call(rbind, strsplit(infections, split = '_'))[,3]
    recurrences = paste(id, inf_no[inf_no > 1], sep = '_')
    # We do a complexity check - need to work out what an OK upper limit should be
    # Aimee: Rather than having two complexity checks, I wonder if we might do this check
    # outside do.par/forloop 
    if(complexity_problem > UpperComplexity){ # Return NAs and skip to next ID
      
      writeLines(sprintf('\nSkipping this problem (ID %s), too complex (%s).', id, complexity_problem))
      Post_probs = data.frame(C=rep(NA,length(recurrences)),
                              L=rep(NA,length(recurrences)),
                              I=rep(NA,length(recurrences)))
      rownames(Post_probs) = recurrences
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
      
      if(verbose){writeLines(sprintf('Run time (ms) over all graphs in graph space: 
                                     %s\nRun time (ms) per graph in graph space: %s', 
                                     round(ms_per_graph_space), round(ms_per_graph,2)))}
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
      
      if(Tn == 2){ # Return a data.frame C, L, I for single recurrence
        Post_probs = data.frame(C = Pr_Rn_yn['C'],
                                L = Pr_Rn_yn['L'],
                                I = Pr_Rn_yn['I'])
        rownames(Post_probs) = recurrences
      }
      if(Tn == 3){ # Return a data.frame C, L, I for recurrences 1 and 2
        Post_probs = data.frame(C = c(sum(Pr_Rn_yn[c('CL','CI','CC')]), sum(Pr_Rn_yn[c('LC','CC','IC')])), 
                                L = c(sum(Pr_Rn_yn[c('LL','LI','LC')]), sum(Pr_Rn_yn[c('LL','CL','IL')])), 
                                I = c(sum(Pr_Rn_yn[c('IL','II','IC')]), sum(Pr_Rn_yn[c('LI','CI','II')])))
        rownames(Post_probs) = recurrences
      }
      
      Post_probs # return vector
    }
  }
  
  # # Store to check relationship between problem complexity and time
  # complexity_time = array(NA, c(N,2), dimnames = list(IDs, c('Number of graphs in viable graph space', 
  #                                                            'Run time (ms) per graph')))
  #
  ## Uncomment for doParallel
  # Post_probs[names(Post_probs)] = Post_probs
  # Aimee I'm saving this so plot can be generated outside of function
  # James: need to think how to incorporate this into the parallel version
  # save(complexity_time, file = '../../RData/complexity_time.RData')
  return(Post_probs)
}