# list of functions for the graph algorithm
# Tue Sept 18th: updated notation throughout: 
# COI -> COI

# This implements the Viable graph brute-force search algorithm described in the appendix
test_correct_graph = function(G){
  
  INCORRECT = FALSE # Presume correct until proven incorrect
  
  # First check that all connected components are cliques
  S = components(graph = G)
  for(k in 1:S$no){
    G_sub = induced_subgraph(graph = G, vids = vertex_attr(G)$label[S$membership==k])
    MC = max_cliques(graph = G_sub, min = vcount(G_sub)) 
    if(length(MC) == 0) INCORRECT = TRUE
  }
  
  # If all connected components are cliques, then iterate through all cliques of size 3
  if(!INCORRECT){ 
    all_cliques = cliques(graph = G, min = 3, max=3)
    K = length(all_cliques)
    if(K > 0){
      i = 1
      while(!INCORRECT & i < (K+1) ){ # stop searching as soon as we find it's incorrect
        G_sub = induced_subgraph(graph = G, vids = all_cliques[[i]])
        Score_sub = sum(edge_attr(G_sub)$w) # the sum of the edges tells us if it's correct
        if(Score_sub == 2.5){ INCORRECT = TRUE } # 2.5 is the only possible incorrect score so we reject
        i = i+1
      }
    }
  }
  if(INCORRECT) return(FALSE) else return(TRUE)
}

generate_random_Vivax_model = function(COI_T1, COI_T2){
  # generate the adjacency matrix within time point t1
  Adj_T1 = matrix( sample( c(0,.5), COI_T1^2, replace=T), ncol = COI_T1)
  # then within t2
  Adj_T2 = matrix( sample( c(0,.5), COI_T2^2, replace=T), ncol = COI_T2)
  # and across timepoints
  Adj_T1_T2 = matrix( sample( c(0,.5, 1), COI_T1*COI_T2, replace=T), ncol = COI_T1)
  # merge matrices
  Adj_all = bdiag(Adj_T1, Adj_T2)
  # add the across component
  Adj_all[(COI_T1+1):(COI_T1+COI_T2), 1:COI_T1] = Adj_T1_T2
  
  # Just take the lower triangular matrix
  Adj_all[upper.tri(Adj_all)] = 0
  Adj_all = as.matrix(Adj_all)
  colnames(Adj_all) = 1:(COI_T1+COI_T2)
  diag(Adj_all) = 0
  # make it symmetrical
  Adj_all = Adj_all + t(Adj_all)
  
  # Make the graph in igraph
  G1 = graph_from_adjacency_matrix(adjmatrix = Adj_all, mode='undirected', 
                                   diag=F, weighted = 'w', add.colnames = T)
  # set the group attributes for G1, this is for plotting
  G1 = set_vertex_attr(G1, "group", value = c(rep(1,COI_T1), rep(2,COI_T2)))
  G1 = set_vertex_attr(G1, "label", value = 1:(COI_T1+COI_T2))
  return(G1)
}

plot_Vivax_model = function(G, vertex_labels = TRUE, 
                            edge_col = c('0' = 'white', '0.5' = 'red', '1' = 'red'), 
                            edge_lty = c('0' = 1, '0.5' = 2, '1' = 1)){
  
  Ms = as.vector(table(vertex_attr(G)$group))
  
  N = sum(Ms,na.rm=TRUE)
  if( ! vcount(G) == N ) stop('COI vector not correct')
  if( sum(is.na(Ms)) > 0 ) stop('Don\'t allow NAs in COI vector')
  
  Ts = length(Ms) # No. of episodes 
  LL = array(dim = c(N,2))
  X = rbind(seq(0,1,length.out = Ts), Ms)
  # AT: jitter within time point
  z = as.numeric(rep(Ms, Ms) > 2)
  jitter = (-1)^(1:length(z))*(.2/Ts)
  LL[,1] = unlist(apply(X, 2, function(x) rep(x[1], x[2]))) + jitter*z #+ rep(c(0, .2/Ts)
  LL[,2] = unlist(sapply(Ms, function(x) {
    if(x==1) return(0.5) else return(seq(0,1,length.out = x))
  }
  ))
  # # The layout vector so that it is legible when plotted
  # l = cbind(c(array(c(0,.1), M1),array(c(.9,1), M2)), 
  #            c(seq(0,1,length.out = M1), seq(0,1,length.out = M2)))
  
  # plot the graph
  if(vertex_labels){
    plot(G, vertex.color = vertex_attr(G)$group, 
         edge.color = edge_col[as.character(edge_attr(G)$w)], 
         edge.width = 2, 
         edge.lty = edge_lty[as.character(edge_attr(G)$w)], 
         layout=LL)
  } else { 
    G_no_label = set_vertex_attr(G, "label", value = NA) 
    plot(G_no_label, vertex.color = vertex_attr(G)$group, 
         edge.color = edge_col[as.character(edge_attr(G)$w)], 
         edge.width = 2, 
         edge.lty = edge_lty[as.character(edge_attr(G)$w)], 
         layout=LL)
  }
  
}


##############################################################################
##############################################################################
generate_all_models_2Ts = function(M1=3, M2=3){
  
  if(sum(c(M1,M2))>6) stop('Whooh there, too complex') 
  
  # Set up all the within timepoint models
  T1models = if(M1 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M1,2), v = c(0,.5), repeats.allowed = TRUE)}
  T2models = if(M2 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M2,2), v = c(0,.5), repeats.allowed = TRUE)}
  
  # and the across timepoint models
  T12models = as.matrix(expand.grid( rep(list(c(0,.5, 1)), M1*M2)))
  
  # compute size of problem 
  KK1 = nrow(T1models)
  KK2 = nrow(T2models)
  KK12= nrow(T12models)
  KK = KK1 * KK2 * KK12
  print(sprintf('No. of graphs compatible with COIs: %s', KK))
  
  res = array(FALSE, dim = c(KK1, KK2, KK12))
  Pbar = txtProgressBar(min=1, max = KK)
  ind = 1
  correct_ind = 1
  
  correct_models = list()
  for(i in 1:KK1){
    for(j in 1:KK2){
      for(k in 1:KK12){
        # plotting progress...
        setTxtProgressBar(Pbar, ind)
        ind = ind+1
        
        # Lets make the sub adjacency matrices
        Adj_T1 = matrix(rep(0, M1^2), ncol = M1)
        Adj_T2 = matrix(rep(0, M2^2), ncol = M2)
        
        Adj_T1[lower.tri(Adj_T1)] = as.vector(T1models[i,])
        Adj_T2[lower.tri(Adj_T2)] = as.vector(T2models[j,])
        
        Adj_T1_T2 = matrix(T12models[k,], ncol = M1)
        
        # now we put them together for a full matrix
        Adj_all = matrix(rep(0,(M1+M2)^2), ncol = M1+M2)
        Adj_all[1:M1, 1:M1] = Adj_T1
        Adj_all[(M1+1):(M1+M2), (M1+1):(M1+M2)] = Adj_T2
        Adj_all[(M1+1):(M1+M2), 1:M1] = Adj_T1_T2
        Adj_all = Adj_all + t(Adj_all)
        colnames(Adj_all) = 1:(M1+M2)
        # turn this into an igraph item
        G1 = graph_from_adjacency_matrix(adjmatrix = Adj_all, mode='undirected', 
                                         diag=F, weighted = 'w', add.colnames = T)
        # set the group attributes for G1
        G1 = set_vertex_attr(G1, "group", value = c(rep(1,M1), rep(2,M2)))
        G1 = set_vertex_attr(G1, "label", value = 1:(M1+M2))
        
        # test correctness
        res[i,j,k] = test_correct_graph(G = G1)
        if(res[i,j,k]){
          correct_models[[correct_ind]] = G1
          correct_ind=correct_ind+1
        }
      }
    }
  }
  return(correct_models)
}

# Check if Pr_yn_G = 0 due to a clonal edge between non-identical vertices
test_cln_incompatible = function(A, vertex_data_matrix){
  
  cln_edges = which(A == 1, arr.ind = TRUE) # Which edges are clones
  num_cln_edges = nrow(cln_edges) # Number of edges that are clones
  
  if(num_cln_edges == 0){ # If there are no clones, return FALSE
    return(FALSE) 
  } else { # For each clonal edge
    for(i in 1:nrow(cln_edges)){ 
      score = !identical(vertex_data_matrix[cln_edges[i,'row'],],vertex_data_matrix[cln_edges[i,'col'],]) 
      if(score){break()} # as soon as different, fail and break
    }
    return(score)
  }
}

# Function to sum over vertex haploid genotype labels given edge kinship labelled graph
# Returns: Pr(yn | Gnb) unnormalised = sum from a = 1 to A over Pr(yn | Gnab)
# Unnormalised because  Pr(yn | Gnb) = sum from a = 1 to A over Pr(yn | Gnab) * (1/A)
Log_Pr_yn_Gnb_unnormalised = function(G, Gabs, cn, Tn, log_Fs, MSs, alpha_terms){
  
  sum_cn = sum(cn)
  
  # Get adjacency matrix for Gnb
  if(is.null(edge_attr(G)$w)){ # If all edges are strangers generate A afresh
    A = array(0, dim = c(sum_cn, sum_cn))
  } else { # Extract from A from graph
    A = as.matrix(get.adjacency(G, attr = 'w'))
  } 
  A[upper.tri(A, diag = TRUE)] = NA # Work with lower triangular only
  
  # Sum over all Gna given adjacency matrix for Gnb
  log_Pr_yn_Gabs = sapply(Gabs, function(x){ # labelled_G_ind inc. compatible labelled graphs only
    Log_Pr_yn_Gab(G, log_Fs, MSs, Z = alpha_terms, A = A, vertex_data_matrix = x) # Pass Gn to Log_Pr_yn_Gab
  })
  
  log_pr_yn_Gb_unnormalised = logSumExp(log_Pr_yn_Gabs, na.rm = TRUE) 
  return(log_pr_yn_Gb_unnormalised)
}


###############################################################################################
# MS where one or more NA are ignored when summing over m
# G is assumed labelled
# Z is a vector of alpha_terms (calculated above since not data dependent)
###############################################################################################
Log_Pr_yn_Gab = function(G, log_Fs, MSs, Z, A, vertex_data_matrix){ 
  
  M = length(MSs) # Number of microsatellites

  # Check if data have zero prob given G due to clones
  if(test_cln_incompatible(A, vertex_data_matrix)){ # If labelled G is incompatible
    
    log_Pr_yn_Gab = NA 
    return(log_Pr_yn_Gab) # Function ends here 
    
  } else { 
    
    # Calculate non-zero Pr( yn | G )
    log_eij = array(NA, dim = dim(A))
    
    for(i in 2:nrow(A)){ # Over all edges 
      for(j in 1:(i-1)){
        
        log_eijm = rep(NA, M) # Store 
        hi = vertex_data_matrix[i,] # Extract h_im 
        hj = vertex_data_matrix[j,] # Extract h_jm
        I_hj = hi == hj # Marker identity indicator
        
        # Stranger edge indicator
        if(A[i,j] == 0){
          for(m in 1:M){ # For each microsatellite m = 1:M 
            log_fs = log_Fs[[ MSs[m] ]][c(hi[m], hj[m])] # Extract the log frequencies for h_im and h_jm
            LX = c(Z['log_alpha'] + log_fs[1], Z['log_1_alpha'] + sum(log_fs)) # log(x) vector for logSumExp
            if(Z['alpha'] == 0){ # Avoid logSumExp as log(alpha = 1) = -Inf will overwhelm
              log_eijm[m] = sum(log_fs) 
            } else {
              log_eijm[m] = ifelse(I_hj[m], logSumExp(LX), sum(Z['log_1_alpha'], log_fs))
            }
          }
        }
        
        # Sibling edge indicator
        if(A[i,j] == .5){
          for(m in 1:M){ # For each microsatellite m = 1:M 
            log_fs = log_Fs[[ MSs[m] ]][c(hi[m], hj[m])] # Extract the log frequencies for h_im and h_jm
            LX = c(Z['log_0.5alpha'] + log_fs[1], Z['log_0.5_alpha'] + sum(log_fs)) # log(x) vector for logSumExp
            log_eijm[m] = ifelse(I_hj[m], logSumExp(LX), sum(Z['log_0.5_alpha'],log_fs))
          }
        }
          
        # Clone edge indicator
        if(A[i,j] == 1){
          for(m in 1:M){ # For each microsatellite m = 1:M 
            log_fs = log_Fs[[ MSs[m] ]][c(hi[m], hj[m])] # Extract the log frequencies for h_im and h_jm
            log_eijm[m] = ifelse(I_hj[m], log_fs[1], NA)
          }
        }
        log_eij[i,j] = sum(log_eijm, na.rm = TRUE) # Calculate Pr_hij_eij 
      }
    }
    log_Pr_yn_Gab = sum(log_eij, na.rm = TRUE)
    return(log_Pr_yn_Gab)
  }
}


################################################################
################################################################

generate_all_models_3Ts = function(M1=2, M2=2, M3=2){
  if(M1>4 | M2>4 | M3>4 | sum(c(M1,M2,M3))>6) stop('Whooh there, too complex')
  
  # # Set up all the within timepoint models
  
  # AT: I think this is a only slightly more elegant way for within time points
  T1models = if(M1 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M1,2), v = c(0,.5), repeats.allowed = TRUE)}
  T2models = if(M2 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M2,2), v = c(0,.5), repeats.allowed = TRUE)}
  T3models = if(M3 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M3,2), v = c(0,.5), repeats.allowed = TRUE)}
  
  # and the across timepoint models
  T12models = as.matrix(expand.grid( rep(list(c(0,.5, 1)), M1*M2) ))
  T23models = as.matrix(expand.grid( rep(list(c(0,.5, 1)), M2*M3) ))
  T13models = as.matrix(expand.grid( rep(list(c(0,.5, 1)), M1*M3) ))
  
  # compute size of problem
  KK1 = nrow(T1models)
  KK2 = nrow(T2models)
  KK3 = nrow(T3models)
  
  KK12= nrow(T12models)
  KK23= nrow(T23models)
  KK13= nrow(T13models)
  
  KK = (KK1 * KK2 * KK3) * (KK12 * KK23 * KK13)
  
  print(paste('complexity of problem is: ', KK))
  
  
  Pbar = txtProgressBar(min=1, max = KK)
  ind = 1
  correct_ind = 1; correct_models = list()
  
  # Now for the iteration with brute force...
  for(i1 in 1:KK1){
    for(i2 in 1:KK2){
      for(i3 in 1:KK3){
        for(i12 in 1:KK12){
          for(i23 in 1:KK23){
            for(i13 in 1:KK13){
              # plotting progress...
              setTxtProgressBar(Pbar, ind); ind = ind+1
              # Lets make the sub adjacency matrices
              Adj_T1 = matrix(rep(0, M1^2), ncol = M1)
              Adj_T2 = matrix(rep(0, M2^2), ncol = M2)
              Adj_T3 = matrix(rep(0, M3^2), ncol = M3)
              
              Adj_T1[lower.tri(Adj_T1)] = as.vector(T1models[i1,])
              Adj_T2[lower.tri(Adj_T2)] = as.vector(T2models[i2,])
              Adj_T3[lower.tri(Adj_T3)] = as.vector(T3models[i3,])
              
              Adj_T1_T2 = matrix(T12models[i12,], ncol = M1) # this has M1 cols and M2 rows
              Adj_T2_T3 = matrix(T23models[i23,], ncol = M2) # this has M2 cols and M3 rows
              Adj_T1_T3 = matrix(T13models[i13,], ncol = M1) # this has M1 cols and M3 rows
              
              Adj_all = matrix(rep(0,(M1+M2+M3)^2), ncol = M1+M2+M3)
              Adj_all[1:M1, 1:M1] = Adj_T1
              Adj_all[(M1+1):(M1+M2), (M1+1):(M1+M2)] = Adj_T2
              Adj_all[(M1+M2+1):(M1+M2+M3), (M1+M2+1):(M1+M2+M3)] = Adj_T3
              
              Adj_all[(M1+1):(M1+M2), 1:M1] = Adj_T1_T2
              Adj_all[(M1+M2+1):(M1+M2+M3), (M1+1):(M1+M2)] = Adj_T2_T3
              Adj_all[(M1+M2+1):(M1+M2+M3), 1:M1] = Adj_T1_T3
              
              Adj_all = Adj_all + t(Adj_all)
              colnames(Adj_all) = 1:(M1+M2+M3)
              
              # turn this into an igraph item
              G1 = graph_from_adjacency_matrix(adjmatrix = Adj_all, mode='undirected', 
                                               diag=F, weighted = 'w', add.colnames = T)
              # set the group attributes for G1
              G1 = set_vertex_attr(G1, "group", value = c(rep(1,M1), rep(2,M2), rep(3,M3)))
              G1 = set_vertex_attr(G1, "label", value = 1:(M1+M2+M3))
              
              # test correctness and store if correct
              if(test_correct_graph(G = G1)){
                correct_models[[correct_ind]] = G1
                correct_ind=correct_ind+1
              }
            }
          }
        }
      }
    }
  }
  return(correct_models)
}

