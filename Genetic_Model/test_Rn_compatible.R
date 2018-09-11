# This is not presently scalable beyond Tn = {1,2}
#+++++++++++++++++++++++
# something wrong with rowSums(G_Rn_comp)
#+++++++++++++++++++++++
test_Rn_compatible = function(G, Rns){ 
  
  cns = table(vertex_attr(G)$group) # Extract COIs
  Tn = length(cns) # Number of infections experienced by the nth individual
  Rn = Rns[[as.character(Tn)]]
  
  # Reconstruct the adjaceny matrix
  if(is.null(edge_attr(G)$w)){ 
    Adj_matrix = array(0, dim = rep(sum(cns),2)) # If there are no edge weights
  } else {
    Adj_matrix = as.matrix(get.adjacency(G, attr = 'w'))      
  } 
  
  # Generate NA mask for all block diagonals regardless of Tn
  blocks = as.matrix(bdiag(lapply(cns, function(x){matrix(NA, ncol = x, nrow = x)})))
  dimnames(blocks) = list(vertex_attr(G)$group,vertex_attr(G)$group) # Name when checking 
  
  # Note that !is.na(blocks) accesses across infections edges only
  all_cln_across_inf = all(Adj_matrix[!is.na(blocks)] == 1) # non-block diagonal elements
  all_str_across_inf = all(Adj_matrix[!is.na(blocks)] == 0) # non-block diagonal elements
  
  if(Tn == 2){ # There is only one past infection to compare to 
    # Generate the compatibility of C, L, I given collection of parasites
    G_mid_Rn = c('C' = all_cln_across_inf, 'L' = TRUE, 'I' = all_str_across_inf)
  }
  
  if(Tn == 3){ # There are either 2 past or 1 past + 1 future infections to compare to 
    
    # For Rn = "LC" and "LI" mask the t=1:t=2 non-block diagonal leaving t=1:t=3 and t=2:t=3
    mask_10 = blocks
    mask_10[(cns[1]+1):sum(cns[1:2]), 1:cns[1]] = NA # Mask 1:2
    mask_10[1:cns[1],(cns[1]+1):sum(cns[1:2])] = NA # Mask 2:1
    all_str_across_10 = all(Adj_matrix[!is.na(mask_10)] == 0) 
    all_cln_across_10 = all(Adj_matrix[!is.na(mask_10)] == 1)
    
    # For Rn = "CL" and "IL" mask the t=1:t=3 and t=2:t=3 non-block diagonal leaving t=1:t=2
    mask_01 = blocks
    mask_01[sum(cns[1:2],1):sum(cns[1:3]), 1:cns[1]] = NA # Mask 1:3
    mask_01[1:cns[1], sum(cns[1:2],1):sum(cns[1:3])] = NA # Mask 3:1
    mask_01[sum(cns[1:2],1):sum(cns[1:3]), (cns[1]+1):sum(cns[1:2])] = NA # Mask 2:3
    mask_01[(cns[1]+1):sum(cns[1:2]), sum(cns[1:2],1):sum(cns[1:3])] = NA # Mask 3:2
    all_str_across_01 = all(Adj_matrix[!is.na(mask_01)] == 0)
    all_cln_across_01 = all(Adj_matrix[!is.na(mask_01)] == 1)
    
    G_mid_Rn = c("CC" = all_cln_across_inf, 
                 "CI" = all_cln_across_01 & all_str_across_10, 
                 "CL" = all_cln_across_01, 
                 "IC" = all_str_across_01 & all_cln_across_10,
                 "II" = all_str_across_inf, 
                 "IL" = all_str_across_01, 
                 "LC" = all_cln_across_10, 
                 "LI" = all_str_across_10, 
                 "LL" = TRUE)
  }
  return(G_mid_Rn)
}

