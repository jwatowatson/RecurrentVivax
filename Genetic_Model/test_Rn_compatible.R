# This is not presently scalable beyond Tn = {1,2}
test_Rn_compatible = function(G, Rns){ 
  
  cns = as.vector(table(vertex_attr(G)$group)) # Extract COIs (as.vector removes names which otherwise get appended to names of G_mid_Rn)
  Tn = length(cns) # Number of infections experienced by the nth individual
  Rn = Rns[[as.character(Tn)]] # Recurrence states give Tn
  
  # Reconstruct the adjacency matrix from G
  if(is.null(edge_attr(G)$w)){ 
    Adj_matrix = array(0, dim = rep(sum(cns),2)) # If there are no edge weights create matriz of zeros
  } else {
    Adj_matrix = as.matrix(get.adjacency(G, attr = 'w'))      
  } 
  
  # Work with lower triangle only (imporant for cln sums as otherwise duplicate counts)
  Adj_matrix[upper.tri(Adj_matrix)] = NA 
  
  # Generate NA mask over all block diagonals (within infection entries) regardless of Tn
  # s.t. !is.na(blocks) accesses across infections edges only
  # blocks is lower triangle only
  blocks = as.matrix(bdiag(lapply(cns, function(x){matrix(NA, ncol = x, nrow = x)})))
  # dimnames(blocks) = list(vertex_attr(G)$group,vertex_attr(G)$group) # Useful when checking but can comment after
  
  # Are all non-block diagonal elements (across infections) strangers? 
  all_str_across_inf = all(Adj_matrix[!is.na(blocks)] == 0, na.rm = T) # na.rm important due to NAs in upper triangle
  
  if(Tn == 2){ # There is only one past infection to compare to
    # True if no sibs and cn[2] clones over non-block diagonal elements t1:t2
    cn_cln_across_inf = !any(Adj_matrix[!is.na(blocks)] == 0.5,  na.rm = T) & sum(Adj_matrix[!is.na(blocks)], na.rm = T) == cns[2] 
    # Generate the compatibility of C, L, I given collection of parasites
    G_mid_Rn = c('C' = cn_cln_across_inf, 'L' = TRUE, 'I' = all_str_across_inf)
  }
  
  if(Tn == 3){ # There are either 2 past or 1 past + 1 future infections to compare to 
    
    # For Rn = "LI" mask the t=1:t=2 non-block diagonal leaving t=1:t=3 and t=2:t=3
    mask_12 = blocks
    mask_12[(cns[1]+1):sum(cns[1:2]), 1:cns[1]] = NA # Mask 1:2
    mask_12[1:cns[1],(cns[1]+1):sum(cns[1:2])] = NA # Mask 2:1
    all_str_across_13_23 = all(Adj_matrix[!is.na(mask_12)] == 0, na.rm = T) 
    
    # For Rn = "CL" and "IL" mask the t=1:t=3 and t=2:t=3 non-block diagonal leaving t=1:t=2
    mask_13_23 = blocks
    mask_13_23[sum(cns[1:2],1):sum(cns[1:3]), 1:cns[1]] = NA # Mask 1:3
    mask_13_23[1:cns[1], sum(cns[1:2],1):sum(cns[1:3])] = NA # Mask 3:1
    mask_13_23[sum(cns[1:2],1):sum(cns[1:3]), (cns[1]+1):sum(cns[1:2])] = NA # Mask 2:3
    mask_13_23[(cns[1]+1):sum(cns[1:2]), sum(cns[1:2],1):sum(cns[1:3])] = NA # Mask 3:2
    all_str_across_12 = all(Adj_matrix[!is.na(mask_13_23)] == 0, na.rm = T)
    cn_cln_across_12 = !any(Adj_matrix[!is.na(mask_13_23)] == 0.5, na.rm = T) & sum(Adj_matrix[!is.na(mask_13_23)], na.rm = T) == cns[2]
    
    # For Rn = "LC" mask the t=1:t=2 and t=1:t=3 non-block diagonals leaving only t=2:t=3
    mask_12_13 = mask_12
    mask_12_13[sum(cns[1:2],1):sum(cns[1:3]), 1:cns[1]] = NA # Mask 1:3
    mask_12_13[1:cns[1], sum(cns[1:2],1):sum(cns[1:3])] = NA # Mask 3:1
    cn_cln_across_23 = !any(Adj_matrix[!is.na(mask_12_13)] == 0.5, na.rm = T) & sum(Adj_matrix[!is.na(mask_12_13)], na.rm = T) == cns[3]
    
    G_mid_Rn = c("CC" = cn_cln_across_12 & cn_cln_across_23, 
                 "CI" = cn_cln_across_12 & all_str_across_13_23, 
                 "CL" = cn_cln_across_12, 
                 "IC" = all_str_across_12 & cn_cln_across_23,
                 "II" = all_str_across_inf, 
                 "IL" = all_str_across_12, 
                 "LC" = cn_cln_across_23, 
                 "LI" = all_str_across_13_23, 
                 "LL" = TRUE)
  }
  return(G_mid_Rn)
}

