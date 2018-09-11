test_Rn_compatible = function(G){  # This is not yet scalable beyond Tn = {2,3}
  
  cn = table(vertex_attr(G)$group) # Extract MOIs
  Tn = length(cn) # Number of infections experienced by the nth individual
  
  # If there are no edge weights (i.e. all strangers), all compatible
  if(is.null(edge_attr(G)$w)){
    if(Tn == 2){G_mid_Rn = c('1' = TRUE, '0' = TRUE)}
    if(Tn == 3){G_mid_Rn = c('11' = TRUE, '10' = TRUE, '01' = TRUE, '00' = TRUE)}
    return(G_mid_Rn)
  } else{ # If there are edge weights, must test compatibility
    
    # Reconstruct the adjaceny matrix
    A = as.matrix(get.adjacency(G, attr = 'w'))
    
    # Generate mask for all block diagonals regardless of Tn
    blocks = as.matrix(bdiag(lapply(cn, function(x){matrix(NA, ncol = x, nrow = x)})))
    # dimnames(blocks) = list(vertex_attr(G)$group,vertex_attr(G)$group) # useful to name for checking purposes
    all_unrelated_acrossInf = all(A[!is.na(blocks)] == 0) # non-block diagonal elements
    G_mid_Rn = c('1' = TRUE, '0' = all_unrelated_acrossInf)
    
    if(Tn == 3){
      
      # For Rn = 01 mask the t=0:t=1 non-block diagonal
      mask_10 = blocks
      mask_10[(cn[1]+1):sum(cn[1:2]), 1:cn[1]] = NA # Mask 1:2
      mask_10[1:cn[1],(cn[1]+1):sum(cn[1:2])] = NA # Mask 2:1
      score_10 = all(A[!is.na(mask_10)] == 0)
      
      # For Rn = 10 mask the t=1:t=2 and t=0:t=1 non-block diagonal
      mask_01 = blocks
      mask_01[sum(cn[1:2],1):sum(cn[1:3]), 1:cn[1]] = NA # Mask 1:3
      mask_01[1:cn[1], sum(cn[1:2],1):sum(cn[1:3])] = NA # Mask 3:1
      mask_01[sum(cn[1:2],1):sum(cn[1:3]), (cn[1]+1):sum(cn[1:2])] = NA # Mask 2:3
      mask_01[(cn[1]+1):sum(cn[1:2]), sum(cn[1:2],1):sum(cn[1:3])] = NA # Mask 3:2
      score_01 = all(A[!is.na(mask_01)] == 0)
      
      G_mid_Rn = c('11' = TRUE, '10' = score_10, '01' = score_01, '00' = all_unrelated_acrossInf)
    }
    return(G_mid_Rn)
  }
}