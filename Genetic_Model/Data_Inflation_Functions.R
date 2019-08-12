Inflate_into_pairs = function(MS_data){
  # get the true IDs (patients)
  ids = unique(MS_data$ID)
  # create an empty matrix
  inflated_data = data.frame(matrix(rep(NA,ncol(MS_data)), ncol = ncol(MS_data)))
  # Some extra information columns
  colnames(inflated_data) = colnames(MS_data)
  inflated_data$ID_True = NA
  inflated_data$First_EpNumber = NA
  inflated_data$Second_EpNumber = NA
  pb = txtProgressBar(min=0, max = length(ids),style = 3)
  for(id in ids){
    setTxtProgressBar(pb, value = which(id==ids))
    sub_data = filter(MS_data, ID==id)
    eps = unique(sub_data$Episode_Identifier)
    real_episodes = unique(sub_data$Episode)
    N_eps = length(eps)
    if(N_eps > 1){ 
      for(i in 1:(N_eps-1)){
        for(j in (i+1):N_eps){
          sub_data_ep_ij = filter(sub_data, Episode_Identifier %in% c(eps[i], eps[j]))
          sub_data_ep_ij$ID_True = sub_data_ep_ij$ID
          # make new ID variable for the model to parse
          sub_data_ep_ij$ID = paste('TID',id, '%Primary%',i,'%Recurrence%', j,sep='')
          
          # This is a weird hack to turn the vector of episodes into 1 & 2 (if it was 4 & 7 for example)
          sub_data_ep_ij$Episode = as.numeric(as.factor(sub_data_ep_ij$Episode)) 
          # update Episode Identifier column
          sub_data_ep_ij$Episode_Identifier = apply(sub_data_ep_ij, 1, function(x) paste(x['ID'],'_',x['Episode'],sep=''))
          
          sub_data_ep_ij$First_EpNumber = real_episodes[i]
          sub_data_ep_ij$Second_EpNumber = real_episodes[j]
          inflated_data = rbind(inflated_data, sub_data_ep_ij)
        }
      }
    }
  }
  
  length(unique(inflated_data$ID))
  inflated_data = inflated_data[-1,]
  MS_data = inflated_data
  
  return(MS_data)
}


## This function returns an inflated MS_data data.frame which has 
# all possible pairwise comparisons between different individuals
Make_All_Pairwise_Comparisons = function(MS_data, ncores=6){
  # takes a while to make the dataset
  require(stringr)
  
  # get rid of duplicated rows for the episode IDs
  MS_summary = MS_data[!duplicated(MS_data$Episode_Identifier),]
  ep_ids = MS_summary$Episode_Identifier
  IDs = unique(MS_data$ID)
  
  # indexing all pairwise comparisons of episodes
  indices = expand.grid(1:(length(ep_ids)-1), 2:length(ep_ids))
  IDs = data.frame(ID1 = MS_summary$ID[indices$Var1],ID2 = MS_summary$ID[indices$Var2])
  ind1 = indices$Var1 < indices$Var2
  ind2 = IDs$ID1 != IDs$ID2
  
  # Select those with different ids and only in one order
  indices = indices[ind1&ind2,]
  writeLines(sprintf('There are %s across individual comparisons possible',nrow(indices)))
  library(doParallel)
  registerDoParallel(cores = ncores)
  inflated_data = foreach(s=1:nrow(indices), .combine = rbind) %dopar% {
    
    i = indices$Var1[s]
    j = indices$Var2[s]
    
    sub_data_ep_ij = dplyr::filter(MS_data, Episode_Identifier %in% c(ep_ids[i], ep_ids[j]))
    # This is a weird hack to turn the vector of episodes into 1..3
    sub_data_ep_ij$Episode[sub_data_ep_ij$Episode_Identifier==ep_ids[i]] = 1
    sub_data_ep_ij$Episode[sub_data_ep_ij$Episode_Identifier==ep_ids[j]] = 2
    sub_data_ep_ij$ID1 = unique(sub_data_ep_ij$ID)[1]
    sub_data_ep_ij$ID2 = unique(sub_data_ep_ij$ID)[2]
    
    sub_data_ep_ij$ID = paste('APC_EP1%',unique(sub_data_ep_ij$Episode_Identifier)[1], '%EP2%',
                              unique(sub_data_ep_ij$Episode_Identifier)[2], 
                              sep='')
    sub_data_ep_ij$Episode_Identifier = apply(sub_data_ep_ij, 1, function(x) paste(x['ID'],'_',x['Episode'],sep=''))
    
    #sub_data_ep_ij = check_MS_validity(sub_data_ep_ij)
    sub_data_ep_ij
    
  }
  return(inflated_data)
}

# this only works on pairs of episodes
# This function reduces the complexity in pairwise comparisons: could impact inference
check_MS_validity = function(MS_data){
  if(length(unique(MS_data$ID))>1) stop('Single ID MS data with 2 episodes only in check_MS_validity function')
  
  ep_1_ind = MS_data$Episode==1
  ep_2_ind = MS_data$Episode==2
  
  # sets to NA when no comparison data
  MS_ind = grep('PV.', colnames(MS_data))
  for(j in MS_ind){
    if(sum(!is.na(MS_data[ep_1_ind,j])) ==0 ){
      MS_data[ep_2_ind,j] = NA
    }
    if(sum(!is.na(MS_data[ep_2_ind,j])) ==0 ){
      MS_data[ep_1_ind,j] = NA
    }
  }
  
  # gets rid of NA only rows in the MS data
  ind_keep = apply(MS_data[,MS_ind], 1, function(x) sum(!is.na(x)) > 0)
  MS_data = MS_data[ind_keep,,drop=FALSE]
  
  # updates the MOI values
  if(nrow(MS_data)>0){
    ep_1_ind = MS_data$Episode==1
    ep_2_ind = MS_data$Episode==2
    MS_data$MOI_id[ep_1_ind] = 1:sum(ep_1_ind)
    MS_data$MOI_id[ep_2_ind] = 1:sum(ep_2_ind)
  }
  return(MS_data)
}