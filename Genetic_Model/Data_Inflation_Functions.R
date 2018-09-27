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
  pb = txtProgressBar(min=1, max = length(ids),style = 3)
  for(id in ids){
    setTxtProgressBar(pb, value = which(id==ids))
    sub_data = filter(MS_data, ID==id)
    eps = unique(sub_data$Episode_Identifier)
    N_eps = length(eps)
    if(N_eps > 1){ 
      for(i in 1:(N_eps-1)){
        for(j in (i+1):N_eps){
          sub_data_ep_ij = filter(sub_data, Episode_Identifier %in% c(eps[i], eps[j]))
          sub_data_ep_ij$ID_True = sub_data_ep_ij$ID
          sub_data_ep_ij$ID = paste('TID',id, '%Primary%',i,'%Recurrence%', j,sep='')
          sub_data_ep_ij$Episode_Identifier = apply(sub_data_ep_ij, 1, function(x) paste(x['ID'],'_',x['Episode'],sep=''))
          
          # This is a weird hack to turn the vector of episodes into 1..3
          sub_data_ep_ij$Episode = as.numeric(as.factor(sub_data_ep_ij$Episode)) 
          
          sub_data_ep_ij$timeSinceLastEpisode[sub_data_ep_ij$Episode==2] = diff(unique(sub_data_ep_ij$timeSinceEnrolment))
          sub_data_ep_ij$First_EpNumber = i
          sub_data_ep_ij$Second_EpNumber = j
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
# all possible pairwise comparisons
Make_All_Pairwise_Comparisons = function(MS_data){
  # takes a while to make the dataset
  require(stringr)
  
  ids = unique(MS_data$Episode_Identifier)
  IDs = MS_data$ID[!duplicated(MS_data$Episode_Identifier)]
  
  indices = expand.grid(1:(length(ids)-1), 2:length(ids))
  indices = indices[indices$Var1<indices$Var2,]
  library(doParallel)
  registerDoParallel(cores = 44)
  inflated_data = foreach(s=1:nrow(indices), .combine = rbind) %dopar% {
    i = indices[s,1]
    j = indices[s,2]
    sub_data_ep_ij = dplyr::filter(MS_data, Episode_Identifier %in% c(ids[i], ids[j]))
    # This is a weird hack to turn the vector of episodes into 1..3
    if(length(unique(sub_data_ep_ij$ID))>1){
      sub_data_ep_ij$Episode[sub_data_ep_ij$Episode_Identifier==ids[i]] = 1
      sub_data_ep_ij$Episode[sub_data_ep_ij$Episode_Identifier==ids[j]] = 2
      sub_data_ep_ij$ID1 = unique(sub_data_ep_ij$ID)[1]
      sub_data_ep_ij$ID2 = unique(sub_data_ep_ij$ID)[2]
      
      
      sub_data_ep_ij$ID = paste('EP1__',unique(sub_data_ep_ij$Episode_Identifier)[1], '___EP2_',
                                unique(sub_data_ep_ij$Episode_Identifier)[2], 
                                sep='')
      sub_data_ep_ij$Episode_Identifier = apply(sub_data_ep_ij, 1, function(x) paste(x['ID'],'__',x['Episode'],sep=''))
      
      sub_data_ep_ij
    } else {
      
      return(NA)
      
    }
    
  }
  ind = apply(inflated_data, 1, function(x) sum(is.na(x))) < 22
  MS_data = inflated_data[ind, ]
  return(MS_data)
}

