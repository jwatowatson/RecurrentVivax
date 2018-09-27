### Taking Aimee's code and attempting to make a function from it ###
reformat_MSdata = function(MSdata, MSs){
  
  library(dplyr) # Needed for filter
  
  #==========================================================================
  # Reformat the data s.t. there are no NA gaps in mixed infections
  #==========================================================================
  MS_data_reformated = MSdata
  
  for(id in unique(MS_data_reformated$ID)){ # IDs over n = 1,..,N individuals
    
    yn = filter(MS_data_reformated, ID == id) # Extract data for the nth individual 
    cn = table(yn$Episode_Identifier) # Extract MOIs for each of the nth individual's infections

    if(any(cn > 1)){ 
      for(inf in names(which(cn > 1))){
        ynt = filter(yn, Episode_Identifier == inf) # Extract data for the tth infection
        for(moi_id in 2:max(ynt$MOI_id)){
          data_moi_minus = ynt[(ynt$MOI_id == moi_id - 1), MSs] # Data from previous MOI_id entry
          data_moi = ynt[ynt$MOI_id == moi_id, MSs] # Data from subsequent MOI_id entry
          na_ind = is.na(data_moi) # Missing markers in MOI_id > 1 entry
          row_ind = (MS_data_reformated$Episode_Identifier == inf) & (MS_data_reformated$MOI_id == moi_id)
          ynt[ynt$MOI_id == moi_id, MSs][na_ind] = data_moi_minus[,na_ind] # Filling the gaps in ynt for moi_id +1 
          MS_data_reformated[row_ind, MSs][na_ind] = unlist(data_moi_minus[,na_ind]) # Filling the gaps in reformatted data
        }
      }
    }
  }
  return(MS_data_reformated)
}