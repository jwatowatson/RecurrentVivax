####################################################################
# This is a script designed to check that the genetic model returns 
# the prior when its given no genetic data
# Note that because we replace all genetic data with NA, there is no
# need to iterate over cardinalities, relationships, different 
# individuals nor numbers of markers 
# (just need to specify to sim data, then replace by missing)
####################################################################

# Set up
rm(list = ls())
load('../RData/RPackages_List.RData')
new.pkg <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(new.pkg) > 0){ # Install using .R script
  stop('Please run the script Install_and_load_required_packages.R before returning to this script. Thanks.')}else{
    sapply(pkgs, require, character.only = TRUE) # Load all packages
  }
source('../Genetic_Model/post_prob_CLI.R')
source('../Genetic_Model/iGraph_functions.R')
source('../Genetic_Model/Likelihood_function.R')
source('../Genetic_Model/Data_functions.R')
source('../Genetic_Model/test_Rn_compatible.R')
source('../Genetic_Model/hap_combinations_functions.R')
source('BuildSimData.R')

RUN_MODELS = T # Set to true to run model

# Specify but irrelevant
M = 6 # Number of markers
N = 1 # Number of individuals
N_alleles = 13 # Number of alleles of markers
relatedness = 'Sibling' # Relationship

# Enumerated settings
COIs_1 = 1:3 # COIs of primary episode (COMMENT)
COIs_2 = 1:3 # COIs of recurrent episode
COIs_3 = 0:3 # COIs of the second recurrent episode
COI_c_max = 6 # Maximum COI we currently consider

settings = expand.grid(COIs_1, COIs_2, COIs_3) # Enumerate all possible combinations of COIs  
names(settings) = c('COI_1', 'COI_2', 'COI_3')
settings$COI_cum = rowSums(settings) # Add cumulative
settings = settings[settings$COI_cum<=COI_c_max,] # Remove any with cumulative greater than that feasible
settings = arrange(settings, COI_cum) # Arrange in order of increasing COI
settings$COI_pattern <- apply(settings, 1, function(x){paste(x[1:3], collapse = "_")}) # Store the COIs as a patter

jobs = nrow(settings) # Number of simulations jobs

if(RUN_MODELS){
  
  thetas_all = foreach(job = 1:jobs, .combine = rbind, # iterate over jobs parameter settings
                       .packages = c('dplyr','Matrix','gtools','igraph','matrixStats','doParallel')
  ) %do% { # parallisation happening inside the function
    
    # Simulated data
    if(settings$COI_3[job] == 0){
      Tn = 2
      COIs = c(settings$COI_1[job], settings$COI_2[job])
    } else {
      Tn = 3
      COIs = c(settings$COI_1[job], settings$COI_2[job],  settings$COI_3[job])
    }
    
    sim_output = BuildSimData(Tn, COIs, M, N, N_alleles, relatedness) # Simulated data 
    MSnames = paste0('MS',1:M,sep = '') # Reconstruct microsatellite names
    sim_output$MS_data_sim[,MSnames] = NA # Replace all data with NA 
    
    # Run the model on the data using default uniform prior over states
    TH = post_prob_CLI(MSdata = sim_output$MS_data_sim, Fs = sim_output$FS)
    TH$COI_pattern = settings$COI_pattern[job]
    
    TH$Episode_ID = rownames(TH) # Record the episode ID
    rownames(TH) = NULL # Replace rownames to prevent odd renaming due to dopar
    
    TH # return results
  }
  
  # Save 
  save(thetas_all, file ='SimulationOutputs/Sim_Genetic_Results/PriorReturnCheck.RData')
}


# Not that recrudescence has zero probability when diversity excceds
# that of previous infection, hence some COI patterns have zero probability
# of being a recrudescence under the model: 
load('SimulationOutputs/Sim_Genetic_Results/PriorReturnCheck.RData')
print(dlply(thetas_all, 'COI_pattern'))




