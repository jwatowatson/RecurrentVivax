##############################################################################
# This script runs the genetic model on simulated data (takes a long time).
# 12-th job for 3 relationships and 2 cardinalities each with N_indivs = 500 
# took approx. 4 hrs with 8 cores
##############################################################################

rm(list = ls())
load('../RData/RPackages_List.RData')
new.pkg <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(new.pkg) > 0){install.packages(new.pkg)} # Install using .R script
sapply(pkgs, require, character.only = TRUE) # Load packages

source('../Genetic_Model/iGraph_functions.R')
source('../Genetic_Model/Likelihood_function.R')
source('../Genetic_Model/Data_functions.R')
source('../Genetic_Model/test_Rn_compatible.R')
source('../Genetic_Model/post_prob_CLI.R')
source('../Genetic_Model/hap_combinations_functions.R')

RUN_MODELS = T
#+++++++++++++++
CORES = parallel::detectCores()-2
#+++++++++++++++

set.seed(1) # For reproducibility
cardinalities = c(4,13) # Marker cardinalities (set to match min and mean of our panel)
Ms = seq(3,12,3) # Number of MS markers
COIs_1 = c(1,2) # COIs of primary episode (COMMENT)
COIs_2 = c(1,2) # COIs of recurrent episode
relationships = c('Sibling','Stranger','Clone') # Relationships of parasite
COI_c_max = 3 # Maximum COI we currently consider for cor

# Enumerate all combinations of COI complexity and numbers of markers
settings = expand.grid(Ms = Ms, COI_1 = COIs_1, COI_2 = COIs_2) # Enumerate all possible combinations of COIs and markers 
settings = settings[settings$COI_1+settings$COI_2<=COI_c_max,] # Remove any with cumulative COI greater than that feasible
settings = rbind(settings, expand.grid(Ms = Ms, COI_1 = 3, COI_2 = 1)) # Add settings for 3_1 example
settings$COI_pattern <- paste(settings$COI_1, settings$COI_2, sep = "_") # Store the COIs as a patter
rownames(settings) = NULL
jobs_cor = which(settings$COI_pattern != '3_1') # Indices of core simulations jobs
jobs_add = which(settings$COI_pattern == '3_1') # Indices of addn simulations jobs

#====================================================================
# Run model on well-specified data: core jobs only
# Cardinalities 4 and 13
#====================================================================
tic()
if(RUN_MODELS){
  for(cardinality in cardinalities){ # iterate over cardinality of markers
    for(relationship in relationships){ # iterative over simulation scenario 
      
      thetas_all = foreach(job = jobs_cor, .combine = rbind, # iterate over jobs parameter settings
                           .packages = c('dplyr','Matrix','gtools','igraph','matrixStats','doParallel')
      ) %do% { # parallisation happening inside the function
        
        # load well-specified data
        load(sprintf('SimulationOutputs/Sim_Genetic_Data/MS_Data_Cardinality%s_M%s_COIs%s_%s.RData', 
                     cardinality, settings$M[job], settings$COI_pattern[job], relationship))
        
        # #+++++++++++++++++++
        # Ind = sim_output$MS_data_sim$ID == "SIM_1"
        # #+++++++++++++++++++
        
        # Run the model
        TH = post_prob_CLI(MSdata = sim_output$MS_data_sim, Fs = sim_output$FS, cores = CORES) 
        TH$setting = job # Add setting number for plotting
        TH # return results
        
      }
      writeLines(sprintf('********** Done for %s, cardinality %s **********', relationship, cardinality))
      save(thetas_all, file = sprintf('SimulationOutputs/Sim_Genetic_Results/%s_%s.RData', relationship, cardinality))
    }
  }
}

#====================================================================
# Run models on erroneous data: cor jobs only
# Cardinality fixed to 13
#====================================================================
tic()
if(RUN_MODELS){
  for(relationship in relationships){ # iterative over simulation scenario 
    
    thetas_all = foreach(job = jobs_cor, .combine = rbind, # iterate over jobs parameter settings
                         .packages = c('dplyr','Matrix','gtools','igraph','matrixStats','doParallel')
    ) %do% { # parallisation happening inside the function
      
      # load miss-specified data
      load(sprintf('./SimulationOutputs/Sim_Genetic_Data/MS_Data_Cardinality13_M%s_COIs%s_%s_Error.RData', 
                   settings$M[job], settings$COI_pattern[job], relationship))
      
      # #+++++++++++++++++++
      # Ind = sim_output$MS_data_sim$ID == "SIM_1"
      # #+++++++++++++++++++
      
      # Run the model 
      TH = post_prob_CLI(MSdata = sim_output$MS_data_sim, Fs = sim_output$FS, cores = CORES) 
      TH$setting = job # Add setting number for plotting
      TH # return results
      
    }
    writeLines(sprintf('********** Done for %s **********', relationship))
    save(thetas_all, file = sprintf('./SimulationOutputs/Sim_Genetic_Results/%s_13_Error.RData', relationship))
  }
}
toc()


#====================================================================
# Run models on well-specified data: add jobs only
# Cardinality fixed to 13
#====================================================================
tic()
if(RUN_MODELS){
  for(relationship in relationships){ # iterative over simulation scenario 
    
    thetas_all = foreach(job = jobs_add, .combine = rbind, # iterate over jobs parameter settings
                         .packages = c('dplyr','Matrix','gtools','igraph','matrixStats','doParallel')
    ) %do% { # parallisation happening inside the function
      
      # load data
      load(sprintf('./SimulationOutputs/Sim_Genetic_Data/MS_Data_Cardinality13_M%s_COIs%s_%s.RData', 
                   settings$M[job], settings$COI_pattern[job], relationship))
      
      # #+++++++++++++++++++
      # Ind = sim_output$MS_data_sim$ID == "SIM_1"
      # #+++++++++++++++++++
      
      # Run the model 
      TH = post_prob_CLI(MSdata = sim_output$MS_data_sim, Fs = sim_output$FS, cores = CORES) 
      TH$setting = job # Add setting number for plotting
      TH # return results
      
    }
    writeLines(sprintf('********** Done for %s **********', relationship))
    save(thetas_all, file = sprintf('./SimulationOutputs/Sim_Genetic_Results/%s_13_COIs3_1.RData', relationship))
  }
}
toc()


