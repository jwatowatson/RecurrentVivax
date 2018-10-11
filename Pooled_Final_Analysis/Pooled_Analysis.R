## ---- echo=FALSE, include=FALSE------------------------------------------
#==========================================================================
# Set up
#==========================================================================
knitr::opts_chunk$set(cache = TRUE, cache.comments = FALSE, 
                      echo = FALSE, include = TRUE, 
                      fig.width = 7, fig.height = 7,
                      fig.pos = 'H', 
                      dev = 'png', dpi = 300)

require(plyr)
require(dplyr) # For filter
require(gtools) # For permutations
require(tictoc) # For bdiag()
require(doParallel) # For post_prob_R
require(Matrix)
require(matrixStats) # for logSumExp
library(igraph) # For igraph
require(RColorBrewer)
load('../RData/TimingModel/MOD3_theta_estimates.RData')
load('../RData/TimingModel/MOD3_Posterior_samples.RData')
load('../RData/TimingModel/Combined_Time_Event.RData')

source("../Genetic_Model/Data_functions.R")
source('../Genetic_Model/iGraph_functions.R')
source("../Genetic_Model/post_prob_CLI.R") 
source("../Genetic_Model/post_prob_CLI_sequential.R") 
source("../Genetic_Model/test_Rn_compatible.R") 
source("../Genetic_Model/Data_Inflation_Functions.R")
source('../Plotting_MS_Data/PlottingFunction.R')

# The pooled MS data from BPD and VHX
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData')


# This is an important parameter: pseudo weight in the Dirichlet prior
D_weight_Prior = 1

# Prior weight for the Dirichlet (setting weight to 0 recovers empirical freq):
MSs_all = c("PV.3.502","PV.3.27","PV.ms8",
            "PV.1.501","PV.ms1","PV.ms5",
            "PV.ms6", "PV.ms7","PV.ms16")

Fs_Combined =  apply(MS_pooled[,MSs_all], 2, function(x, Ind_Primary){
  # Extract xmax 
  xmax = max(x,na.rm=T)
  # prior parameter vector (iterpolates unobserved repeat lengths < xmax)
  param_vector = array(D_weight_Prior, dim = xmax, dimnames = list(1:xmax)) 
  # observed data summarised as counts
  obs_counts = table(x[Ind_Primary]) 
  # posterior parameter vector
  param_vector[names(obs_counts)] = param_vector[names(obs_counts)] + obs_counts
  # posterior mean
  posterior_mean = param_vector/sum(param_vector)
  return(posterior_mean)
})


# The pooled MS data from BPD and VHX
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData')
tic()
APC_MS_data = Make_All_Pairwise_Comparisons(MS_data = MS_pooled, ncores=42)
save(APC_MS_data, file = 'APC_MS_data.bigRData')
toc()
Inflated_Results = post_prob_CLI(MS_data = APC_MS_data, 
                                 Fs = Fs_Combined, 
                                 UpperComplexity = 10^6, 
                                 verbose = F,
                                 cores = 42)
save(Inflated_Results, file = 'Inflated_Results.bigRData')

