##############################################################################
# To explore solutions to sim recrudescence frailty
# Timing based on real data set in 
# Generated_analyse_probabilistic_v_deterministic_inc.NA.IDs.R
##############################################################################
rm(list = ls())
load('../RData/RPackages_List.RData')
new.pkg <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(new.pkg) > 0){ # Install using .R script
  stop('Please run the script Install_and_load_required_packages.R before returning to this script. Thanks.')}else{
    sapply(pkgs, require, character.only = TRUE) # Load all packages
  }

source('../Genetic_Model/iGraph_functions.R')
source('../Genetic_Model/Likelihood_function.R')
source('../Genetic_Model/Data_functions.R')
source('../Genetic_Model/test_Rn_compatible.R')
source('../Genetic_Model/post_prob_CLI.R')
source('../Genetic_Model/hap_combinations_functions.R')
set.seed(1)


#====================================================================
# Does increasing Max_Hap_comb recover non-zero recrudescence: no...
# Recovered non-zero for one case only with M = 9; for zero with M = 12
# Meanwhile computation time increases linearly
# A.s. increasing Max_Hap_comb is not a solution to overcoming frailty 
# in highly complex simulated data
#====================================================================
load('SimulationOutputs/Sim_Genetic_Data/MS_Data_Cardinality13_M3_COIs1_1_Clone.RData') # load data
inds = sim_output$MS_data_sim$ID == 'SIM_1' # All probablistic because complicated data
max_hap_combs = rev(c(100, 300, 500, 800, 1000)) # When we allow more than 1000 with sim data, start to skip problems as too complex
tic.clearlog()
Results = t(sapply(max_hap_combs, function(x){
  tic(sprintf('prob_%s',x))
  TH = post_prob_CLI(MSdata = sim_output$MS_data_sim[inds,], Fs = sim_output$FS, 
                     Max_Hap_comb = x, Max_Hap_genotypes = 0)
  toc(log = T)
  TH 
}))

Time_msg = unlist(tic.log(format = T))
Times = as.numeric(do.call(rbind, strsplit(Time_msg, split = ' '))[,1])
plot(x = max_hap_combs, y = Times)




#====================================================================
# Does increasing Max_Hap_genotypes recover non-zero recrudescence: no
# Problem sim cases are all in excess of any reasonable Max_Hap_genotypes
# Moreover, setting Max_Hap_genotypes above 300 leads to problems with combinations
# It's actually faster when deterministic, regardless of Max_Hap_comb = 300 or 800 
# (since both P1)
#====================================================================
load('SimulationOutputs/Sim_Genetic_Data/MS_Data_Cardinality13_M6_COIs2_1_Clone.RData') # load data
inds = sim_output$MS_data_sim$ID == 'SIM_1' # Need a not so complicated example otherwise Det not poss
max_genotypes = c(0, 10, 100, 200, 300) # Cannot go over 300
tic.clearlog()
Results = t(sapply(max_genotypes, function(x){
  tic()
  TH = post_prob_CLI(MSdata = sim_output$MS_data_sim[inds,], Fs = sim_output$FS, verbose = T, 
                     Max_Hap_genotypes = x, Max_Hap_comb = 300) 
  toc(log = T)
  TH 
}))

Time_msg = unlist(tic.log(format = T))
Times = as.numeric(do.call(rbind, strsplit(Time_msg, split = ' '))[,1])
plot(x = max_genotypes, y = Times)



#=============================================
# Can we do a potential for clone check? This has 
# been aborted because too hard 
# 
# Nevertheless I keep notes here in case of future
# editions, esp. with different data types. 
# 
# This is hard because it requires infections to 
# be processed as pairs whereas currently
# we phase each infection in turn
#
# Currently we have for a pair of infections a way to 
# extract intersecting allels and check if haplotypes using 
# those allels have been sampled probabilistically
#
# This is not yet automated over multiple infections 
# Moveover, I'm not sure how we might make haplotype combinations 
# from intersecting haplotypes were non already sampled 
# e.g. 
# could sample from the intesecting haplotypes then 
# add additional markers so compatible with data observed
#
#=============================================
load('SimulationOutputs/Sim_Genetic_Data/MS_Data_Cardinality13_M9_COIs3_1_Clone.RData') # load data
MSs = names(sim_output$FS) # extract microsatellite names 
yns = dlply(sim_output$MS_data_sim, 'ID') # list of per-ID data
ynts = lapply(yns, function(yn) dlply(yn, 'Episode_Identifier')) # list of list of per-episode data



# Make a fake challenging data set by adding some extra 
# rows to a random example to stress test
ynts[["SIM_10"]]$SIM_10_2 = rbind(ynts[["SIM_10"]]$SIM_10_2, 
                                  ynts[["SIM_10"]]$SIM_10_2,  
                                  ynts[["SIM_10"]]$SIM_10_2)
ynts[["SIM_10"]]$SIM_10_2$MOI_id[2:3] = 2:3
ynts[["SIM_10"]]$SIM_10_2$MS1[2:3] = c(9,12)
ynts[["SIM_10"]]$SIM_10_2$MS2[2:3] = c(3,11)

ID = "SIM_10"

# Specifiy data set
MSdata = do.call(rbind, ynts[[ID]])
Fs = sim_output$FS

t = 2 # Specify a recurrence (in this case there is only one) 
ynts_pair = ynts[[ID]][c(t-1,t)] # Extract data across recurrence and preceding infection

# Define a function capable of extracting alleles (as) that are shared ACROSS infections 
extract_intersecting_as = function(z){ 
  sapply(MSs, function(MS){ # For each MS, 
    ynm = lapply(z, function(ynt){ynt[,MS]}) # Extract alleles ACROSS infection pairs
    # For each recurrence, check if any alleles intersect with the previous infection
    # Sequential comparisons are valid here since we're checking clones for evidence of recrudescence
    shared_as = lapply(2:length(ynm), function(t){
      intersect(x = ynm[[t-1]], y = ynm[[t]])
    })
  })
}

# Extract intersecting alleles
Clonal_intersect = extract_intersecting_as(z = ynts_pair)

# Expand intersecting alleles into haplotypes
Hnt_clonal_intersect = expand.grid(Clonal_intersect)

# This next bit relies on Hap_combinations which have to be 
# generated by hand within post_prob_CLI
if(!exists(Hap_combinations)){stop("Need to generate Hap_combinations inside post_prob_CLI")}

# Convert Hap_combinations into list of strings
All_hap_comb_strings = lapply(Hap_combinations, function(x){
  apply(do.call(rbind, x), 1, paste, collapse = '')})

# Extract intersecting clones as strings
h_intersects = apply(Hnt_clonal_intersect[,MSs], 1, paste, collapse = '')

# For each intersecting clone, check if its in any of the hap combinations
sapply(h_intersects, function(h_intersect){ 
  sapply(All_hap_comb_strings, function(x){h_intersect %in% x})
})


