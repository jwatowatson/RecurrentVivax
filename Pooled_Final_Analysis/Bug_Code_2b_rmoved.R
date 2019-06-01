
rm(list = ls())
load('../RData/RPackages_List.RData')
new.pkg <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if(length(new.pkg) > 0){ # Install using .R script
  stop('Please run the script Install_and_load_required_packages.R before returning to this script. Thanks.')}else{
    sapply(pkgs, require, character.only = TRUE) # Load all packages
  }


load('../RData/TimingModel/MOD2_theta_estimates.RData')
load('../RData/TimingModel/MOD2_Posterior_samples.RData')
load('../RData/TimingModel/Combined_Time_Event.RData')

source("../Genetic_Model/Data_functions.R")
source('../Genetic_Model/iGraph_functions.R')
source("../Genetic_Model/post_prob_CLI.R") 
source("../Genetic_Model/test_Rn_compatible.R") 
source("../Genetic_Model/Data_Inflation_Functions.R")
source("../Genetic_Model/hap_combinations.R")

# The pooled MS data from BPD and VHX
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData')



MSs_all = c("PV.3.502","PV.3.27","PV.ms8",
            "PV.1.501","PV.ms1","PV.ms5",
            "PV.ms6","PV.ms7","PV.ms16")


# Prior weight for the Dirichlet (setting weight to 0 recovers empirical freq):
D_weight_Prior = 1

# These are motif lengths: for the plotting
MSs_Motifs = list("PV.3.502"=8,'PV.3.27' = 4, 
                  "PV.ms8" = 3, "PV.1.501"= 7, 
                  "PV.ms1" = 3, "PV.ms5" = 3,
                  "PV.ms6" = 3, "PV.ms16" =3,
                  "PV.ms7" = 3)

writeLines(paste('Number of episodes used to compute frequencies:',
                 sum(MS_pooled$Episode==1 & MS_pooled$MOI_id==1)))
Ind_Primary = which(MS_pooled$Episode==1)

# Posterior Dirichlet parameter vector
Alpha_Posteriors = apply(MS_pooled[,MSs_all], 2, function(x, Ind_Primary){
  # Extract xmax 
  xmax = max(x,na.rm=T)
  # prior parameter vector (iterpolates unobserved repeat lengths < xmax)
  param_vector = array(D_weight_Prior, dim = xmax, dimnames = list(1:xmax)) 
  # observed data summarised as counts
  obs_counts = table(x[Ind_Primary]) 
  # posterior parameter vector
  param_vector[names(obs_counts)] = param_vector[names(obs_counts)] + obs_counts
  return(param_vector)
})

# Calculate posterior mean 
Fs_Combined = sapply(Alpha_Posteriors, function(x){x/sum(x)})

# First we remove MS data for which there are no recurrent data
N_episodes_typed = table(MS_pooled$ID[!duplicated(MS_pooled$Episode_Identifier)])
MS_pooled = filter(MS_pooled, ID %in% names(N_episodes_typed[N_episodes_typed>1]))
load('../RData/GeneticModel/thetas_9MS.RData')

MS_pooled_summary = MS_pooled[!duplicated(MS_pooled$Episode_Identifier),] # Collapse rows due to COI > 2

ind_calculated = which(MS_pooled_summary$Episode_Identifier %in% thetas_9MS$Episode_Identifier)
IDs_calculated = unique(MS_pooled_summary$ID[ind_calculated])
IDs_remaining = unique(MS_pooled_summary$ID[! MS_pooled_summary$ID %in% IDs_calculated])
writeLines(sprintf('individuals with more than two recurrences: %s',length(IDs_remaining)))

MS_inflate = reformat_MSdata(filter(MS_pooled, ID %in% IDs_remaining), MSs = MSs_all)
MS_inflated = Inflate_into_pairs(MS_data = MS_inflate)

all_rec_eps_ind = !duplicated(MS_inflated$Episode_Identifier) & MS_inflated$Episode==2
all_rec_eps = MS_inflated$Episode_Identifier[all_rec_eps_ind]
P_matrix = data.frame(array(dim = c(length(all_rec_eps),4)))
colnames(P_matrix) = c('Episode_Identifier','C','I','L')
P_matrix$Episode_Identifier = all_rec_eps

# AT: not sure what to do with this line but Ksamples doesn't appear to feature elsewhere
# Ksamples = min(Ksamples, length(grep('C',colnames(Post_samples_matrix)))

K_results = sum(!duplicated(MS_inflated$Episode_Identifier[MS_inflated$Episode>1]))


# draw a random distribution over the allele frequencies from posterior
Fs_random = lapply(Alpha_Posteriors, rdirichlet, n=1)

# Name allele: important since we use the haplotype labels to 
# extract allele frequencies in Log_Pr_yn_Gab
for(i in 1:length(Fs_random)) {
  names(Fs_random[[i]]) = 1:length(Fs_random[[i]])
}

id = unique(MS_inflated$ID)[131]
id
Res = post_prob_CLI(MSdata = MS_inflated[MS_inflated$ID==id,], 
                    Fs = Fs_random, 
                    UpperComplexity = 10^6, 
                    verbose = T,
                    cores = 1)