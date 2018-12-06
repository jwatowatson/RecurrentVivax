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
require(knitr) # for kable()
load('../RData/TimingModel/MOD2_theta_estimates.RData')
load('../RData/TimingModel/MOD2_Posterior_samples.RData')
load('../RData/TimingModel/Combined_Time_Event.RData')

source("../Genetic_Model/Data_functions.R")
source('../Genetic_Model/iGraph_functions.R')
source("../Genetic_Model/post_prob_CLI.R") 
source("../Genetic_Model/post_prob_CLI_sequential.R") 
source("../Genetic_Model/test_Rn_compatible.R") 
source("../Genetic_Model/Data_Inflation_Functions.R")
#source('../Plotting_MS_data/PlottingFunction.R')

# The pooled MS data from BPD and VHX
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData')

# Booleans that describe generation of Markdown file. 
# The large jobs should be done on the cluster

# runs the once model on those that can be directly computed
RUN_MODELS_SINGLE_SIMPLE = F           # can do on laptop (5 minutes)
# runs the model for 100 random draws from the Time model posterior on those that can be directly computed
RUN_MODELS_FULL_POSTERIOR_SIMPLE = F    # can do on laptop (2 hours)
# runs the model for 100 random draws from the Time model posterior on those that cannot be directly computed
RUN_MODELS_FULL_POSTERIOR_INFLATED = F  # cluster is better (10 hours)
# Does the false positive discovery estimation
RUN_MODELS_FALSE_POSITIVE = T           # only cluster!! (1 week)
# generates plots
CREATE_PLOTS = F

# Colour scheme
# Previous Set1 not colourblind friendly: display.brewer.all(colorblindFriendly = T)
Dark2 = brewer.pal(8, 'Dark2')
Set2 = brewer.pal(8, 'Set2')

# This needs moving to time analysis: 
drug_cols3 = array(Dark2[c(4,6,1)], dim = 3, dimnames = list(c('AS','CHQ','CHQ/PMQ'))) 
drug_cols2 = array(Dark2[c(2,2,1)], dim = 3, dimnames = list(c('AS','CHQ','CHQ/PMQ')))
drug_cols_light2 = sapply(drug_cols2, adjustcolor, alpha.f = 0.5)

# Vector of states
states = c(relapse = 'L', reinfection = 'I', recrudescence = 'C')

# James, there are some odd bugs something to do with !duplicated() and the tibble...
# https://stackoverflow.com/questions/39041115/fixing-a-multiple-warning-unknown-column
# Also, Thetas_full_post has redundant column names

## ---- echo=FALSE---------------------------------------------------------
writeLines(sprintf('Number of individuals with at least one episode typed: %s',
                   length(unique(MS_pooled$ID))))
writeLines(sprintf('Number of episodes typed: %s',
                   length(unique(MS_pooled$Episode_Identifier))))
writeLines(sprintf('Number of recurrences typed: %s',
                   length(unique(MS_pooled$Episode_Identifier[MS_pooled$Episode>1]))))

# Which drug arms have been typed:
writeLines('\nOverall in the dataset: breakdown by treatment group:')
table(MS_pooled$Treatment[!duplicated(MS_pooled$ID)])
writeLines('\nWithin VHX: breakdown by treatment group:')
table(MS_pooled$Treatment[!duplicated(MS_pooled$ID) & 
                            1:nrow(MS_pooled)%in% grep('VHX',MS_pooled$ID)])

writeLines(sprintf('\nFrom BPD trial there are %s individuals with total of %s episodes typed',
                   length(unique(MS_pooled$ID[grep('BPD',MS_pooled$ID)])),
                   length(unique(MS_pooled$Episode_Identifier[grep('BPD',MS_pooled$ID)]))))
writeLines(sprintf('From VHX trial there are %s individuals with total of %s episodes typed',
                   length(unique(MS_pooled$ID[grep('VHX',MS_pooled$ID)])),
                   length(unique(MS_pooled$Episode_Identifier[grep('VHX',MS_pooled$ID)]))))

## ------------------------------------------------------------------------
# Were all episodes typed if person was selected for genotyping? 
MS_pooled_summary = MS_pooled[!duplicated(MS_pooled$Episode_Identifier),] # Collapse rows due to COI > 2
All_VHX_epi_count = table(Combined_Time_Data$patientid[grep('VHX',Combined_Time_Data$patientid)])
All_BPD_epi_count = table(Combined_Time_Data$patientid[grep('BPD',Combined_Time_Data$patientid)])
Typ_VHX_epi_count =  table(MS_pooled_summary$ID[grep('VHX',MS_pooled_summary$Episode_Identifier)]) # Typed_VHX_epi_count
Typ_BPD_epi_count =  table(MS_pooled_summary$ID[grepl('BPD_',MS_pooled_summary$Episode_Identifier)]) # Typed_BPD_epi_count

# VHX data set
no_person_typed = All_VHX_epi_count[names(All_VHX_epi_count) %in% names(Typ_VHX_epi_count)][names(Typ_VHX_epi_count)]
no_typed_person_typed = Typ_VHX_epi_count 
X0 = sum(no_person_typed == no_typed_person_typed) # 90% completely typed
X1 = length(no_typed_person_typed)
X2 = range((no_person_typed - no_typed_person_typed)[no_person_typed != no_typed_person_typed]) 
writeLines(sprintf('Of %s of %s VHX individual/s selected for genotyping, %s to %s of their episodes were not typed',
                   X0,X1,X2[1],X2[2]))

# BPD data set
no_person_typed = All_BPD_epi_count[names(All_BPD_epi_count) %in% names(Typ_BPD_epi_count)][names(Typ_BPD_epi_count)]
no_typed_of_typed = Typ_BPD_epi_count
X0 = sum(no_person_typed == no_typed_of_typed)
X1 = length(no_typed_of_typed)
X2 = range((no_person_typed - no_typed_of_typed)[no_person_typed != no_typed_of_typed]) # missing 1 if missing
writeLines(sprintf('Of %s of %s BPD individual/s selected for genotyping, %s to %s of their episodes were not typed',
                   X0,X1,X2[1],X2[2]))

## ----VHX_typed_untyped---------------------------------------------------
# ****** Bug in recent version ******
# 
# 
# # Is there a significant difference between VHX that where genetically typed and not?
# All_VHX_rec_count = All_VHX_epi_count[All_VHX_epi_count > 1] # Condition on those that have one or more recurrence
# x1 = Typ_VHX_epi_count
# x2 = All_VHX_rec_count[!names(All_VHX_rec_count) %in% names(x1)] # Untyped_VHX_epi_count
# max_rec = max(All_VHX_rec_count)
# 
# plot_store = array(dim = c(2, max_rec-1), dimnames = list(c('Untyped','Typed'), 2:max_rec))
# plot_store['Untyped', names(table(x2))] = table(x2)
# plot_store['Typed', names(table(x1))] = table(x1)
# 
# barplot(plot_store, col = c('lightgray','red'), beside = T, density = c(100,25),
#         main='VHX subset', xlab = 'Number of episodes')
# legend('topright', legend = c('Untyped', 'Typed'), fill =c('lightgray', 'red'),
#        density = c(100,25), bty = 'n')
# 
# # Bootstap test of difference in median
# theta_hat = median(x1) - median(x2) # Median diff 2
# theta_hat_boot = sapply(1:5000, function(x){
#   x1boot = sample(x1, size = length(x1), replace = T)
#   x2boot = sample(x2, size = length(x2), replace = T)
#   median(x1boot) - median(x2boot)})
# se = sd(theta_hat_boot)
# Normal_95CI = c(theta_hat - 1.96*se, theta_hat + 1.96*se)  # Significant
# percentile_95CI = quantile(theta_hat_boot, probs = c(0.025, 0.975)) # Not significant
# pivotal_95CI = 2*theta_hat - rev(percentile_95CI) # Significant

## ----COIs_VHX_BPD--------------------------------------------------------
COIs = sapply(unique(MS_pooled$Episode_Identifier), function(x){
  ind = MS_pooled$Episode_Identifier == x
  max(MS_pooled$MOI_id[ind])
})
x1 = COIs[grepl('VHX', names(COIs))]
x2 = COIs[grepl('BPD', names(COIs))]
theta_hat = median(x1) - median(x2) # Median diff 0

COIs_store = array(0, dim = c(2, max(COIs)), dimnames = list(c('VHX','BPD'), 1:max(COIs)))
COIs_store['VHX', names(table(x1))] = table(x1)
COIs_store['BPD', names(table(x2))] = table(x2)
barplot(COIs_store, col = c('lightgray','red'), beside = T, density = c(100,25), 
        main='COIs', xlab = 'COIs')
legend('topright', legend = c('VHX', 'BPD'), fill =c('lightgray', 'red'), 
       density = c(100,25), bty = 'n')

# How many could have have sibs? 
X0 = sum(COIs >= 3)
X1 = round(sum(COIs >= 3)/length(COIs)*100, 2)

writeLines(sprintf('Median COI in VHX and BPD: %s and %s, respectively', median(x1), median(x2)))
writeLines(sprintf('Number of episodes with COI >= 3: %s of %s (%s percent)', X0, length(COIs), X1))

## ------------------------------------------------------------------------
MSs_all = c("PV.3.502","PV.3.27","PV.ms8","PV.1.501","PV.ms1","PV.ms5","PV.ms6",
            "PV.ms7","PV.ms16")
MSs_Main = c('PV.3.27', 'PV.3.502', 'PV.ms8') # These are typed for all episodes (the core group)

## ---- echo=F-------------------------------------------------------------
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

## ----AlleleFrequencies, echo=F-------------------------------------------
# Changed s.t. frequences rather than percentages plotted
# No need for colour

if(CREATE_PLOTS){
  par(mfrow=c(3,3), las=1, bty='n', cex.axis=1.2, family = 'serif')
  for(ms in MSs_all){ # As bars
    K = length(Fs_Combined[[ms]]) # Cardinality of ms 
    xs = rdirichlet(n = 1000, alpha = Alpha_Posteriors[[ms]]) # Sample from posterior
    # MC approximation of 0.975 percentile expressed as percentage
    YMAX = max(apply(xs, 2, quantile, probs = .975)) # Round up to nearest decimal
    # Number of observations
    N_MS = length(unique(MS_pooled$ID[MS_pooled$Episode==1 & !is.na(MS_pooled[,ms])]))
    plot(NULL, panel.first = grid(),
         main = sprintf('%s', ms), font.main = 3, 
         pch = 18, ylim=c(0,YMAX), xlim = c(1,ncol(xs)), 
         ylab='', xlab='', yaxt='n')
    abline(h=1/K, lty = 'dashed', lwd = 0.5)
    title(ylab='Frequency', line = 3.25)
    title(xlab= sprintf('Repeat length (motif length = %s)', MSs_Motifs[[ms]]), line = 2.5)
    mtext(side = 3, text = sprintf('(n = %s)',N_MS), line = 0.5, cex = 0.65, adj = 1)
    for(k in 1:ncol(xs)){
      lines(rep(k,2), quantile(xs[,k],probs = c(0.025,0.975)))
    }
    points(1:length(Fs_Combined[[ms]]),Fs_Combined[[ms]],pch=18) # Point estimates: mean 
    axis(2, round(seq(0, round(YMAX,2), length.out = 3),2))
  }
}

## ------------------------------------------------------------------------

# # We also remove MS data for which there are no recurrent data
N_episodes_typed = table(MS_pooled$ID[!duplicated(MS_pooled$Episode_Identifier)])
MS_pooled = filter(MS_pooled, ID %in% names(N_episodes_typed[N_episodes_typed>1]))
# recreate pooled summary dataset
MS_pooled_summary = MS_pooled[!duplicated(MS_pooled$Episode_Identifier),] 

writeLines(sprintf('Number of individuals with at least two episodes typed: %s',
                   length(unique(MS_pooled$ID))))
writeLines(sprintf('Number of episodes in individuals with at least two episodes: %s',
                   length(unique(MS_pooled$Episode_Identifier))))
writeLines(sprintf('Number of recurrences typed: %s',
                   length(unique(MS_pooled$Episode_Identifier[MS_pooled$Episode>1]))))

## ------------------------------------------------------------------------
inds = grepl('mean_theta', colnames(Mod2_ThetaEstimates)) # Extract mean
p = data.frame(Episode_Identifier = Mod2_ThetaEstimates$Episode_Identifier, 
               Mod2_ThetaEstimates[,inds],
               stringsAsFactors = F) # Reformat
colnames(p) = gsub(pattern = 'Recrudescence_mean_theta',replacement = 'C',x = colnames(p))
colnames(p) = gsub(pattern = 'Relapse_mean_theta',replacement = 'L',x = colnames(p))
colnames(p) = gsub(pattern = 'ReInfection_mean_theta',replacement = 'I',x = colnames(p))

genetic_AND_time_data_eps = intersect(p$Episode_Identifier, MS_pooled$Episode_Identifier)
p = p[p$Episode_Identifier %in% genetic_AND_time_data_eps,] 
# Only need priors for those with genetic data
#Post_samples_matrix = Post_samples_matrix[Post_samples_matrix$Episode_Identifier %in% genetic_AND_time_data_eps,]

## ---- include=FALSE------------------------------------------------------
if(RUN_MODELS_SINGLE_SIMPLE){
  #===============================================
  # Run (with time-to-event)
  #===============================================
  # Approx 100 secs per full model run
  tic()
  thetas_9MS = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, 
                             p = p, cores = 42, verbose = F) 
  thetas_9MS$Episode_Identifier = rownames(thetas_9MS)
  save(thetas_9MS, file = '../RData/GeneticModel/thetas_9MS.RData')
  toc()
  
  #===============================================
  # Run (without time-to-event)
  #===============================================
  tic()
  thetas_9MS_Tagnostic = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_Combined, 
                                       cores = 42, verbose = F)
  thetas_9MS_Tagnostic$Episode_Identifier = rownames(thetas_9MS_Tagnostic)
  save(thetas_9MS_Tagnostic, file = '../RData/GeneticModel/thetas_9MS_Tagnostic.RData')
  toc()
} else {
  load('../RData/GeneticModel/thetas_9MS.RData')
  load('../RData/GeneticModel/thetas_9MS_Tagnostic.RData')
}

## ---- include=FALSE------------------------------------------------------
if(RUN_MODELS_FULL_POSTERIOR_SIMPLE){
  #===============================================
  # Run full Bayesian with sampling at random from time prior (treats time posterior as discrete uniform prior)
  # takes about 1.1 hour on 6 cores
  # The outer loop is not parallelised - suboptimal coding
  #===============================================
  
  Ksamples = max(length(grep('C',colnames(Post_samples_matrix))), 100)
  tic()
  Thetas_full_post = foreach(ss = 1:Ksamples, .combine = cbind) %do% {
    
    # draw a random distribution over the allele frequencies from posterior
    Fs_random = lapply(Alpha_Posteriors, rdirichlet, n=1)
    # I'm not sure if this is needed - how are these called deep inside?
    for(i in 1:length(Fs_random)){
      names(Fs_random[[i]]) = 1:length(Fs_random[[i]])
    }
    # take the ss sample from the time prior
    indices = c(grep('C',colnames(Post_samples_matrix))[ss],
                grep('L',colnames(Post_samples_matrix))[ss],
                grep('I',colnames(Post_samples_matrix))[ss],
                grep('Episode_Identifier',colnames(Post_samples_matrix)))
    p = Post_samples_matrix[,indices]
    # run the model
    thetas_9MS = post_prob_CLI(MSdata = MS_pooled, Fs = Fs_random, 
                               p = p, cores = 42, verbose = F) 
    thetas_9MS$Episode_Identifier = rownames(thetas_9MS)
    thetas_9MS
  }
  
  save(Thetas_full_post, file = '../RData/GeneticModel/Full_Posterior_Model_samples.RData')
  toc()
} else {
  load(file = '../RData/GeneticModel/Full_Posterior_Model_samples.RData')
}

## ------------------------------------------------------------------------
# Output of time-to-event model (sorted by episode number s.t. columns correspond)
Time_Estimates_1 = filter(Mod2_ThetaEstimates, Episode_Identifier %in% thetas_9MS$Episode_Identifier)
Time_Estimates_1 = arrange(Time_Estimates_1, Episode_Identifier)

# Outputs of genetic model w/wo time prior 
# sorted by episode number s.t. columns correspond and drug added
thetas_9MS = arrange(thetas_9MS, Episode_Identifier)
thetas_9MS_Tagnostic = arrange(thetas_9MS_Tagnostic, Episode_Identifier)
thetas_9MS$drug = Time_Estimates_1$arm_num # Add drug
thetas_9MS_Tagnostic$drug = Time_Estimates_1$arm_num # Add drug


# Extract BPD only for BPD only plots
BPD_data = Thetas_full_post[grep('BPD',rownames(Thetas_full_post)),]
Thetas_BPD = thetas_9MS[grep('BPD', thetas_9MS$Episode_Identifier),]

# Extract prior used in absence of time-to-event
Time_agnostic_p = as.list(formals(post_prob_CLI)$p)

## ------------------------------------------------------------------------
# Highlight some specific examples (useful for slides when presenting work)
Hightlight_BPD_examples = T 
if(Hightlight_BPD_examples){
  example_inds = grepl('BPD_644_', Thetas_BPD$Episode_Identifier) | 
    grepl('BPD_598_', Thetas_BPD$Episode_Identifier) | 
    grepl('BPD_562_', Thetas_BPD$Episode_Identifier) |
    grepl('BPD_53_', Thetas_BPD$Episode_Identifier) 
  example_ids = Thetas_BPD$Episode_Identifier[example_inds]
  example_inds_times = MS_pooled_summary$timeSinceLastEpisode[MS_pooled_summary$Episode_Identifier %in% example_ids]
  Tagnostic_example_inds = thetas_9MS_Tagnostic$Episode_Identifier %in% example_ids
  thetas9MS_example_inds = thetas_9MS$Episode_Identifier %in% example_ids
  Time1_example_inds = Time_Estimates_1$Episode_Identifier %in% example_ids
}


## ------------------------------------------------------------------------
if(RUN_MODELS_FALSE_POSITIVE){
  writeLines('Starting the FDR computations')
  # check if the massive pairwise dataset has been made, if not make it 
  # (takes a long time ~20hours)
  if(!"APC_MSdata.bigRData"%in%list.files()){
    # The pooled MS data from BPD and VHX
    writeLines('Making the Inflated FDR dataset')
    load('../RData/GeneticModel/MS_data_PooledAnalysis.RData')
    tic()
    APC_MSdata = Make_All_Pairwise_Comparisons(MS_data = MS_pooled, ncores=42)
    save(APC_MSdata, file = 'APC_MSdata.bigRData')
    toc()
  } 
  load('APC_MSdata.bigRData')
  print('The inflated pairwise dataset is available, now running the analysis...')
  # Run the genetic model on the pairwise data
  tic()
  writeLines('Running the model on inflated dataset')
  Inflated_Results = post_prob_CLI(MSdata = APC_MSdata, 
                                   Fs = Fs_Combined, 
                                   UpperComplexity = 10^6, 
                                   verbose = F,
                                   cores = 42)
  toc()
  save(Inflated_Results, file = 'Inflated_Results.bigRData')
} else {
  load('Inflated_Results.bigRData')
  Inflated_Results = Inflated_Results[!is.na(Inflated_Results$L),]
  load('APC_MSdata.bigRData')
}

## ------------------------------------------------------------------------
round(100*sum(Inflated_Results$L>0.5)/nrow(Inflated_Results),1)

## ------------------------------------------------------------------------
#BPD_data = Combined_Time_Data[grep('BPD', Combined_Time_Data$patientid),]
# hist(table(BPD_data$patientid),breaks = seq(0.5,3.5,by=1),main='',
#      xlab = '', xaxt='n')
# axis(1, at = 1:3, labels = c('Enrollment\nonly','1\nrecurrence','2\nrecurrences'),
#      tick = FALSE)

Combined_Time_Data$Episode_Identifier = apply(Combined_Time_Data,1,
                                              function(x){
                                                paste(x['patientid'],as.integer(x['episode']),
                                                      sep='_')} )

# iterate over every episode and use either the joint posterior 
# or if missing the time probability (this could be time censored probability)
Combined_Time_Data$Reinfection_Probability=
  Combined_Time_Data$Reinfection_Probability_LL=
  Combined_Time_Data$Reinfection_Probability_UL = NA

Mod2_ThetaEstimates$Failure_Identifier = 
  apply(Mod2_ThetaEstimates, 1, 
        function(x) paste(x['patientid'], as.integer(x['episode'])-1,sep='_'))

sss=0
for(i in 1:nrow(Combined_Time_Data)){
  ep_id = Combined_Time_Data$Episode_Identifier[i]
  MS_id = paste(Combined_Time_Data$patientid[i],
                as.integer(Combined_Time_Data$episode[i])+1, sep='_')
  # If in MS_final then use full probability
  if(MS_id %in% MS_final$Episode_Identifier){
    Combined_Time_Data$Reinfection_Probability[i] =
      MS_final$I_mean[MS_final$Episode_Identifier==MS_id]
    Combined_Time_Data$Reinfection_Probability_UL[i] =
      MS_final$I_upper[MS_final$Episode_Identifier==MS_id]
    Combined_Time_Data$Reinfection_Probability_LL[i] =
      MS_final$I_lower[MS_final$Episode_Identifier==MS_id]
  } else { # use the time to event model
    ind = which(Mod2_ThetaEstimates$Failure_Identifier==ep_id)
    if(length(ind)>0){
      Combined_Time_Data$Reinfection_Probability[i] =
        Mod2_ThetaEstimates$ReInfection_mean_theta[ind]
      Combined_Time_Data$Reinfection_Probability_UL[i] =
        Mod2_ThetaEstimates$ReInfection_975_theta[ind]
      Combined_Time_Data$Reinfection_Probability_LL[i] =
        Mod2_ThetaEstimates$ReInfection_025_theta[ind]
      sss=sss+1
    }
  }
}

## ------------------------------------------------------------------------
Combined_Time_Data = arrange(Combined_Time_Data, patientid, episode)
load('../RData/PK_data/BPD_pk.RData')
BPD_pk = filter(BPD_pk, !is.na(Episode))
Combined_Time_Data$log10_carboxyPMQ = NA
Combined_Time_Data$log10_PMQ = NA
Combined_Time_Data$NumberDaysPMQ = 14
for(i in 1:nrow(Combined_Time_Data)){
  id = Combined_Time_Data$patientid[i]
  ep_i = Combined_Time_Data$episode[i]
  all_id_eps = Combined_Time_Data$episode[Combined_Time_Data$patientid==id]
  pk_ind = which(BPD_pk$ID == id & BPD_pk$Episode==ep_i)
  if(length(pk_ind)>0){
    if(length(pk_ind)>1) print(id)
    Combined_Time_Data$log10_carboxyPMQ[i] = mean(BPD_pk$log10_carboxyPQ_PK[pk_ind])
    Combined_Time_Data$log10_PMQ[i] = mean(BPD_pk$log10_PQ_PK[pk_ind])
    Combined_Time_Data$NumberDaysPMQ[i] = BPD_pk$NumberofPKDays[pk_ind[1]]
  }
}

## ----CarboxyPredictionFailure--------------------------------------------
# These are two outliers - have discussed with Cindy
BPD444_recurrences = Combined_Time_Data$patientid=='BPD_44' & Combined_Time_Data$episode>1
BPD_598 = which(Combined_Time_Data$patientid=='BPD_598')
ind_keep = !BPD444_recurrences #& !BPD_598
require(lme4)
Combined_Time_Data$Failure_YN = Combined_Time_Data$Reinfection_Probability < 0.5
mod = glmer(Failure_YN ~ log10_carboxyPMQ + NumberDaysPMQ + 
              (1 | patientid), 
            family = 'binomial', data=Combined_Time_Data[ind_keep,])
summary(mod)
# Plot the data and model
xs=seq(0,4,by=.01)
par(las = 1, bty='n')
regimen_colors = brewer.pal(8, 'Dark2')[c(1,5)]
plot(Combined_Time_Data$log10_carboxyPMQ[ind_keep]*Combined_Time_Data$NumberDaysPMQ[ind_keep],
     jitter(as.numeric(Combined_Time_Data$Failure_YN[ind_keep]), factor = 0.02), 
     col = regimen_colors[mapvalues(Combined_Time_Data$NumberDaysPMQ[ind_keep],c(7,14),c(1,2))],
     pch = mapvalues(Combined_Time_Data$NumberDaysPMQ[ind_keep],c(7,14),c(3,4)),
     xlab = 'Carboxy-primaquine exposure: days * log(ng/mL)',
     panel.first = grid(), ylab = 'Probability of failure')
legend(x = 40, y = 0.3, bty='n', col =regimen_colors, 
       pch=c(3,4), legend = c('7 days','14 days'))
lines(xs*7, predict(mod, newdata=data.frame(log10_carboxyPMQ=xs, NumberDaysPMQ=7,
                                            patientid='new'),allow.new.levels=T,
                    type='response'), lwd=2, col= regimen_colors[1], lty=2)
lines(xs*14, predict(mod, newdata=data.frame(log10_carboxyPMQ=xs, NumberDaysPMQ=14,
                                             patientid='new'),allow.new.levels=T,
                     type='response'), lwd=2, col= regimen_colors[2], lty = 2)
points(Combined_Time_Data$log10_carboxyPMQ[BPD_598]*Combined_Time_Data$NumberDaysPMQ[BPD_598],
       Combined_Time_Data$Failure_YN[BPD_598], cex=2)
text(Combined_Time_Data$log10_carboxyPMQ[BPD_598[1]]*
       Combined_Time_Data$NumberDaysPMQ[BPD_598[1]]+3,
     Combined_Time_Data$Failure_YN[BPD_598[1]]-0.05, labels = 'BPD_598')

## ------------------------------------------------------------------------
mu_hat_7 = mean(Combined_Time_Data$log10_carboxyPMQ[Combined_Time_Data$NumberDaysPMQ==7])
sd_hat_7 = sd(Combined_Time_Data$log10_carboxyPMQ[Combined_Time_Data$NumberDaysPMQ==7])
mu_hat_14 = mean(Combined_Time_Data$log10_carboxyPMQ[Combined_Time_Data$NumberDaysPMQ==14],na.rm=T)
sd_hat_14 = sd(Combined_Time_Data$log10_carboxyPMQ[Combined_Time_Data$NumberDaysPMQ==14],na.rm=T)

par(mfrow=c(1,2))
hist(Combined_Time_Data$log10_carboxyPMQ[Combined_Time_Data$NumberDaysPMQ==7],
     main='', xlab='carboxy primaquine day 7 (ng/mL)')
abline(v = mu_hat_7 - sd_hat_7*3)

hist(Combined_Time_Data$log10_carboxyPMQ[Combined_Time_Data$NumberDaysPMQ==14],
     main = '', xlab='carboxy primaquine day 7 (ng/mL)')
abline(v = mu_hat_14 - sd_hat_14*3)

outliers7 = Combined_Time_Data$NumberDaysPMQ==7 & Combined_Time_Data$log10_carboxyPMQ<mu_hat_7 - sd_hat_7*3
outliers14 = Combined_Time_Data$NumberDaysPMQ==14 & Combined_Time_Data$log10_carboxyPMQ<mu_hat_14 - sd_hat_14*3

mod_No_Outliers = glmer(Failure_YN ~ log10_carboxyPMQ + NumberDaysPMQ + 
                          (1 | patientid), 
                        family = 'binomial', data=Combined_Time_Data[ind_keep & !outliers14 & !outliers7,])
summary(mod_No_Outliers)

## ----CarboxyPredictionFailure_NoOutliers---------------------------------
par(las = 1, bty='n')
plot(Combined_Time_Data$log10_carboxyPMQ[ind_keep]*Combined_Time_Data$NumberDaysPMQ[ind_keep],
     jitter(as.numeric(Combined_Time_Data$Failure_YN[ind_keep]),factor = 0.03), 
     col = regimen_colors[mapvalues(Combined_Time_Data$NumberDaysPMQ[ind_keep],c(7,14),1:2)],
     pch = mapvalues(Combined_Time_Data$NumberDaysPMQ[ind_keep],c(7,14),c(3,4)),
     xlab = 'Carboxy-primaquine trough exposure: days * log(ng/mL)',
     ylab = 'Probability of failure', panel.first = grid())
legend(x = 25, y = 0.7, bty='n', col = c(1, regimen_colors, 1,1),
       pch=c(1,3,4,NA,NA), lty = c(NA,NA,NA,1,2),lwd=c(NA,NA,NA,2,2),
       legend = c('Outlier','7 days','14 days',
                  'Regression fit: all data', 'Regression fit: outliers removed'))
lines(xs*7, predict(mod, newdata=data.frame(log10_carboxyPMQ=xs, NumberDaysPMQ=7,
                                            patientid='new'),allow.new.levels=T,
                    type='response'), lwd=2, col= regimen_colors[1], lty=1)
lines(xs*14, predict(mod, newdata=data.frame(log10_carboxyPMQ=xs, NumberDaysPMQ=14,
                                             patientid='new'),allow.new.levels=T,
                     type='response'), lwd=2, col= regimen_colors[2], lty = 1)

lines(xs*7, predict(mod_No_Outliers, newdata=data.frame(log10_carboxyPMQ=xs, NumberDaysPMQ=7,
                                                        patientid='new'),allow.new.levels=T,
                    type='response'), lwd=2, col= regimen_colors[1], lty=2)
lines(xs*14, predict(mod_No_Outliers, newdata=data.frame(log10_carboxyPMQ=xs, NumberDaysPMQ=14,
                                                         patientid='new'),allow.new.levels=T,
                     type='response'), lwd=2, col= regimen_colors[2], lty = 2)

# outline the outliers
points(Combined_Time_Data$log10_carboxyPMQ[outliers14|outliers7]*Combined_Time_Data$NumberDaysPMQ[outliers14|outliers7],
       Combined_Time_Data$Failure_YN[outliers14|outliers7], cex=2)

## ------------------------------------------------------------------------
# now we calculate the primaquine failure rate
# For individuals with two episodes: P(failure) = 1 - P(Rec 1 = I)*P(Rec 2 = I)
Summary_data = Combined_Time_Data[!duplicated(Combined_Time_Data$patientid),]
Summary_data$Failure_UL = Summary_data$Failure_LL = 
  Summary_data$Failure = Summary_data$CPMQ = 
  Summary_data$CPMQ = NA
for(i in 1:nrow(Summary_data)){
  ind = which(Combined_Time_Data$patientid==Summary_data$patientid[i])
  Summary_data$Failure[i] = 1-prod(Combined_Time_Data$Reinfection_Probability[ind],na.rm=T)
  Summary_data$Failure_UL[i] = 1-prod(Combined_Time_Data$Reinfection_Probability_UL[ind],na.rm=T)
  Summary_data$Failure_LL[i] = 1-prod(Combined_Time_Data$Reinfection_Probability_LL[ind],na.rm=T)
  Summary_data$CPMQ[i] = median(Combined_Time_Data$log10_carboxyPMQ[ind],na.rm=T)
}
BPD_data = Summary_data[grep('BPD', Summary_data$patientid),]

P_Failure=100*sum(BPD_data$Failure)/nrow(BPD_data)
# invert the intervals here - optimistic for not failure = pessimistic for failure
P_Failure_UL = 100*sum(BPD_data$Failure_LL)/nrow(BPD_data)
P_Failure_LL = 100*sum(BPD_data$Failure_UL)/nrow(BPD_data)

writeLines(sprintf('The primaquine failure rate in the %s individuals is %s%% (%s-%s) over the course of %s years total follow-up.',
                   nrow(BPD_data), round(P_Failure,2),
                   round(P_Failure_LL,2),
                   round(P_Failure_UL,2), round(sum(BPD_data$FU_time)/365)))

## ------------------------------------------------------------------------
TwoD6_dat = read.csv('../RData/PK_data/TwoD6&Vivax Genotyping_ASscore.csv')
TwoD6_dat$ID = apply(TwoD6_dat, 1, function(x) paste(x['Study'],
                                                     as.integer(x['Patient.ID']),
                                                     sep = '_'))
TwoD6_dat$Phenotype = mapvalues(TwoD6_dat$X2D6.Phenotype, 
                                from = c('PM','IM', 'EM'), to = 1:3)
Combined_Time_Data$Phenotype = Combined_Time_Data$ASscore = NA
for(i in 1:nrow(Combined_Time_Data)){
  id = Combined_Time_Data$patientid[i]
  if(sum(TwoD6_dat$ID==id)>0){
    Combined_Time_Data$ASscore[i] = TwoD6_dat$AS.score[TwoD6_dat$ID==id]
    Combined_Time_Data$Phenotype[i] = TwoD6_dat$Phenotype[TwoD6_dat$ID==id]
  }
}

mod_2D6 = lmer(log10_carboxyPMQ ~ ASscore + NumberDaysPMQ + (1 | patientid), 
               data = Combined_Time_Data)
summary(mod_2D6)
plot(Combined_Time_Data$ASscore, Combined_Time_Data$log10_carboxyPMQ, pch=20)
lines(xs, predict(mod_2D6, data.frame(ASscore=xs,NumberDaysPMQ=7,patientid='new'),allow.new.levels=T),lwd=3)

## ------------------------------------------------------------------------
Combined_2D6data = filter(Combined_Time_Data, !is.na(ASscore), !Censored)
for(id in unique(Combined_2D6data$patientid)){
  ind = Combined_2D6data$patientid==id
  Combined_2D6data$Failure_YN[ind] = max(Combined_2D6data$Failure_YN[ind])
}
Combined_2D6data = Combined_2D6data[!duplicated(Combined_2D6data$patientid),]
mod_Failure = glm(Failure_YN ~ ASscore, 
                  data = Combined_2D6data, family = 'binomial')
summary(mod_Failure)

