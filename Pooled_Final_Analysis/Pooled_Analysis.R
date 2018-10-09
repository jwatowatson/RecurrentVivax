## ---- echo=FALSE, include=FALSE------------------------------------------
#==========================================================================
# Set up
#==========================================================================
knitr::opts_chunk$set(cache = TRUE, cache.comments = FALSE, 
                      echo = TRUE, include = TRUE, 
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

# The pooled MS data from BPD and VHX
load('../RData/GeneticModel/MS_data_PooledAnalysis.RData')

RUN_MODELS = T
RUN_MODELS_CLUSTER = T
CREATE_PLOTS = T

## ------------------------------------------------------------------------
writeLines(sprintf('Number of individuals with at least one episode typed: %s',
                   length(unique(MS_pooled$ID))))
writeLines(sprintf('Number of episodes typed: %s',
                   length(unique(MS_pooled$Episode_Identifier))))
writeLines(sprintf('Number of recurrences typed: %s',
                   length(unique(MS_pooled$Episode_Identifier[MS_pooled$Episode>1]))))
# We also remove MS data for which there are no recurrent data
N_episodes_typed = table(MS_pooled$ID[!duplicated(MS_pooled$Episode_Identifier)])
MS_pooled = filter(MS_pooled, ID %in% names(N_episodes_typed[N_episodes_typed>1]))

writeLines(sprintf('Number of individuals with at least two episodes typed: %s',
                   length(unique(MS_pooled$ID))))
writeLines(sprintf('Number of recurrences typed: %s',
                   length(unique(MS_pooled$Episode_Identifier[MS_pooled$Episode>1]))))

## ------------------------------------------------------------------------
MSs_VHX = c("PV.3.502","PV.3.27","PV.ms8","PV.1.501","PV.ms1","PV.ms5","PV.ms6")
MSs_all = c("PV.3.502","PV.3.27","PV.ms8","PV.1.501","PV.ms1","PV.ms5","PV.ms6",
            "PV.ms7","PV.ms16")
MSs_BPD = MSs_all
MSs_Main = c('PV.3.27', 'PV.3.502', 'PV.ms8') # These are typed for all episodes (the core group)

## ---- echo=F-------------------------------------------------------------
#--------------------------------------------------------------------------
# Estimate allele frequencies from reformated data ignoring NAs.

# Prior weight for the Dirichlet (setting weight to 0 recovers empirical freq):
MSs_all = c("PV.3.502","PV.3.27","PV.ms8",
            "PV.1.501","PV.ms1","PV.ms5",
            "PV.ms6", "PV.ms7","PV.ms16")

# These are motif lengths: for the plotting
MSs_Motifs = list("PV.3.502"=8,'PV.3.27' = 4, 
                  "PV.ms8" = 3, "PV.1.501"= 7, 
                  "PV.ms1" = 3, "PV.ms5" = 3,
                  "PV.ms6" = 3, "PV.ms16" =3,
                  "PV.ms7" = 3)

# This is an important parameter: pseudo weight in the Dirichlet prior
D_weight_Prior = 1

writeLines(paste('Number of episodes used to compute frequencies:',
                 sum(MS_pooled$Episode==1 & MS_pooled$MOI_id==1)))
Ind_Primary = which(MS_pooled$Episode==1)

# I wrote the following to check I understood - suggest as a more-readable alternative (agrees with ms text)
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

## ----AlleleFrequencies, echo=F-------------------------------------------
if(CREATE_PLOTS){
  par(mfrow=c(3,3), las=1, bty='n', cex.axis=1.2)
  
  for(ms in MSs_all){ # As bars
    K = length(Fs_Combined[[ms]]) # Cardinality of ms 
    xs = rdirichlet(n = 1000, alpha = Alpha_Posteriors[[ms]]) # Sample from posterior
    # MC approximation of 0.9 percentile expressed as percentage
    YMAX = max(apply(100*xs, 2, quantile, probs = .9)) 
    # Number of observations
    N_MS = length(unique(MS_pooled$ID[MS_pooled$Episode==1 & !is.na(MS_pooled[,ms])]))
    plot(1:ncol(xs), rep(NA,K), 
         main = sprintf('%s (n = %s)', ms, N_MS),
         pch = 18, ylim=c(0,YMAX),col='red',
         ylab='%', xlab='', yaxt='n')
    abline(h=100/K)
    mtext(text = paste('Motif length:',MSs_Motifs[[ms]]),side = 1,line=3)
    for(k in 1:ncol(xs)){
      lines(rep(k,2), 100*quantile(xs[,k],probs = c(0.025,0.975)), col='blue')
    }
    points(1:length(Fs_Combined[[ms]]),
           100*Fs_Combined[[ms]],col='red',pch=18)
    
    axis(2, seq(0, round(YMAX), length.out = 3))
  }
}

## ------------------------------------------------------------------------
inds = grepl('mean_theta', colnames(Mod3_ThetaEstimates)) # Extract mean
Episode_Identifier = Mod3_ThetaEstimates$Episode_Identifier
p = data.frame(Episode_Identifier = Episode_Identifier, Mod3_ThetaEstimates[,inds],
               stringsAsFactors = F) # Reformat
colnames(p) = c('Episode_Identifier', 'C', 'L', 'I')

genetic_AND_time_data_eps = intersect(p$Episode_Identifier, MS_pooled$Episode_Identifier)
p = p[p$Episode_Identifier %in% genetic_AND_time_data_eps,]
Post_samples_matrix = Post_samples_matrix[Post_samples_matrix$Episode_Identifier %in% genetic_AND_time_data_eps,]

## ---- include=FALSE------------------------------------------------------
if(RUN_MODELS){
  #===============================================
  # Run new version (with time-to-event)
  #===============================================
  tic()
  thetas_9MS = post_prob_CLI(MS_data = MS_pooled, Fs = Fs_Combined, 
                             p = p, cores = 6, verbose = F) 
  thetas_9MS$Episode_Identifier = rownames(thetas_9MS)
  save(thetas_9MS, file = '../RData/GeneticModel/thetas_9MS.RData')
  toc()
  
  #===============================================
  # Run new version (without time-to-event)
  #===============================================
  tic()
  thetas_9MS_Tagnostic = post_prob_CLI(MS_data = MS_pooled, Fs = Fs_Combined, 
                                       cores = 6, verbose = F)
  thetas_9MS_Tagnostic$Episode_Identifier = rownames(thetas_9MS_Tagnostic)
  save(thetas_9MS_Tagnostic, file = '../RData/GeneticModel/thetas_9MS_Tagnostic.RData')
  toc()
} else {
  load('../RData/GeneticModel/thetas_9MS.RData')
  load('../RData/GeneticModel/thetas_9MS_Tagnostic.RData')
}

## ---- include=FALSE------------------------------------------------------
if(RUN_MODELS_CLUSTER){
  #===============================================
  # Run full Bayesian with sampling at random from time prior
  # takes about 220 seconds on 6 cores
  #===============================================
  registerDoParallel(cores = 42)
  Ksamples = length(grep('C',colnames(Post_samples_matrix)))
  tic()
  Thetas_full_post = foreach(ss = 1:Ksamples, .combine = cbind) %dopar% {
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
    thetas_9MS = post_prob_CLI(MS_data = MS_pooled, Fs = Fs_random, 
                               p = p, cores = 1, verbose = F) 
    thetas_9MS$Episode_Identifier = rownames(thetas_9MS)
    thetas_9MS
  }
  save(Thetas_full_post, file = '../RData/GeneticModel/Full_Posterior_Model_samples.RData')
  toc()
} else {
  load(file = '../RData/GeneticModel/Full_Posterior_Model_samples.RData')
}

## ------------------------------------------------------------------------
mycols = brewer.pal(n=3, name = 'Set1')
thetas_9MS = arrange(thetas_9MS, Episode_Identifier)
thetas_9MS_Tagnostic = arrange(thetas_9MS_Tagnostic, Episode_Identifier)

Time_Estimates_1 = filter(Mod3_ThetaEstimates, 
                          Episode_Identifier %in% thetas_9MS$Episode_Identifier)
Time_Estimates_1 = arrange(Time_Estimates_1, Episode_Identifier)

thetas_9MS$drug = Time_Estimates_1$arm_num
thetas_9MS_Tagnostic$drug = Time_Estimates_1$arm_num

# for plotting
thetas_9MS$drug_col = mapvalues(x = thetas_9MS$drug, 
                                c('AS','CHQ','CHQ/PMQ'), mycols)
thetas_9MS_Tagnostic$drug_col = mapvalues(x = thetas_9MS_Tagnostic$drug, 
                                          c('AS','CHQ','CHQ/PMQ'), mycols)


## ------------------------------------------------------------------------
BPD_data = Thetas_full_post[grep('BPD',rownames(Thetas_full_post)),]
Thetas_BPD = thetas_9MS[grep('BPD', thetas_9MS$Episode_Identifier),]


# Added by Aimee: some examples
# Colour some specific examples 
example_inds = grepl('_644_', Thetas_BPD$Episode_Identifier) | 
  grepl('BPD_598_', Thetas_BPD$Episode_Identifier) | 
  grepl('BPD_562_', Thetas_BPD$Episode_Identifier) |
  grepl('BPD_53_', Thetas_BPD$Episode_Identifier) 
example_ids = Thetas_BPD$Episode_Identifier[example_inds]
example_inds_times = MS_pooled$timeSinceLastEpisode[MS_pooled$Episode_Identifier %in% example_ids]
Tagnostic_example_inds = thetas_9MS_Tagnostic$Episode_Identifier %in% example_ids
thetas9MS_example_inds = thetas_9MS$Episode_Identifier %in% example_ids
Time1_example_inds = Time_Estimates_1$Episode_Identifier %in% example_ids


## ------------------------------------------------------------------------
if(CREATE_PLOTS){
  par(mfrow=c(2,2),las=1, bty='n')
  # Time agnostic versus full posterior 
  plot(log10(thetas_9MS_Tagnostic$L), log10(thetas_9MS$L), 
       col = thetas_9MS$drug_col, main = 'Relapse',pch=20,
       xlab = 'Time agnostic', ylab = 'Time included')
  lines(-10:0,-10:0)
  
  # Annotate by examples
  points(x = log10(thetas_9MS_Tagnostic$L[Tagnostic_example_inds]), 
         y = log10(thetas_9MS$L[thetas9MS_example_inds]), 
         pch=1, cex = 1.5, col='black')
  text(x = log10(thetas_9MS_Tagnostic$L[Tagnostic_example_inds]), 
       y = log10(thetas_9MS$L[thetas9MS_example_inds]), 
       labels = example_ids, pos = 2, cex = 0.7)
  
  
  plot(log10(thetas_9MS_Tagnostic$I), log10(thetas_9MS$I), 
       col=thetas_9MS$drug_col, main = 'Reinfection',pch=20,
       xlab = 'Time agnostic', ylab = 'Time included')
  lines(-20:0,-20:0)
  
  
  
  ##### Prior versus full posterior
  plot(log10(Time_Estimates_1$Relapse_mean_theta),
       log10(thetas_9MS$L),main = 'Relapse',
       col=thetas_9MS$drug_col,pch=20,
       xlab = 'Time based prior', ylab = 'Full posterior')
  lines(-10:10,-10:10)
  
  # Annotate by examples
  points(x = log10(Time_Estimates_1$Relapse_mean_theta[Time1_example_inds]), 
         y = log10(thetas_9MS$L[thetas9MS_example_inds]), 
         pch=1, cex = 1.5, col='black')
  text(x = log10(Time_Estimates_1$Relapse_mean_theta[Time1_example_inds]), 
       y = log10(thetas_9MS$L[thetas9MS_example_inds]), 
       labels = example_ids, pos = 2, cex = 0.7)
  
  plot(log10(Time_Estimates_1$ReInfection_mean_theta),
       log10(thetas_9MS$I),main = 'Reinfection',
       col=thetas_9MS$drug_col,pch=20,
       xlab = 'Time based prior', ylab = 'Full posterior')
  lines(-10:10,-10:10)
}

## ------------------------------------------------------------------------
if(CREATE_PLOTS){
  
  par(mfrow=c(1,2),las=1, bty='n')
  
  reLapse_ordered = sort.int(thetas_9MS$L, decreasing = TRUE, index.return = TRUE)
  plot(reLapse_ordered$x, pch=20, col = thetas_9MS$drug_col[reLapse_ordered$ix],
       xlab = 'Recurrence index', ylab = 'Probability of relapse state',
       main = 'Full posterior: reLapse')
  CI = cbind(apply(
    Thetas_full_post[reLapse_ordered$i,grep('L',colnames(Thetas_full_post)),],
    1, quantile, probs = 0.025), 
    apply(Thetas_full_post[reLapse_ordered$i,grep('L',colnames(Thetas_full_post)),],
          1, quantile, probs = 0.975))
  for(i in 1:length(reLapse_ordered$x)){
    if(diff(CI[i,]) > 0.005) arrows(i,CI[i,1],i,CI[i,2], 
                                    length = 0.02,angle = 90, 
                                    code = 3,
                                    col=thetas_9MS$drug_col[reLapse_ordered$ix[i]])
  }
  
  legend('topright',col = mycols[2:3], legend = c('No radical cure','Radical cure'),pch=20)
  
  reLapse_ordered_Tagn = sort.int(thetas_9MS_Tagnostic$L, 
                                  decreasing = TRUE, index.return = TRUE)
  plot(reLapse_ordered_Tagn$x, pch=20, cex=.8,
       col = thetas_9MS_Tagnostic$drug_col[reLapse_ordered_Tagn$ix],
       xlab = 'Recurrence index', ylab = 'Probability of relapse state',
       main = 'Time agnostic posterior: reLapse')
  legend('topright',col = mycols[2:3], 
         legend = c('No radical cure','Radical cure'),pch=20)
}

## ------------------------------------------------------------------------
if(CREATE_PLOTS){
  
  par(mfrow=c(1,2),las=1, bty='n')
  reinfection_ordered = sort.int(thetas_9MS$I, decreasing = TRUE, index.return = TRUE)
  plot(reinfection_ordered$x, pch=20, col = thetas_9MS$drug_col[reinfection_ordered$ix],
       xlab = 'Recurrence index', ylab = 'Probability of reinfection state',
       main = 'Full posterior: reInfection')
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
  
  reinfection_ordered_Tagn = sort.int(thetas_9MS_Tagnostic$I, decreasing = TRUE, index.return = TRUE)
  plot(reinfection_ordered_Tagn$x, pch=20, cex=.8,
       col = thetas_9MS_Tagnostic$drug_col[reinfection_ordered_Tagn$ix],
       xlab = 'Recurrence index', ylab = 'Probability of reinfection state',
       main = 'Time agnostic posterior: reInfection')
  legend('topright',col = mycols[2:3], legend = c('No radical cure','Radical cure'),pch=20)
}

## ------------------------------------------------------------------------
if(CREATE_PLOTS){
  par(mfrow=c(1,2),las=1, bty='n')
  recrud_ordered = sort.int(thetas_9MS$C, decreasing = TRUE, index.return = TRUE)
  plot(recrud_ordered$x, pch=20, col = thetas_9MS$drug_col[recrud_ordered$ix],
       xlab = 'Recurrence index', ylab = 'Probability of recrudescence state',
       main = 'Full posterior: reCrudescence')
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
  
  recrud_ordered_Tagn = sort.int(thetas_9MS_Tagnostic$C, decreasing = TRUE, index.return = TRUE)
  plot(recrud_ordered_Tagn$x, pch=20, cex=.8,
       col = thetas_9MS_Tagnostic$drug_col[recrud_ordered_Tagn$ix],
       xlab = 'Recurrence index', ylab = 'Probability of recrudescence state',
       main = 'Time agnostic posterior: reCrudescence')
  legend('topright',col = mycols[2:3], legend = c('No radical cure','Radical cure'),pch=20)
}

## ----BPD_efficacy--------------------------------------------------------
if(CREATE_PLOTS){
  MS_pooled = MS_pooled[!duplicated(MS_pooled$Episode_Identifier),]
  par(mfrow=c(1,2),las=1, bty='n')
  reLapse_ordered = sort.int(Thetas_BPD$L, decreasing = TRUE, index.return = TRUE)
  plot(reLapse_ordered$x, pch=20, col = Thetas_BPD$drug_col[reLapse_ordered$ix],
       xlab = 'Recurrence index', ylab = 'Probability of relapse',
       main = '')
  CI = cbind(apply(
    BPD_data[reLapse_ordered$ix,grep('L',colnames(BPD_data)),],
    1, quantile, probs = 0.025), 
    apply(BPD_data[reLapse_ordered$ix,grep('L',colnames(BPD_data)),],
          1, quantile, probs = 0.975))
  for(i in 1:length(reLapse_ordered$x)){
    if(diff(CI[i,]) > 0.005) arrows(i,CI[i,1],i,CI[i,2], 
                                    length = 0.02,angle = 90, 
                                    code = 3,
                                    col=Thetas_BPD$drug_col[reLapse_ordered$ix[i]])
  }
  # # Annotate by examples
  # points(reLapse_ordered$ix[example_inds],
  #        Thetas_BPD$L[example_inds], pch=1, cex = 1.5, col='black')
  # text(x = example_inds_times, y = Thetas_BPD$L[example_inds], labels = example_ids, pos = 3)
  
  
  writeLines(sprintf('The mean percentage of recurrences which are estimated to be relapses is %s%%',
                     round(100*sum(Thetas_BPD$L + Thetas_BPD$C)/nrow(Thetas_BPD))))
  
  
  plot(NA,NA,xlim=c(0,max(MS_pooled$timeSinceLastEpisode,na.rm=T)), ylim=c(0,1),
       ylab = 'Probability of relapse', xlab = 'Time since last episode')
  for(i in 1:length(reLapse_ordered$x)){
    kk = reLapse_ordered$ix[i]
    x_time = MS_pooled$timeSinceLastEpisode[Thetas_BPD$Episode_Identifier[kk]==
                                              MS_pooled$Episode_Identifier]
    points(x_time,
           Thetas_BPD$L[kk], pch=20, col=mycols[3])
    if(diff(CI[i,]) > 0.005) arrows(x_time,CI[i,1],x_time,CI[i,2], 
                                    length = 0.02,angle = 90, 
                                    code = 3,
                                    col=Thetas_BPD$drug_col[reLapse_ordered$ix[i]])
  }
  
  # # Annotate by examples
  # points(example_inds_times,Thetas_BPD$L[example_inds], pch=1, cex = 1.5, col='black')
  # text(x = example_inds_times, y = Thetas_BPD$L[example_inds], labels = example_ids, pos = 3)
}

## ------------------------------------------------------------------------
ind_calculated = which(MS_pooled$Episode_Identifier %in% thetas_9MS$Episode_Identifier)
IDs_calculated = unique(MS_pooled$ID[ind_calculated])
IDs_remaining = unique(MS_pooled$ID[! MS_pooled$ID %in% IDs_calculated])

## ---- include=FALSE------------------------------------------------------
MS_inflate = filter(MS_pooled, MS_pooled$ID %in% IDs_remaining)
MS_inflated = Inflate_into_pairs(MS_data = MS_inflate)

if(RUN_MODELS_CLUSTER){  
  
  all_rec_eps = unique(MS_inflated$Episode_Identifier[MS_inflated$Episode==2])
  P_matrix = data.frame(array(dim = c(length(all_rec_eps),4)))
  colnames(P_matrix) = c('Episode_Identifier','C','I','L')
  P_matrix$Episode_Identifier = all_rec_eps
  Ksamples = length(grep('C',colnames(Post_samples_matrix)))
  
  K_results = sum(!duplicated(MS_inflated$Episode_Identifier[MS_inflated$Episode>1]))
  Res_total = array(NA, dim = c(K_results, 3*Ksamples))
  
  tic()
  Res_total=foreach(ss = 1:Ksamples, .combine = cbind) %do% {
    
    # draw a random distribution over the allele frequencies from posterior
    Fs_random = lapply(Alpha_Posteriors, rdirichlet, n=1)
    # I'm not sure if this is needed - how are these called deep inside?
    for(i in 1:length(Fs_random)){
      names(Fs_random[[i]]) = 1:length(Fs_random[[i]])
    }
    # take the ss sample from the time prior
    rec_states = c('C','L','I')
    indices = c(grep(rec_states[1],colnames(Post_samples_matrix))[ss],
                grep(rec_states[2],colnames(Post_samples_matrix))[ss],
                grep(rec_states[3],colnames(Post_samples_matrix))[ss])
    for(ep in all_rec_eps){
      i = which(P_matrix$Episode_Identifier==ep)
      j = which(MS_inflated$Episode_Identifier==ep)[1]
      k = which(Post_samples_matrix$Episode_Identifier == paste(MS_inflated$ID_True[j],
                                                                MS_inflated$Second_EpNumber[j],
                                                                sep='_'))
      P_matrix[i,rec_states] = Post_samples_matrix[k,indices]
    }
    
    Res = post_prob_CLI(MS_data = MS_inflated, 
                        Fs = Fs_random, 
                        p = P_matrix,
                        UpperComplexity = 10^8, 
                        verbose = F,
                        cores = 42)
    Res
  }
  rownames(Res_total) = rownames(Res)
  save(Res_total, file = 'FullPosterior_INF.bigRData')
  toc()
} else {
  load('FullPosterior_INF.bigRData')
}
