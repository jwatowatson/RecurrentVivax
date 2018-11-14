---
title: "Time-to-event model of vivax recurrence"
author: James Watson & Aimee Taylor
output:
  html_document:
    fig_caption: yes
    keep_md: yes
---


# Preliminaries

Load R packages.


```
## Warning: package 'dplyr' was built under R version 3.4.4
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:plyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```
## Warning: package 'rstan' was built under R version 3.4.4
```

```
## Loading required package: ggplot2
```

```
## Warning: package 'ggplot2' was built under R version 3.4.4
```

```
## Loading required package: StanHeaders
```

```
## Warning: package 'StanHeaders' was built under R version 3.4.4
```

```
## rstan (Version 2.18.1, GitRev: 2e1f913d3ca3)
```

```
## For execution on a local, multicore CPU with excess RAM we recommend calling
## options(mc.cores = parallel::detectCores()).
## To avoid recompilation of unchanged Stan programs, we recommend calling
## rstan_options(auto_write = TRUE)
```

```
## Loading required package: gdata
```

```
## gdata: read.xls support for 'XLS' (Excel 97-2004) files ENABLED.
```

```
## 
```

```
## gdata: read.xls support for 'XLSX' (Excel 2007+) files ENABLED.
```

```
## 
## Attaching package: 'gdata'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, first, last
```

```
## The following object is masked from 'package:stats':
## 
##     nobs
```

```
## The following object is masked from 'package:utils':
## 
##     object.size
```

```
## The following object is masked from 'package:base':
## 
##     startsWith
```


# Load Stan models


## Model 1

* ReInfections occur at 'random' (exponential distribution). A random effects term is used to adjust for inter-individual variability in propensity to be reinfected.

* Recrudescences can happen in the first couple of weeks. Relapses happen as described by a Weibull with blood-stage drug dependent parameters.

* Primaquine is assumed to have 100% efficacy.






## Model 2: Mixture of three components


* ReInfections occur at 'random' (exponential distribution).
A random effects term is used to adjust for inter-individual variability in propensity to be reinfected.

* Recrudescences can happen in the first couple of weeks. 

* Relapses are broken into two components. The fast/periodic relapse component happens as described by a Weibull with blood-stage drug dependent parameters (same as in Model 1). The slow/random relapse component is described by an exponential disrtibution.

* Primaquine is assumed to have 100% efficacy.





## Model 3: Mixture of three components and reLapses after PMQ

* ReInfections occur at 'random' (exponential distribution).
A random effects term is used to adjust for inter-individual variability in propensity to be reinfected.

* Recrudescences can happen in the first couple of weeks. 

* Relapses are broken into two components (same as in Model 2).

* Primaquine is not assumed to have 100% efficacy and a random effects term describes the propensity to relapse after primaquine.





# VHX and BPD combined dataset



```r
load('../RData/TimingModel/Combined_Time_Event.RData')
# Get rid of the very short durations 
Combined_Time_Data = filter(Combined_Time_Data, Time_to_event > 5)
```

```
## Warning: package 'bindrcpp' was built under R version 3.4.4
```

```r
# deal with non-integer ids
Combined_Time_Data$patientid = factor(Combined_Time_Data$patientid)
# create a mapping from the factor to integer
map = mapLevels(x=Combined_Time_Data$patientid)
# Convert ID to integer
Combined_Time_Data$ID = as.integer(Combined_Time_Data$patientid)
Combined_Time_Data = arrange(Combined_Time_Data, ID, episode)
# We sort the ids so that ids go from 1...N
ids = unique(Combined_Time_Data$ID)
N = as.integer(length(ids))

# Create a vector that maps the kth episode to the correct individual (for storing in the stan code)
rank_ids = 1:length(ids)
Combined_Time_Data$ordered_ids = unlist(sapply(1:length(ids),
                                               function(x) rep(rank_ids[x], 
                                                               sum(Combined_Time_Data$ID==ids[x])))
)
Combined_Time_Data$N_index = as.integer(Combined_Time_Data$ordered_ids)

# Create a vector that maps the nth person to 0 (always recived Primaquine) or
noPMQ_ind = which(Combined_Time_Data$arm_num != "CHQ/PMQ")
N_noPMQ = length(unique(Combined_Time_Data$ID[noPMQ_ind]))
Combined_Time_Data$N_noPMQ_index = 0
Combined_Time_Data$N_noPMQ_index[noPMQ_ind] = as.integer(as.factor(Combined_Time_Data$ID[noPMQ_ind]))
N_noPMQ_index = Combined_Time_Data$N_noPMQ_index[!duplicated(Combined_Time_Data$ID)]

# Create a vector that maps the nth person to 0 or index
PMQ_ind = which(Combined_Time_Data$arm_num == "CHQ/PMQ")
N_PMQ = length(unique(Combined_Time_Data$ID[PMQ_ind]))
Combined_Time_Data$N_PMQ_index = 0
Combined_Time_Data$N_PMQ_index[PMQ_ind] = as.integer(as.factor(Combined_Time_Data$ID[PMQ_ind]))
N_PMQ_index = Combined_Time_Data$N_PMQ_index[!duplicated(Combined_Time_Data$ID)]

# Turn drug into a numeric vector
numeric_drug = as.integer(revalue(Combined_Time_Data$arm_num,
                                  c('AS'='0','CHQ'='1','CHQ/PMQ'='2')))
ind = !duplicated(Combined_Time_Data$patientid)
drug_received = numeric_drug[ind]
```

# Prior specification


```r
# The hierachical parameters defining the prior distributions for model 1
Prior_params_M1 = list(mu_inv_lambda = 400,
                       sigma_inv_lambda = 50,
                       mu_AS_shape = 2,
                       sigma_AS_shape = 1,
                       mu_AS_scale = 25,
                       sigma_AS_scale = 5,
                       mu_CQ_shape = 2,
                       sigma_CQ_shape = 1,
                       mu_CQ_scale = 42,
                       sigma_CQ_scale = 5,
                       Hyper_logit_mean_p = -3,
                       Hyper_logit_sd_p = .2,
                       Hyper_logit_c1_mean = -3,
                       Hyper_logit_c1_sd = .25,
                       Hyper_logit_exp_p = 1)
# Model 2 has the same parameters with a few extra
Prior_params_M2 = c(Prior_params_M1, 
                    Early_L_logit_mean = 0,
                    Early_L_logit_sd = 1,
                    mu_inv_gamma = 80,
                    sigma_inv_gamma = 10)
# Model 3: extra parameters
Prior_params_M3 = c(Prior_params_M2, 
                    Hyper_logit_mean_p_PMQ = 3,
                    Hyper_logit_sd_p_PMQ = .25)
```

# Run Model

Set up parameters for the MCMC runs.

```r
# Choose as many chains as available cores
Chains = 6
options(mc.cores = Chains)
IT = 4*10^5
WarmUp = .5*IT
thin = 10^4

# put the data into stan format
# For models 1 & 2
Time_data_1_2 = list(N         = N,
                     #Number of individuals
                     Neps      = as.integer(nrow(Combined_Time_Data)),
                     #Number of durations
                     N_noPMQ   = N_noPMQ,
                     # Number of individuals who do not recieve PMQ
                     Durations = as.double(Combined_Time_Data$Time_to_event),
                     #Time to reinfection or time to censoring
                     Censored  = as.integer(Combined_Time_Data$Censored),
                     #If the duration is right censored or not
                     Drug      = numeric_drug,
                     # drug coded as an integer
                     N_index   = Combined_Time_Data$N_index,
                     # the index of the individual for each time interval
                     N_noPMQ_index = N_noPMQ_index
                     # the index mapping no PMQ individuals to their rank
)
# For model 3
Time_data_3 =list(N         = N,
                  #Number of individuals
                  Neps      = as.integer(nrow(Combined_Time_Data)),
                  #Number of durations
                  N_noPMQ   = N_noPMQ,
                  # Number of individuals who do not recieve PMQ
                  N_PMQ     = N_PMQ,
                  # Number of individuals who do not recieve PMQ
                  Durations = as.double(Combined_Time_Data$Time_to_event),
                  #Time to reinfection or time to censoring
                  Censored  = as.integer(Combined_Time_Data$Censored),
                  #If the duration is right censored or not
                  Drug      = numeric_drug,
                  # drug coded as an integer
                  N_index   = Combined_Time_Data$N_index,
                  # the index of the individual for each time interval
                  N_noPMQ_index = N_noPMQ_index,
                  # the index mapping no PMQ individuals to their rank
                  N_PMQ_index = N_PMQ_index
                  # the index mapping PMQ individuals to their rank
)
```

Run with X parallel chains. This depends on local computing power.


```r
if(RUN_MODELS){
  mod1_RE = sampling(Timing_Model1_RE,
                     data = c(Time_data_1_2, Prior_params_M1),
                     iter = IT, warmup = WarmUp,
                     chains=Chains, thin = thin)
  save(mod1_RE, file = 'OutputResults/StanModels_mod1.RData')
  
  mod2_RE = sampling(Timing_Model2_RE,
                     data = c(Time_data_1_2, Prior_params_M2),
                     iter = IT, warmup = WarmUp,
                     chains=Chains, thin = thin)
  save(mod2_RE, file = 'OutputResults/StanModels_mod2.RData')
  
  mod3_RE = sampling(Timing_Model3_RE,
                     data = c(Time_data_3, Prior_params_M3),
                     iter = IT, warmup = WarmUp,
                     chains=Chains, thin = thin)
  save(mod3_RE, file = 'OutputResults/StanModels_mod3.RData')
}
load('OutputResults/StanModels_mod1.RData')
load('OutputResults/StanModels_mod2.RData')
load('OutputResults/StanModels_mod3.RData')
```

# Plot output


Let's make some nice colors for the plotting

```r
# Colour scheme
# Previous Set1 not colourblind friendly: display.brewer.all(colorblindFriendly = T)
Dark2 = brewer.pal(8, 'Dark2')
Set2 = brewer.pal(8, 'Set2')

drug_cols3 = array(Dark2[c(4,6,1)], dim = 3, dimnames = list(c('AS','CHQ','CHQ/PMQ'))) 
drug_cols2 = array(Dark2[c(2,2,1)], dim = 3, dimnames = list(c('AS','CHQ','CHQ/PMQ')))
drug_cols_light2 = array(Set2[c(2,2,1)], dim = 3, dimnames = list(c('AS','CHQ','CHQ/PMQ')))

# Vector of states
states = c(relapse = 'L', reinfection = 'I', recrudescence = 'C')

#mycols = brewer.pal(n=3, name = 'Set1')
# Do we include censored time intervals in the plots:
PLOT_Censored_Obs = F

if(PLOT_Censored_Obs){
  ind_plotting = 1:nrow(Combined_Time_Data)
} else {
  ind_plotting = !Combined_Time_Data$Censored
}

# We jitter the time to event (days) to improve plotting
# This alters the saved output used in the genetic model... 
# Combined_Time_Data$Time_to_event = jitter(Combined_Time_Data$Time_to_event, amount = 1)
```


## Model 1


```r
# # Traceplots
# traceplot(mod1_RE,c('AS_shape', 'CQ_shape', 'AS_scale', 'CQ_scale'))
# traceplot(mod1_RE, c('inv_lambda','logit_c1','Recrud_shape','Recrud_scale'))
# traceplot(mod1_RE, c('logit_mean_p','logit_sd_p'))
# 
# # Extract samples
# thetas = extract(mod1_RE)
# 
# par(las=1, mfrow=c(1,2))
# hist(thetas$inv_lambda, main='Mean time to reinfection',
#      xlab='1/lambda (days)')
# hist(inv.logit(apply(thetas$logit_p,2,mean)),
#      main = 'proportion of reInfections', xlab='p')
# 
# par(las=1, mfrow=c(1,1))
# # Plot the outcome of the predicted labels
# #****************************** Reinfection ******************************
# labels = extract(mod1_RE, 'prob_labels')$prob_labels
# mean_labels_Reinfection = apply(labels[,ind_plotting,1,drop=T], 2, mean)
# plot(Combined_Time_Data$Time_to_event[ind_plotting], log10(mean_labels_Reinfection),
#      col = mycols[numeric_drug[ind_plotting]+1],
#      pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
#      ylab='Probability of ReInfection', yaxt='n',xaxt='n',
#      xlab='Months from last episode', xlim=c(0,400))
# axis(1, at = seq(0, 420, by=60), labels = seq(0, 420, by=60)/30)
# axis(2, at = -2:0, labels= 10^(-2:0))
# axis(2, at = log10(seq(.1,1,by=.1)), labels = NA)
# axis(2, at = log10(seq(.01,.1,by=.01)), labels = NA)
# axis(2, at = log10(seq(.001,.01,by=.001)), labels = NA)
# legend('bottomright',legend = c('AS','CQ','CQ+PMQ'),
#        col=c(mycols),pch = c(rep(1,3)), bty='n',lwd=2,lty=NA)
# 
# 
# #****************************** Recrudescence ****************************
# mean_labels_ReCrud = apply(labels[,ind_plotting,3,drop=T], 2, mean)
# plot(Combined_Time_Data$Time_to_event[ind_plotting], mean_labels_ReCrud,
#      col = mycols[numeric_drug[ind_plotting]+1], xaxt='n',
#      pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
#      ylab='ReCrudescence',
#      xlab='Weeks from last episode', xlim=c(0,30))
# axis(1, at = seq(0,28,by=7), labels = seq(0,28,by=7)/7)
# 
# 
# #****************************** Relapse ****************************
# mean_labels_ReLap = apply(labels[,ind_plotting,2,drop=T], 2, mean)
# plot(Combined_Time_Data$Time_to_event[ind_plotting], log10(mean_labels_ReLap),
#      col = mycols[numeric_drug[ind_plotting]+1], xaxt='n',
#      pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
#      ylab='log10 Probability of ReLapse',yaxt='n',
#      xlab='Months from last episode', xlim=c(0,100),
#      ylim = c(-10,0))
# axis(1, at = seq(0, 90, by= 30), labels = seq(0, 90, by=30)/30)
# axis(2, at = -2:0, labels= 10^(-2:0))
```

## Model 2


```r
# par(las=1)
# traceplot(mod2_RE,c('AS_shape', 'CQ_shape', 'AS_scale', 'CQ_scale'))
# traceplot(mod2_RE, c('inv_lambda','logit_c1','Recrud_shape','Recrud_scale'))
# traceplot(mod2_RE, c('logit_mean_p','logit_sd_p','logit_EarlyL'))
# par(mfrow=c(1,2))
# thetas = extract(mod2_RE)
# hist(thetas$inv_lambda, main='Mean time to reinfection', xlab='1/lambda (days)')
# hist(thetas$inv_gamma, main='Mean time to late reLapse', xlab='1/gamma (days)')
# 
# par(las=1, mfrow=c(1,3))
# 
# hist(inv.logit(apply(thetas$logit_p,2,mean)), xlab = 'Probability reinfection', main = '')
# hist(inv.logit(thetas$logit_c1), xlab = 'proportion recrudescence', main='')
# hist(inv.logit(thetas$logit_EarlyL), xlab = 'proportion early Relapse', main='')
# 
# par(las=1, mfrow=c(1,1))
# # Plot the outcome of the predicted labels
# #****************************** Reinfection ******************************
# labels = extract(mod2_RE, 'prob_labels')$prob_labels
# mean_labels_Reinfection = apply(labels[,ind_plotting,1,drop=T], 2, mean)
# plot(Combined_Time_Data$Time_to_event[ind_plotting], log10(mean_labels_Reinfection),
#      col = mycols[numeric_drug[ind_plotting]+1],
#      pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
#      ylab='Probability of ReInfection', yaxt='n',xaxt='n',
#      xlab='Months from last episode', xlim=c(0,400))
# axis(1, at = seq(0, 420, by=60), labels = seq(0, 420, by=60)/30)
# axis(2, at = -2:0, labels= 10^(-2:0))
# axis(2, at = log10(seq(.1,1,by=.1)), labels = NA)
# axis(2, at = log10(seq(.01,.1,by=.01)), labels = NA)
# axis(2, at = log10(seq(.001,.01,by=.001)), labels = NA)
# legend('bottomright',legend = c('Artesunate','Chloroquine','Chloroquine+\nPrimaquine'),
#        col=c(mycols),pch = rep(1,3), bty='n',lwd=2,lty=NA)
# 
# 
# #****************************** Recrudescence ****************************
# mean_labels_ReCrud = apply(labels[,ind_plotting,4,drop=T], 2, mean)
# plot(Combined_Time_Data$Time_to_event[ind_plotting], mean_labels_ReCrud,
#      col = mycols[numeric_drug[ind_plotting]+1], xaxt='n',
#      pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
#      ylab='ReCrudescence',
#      xlab='Weeks from last episode', xlim=c(0,60))
# axis(1, at = seq(0,54,by=7), labels = seq(0,54,by=7)/7)
# 
# 
# #****************************** Relapse ****************************
# mean_labels_ReLap1 = apply(labels[,ind_plotting,2,drop=T], 2, mean)
# mean_labels_ReLap2 = apply(labels[,ind_plotting,3,drop=T], 2, mean)
# mean_labels_ReLap = mean_labels_ReLap1 + mean_labels_ReLap2
# plot(Combined_Time_Data$Time_to_event[ind_plotting], mean_labels_ReLap,
#      col = mycols[numeric_drug[ind_plotting]+1], xaxt='n',
#      pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
#      ylab='log10 Probability of ReLapse',
#      xlab='Months from last episode')
# axis(1, at = seq(0, 420, by= 60), labels = seq(0, 420, by=60)/30)
```

## Model 3


```r
traceplot(mod3_RE,c('AS_shape', 'CQ_shape', 'AS_scale', 'CQ_scale'))
```

![](TimingModelStan_files/figure-html/plotModel3-1.png)<!-- -->

```r
traceplot(mod3_RE, c('inv_lambda','Recrud_shape','Recrud_scale'))
```

![](TimingModelStan_files/figure-html/plotModel3-2.png)<!-- -->

```r
traceplot(mod3_RE, c('logit_c1_AS','logit_c1_CQ','logit_c1_CQ_PMQ'))
```

![](TimingModelStan_files/figure-html/plotModel3-3.png)<!-- -->

```r
traceplot(mod3_RE, c('logit_mean_p','logit_sd_p','logit_EarlyL'))
```

![](TimingModelStan_files/figure-html/plotModel3-4.png)<!-- -->

```r
par(mfrow=c(1,2))
thetas = extract(mod3_RE)
hist(thetas$inv_lambda, xlab='Mean time to reinfection (days)', main='')
hist(thetas$inv_gamma, xlab='Mean time to late reLapse (days)', main='')
```

![](TimingModelStan_files/figure-html/plotModel3-5.png)<!-- -->

```r
hist(inv.logit(apply(thetas$logit_p,2,mean)),
     xlab = 'Reinfection: no PMQ', main = '',
     yaxt='n',ylab='')

hist(inv.logit(apply(thetas$logit_p_PMQ,2,mean)),
     xlab = 'Reinfection: PMQ', main = '',
     yaxt='n',ylab='')
```

![](TimingModelStan_files/figure-html/plotModel3-6.png)<!-- -->

```r
# Recrudescence weights
plot(density(100*inv.logit(thetas$logit_c1_AS)), col = drug_cols3['AS'],lwd=3,
     xlab = 'Recrudescence (%)', main='',yaxt='n',ylab='', xlim=c(0,10))
lines(density(100*inv.logit(thetas$logit_c1_CQ)), col = drug_cols3['CHQ'],lwd=3)
lines(density(100*inv.logit(thetas$logit_c1_CQ_PMQ)), col = drug_cols3['CHQ/PMQ'], lwd=3)
legend('topright',col=drug_cols3, legend = c('Artesunate monotherapy','Chloroquine monotherapy','Chloroquine+Primaquine'),lwd=3)

hist(inv.logit(thetas$logit_EarlyL), xlab = 'Early Relapse',
     main='',yaxt='n',ylab='')
```

![](TimingModelStan_files/figure-html/plotModel3-7.png)<!-- -->

```r
par(las=1, mfrow=c(1,1))
# Plot the outcome of the predicted labels
#****************************** Reinfection ******************************
labels = extract(mod3_RE, 'prob_labels')$prob_labels
mean_labels_Reinfection = apply(labels[,ind_plotting,1,drop=T], 2, mean)
plot(Combined_Time_Data$Time_to_event[ind_plotting], log10(mean_labels_Reinfection),
     col = drug_cols3[numeric_drug[ind_plotting]+1],
     pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
     ylab='Probability of ReInfection', yaxt='n',xaxt='n',
     xlab='Months from last episode', xlim=c(0,400))
axis(1, at = seq(0, 420, by=60), labels = seq(0, 420, by=60)/30)
axis(2, at = -2:0, labels= 10^(-2:0))
axis(2, at = log10(seq(.1,1,by=.1)), labels = NA)
axis(2, at = log10(seq(.01,.1,by=.01)), labels = NA)
axis(2, at = log10(seq(.001,.01,by=.001)), labels = NA)
legend('bottomright',legend = c('Artesunate','Chloroquine','Chloroquine\nPrimaquine'),
       col=c(drug_cols3),pch = rep(1,3), bty='n',lwd=2,lty=NA)
```

![](TimingModelStan_files/figure-html/plotModel3-8.png)<!-- -->

```r
#****************************** Recrudescence ****************************
mean_labels_ReCrud = apply(labels[,ind_plotting,4,drop=T], 2, mean)
plot(Combined_Time_Data$Time_to_event[ind_plotting], mean_labels_ReCrud,
     col = drug_cols3[numeric_drug[ind_plotting]+1], xaxt='n',
     pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
     ylab='ReCrudescence',
     xlab='Weeks from last episode', xlim=c(0,60))
axis(1, at = seq(0,54,by=7), labels = seq(0,54,by=7)/7)
```

![](TimingModelStan_files/figure-html/plotModel3-9.png)<!-- -->

```r
#****************************** Relapse ****************************
mean_labels_ReLap1 = apply(labels[,ind_plotting,2,drop=T], 2, mean)
mean_labels_ReLap2 = apply(labels[,ind_plotting,3,drop=T], 2, mean)
mean_labels_ReLap = mean_labels_ReLap1 + mean_labels_ReLap2
plot(Combined_Time_Data$Time_to_event[ind_plotting], mean_labels_ReLap,
     col = drug_cols3[numeric_drug[ind_plotting]+1], xaxt='n',
     pch = as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
     ylab='log10 Probability of ReLapse',
     xlab='Months from last episode')
axis(1, at = seq(0, 420, by= 60), labels = seq(0, 420, by=60)/30)
```

![](TimingModelStan_files/figure-html/plotModel3-10.png)<!-- -->



# Model Evaluation


```r
library(loo)
```

```
## Warning: package 'loo' was built under R version 3.4.4
```

```
## This is loo version 2.0.0.
## **NOTE: As of version 2.0.0 loo defaults to 1 core but we recommend using as many as possible. Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session. Visit mc-stan.org/loo/news for details on other changes.
```

```
## 
## Attaching package: 'loo'
```

```
## The following object is masked from 'package:rstan':
## 
##     loo
```

```r
log_lik3 <- extract_log_lik(mod3_RE)

#waic3 = waic(log_lik3)
loo3 = loo(log_lik3)
```

```
## Warning: Relative effective sample sizes ('r_eff' argument) not specified.
## For models fit with MCMC, the reported PSIS effective sample sizes and 
## MCSE estimates will be over-optimistic.
```

```
## Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
```

```r
print(loo3)
```

```
## 
## Computed from 120 by 2707 log-likelihood matrix
## 
##          Estimate    SE
## elpd_loo  -7765.3 123.4
## p_loo       166.1   8.9
## looic     15530.6 246.8
## ------
## Monte Carlo SE of elpd_loo is NA.
## 
## Pareto k diagnostic values:
##                          Count Pct.    Min. n_eff
## (-Inf, 0.5]   (good)     2398  88.6%   22        
##  (0.5, 0.7]   (ok)        246   9.1%   20        
##    (0.7, 1]   (bad)        60   2.2%   21        
##    (1, Inf)   (very bad)    3   0.1%   14        
## See help('pareto-k-diagnostic') for details.
```

# Interpretation of results

The cumulative probability of time to relapse. We draw from the posterior distribution to get a predicted time to relapse (in those who will relapse before they are reinfected).


```r
thetas = extract(mod3_RE)
K = 500000
Early_relapse = sample(x = inv.logit(thetas$logit_EarlyL),replace = T, size = K) 
Mixture = sapply(Early_relapse, FUN = function(x){
  sample(x = 1:2, replace = T, size = 1 , prob = c(x,1-x))
})
# Artesunate mono-therapy
T1s = rweibull(n = K, shape = sample(x = thetas$AS_shape, size = K, replace = T),
               scale = sample(x = thetas$AS_scale, size = K, replace = T))
T2s = rexp(n = K, rate =  sample(x = 1/thetas$inv_gamma, size = K, replace = T))
Times_RelapseAS = c(T1s[Mixture==1], T2s[Mixture==2])
```


```r
par(mfrow=c(1,2), las=1)

# Chloroquine mono-therapy
T1s = rweibull(n = K, shape = sample(x = thetas$CQ_shape, size = K, replace = T),
               scale = sample(x = thetas$CQ_scale, size = K, replace = T))
T2s = rexp(n = K, rate =  sample(x = 1/thetas$inv_gamma, size = K, replace = T))
Times_RelapseCQ = c(T1s[Mixture==1], T2s[Mixture==2])

hist(Times_RelapseAS, breaks = seq(0,7000, by=7), xlim=c(0,100),
     main='', col=drug_cols3['AS'], density = 15, yaxt='n',ylab = '',
     xlab='Months after treatment', xaxt='n')
hist(Times_RelapseCQ, breaks = seq(0,7000, by=7), xlim=c(0,100), 
     main='', add=T, density=15, col=drug_cols3['CHQ'], angle= -30)
axis(1, at = seq(0,90, by=30), labels =  seq(0,90, by=30)/30)
legend('topright', col=drug_cols3[1:2], lwd=3,bty='n',
       legend = c('Artesunate','Chloroquine'))

plot(ecdf(Times_RelapseAS[Times_RelapseAS<360]), 
     col=drug_cols3['AS'], xaxt='n',main='',
     xlab='Months after treatment', bty='n', yaxt='n',
     ylab = 'Relapses observed (%)', lwd=3, lty=2)
axis(2, at = c(0,.25,.5,.75,.9,1), 
     labels = 100*c(0,.25,.5,.75,.9,1))
axis(1, at = seq(0,360, by=60), labels =  seq(0,360, by=60)/30)
lines(ecdf(Times_RelapseCQ[Times_RelapseCQ<360]), 
      col=drug_cols3['CHQ'], lwd=3, lty=2)
axis(1, at = seq(0,360, by=60), labels =  seq(0,360, by=60)/30)
abline(h=.9, lwd=2, lty=2)
abline(v=4*30, lwd=2, lty=2)
```

![](TimingModelStan_files/figure-html/PredictedRelapse-1.png)<!-- -->


```r
# Now we compute the fixed effect probabilities over time
Ts = 1:360
# AS
prob_labels_raw_AS = array(dim = c(4,length(Ts)))
# Reinfection
prob_labels_raw_AS[1,] = exp(mean(thetas$log_p))*dexp(Ts, rate = 1/mean(thetas$inv_lambda));
# Early Relapse
prob_labels_raw_AS[2,] = exp(mean(thetas$log_1m_p))*exp(mean(thetas$log_1m_c1_AS))*exp(mean(thetas$log_EarlyL))*dweibull(Ts , shape = mean(thetas$AS_shape), scale = mean(thetas$AS_scale))
# Late Relapse
prob_labels_raw_AS[3,] = exp(mean(thetas$log_1m_p))*exp(mean(thetas$log_1m_c1_AS))*exp(mean(thetas$log_1m_EarlyL))*dexp(Ts, rate = 1/mean(thetas$inv_gamma))
# Recrudescence
prob_labels_raw_AS[4,] = exp(mean(thetas$log_1m_p))*exp(mean(thetas$log_c1_AS))*dweibull(Ts, shape =  mean(thetas$Recrud_shape), scale = mean(thetas$Recrud_scale))

# CQ
prob_labels_raw_CQ = array(dim = c(4,length(Ts)))
# Reinfection
prob_labels_raw_CQ[1,] = exp(mean(thetas$log_p))*dexp(Ts, rate = 1/mean(thetas$inv_lambda));
# Early Relapse
prob_labels_raw_CQ[2,] = exp(mean(thetas$log_1m_p))*exp(mean(thetas$log_1m_c1_CQ))*exp(mean(thetas$log_EarlyL))*dweibull(Ts , shape = mean(thetas$CQ_shape), scale = mean(thetas$CQ_scale))
# Late Relapse
prob_labels_raw_CQ[3,] = exp(mean(thetas$log_1m_p))*exp(mean(thetas$log_1m_c1_CQ))*exp(mean(thetas$log_1m_EarlyL))*dexp(Ts, rate = 1/mean(thetas$inv_gamma))
# Recrudescence
prob_labels_raw_CQ[4,] = exp(mean(thetas$log_1m_p))*exp(mean(thetas$log_c1_CQ))*dweibull(Ts, shape =  mean(thetas$Recrud_shape), scale = mean(thetas$Recrud_scale))

# CQ+PMQ
prob_labels_raw_CQPMQ = array(dim = c(4,length(Ts)))
# Reinfection
prob_labels_raw_CQPMQ[1,] = exp(mean(thetas$log_p_PMQ))*dexp(Ts, rate = 1/mean(thetas$inv_lambda));
# Early Relapse
prob_labels_raw_CQPMQ[2,] = exp(mean(thetas$log_1m_p_PMQ))*exp(mean(thetas$log_1m_c1_CQ_PMQ))*exp(mean(thetas$log_EarlyL))*dweibull(Ts , shape = mean(thetas$CQ_shape), scale = mean(thetas$CQ_scale))
# Late Relapse
prob_labels_raw_CQPMQ[3,] = exp(mean(thetas$log_1m_p_PMQ))*exp(mean(thetas$log_1m_c1_CQ_PMQ))*exp(mean(thetas$log_1m_EarlyL))*dexp(Ts, rate = 1/mean(thetas$inv_gamma))
# Recrudescence
prob_labels_raw_CQPMQ[4,] = exp(mean(thetas$log_1m_p_PMQ))*exp(mean(thetas$log_c1_CQ_PMQ))*dweibull(Ts, shape =  mean(thetas$Recrud_shape), scale = mean(thetas$Recrud_scale))

for(i in 1:length(Ts)){
  prob_labels_raw_AS[,i] = prob_labels_raw_AS[,i]/sum(prob_labels_raw_AS[,i])
  prob_labels_raw_CQ[,i] = prob_labels_raw_CQ[,i]/sum(prob_labels_raw_CQ[,i])
  prob_labels_raw_CQPMQ[,i] = prob_labels_raw_CQPMQ[,i]/sum(prob_labels_raw_CQPMQ[,i])
}
```



```r
layout(mat = matrix(data = c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7),byrow = T,nrow = 4))
par(las=1,bty='n', cex.lab=1.3, cex.axis=1.3, mar=c(3,4,1,1))
#****************************** Reinfection ******************************
labels = extract(mod3_RE, 'prob_labels')$prob_labels
mean_labels_Reinfection = apply(labels[,ind_plotting,1,drop=T], 2, mean)
plot(Combined_Time_Data$Time_to_event[ind_plotting], log10(mean_labels_Reinfection),
     col = drug_cols3[numeric_drug[ind_plotting]+1],
     pch = 20,
     cex=1,panel.first = grid(),
     ylab='', yaxt='n',xaxt='n',
     xlab='', xlim=c(0,360))
mtext(text = 'Probability of Reinfection',side = 2, las=3,line=3,cex=.8)
mtext(text = 'Months from last episode',side = 1,line=2,cex=1)
axis(1, at = seq(0, 360, by=60), labels = seq(0, 360, by=60)/30)
axis(2, at = -2:0, labels= 10^(-2:0))
axis(2, at = log10(seq(.1,1,by=.1)), labels = NA)
axis(2, at = log10(seq(.01,.1,by=.01)), labels = NA)
axis(2, at = log10(seq(.001,.01,by=.001)), labels = NA)
legend('bottomright',legend = c('Artesunate',
                                'Chloroquine',
                                'Primaquine+'),
       col=drug_cols3,pch = 20, bty='n',lwd=2,lty=NA, cex=1.35)

lines(Ts, log10(prob_labels_raw_AS[1,]), col=drug_cols3[1], lwd=2, lty=2)
lines(Ts, log10(prob_labels_raw_CQ[1,]), col=drug_cols3[2], lwd=2, lty=2)
lines(Ts, log10(prob_labels_raw_CQPMQ[1,]), col=drug_cols3[3], lwd=2, lty=2)


#****************************** Recrudescence ****************************
mean_labels_ReCrud = apply(labels[,ind_plotting,4,drop=T], 2, mean)
plot(Combined_Time_Data$Time_to_event[ind_plotting], (mean_labels_ReCrud),
     col = drug_cols3[numeric_drug[ind_plotting]+1], xaxt='n',
     pch = 20,
     cex=1,panel.first = grid(),
     ylab='Probability of Recrudescence',yaxt='n',
     xlab='', xlim=c(0,30))
mtext(text = 'Months from last episode',side = 1,line=2,cex=1)
axis(1, at = c(0,15,30), labels = c(0,0.5,1))
axis(2, at = c(0,.1,.2))

lines(Ts, (prob_labels_raw_AS[4,]), col=drug_cols3[1], lwd=2, lty=2)
lines(Ts, (prob_labels_raw_CQ[4,]), col=drug_cols3[2], lwd=2, lty=2)
lines(Ts, (prob_labels_raw_CQPMQ[4,]), col=drug_cols3[3], lwd=2, lty=2)

#****************************** Relapse **********************************
mean_labels_ReLap1= apply(labels[,ind_plotting,2,drop=T], 2, mean)
mean_labels_ReLap2= apply(labels[,ind_plotting,3,drop=T], 2, mean)
mean_labels_ReLap = mean_labels_ReLap1 + mean_labels_ReLap2
plot(Combined_Time_Data$Time_to_event[ind_plotting], log10(mean_labels_ReLap),
     col = drug_cols3[numeric_drug[ind_plotting]+1], xaxt='n',
     pch = 20,#as.numeric(Combined_Time_Data$Censored[ind_plotting])+16,
     cex=1,panel.first = grid(),
     ylab='Probability of Relapse',yaxt='n',
     xlab='', xlim=c(0,360))
axis(1, at = seq(0, 360, by=60), labels = seq(0, 360, by=60)/30)
axis(2, at = -2:0, labels= 10^(-2:0))
axis(2, at = log10(seq(.1,1,by=.1)), labels = NA)
axis(2, at = log10(seq(.01,.1,by=.01)), labels = NA)
axis(2, at = log10(seq(.001,.01,by=.001)), labels = NA)
mtext(text = 'Months from last episode',side = 1,line=2,cex=1)

lines(Ts, log10(prob_labels_raw_AS[2,]+prob_labels_raw_AS[3,]), col=drug_cols3[1], lwd=2, lty=2)
lines(Ts, log10(prob_labels_raw_CQ[2,]+prob_labels_raw_CQ[3,]), col=drug_cols3[2], lwd=2, lty=2)
lines(Ts, log10(prob_labels_raw_CQPMQ[2,]+prob_labels_raw_CQPMQ[3,]), col=drug_cols3[3], lwd=2, lty=2)

#****************************** Histograms of key model parameters ***********************

par(cex.lab=.8, cex.axis=1)
indAS = drug_received[which(N_noPMQ_index>0)] == 0
indCQ = drug_received[which(N_noPMQ_index>0)] == 1
plot(density((1-inv.logit(apply(thetas$logit_p,2,mean)[indAS]))*(1-inv.logit(mean(thetas$logit_c1_AS)))),
     xaxt='n',
     xlab = '', main = '',
     col = drug_cols3['AS'], lwd=3,
     yaxt='n',ylab='', xlim=c(0,1))
lines(density(1-inv.logit(apply(thetas$logit_p,2,mean)[indCQ])),
      col = drug_cols3['CHQ'], lwd=3)
mtext(text = 'Relapse: no PMQ',side = 1,line=2,cex=.8)
axis(1, at=c(0,.5,1))

plot(density(1-inv.logit(apply(thetas$logit_p_PMQ,2,mean))),
     xlab = '', main = '', col=drug_cols3['CHQ/PMQ'],lwd=3,
     yaxt='n',ylab='')
mtext(text = 'Relapse: PMQ',side = 1,line=2,cex=.8)

plot(density(inv.logit(mean(thetas$logit_c1_AS))*(1-inv.logit(apply(thetas$logit_p,2,mean)[indAS]))), 
     col = drug_cols3['AS'],lwd=3,xlim = c(0, 0.03),
     xlab = 'Recrudescence', main='',yaxt='n',ylab='')
lines(density(inv.logit(mean(thetas$logit_c1_CQ))*(1-inv.logit(apply(thetas$logit_p,2,mean)[indCQ]))), 
      col = drug_cols3['CHQ'],lwd=3)
#lines(density(inv.logit(mean(thetas$logit_c1_CQ_PMQ))*(1-inv.logit(apply(thetas$logit_p_PMQ,2,mean)))), 
#       col = drug_cols3['CHQ/PMQ'], lwd=3)
mtext(text = 'Recrudescence',side = 1,line=2,cex=.8)
plot(density(inv.logit(thetas$logit_EarlyL)), xlab = '',
     main='',yaxt='n',ylab='', col = 'black', lwd=3)
mtext(text = 'Early Relapse',side = 1,line=2,cex=.8)
```

![](TimingModelStan_files/figure-html/Model3FinalPlot-1.png)<!-- -->



# Extract probabilistic information


```r
labels = extract(mod3_RE, 'prob_labels')$prob_labels
# recrudescence
Combined_Time_Data$Recrudescence_mean_theta = apply(labels[,,4,drop=T],2,mean)
Combined_Time_Data$Recrudescence_975_theta = apply(labels[,,4,drop=T],2,quantile, 
                                                   probs = 0.975)
Combined_Time_Data$Recrudescence_025_theta = apply(labels[,,4,drop=T],2,quantile, 
                                                   probs = 0.025)

# relapse
Relapse_labels = labels[,,2,drop=T] + labels[,,3,drop=T]
Combined_Time_Data$Relapse_mean_theta = apply(Relapse_labels,2,mean)
Combined_Time_Data$Relapse_975_theta = apply(Relapse_labels,2,quantile, probs = 0.975)
Combined_Time_Data$Relapse_025_theta = apply(Relapse_labels,2,quantile, probs = 0.025)

# reinfection
Combined_Time_Data$ReInfection_mean_theta = apply(labels[,,1,drop=T],2,mean)
Combined_Time_Data$ReInfection_975_theta = apply(labels[,,1,drop=T],2,quantile, 
                                                 probs = 0.975)
Combined_Time_Data$ReInfection_025_theta = apply(labels[,,1,drop=T],2,quantile, 
                                                 probs = 0.025)

# Save this for use by the genetic model
Mod3_ThetaEstimates = Combined_Time_Data
Mod3_ThetaEstimates$episode = as.integer(Mod3_ThetaEstimates$episode)
Mod3_ThetaEstimates$episode = Mod3_ThetaEstimates$episode+1
Mod3_ThetaEstimates$patientid = as.character(Mod3_ThetaEstimates$patientid)
Mod3_ThetaEstimates$Episode_Identifier = 
  apply(Mod3_ThetaEstimates,
        1, 
        function(x) {
          paste(x['patientid'], as.integer(x['episode']), sep = '_')
        } )
save(Mod3_ThetaEstimates, file = '../RData/TimingModel/MOD3_theta_estimates.RData')
```

We also save a matrix of posterior samples to be used by the genetic model for full computation of the posterior:

```r
labels = extract(mod3_RE, 'prob_labels')$prob_labels
K_samples = min(100, dim(labels)[1])
random_ind = sample(1:dim(labels)[1], K_samples)
Relapse_labels = labels[random_ind,,2,drop=T] + labels[random_ind,,3,drop=T]
Post_samples_matrix = data.frame(cbind(t(labels[random_ind, ,4,drop=T]),
                                       t(Relapse_labels[, ]),
                                       t(labels[random_ind, ,1,drop=T])))
colnames(Post_samples_matrix) = c(sapply(c('C','L','I'), rep, K_samples))
Post_samples_matrix$Episode_Identifier = apply(Combined_Time_Data, 1, 
                                               function(x) {
                                                 paste(x['patientid'], as.integer(x['episode'])+1, 
                                                       sep = '_')} )
save(Post_samples_matrix, file = '../RData/TimingModel/MOD3_Posterior_samples.RData')
```




```r
addalpha <- function(colors, alpha=0.5) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}
transp_grey = addalpha('grey')
# Order by difference in posterior interval
Combined_Time_Data = Combined_Time_Data[order(-log10(Combined_Time_Data$Relapse_mean_theta)),]

par(las = 1, mfcol=c(2,2), bty='n')
## No PMQ group
ind = Combined_Time_Data$arm_num!='CHQ/PMQ' & !Combined_Time_Data$Censored
plot(log10(Combined_Time_Data$Relapse_mean_theta[ind]), 
     ylim = c(-2.2, 0),type='l',yaxt='n', main ='No primaquine',
     ylab = 'Probability of relapse', xlab = '',panel.first = grid(),
     lwd=2, col = drug_cols2[2], xaxt='n')
axis(1, at = c(1, 400, 800, 1200))
abline(h = log10(0.5))
mtext(text = 'Recurrence rank',side = 1, line = 3)
axis(2, at = c(0:(-2),log10(0.5)), labels = c(10^(0:(-2)),0.5))
polygon(c(1:sum(ind), rev(1:sum(ind))), 
        y = c(log10(Combined_Time_Data$Relapse_025_theta[ind]),
              rev(log10(Combined_Time_Data$Relapse_975_theta[ind]))), 
        border = NA, col = transp_grey)
lines(log10(Combined_Time_Data$Relapse_mean_theta[ind]), 
      lwd=2, col = drug_cols2[2])
# Time of event versus uncertainty in the interval
df = data.frame(time=Combined_Time_Data$Time_to_event[ind],
                uncertainty=log10(Combined_Time_Data$Relapse_975_theta[ind]) -
       log10(Combined_Time_Data$Relapse_025_theta[ind]))
plot(df$time,df$uncertainty,  xaxt='n',
     ylab = 'Uncertainty (log units)', panel.first = grid(),
     col = drug_cols2[as.integer(Combined_Time_Data$arm_num[ind] == 'CHQ/PMQ') + 2],
     pch=20, xlab='', main = 'No primaquine')
axis(1, at = seq(0,360, by = 60), labels = seq(0,360,by=60)/30)
f = loess(uncertainty ~ time, df)
lines(0:400, predict(f, data.frame(time=0:400)),lwd=2)
mtext(text = 'Time-to-event (months)',side = 1, line = 3)

#PMQ group
ind = Combined_Time_Data$arm_num=='CHQ/PMQ' & !Combined_Time_Data$Censored
plot(log10(Combined_Time_Data$Relapse_mean_theta[ind]), 
     ylim = c(-3, 0),type='l',yaxt='n', xlab = '',lwd=2,col=drug_cols2[3],
     ylab = 'Probability of relapse',panel.first = grid(),
     main = 'High-dose primaquine')
mtext(text = 'Recurrence rank',side = 1, line = 3)
abline(h = log10(0.5))
axis(2, at = c(0:(-3),log10(0.5)), labels = c(10^(0:(-3)),0.5))
polygon(c(1:sum(ind), rev(1:sum(ind))), 
        y = c(log10(Combined_Time_Data$Relapse_025_theta[ind]),
              rev(log10(Combined_Time_Data$Relapse_975_theta[ind]))), 
        border = NA, col = transp_grey)
lines(log10(Combined_Time_Data$Relapse_mean_theta[ind]), 
      lwd=2,col=drug_cols2[3])
# Time of event versus uncertainty in the interval
df = data.frame(time=Combined_Time_Data$Time_to_event[ind],
                uncertainty=log10(Combined_Time_Data$Relapse_975_theta[ind]) -
       log10(Combined_Time_Data$Relapse_025_theta[ind]))
plot(df$time, df$uncertainty,  
     ylab = 'Uncertainty (log units)', panel.first = grid(),
     col = drug_cols2[3], xaxt='n',
     pch=20, xlab='', main = 'High-dose primaquine')
axis(1, at = seq(0,360, by = 60), labels = seq(0,360,by=60)/30)
mtext(text = 'Time-to-event (months)',side = 1, line = 3)
f = loess(uncertainty ~ time, df)
lines((0:400), predict(f, data.frame(time=0:400)),lwd=2)
```

![](TimingModelStan_files/figure-html/UncertaintyOutputsModel3-1.png)<!-- -->

Some rough calculations to make sure we're not completely off track with the model output


```
## No radical cure: episodes per year:  3.97
```

```
## Radical cure: episodes per year:  0.2
```


We look at weighted averages of relapses, recrudescences and reinfections.
We pick 500 random draws from the posterior to calculate credible intervals.



The labels on the observed recurrences along with 95% credible intervals:


```
## Relapses are approximately  95.41 ( 93.6 - 96.7 ) % of recurrences after AS
```

```
## Recrudescences are approximately  0.69 ( 0.1 - 1.3 ) % of recurrences after AS
```

```
## Reinfections are approximately  3.9 ( 2.6 - 5.7 ) % of recurrences after AS
```

```
## Relapses are approximately  94.2 ( 92.3 - 95.8 ) % of recurrences after CQ
```

```
## Recrudescences are approximately  0.66 ( 0.1 - 1.1 ) % of recurrences after CQ
```

```
## Reinfections are approximately  5.17 ( 3.7 - 6.8 ) % of recurrences after CQ
```

```
## Relapses are approximately  6.45 ( 4.5 - 8.3 ) % of recurrences after CQ+PMQ
```

```
## Recrudescences are approximately  0.06 ( 0 - 0.1 ) % of recurrences after CQ+PMQ
```

```
## Reinfections are approximately  93.49 ( 91.7 - 95.5 ) % of recurrences after CQ+PMQ
```

