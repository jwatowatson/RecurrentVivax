## ---- echo=FALSE---------------------------------------------------------

.libPaths("/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(plyr)
library(dplyr)
library(rstan)

## ---- include=FALSE------------------------------------------------------
# The code is in a different file as it's very long
writeLines('Compiling model 0....')
source('StanModel0.R')

## ---- include=FALSE------------------------------------------------------
writeLines('Compiling model 1....')
source('StanModel1.R')

## ----include=FALSE-------------------------------------------------------
source('StanModel2.R')

## ---- include=FALSE------------------------------------------------------
source('StanModel3.R')

## ------------------------------------------------------------------------
# Load the dataset - already been partially cleaned
load('../../RData/VHX_ClinicalData.RData')
# Get rid of the very short durations 
vhx_dat = filter(vhx_dat, Time_to_event>10)
#vhx_dat = filter(vhx_dat, arm_num=='CHQ')
# We sort the ids so that ids go from 1...Nmax
ids = unique(vhx_dat$patientid)
rank_ids = 1:length(ids)
vhx_dat$ordered_ids = unlist(sapply(1:length(ids),
                                    function(x) rep(rank_ids[x], 
                                                    sum(vhx_dat$patientid==ids[x])))
)
N_index = as.integer(vhx_dat$ordered_ids)
# Turn drug into a numeric vector
numeric_drug = as.integer(revalue(vhx_dat$arm_num, c('AS'='0','CHQ'='1','CHQ/PMQ'='2')))

## ------------------------------------------------------------------------
# The hierachical parameters defining the prior distributions for model 1
Prior_params_M0 = list(mu_inv_lambda = 600,
                       sigma_inv_lambda = 25,
                       mu_inv_gamma = 40,
                       sigma_inv_gamma = 2,
                       Hyper_p_a1 = 1/2,
                       Hyper_p_a2 = 1/20)
# The hierachical parameters defining the prior distributions for model 1
Prior_params_M1 = list(mu_inv_lambda = 600,
                       sigma_inv_lambda = 50,
                       mu_AS_shape = 2,
                       sigma_AS_shape = 1,
                       mu_AS_scale = 25,
                       sigma_AS_scale = 5,
                       mu_CQ_shape = 2,
                       sigma_CQ_shape = 1,
                       mu_CQ_scale = 42,
                       sigma_CQ_scale = 5,
                       Hyper_p_a1 = 1/2,
                       Hyper_p_a2 = 1/20)
# Model 2 has the same parameters with a few extra
Prior_params_M2 = c(Prior_params_M1, 
                    Hyper_Early_L_a1 = 1,
                    Hyper_Early_L_a2 = 1,
                    mu_inv_gamma = 100,
                    sigma_inv_gamma = 25)
# Model 3: extra parameters
Prior_params_M3 = c(Prior_params_M2, 
                    Hyper_p_PMQ_a1 = 1/9,
                    Hyper_p_PMQ_a2 = 1)

## ----runModels, include=FALSE--------------------------------------------
# Choose as many chains as available cores
Chains = parallel::detectCores()-1
options(mc.cores = Chains)
IT = 10000
WarmUp = .5*IT
thin=100

VHXdata = list(N         = as.integer(length(ids)),#Number of individuals
               Neps      = as.integer(nrow(vhx_dat)),#Number of durations
               Durations = as.double(vhx_dat$Time_to_event),#Time to reinfection or time to censoring
               Censored  = as.integer(vhx_dat$Right_Truncated),#If the duration is right censored or not
               Drug      = numeric_drug,
               N_index   = N_index)


mod0_RE = sampling(Timing_Model0_RE,
                   data = c(VHXdata,Prior_params_M0),                     
                   iter = IT, warmup = WarmUp, 
                   chains=Chains, thin = thin)

mod1_RE = sampling(Timing_Model1_RE,
              data = c(VHXdata, Prior_params_M1),
              iter = IT, warmup = WarmUp, chains=Chains, thin=thin)

mod2_RE = sampling(Timing_Model2_RE,
              data = c(VHXdata, Prior_params_M2),
              iter = IT, warmup = WarmUp, chains=Chains, thin=thin)

mod3_RE = sampling(Timing_Model3_RE,
              data = c(VHXdata, Prior_params_M3),
              iter = IT, warmup = WarmUp, chains=Chains,thin=thin)

save(mod0_RE, file = 'OutputResults/StanModels_mod0.RData')
save(mod1_RE, file = 'OutputResults/StanModels_mod1.RData')
save(mod2_RE, file = 'OutputResults/StanModels_mod2.RData')
save(mod3_RE, file = 'OutputResults/StanModels_mod3.RData')

## ----plotModel0----------------------------------------------------------
par(las=1)
traceplot(mod0_RE, c('inv_lambda','p_a1', 'p_a2','inv_gamma'))
thetas = extract(mod0_RE)
hist(thetas$inv_lambda, main='Mean time to reinfection',xlab='1/lambda (days)')

Ps = thetas$p
hist(apply(Ps,1,median), main = 'proportion of reInfections', xlab='p')

labels = extract(mod0_RE, 'prob_labels')$prob_labels
mean_labels = apply(labels, c(2,3), mean)
plot(vhx_dat$Time_to_event, log10(mean_labels[,1]), 
     col = numeric_drug+1, 
     pch = as.numeric(vhx_dat$Right_Truncated)+1, 
     ylab='Probability of ReInfection', 
     xlab='days from last infection',yaxt='n')
axis(2, at= -2:0, labels= 10^(-2:0))
legend('bottomright',
       legend = c('AS','CQ','CQ/PMQ',
                  'Not censored','Censored observation'), 
       col=c(1:3,1,1),pch = c(rep(1,3),1:2))

## ----plotModel1----------------------------------------------------------
pdf('AllResults.pdf')
par(las=1)
# traceplot(mod1_RE, c('inv_lambda','p_a1', 'p_a2', 
#                      'AS_shape', 'CQ_shape', 'AS_scale', 'CQ_scale'))
thetas = extract(mod1_RE)
hist(thetas$inv_lambda, main='Mean time to reinfection',xlab='1/lambda (days)')

hist(apply(thetas$p,1,median), main = 'proportion of reInfections without RC', xlab='p')

labels = extract(mod1_RE, 'prob_labels')$prob_labels
mean_labels = apply(labels, c(2,3), mean)
plot(vhx_dat$Time_to_event, log10(mean_labels[,1]), col = numeric_drug+1, pch = as.numeric(vhx_dat$Right_Truncated)+1, ylab='Probability of ReInfection', xlab='days from last infection',yaxt='n')
axis(2, at= -2:0, labels= 10^(-2:0))
legend('bottomright',legend = c('AS','CQ','CQ/PMQ', 'Not censored','Censored observation'), col=c(1:3,1,1),pch = c(rep(1,3),1:2))

## ----plotModel2----------------------------------------------------------
par(las=1)
# traceplot(mod2_RE, c('inv_lambda','Early_L','p_a1','p_a2',
#                      'inv_gamma','Early_L_a1','Early_L_a2',
#                      'AS_shape', 'CQ_shape', 'AS_scale', 'CQ_scale'))
thetas = extract(mod2_RE)
hist(thetas$inv_lambda, main='Mean time to reinfection', xlab='1/lambda (days)')
hist(thetas$inv_gamma, main='Mean time to late reLapse', xlab='1/gamma (days)')

hist(apply(thetas$p,1,median), main = 'proportion of reInfections', xlab = 'p')
hist(apply(thetas$Early_L,1,median), main = 'proportion of early reLapses', xlab='q')

labels = extract(mod2_RE, 'prob_labels')$prob_labels
mean_labels = apply(labels, c(2,3), mean)
plot(vhx_dat$Time_to_event, log10(mean_labels[,1]), col = numeric_drug+1, pch = as.numeric(vhx_dat$Right_Truncated)+1, ylab='Probability of ReInfection', xlab='days from last infection',yaxt='n')
axis(2, at= -2:0, labels= 10^(-2:0))
legend('bottomright',legend = c('AS','CQ','CQ/PMQ', 'Not censored','Censored observation'), col=c(1:3,1,1),pch = c(rep(1,3),1:2))

## ----plotModel3----------------------------------------------------------
par(las=1)
traceplot(mod3, c('inv_lambda','Early_L','p','inv_gamma'))
traceplot(mod3, c('p_PMQ','AS_scale','CQ_scale'))

thetas = extract(mod3)
hist(thetas$inv_lambda, main='Mean time to reinfection')
hist(thetas$inv_gamma, main='Mean time to late reLapse')

hist(thetas$p, main = 'proportion of reInfections without RC')
hist(thetas$Early_L, main = 'proportion of early reLapses')
hist(thetas$p_PMQ, main = 'proportion of reInfections with RC')

labels = extract(mod3, 'prob_labels')$prob_labels
mean_labels = apply(labels, c(2,3), mean)
plot(vhx_dat$Time_to_event, log10(mean_labels[,1]), col = numeric_drug+1, pch = as.numeric(vhx_dat$Right_Truncated)+1, ylab='Probability of ReInfection', xlab='days from last infection',yaxt='n')
axis(2, at= -2:0, labels= 10^(-2:0))
legend('bottomright',legend = c('AS','CQ','CQ/PMQ', 'Not censored','Censored observation'), col=c(1:3,1,1),pch = c(rep(1,3),1:2))

## ------------------------------------------------------------------------
plot(vhx_dat$Time_to_event, log10(mean_labels[,1]), 
     col = numeric_drug+1, xlim=c(0,60),
     pch = as.numeric(vhx_dat$Right_Truncated)+1, 
     ylab='Probability of ReInfection', 
     xlab='days from last infection',yaxt='n')
axis(2, at= -2:0, labels= 10^(-2:0))
legend('topleft',legend = c('AS','CQ','CQ/PMQ', 'Not censored','Censored observation'), 
       col=c(1:3,1,1),pch = c(rep(1,3),1:2))

dev.off()
## ------------------------------------------------------------------------
library(loo)
log_lik0 <- extract_log_lik(mod0)
log_lik1 <- extract_log_lik(mod1)
log_lik2 <- extract_log_lik(mod2)
log_lik3 <- extract_log_lik(mod3)
waic1 = waic(log_lik1)
loo1 = loo(log_lik1)
waic2 = waic(log_lik2)
loo2 = loo(log_lik2)
waic3 = waic(log_lik3)
loo3 = loo(log_lik3)
compare(waic1, waic2, waic3)
compare(loo1, loo2, loo3)

