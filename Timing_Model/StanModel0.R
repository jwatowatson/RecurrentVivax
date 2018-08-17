mod0 = 
  '
data {
  // Things that are data structure and strictly data items:
  int<lower=0>          N;                // The number of individuals
  int<lower=0>          Neps;             // The total number of observed durations
  real<lower=0>         Durations[Neps];  // Vector of durations
  int<lower=0,upper=1>  Censored[Neps];   // Logical Vector: censored observation or nor
  int<lower=0,upper=2>  Drug[Neps];       // Treatment drug corresponding to each duration
  int<lower=0>          N_index[Neps];    // Index of individual groupings
  
  // Hyper parameters that are part of the `data` section
  real<lower=0>   mu_inv_lambda;          // Mean of the prior on time to reinfection
  real<lower=0>   sigma_inv_lambda;       // SD of the prior on time to reinfection
  real<lower=0>   mu_inv_gamma;           // Mean of the prior on time to reinfection
  real<lower=0>   sigma_inv_gamma;        // SD of the prior on time to reinfection
  real            Hyper_logit_mean_p;     // mean of prior on mean of p
  real<lower=0>   Hyper_logit_sd_p;       // sd of prior on mean of p
}

parameters {
  real                  logit_p[N];     // logit proportion of reInfections for each individual
  real                  logit_c1;       // proportion of recrudescences when not reinfection - no radical cure
  real<lower=0>         inv_lambda;     // reInfection mean time
  real<lower=0>         inv_gamma;      // relapse mean time
  real                  logit_mean_p;   // logit mean proportion of reInfection
  real<lower=0>         logit_sd_p;     // logit sd of proportion of reInfection
  real<lower=0>         Recrud_mean;    // Mean time to recrudescence
  real<lower=0>         Recrud_sd;      // SD of time to recrudescence

}

transformed parameters{
  // Turn the inverse rates into rate parameters for the likelihood calculation
  real lambda = 1/inv_lambda;   
  real gamma = 1/inv_gamma;
  // Compute the reInfection rates for each individual
  real log_p[N];
  real log_c1;
  real log_1m_p[N];
  real log_1m_c1;
  log_p = log_inv_logit(logit_p);
  log_c1 = log_inv_logit(logit_c1);
  log_1m_p = log1m_inv_logit(logit_p);
  log_1m_c1 = log1m_inv_logit(logit_c1);
}

model {
  
  // ********* Prior *********
  inv_lambda ~ normal(mu_inv_lambda,sigma_inv_lambda);
  inv_gamma ~ normal(mu_inv_gamma,sigma_inv_gamma);
  logit_mean_p ~ normal(Hyper_logit_mean_p, Hyper_logit_sd_p);
  logit_sd_p ~ normal(1,1);
  Recrud_mean ~ normal(14, 2);
  Recrud_sd ~ normal(5, 1);

  // The probability of reinfection versus relapse when no radical cure
  logit_p ~ normal(logit_mean_p,logit_sd_p);
  // The probability of recrudescence if not reinfection
  logit_c1 ~ normal(-3,.5);

  // ********* Likelihood *********
  for(i in 1:Neps){
    // We iterate through each time to infection and add to log likelihood
    real log_probs[3];
    if(Censored[i] == 0){ // this is an observed time to new infection
      if(Drug[i] == 0){ 
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lpdf(Durations[i] | lambda);          
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd);

        target += log_sum_exp(log_probs);
      } 
      if(Drug[i] == 1){ 
        // Chloroquine Monotherapy
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lpdf(Durations[i] | lambda);
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd);

        target += log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){ 
        
        // Chloroquine + Primaquine: no relapses possible in this model
        // reInfection
        //log_probs[1] = log_1m_c1 + exponential_lpdf(Durations[i] | lambda);
        // Relapse: does not happen so put in a constant value
        //log_probs[2] = log_c1 + weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd);

        target += exponential_lpdf(Durations[i] | lambda);
        //log_sum_exp(log_probs[1], log_probs[2]);
      }
    } else { 
      if(Drug[i] == 0){ 
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lccdf(Durations[i] | lambda);          
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lccdf(Durations[i] | Recrud_mean, Recrud_sd);
        
        target += log_sum_exp(log_probs);

      } 
      if(Drug[i] == 1){ 
        // This is the Chloroquine monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lccdf(Durations[i] | lambda);          
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lccdf(Durations[i] | Recrud_mean, Recrud_sd);
        
        target += log_sum_exp(log_probs);
      } 
      if(Drug[i] == 2){ 
        target += exponential_lccdf(Durations[i] | lambda);
        //log_sum_exp(log_probs[1], log_probs[2]);
      }
    }
  }
} 

generated quantities {

  vector[Neps] log_lik;
  matrix[Neps,3] prob_labels;
  
  // This computes respective densities of each mixture component
  for(i in 1:Neps){
    vector[3] prob_labels_raw;
    if(Censored[i]==0){
      if(Drug[i] == 0){ // Artesunate monotherapy
        // Reinfection
        prob_labels_raw[1] = exp(log_p[N_index[i]])*exp(exponential_lpdf(Durations[i] | lambda));
        // Relapse
        prob_labels_raw[2] = exp(log_1m_p[N_index[i]])*exp(log_1m_c1)*exp(exponential_lpdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[3] = exp(log_1m_p[N_index[i]])*exp(log_c1)*exp(weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd));
      }
      if(Drug[i] == 1){ // Chloroquine Monotherapy
        // Reinfection
        prob_labels_raw[1] = exp(log_p[N_index[i]])*exp(exponential_lpdf(Durations[i] | lambda));
        // Relapse
        prob_labels_raw[2] = exp(log_1m_p[N_index[i]])*exp(log_1m_c1)*exp(exponential_lpdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[3] = exp(log_1m_p[N_index[i]])*exp(log_c1)*exp(weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd));
      }
      if(Drug[i] == 2){ // Chloroquine + Primaquine
        // Reinfection
        prob_labels_raw[1] = 0;
        // Relapse
        prob_labels_raw[2] = 0;
        // Recrudescence
        prob_labels_raw[3] = 1;
      }
    } else {
      if(Drug[i] == 0){ // Artesunate monotherapy:
        // Reinfection
        prob_labels_raw[1] = exp(log_p[N_index[i]])*exp(exponential_lccdf(Durations[i] | lambda));
        // Relapse
        prob_labels_raw[2] = exp(log_1m_p[N_index[i]])*exp(log_1m_c1)*exp(exponential_lccdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[3] = exp(log_1m_p[N_index[i]])*exp(log_c1)*exp(weibull_lccdf(Durations[i] | Recrud_mean, Recrud_sd));
      }
      if(Drug[i] == 1){ // Chloroquine Monotherapy
        // Reinfection
        prob_labels_raw[1] = exp(log_p[N_index[i]])*exp(exponential_lccdf(Durations[i] | lambda));
        // Relapse
        prob_labels_raw[2] = exp(log_1m_p[N_index[i]])*exp(log_1m_c1)*exp(exponential_lccdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[3] = exp(log_1m_p[N_index[i]])*exp(log_c1)*exp(weibull_lccdf(Durations[i] | Recrud_mean, Recrud_sd));
      }
      if(Drug[i] == 2){ // Chloroquine + Primaquine
        // Reinfection
        prob_labels_raw[1] = 0;
        // Relapse
        prob_labels_raw[2] = 0;
        // Recrudescence
        prob_labels_raw[3] = 1;
      }
    }
    // normalise to 1
    for(k in 1:3){
      prob_labels[i,k] = prob_labels_raw[k]/sum(prob_labels_raw);
    }
  }

  // This computes the log likelihood (for model comparison at the end)
  for(i in 1:Neps){
    // We iterate through each time to infection and add to log likelihood
    real log_probs[3];
    
    if(Censored[i] == 0){ // this is an observed time to new infection
      if(Drug[i] == 0){
        // This is the Artesunate monotherapy
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lpdf(Durations[i] | lambda);          
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 1){
        // Chloroquine Monotherapy
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lpdf(Durations[i] | lambda);          
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){
        // Chloroquine + Primaquine: no relapses/recrudescences possible in this model
        // reInfection
        log_lik[i] = weibull_lpdf(Durations[i] | Recrud_mean, Recrud_sd);
      }
    } else {
      // This is the unobserved case (data are right censored)
      if(Drug[i] == 0){
        // This is the Artesunate monotherapy
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lccdf(Durations[i] | lambda);          
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lccdf(Durations[i] | Recrud_mean, Recrud_sd);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 1){
        // Chloroquine Monotherapy
        // reInfection
        log_probs[1] = log_p[N_index[i]] + exponential_lccdf(Durations[i] | lambda);          
        // reLapse
        log_probs[2] = log_1m_p[N_index[i]] + log_1m_c1 + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[3] = log_1m_p[N_index[i]] + log_c1 + weibull_lccdf(Durations[i] | Recrud_mean, Recrud_sd);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){
        // Chloroquine + Primaquine: no relapses possible in this model
        // reInfection
        log_lik[i] = weibull_lccdf(Durations[i] | Recrud_mean, Recrud_sd);
      }
    }
  }
}

'

Timing_Model0_RE = stan_model(model_code = mod0)