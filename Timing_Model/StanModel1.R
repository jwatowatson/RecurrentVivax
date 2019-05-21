mod1 = 
  '
data {
  // Things that are data structure and strictly data items:
  int<lower=0>                N;                          // The number of individuals
  int<lower=0>                Neps;                       // The total number of observed durations
  int<lower=0>                N_noPMQ;                    // The total number of individuals who did not receive PMQ for at least one episode
  int<lower=0>                N_PMQ;                      // The total number of individuals who received PMQ for at least one episode
  real<lower=0>               Durations[Neps];            // Vector of durations
  int<lower=-1,upper=1>       Censored[Neps];             // -1 is left censored; 0: observed; 1: right censored
  int<lower=0,upper=2>        Drug[Neps];                 // Treatment drug corresponding to each duration
  int<lower=0,upper=N>        ID_of_Patient[Neps];        // The patient ID number for each episode
  int<lower=0,upper=N_noPMQ>  ID_mapped_to_noPMQ_rank[N]; // Vector mapping the the individual ID (from 1 to N) to the random effects index: logit_p

  // Hyper parameters that are part of the `data` section
  real<lower=0>   Hyper_lambda_shape;         // Mean of the prior on time to reinfection
  real<lower=0>   Hyper_lambda_rate;          // SD of the prior on time to reinfection
  real<lower=0>   Hyper_gamma_shape;          // Mean of the prior on time to random reLapse
  real<lower=0>   Hyper_gamma_rate;           // SD of the prior on time to random reLapse
  real<lower=0>   Hyper_lambda_recrud_shape;
  real<lower=0>   Hyper_lambda_recrud_rate;
  real            Early_L_logit_mean;           // mean of prior on mean of logit Early_L
  real<lower=0>   Early_L_logit_sd;             // sd of prior on mean of logit Early_L
  real            Hyper_logit_mean_p_mean;      // mean of prior on mean of logit p
  real<lower=0>   Hyper_logit_mean_p_sd;        // sd of prior on mean of logit p
  real<lower=0>   Hyper_logit_sd_p_lambda;      // the prior rate parameter for the sd of of p
  real            Hyper_logit_c1_mean;
  real<lower=0>   Hyper_logit_c1_sd;
  real<lower=0>   Hyper_AS_shape_mean;      
  real<lower=0>   Hyper_AS_shape_sd;      
  real<lower=0>   Hyper_AS_scale_mean;
  real<lower=0>   Hyper_AS_scale_sd;
  real<lower=0>   Hyper_CQ_shape_mean;
  real<lower=0>   Hyper_CQ_shape_sd;
  real<lower=0>   Hyper_CQ_scale_mean;
  real<lower=0>   Hyper_CQ_scale_sd;
}

parameters {
  real                  logit_p[N_noPMQ];   // Logit proportion of reInfections no radical cure
  real                  logit_c1_AS;        // Proportion of recrudescences when not reinfection - AS monotherapy
  real                  logit_c1_CQ;        // Proportion of recrudescences when not reinfection - CQ monotherapy
  real                  logit_EarlyL;       // Proportion of reLapses that are early/periodic
  real<lower=0>         lambda;             // reInfection rate
  real<lower=0>         gamma;              // Late reLapse rate
  real<lower=0>         lambda_recrud;      // Recrudescence rate
  real                  logit_mean_p;       // Logit mean proportion of reInfection (no rad cure)
  real<lower=0>         logit_sd_p;         // Logit sd of proportion of reInfection (no rad cure)
  real<lower=0>         AS_shape;           // Weibull inverse shape parameter: Artesunate
  real<lower=0>         AS_scale;           // Weibull scale parameter: Artesunate
  real<lower=0>         CQ_shape;           // Weibull inverse shape parameter: Chloroquine
  real<lower=0>         CQ_scale;           // Weibull scale parameter: Chloroquine

}

transformed parameters{
  // Compute the reInfection and recrudescence mixing proportions
  // We define log scale parameters from inverse logit transformation
  real log_p[N_noPMQ];
  real log_c1_AS      = log_inv_logit(logit_c1_AS);
  real log_c1_CQ      = log_inv_logit(logit_c1_CQ);
  real log_EarlyL     = log_inv_logit(logit_EarlyL);
  
  // The 1 minus log proportions
  real log_1m_p[N_noPMQ];
  real log_1m_c1_AS     = log1m_inv_logit(logit_c1_AS);
  real log_1m_c1_CQ     = log1m_inv_logit(logit_c1_CQ);
  real log_1m_EarlyL = log1m_inv_logit(logit_EarlyL);

  log_p = log_inv_logit(logit_p);
  log_1m_p = log1m_inv_logit(logit_p);

}

model {
  
  // ********* Prior *********
  //inv_lambda ~ normal(Hyper_inv_lambda_mean,Hyper_inv_lambda_sd);
  lambda ~ gamma(Hyper_lambda_shape, Hyper_lambda_rate);
  gamma ~ gamma(Hyper_gamma_shape,Hyper_gamma_rate);
  lambda_recrud ~ gamma(Hyper_lambda_recrud_shape,Hyper_lambda_recrud_rate);

  logit_mean_p ~ normal(Hyper_logit_mean_p_mean, Hyper_logit_mean_p_sd);
  logit_sd_p ~ exponential(Hyper_logit_sd_p_lambda);

  AS_shape ~ normal(Hyper_AS_shape_mean,Hyper_AS_shape_sd);
  AS_scale ~ normal(Hyper_AS_scale_mean,Hyper_AS_scale_sd);

  CQ_shape ~ normal(Hyper_CQ_shape_mean,Hyper_CQ_shape_sd);
  CQ_scale ~ normal(Hyper_CQ_scale_mean,Hyper_CQ_scale_sd);

  // Random effects term for the reinfection after no radical cure
  // The probability of reinfection versus relapse when no radical cure
  logit_p ~ normal(logit_mean_p,logit_sd_p);

  // The probability of recrudescence if not reinfection
  // Same hyper prior for the three mixing weights
  logit_c1_AS ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);
  logit_c1_CQ ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);

  // The probability of early relapse
  logit_EarlyL ~ normal(Early_L_logit_mean,Early_L_logit_sd);

  // ********* Likelihood *********
  for(i in 1:Neps){
    int Ind;
    // We iterate through each time to infection and add to log likelihood
    real log_probs[4];
    if(Censored[i] == 0){ // this is an observed time to new infection
      if(Drug[i] == 0){ 
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses, reCrudescences, or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | lambda);          
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_AS + log_EarlyL + weibull_lpdf(Durations[i] | AS_shape, AS_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_AS + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_AS + exponential_lpdf(Durations[i] | lambda_recrud);

        target += log_sum_exp(log_probs);
      } 
      if(Drug[i] == 1){ 
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Chloroquine Monotherapy: reLapses, reCrudescences, or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | lambda);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lpdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lpdf(Durations[i] | lambda_recrud);

        target += log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){ 
        // Chloroquine + Primaquine: reInfection is only option
        target += exponential_lpdf(Durations[i] | lambda);
      }
    } 
    if(Censored[i] == 1){ // this is a right censored time to new infection
      if(Drug[i] == 0){ 
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses, reCrudescences, or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | lambda);          
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_AS + log_EarlyL + weibull_lccdf(Durations[i] | AS_shape, AS_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_AS + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_AS + exponential_lccdf(Durations[i] | lambda_recrud);
        
        target += log_sum_exp(log_probs);

      } 
      if(Drug[i] == 1){ 
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Chloroquine monotherapy: reLapses, reCrudescences, or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | lambda);          
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lccdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lccdf(Durations[i] | lambda_recrud);
        
        target += log_sum_exp(log_probs);
      } 
      if(Drug[i] == 2){ 
        // This is the Chloroquine + Primaquine:
        // reInfection is only option
        target += exponential_lccdf(Durations[i] | lambda);
      }
    }
  }
} 
generated quantities {

  vector[Neps] log_lik;
  matrix[Neps,4] prob_labels;

  // This computes respective densities of each mixture component
  for(i in 1:Neps){
    int Ind;
    vector[4] prob_labels_raw;
    if(Censored[i]==0){
      if(Drug[i] == 0){ // Artesunate monotherapy
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lpdf(Durations[i] | lambda));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p[Ind])*exp(log_1m_c1_AS)*exp(log_EarlyL)*exp(weibull_lpdf(Durations[i] | AS_shape, AS_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p[Ind])*exp(log_1m_c1_AS)*exp(log_1m_EarlyL)*exp(exponential_lpdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p[Ind])*exp(log_c1_AS)*exp(exponential_lpdf(Durations[i] | lambda_recrud));
      }
      if(Drug[i] == 1){ // Chloroquine Monotherapy
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lpdf(Durations[i] | lambda));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_EarlyL)*exp(weibull_lpdf(Durations[i] | CQ_shape, CQ_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_1m_EarlyL)*exp(exponential_lpdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p[Ind])*exp(log_c1_CQ)*exp(exponential_lpdf(Durations[i] | lambda_recrud));
      }
      if(Drug[i] == 2){ // Chloroquine + Primaquine
        // Reinfection
        prob_labels_raw[1] = 1;
        // Early Relapse
        prob_labels_raw[2] = 0;
        // Late Relapse
        prob_labels_raw[3] = 0;
        // Recrudescence
        prob_labels_raw[4] = 0;
      }
    }
    if(Censored[i] == 1){ // this is a right censored time to new infection
      if(Drug[i] == 0){ // Artesunate monotherapy:
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lccdf(Durations[i] | lambda));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p[Ind])*exp(log_1m_c1_AS)*exp(log_EarlyL)*exp(weibull_lccdf(Durations[i] | AS_shape, AS_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p[Ind])*exp(log_1m_c1_AS)*exp(log_1m_EarlyL)*exp(exponential_lccdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p[Ind])*exp(log_c1_AS)*exp(exponential_lccdf(Durations[i] | lambda_recrud));
      }
      if(Drug[i] == 1){ // Chloroquine Monotherapy
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lccdf(Durations[i] | lambda));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_EarlyL)*exp(weibull_lccdf(Durations[i] | CQ_shape, CQ_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_1m_EarlyL)*exp(exponential_lccdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p[Ind])*exp(log_c1_CQ)*exp(exponential_lccdf(Durations[i] | lambda_recrud));
      }
      if(Drug[i] == 2){ // Chloroquine + Primaquine
        // Reinfection
        prob_labels_raw[1] = 1;
        // Early Relapse
        prob_labels_raw[2] = 0;
        // Late Relapse
        prob_labels_raw[3] = 0;
        // Recrudescence
        prob_labels_raw[4] = 0;

      }
    }
    // normalise to 1
    for(k in 1:4){
      prob_labels[i,k] = prob_labels_raw[k]/sum(prob_labels_raw);
    }
  }

  // This computes the log likelihood (for model comparison at the end)
  for(i in 1:Neps){
    // We iterate through each time to infection and add to log likelihood
    real log_probs[4];
    int Ind;
    if(Censored[i] == 0){ // this is an observed time to new infection
      if(Drug[i] == 0){
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | lambda);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_AS + log_EarlyL + weibull_lpdf(Durations[i] | AS_shape, AS_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_AS + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_AS + exponential_lpdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 1){
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Chloroquine Monotherapy
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | lambda);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lpdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lpdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){
        // Chloroquine + Primaquine
        // reInfection: only option
        log_lik[i] = exponential_lpdf(Durations[i] | lambda);

      }
    }
    if(Censored[i] == 1){ // this is a right censored time to new infection
      // This is the unobserved case (data are right censored)
      if(Drug[i] == 0){
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | lambda);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_AS + log_EarlyL + weibull_lccdf(Durations[i] | AS_shape, AS_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_AS + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_AS + exponential_lccdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 1){
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Chloroquine monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | lambda);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lccdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lccdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){
        // Chloroquine + Primaquine:
        // reInfection
        log_lik[i] = exponential_lccdf(Durations[i] | lambda);

      }
    }
  }
}
'

Timing_Model1 = stan_model(model_code = mod1)
