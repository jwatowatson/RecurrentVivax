mod2 = 
  '
data {
  // Things that are data structure and strictly data items:
  int<lower=0>                N;                        // The number of individuals
  int<lower=0>                Neps;                     // The total number of observed durations
  int<lower=0>                N_noPMQ;                  // The total number of individuals who did not receive PMQ for at least one episode
  int<lower=0>                N_PMQ;                    // The total number of individuals who received PMQ for at least one episode
  real<lower=0>               Durations[Neps];          // Vector of durations
  int<lower=-1,upper=1>       Censored[Neps];           // -1 is left censored; 0: observed; 1: right censored
  int<lower=0,upper=2>        Drug[Neps];               // Treatment drug corresponding to each duration
  int<lower=0,upper=N>        ID_of_Patient[Neps];            // Index of individual groupings
  int<lower=0,upper=N_noPMQ>  ID_mapped_to_noPMQ_rank[N];     // Vector mapping the the individual ID (from 1 to N) to the random effects index: logit_p
  int<lower=0,upper=N_PMQ>    ID_mapped_to_PMQ_rank[N];       // Vector mapping the the individual ID (from 1 to N) to the random effects index: logit_p_PMQ
  int<lower=1,upper=2>        Study_Period[Neps];       // 1: VHX; 2: BPD

  // Hyper parameters that are part of the `data` section
  real<lower=0>   Hyper_lambda_shape;     // Mean of the prior on time to reinfection
  real<lower=0>   Hyper_lambda_rate;      // SD of the prior on time to reinfection
  real<lower=0>   Hyper_gamma_shape;      // Mean of the prior on time to random reLapse
  real<lower=0>   Hyper_gamma_rate;       // SD of the prior on time to random reLapse
  real<lower=0>   Hyper_lambda_recrud_shape;
  real<lower=0>   Hyper_lambda_recrud_rate;
  real            Early_L_logit_mean;     // mean of prior on mean of logit Early_L
  real<lower=0>   Early_L_logit_sd;       // sd of prior on mean of logit Early_L
  real            Hyper_logit_mean_pPMQ_mean; //
  real<lower=0>   Hyper_logit_mean_pPMQ_sd;   //
  real<lower=0>   Hyper_logit_sd_pPMQ_lambda;
  real            Hyper_logit_mean_p_mean;     // mean of prior on mean of logit p
  real<lower=0>   Hyper_logit_mean_p_sd;       // sd of prior on mean of logit p
  real<lower=0>   Hyper_logit_sd_p_lambda;
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
  real<lower=0>   Hyper_mean_rate_decrease;
  real<lower=0>   Hyper_sd_rate_decrease;
}

parameters {
  real                  logit_p[N_noPMQ];     // logit proportion of reInfections no radical cure
  real                  logit_p_PMQ[N_PMQ];   // The proportion of reInfections after radical cure
  real                  logit_c1_AS;          // proportion of recrudescences when not reinfection - AS monotherapy
  real                  logit_c1_CQ;          // proportion of recrudescences when not reinfection - CQ monotherapy
  real                  logit_c1_CQ_PMQ;      // proportion of recrudescences when not reinfection - CQ + PMQ
  real                  logit_EarlyL;         // Proportion of reLapses that are early/periodic
  real<lower=0>         lambda;               // reInfection rate in VHX study
  real<lower=0>         gamma;                // Late reLapse rate
  real<lower=0>         lambda_recrud;        // Recrudescence rate
  real                  logit_mean_p;         // logit mean proportion of reInfection (no rad cure)
  real<lower=0>         logit_sd_p;           // logit sd of proportion of reInfection (no rad cure)
  real                  logit_mean_p_PMQ;     // logit mean proportion of reInfection (after rad cure)
  real<lower=0>         logit_sd_p_PMQ;       // logit sd of proportion of reInfection (after rad cure)
  real<lower=0>         AS_shape;             // Weibull shape parameter: Artesunate
  real<lower=0>         AS_scale;             // Weibull scale parameter: Artesunate
  real<lower=0>         CQ_shape;             // Weibull shape parameter: Chloroquine
  real<lower=0>         CQ_scale;             // Weibull scale parameter: Chloroquine
  real<lower=0>         rate_decrease;        // Decrease in reinfection rate between VHX and BPD studies
}

transformed parameters{

  // Compute the reInfection and recrudescence mixing proportions
  // We define log scale parameters from inverse logit transformation
  real log_p[N_noPMQ];
  real log_p_PMQ[N_PMQ];
  real log_c1_AS      = log_inv_logit(logit_c1_AS);
  real log_c1_CQ      = log_inv_logit(logit_c1_CQ);
  real log_c1_CQ_PMQ  = log_inv_logit(logit_c1_CQ_PMQ);
  real log_EarlyL     = log_inv_logit(logit_EarlyL);
  
  // The 1 minus log proportions
  real log_1m_p[N_noPMQ];
  real log_1m_p_PMQ[N_PMQ];
  real log_1m_c1_AS     = log1m_inv_logit(logit_c1_AS);
  real log_1m_c1_CQ     = log1m_inv_logit(logit_c1_CQ);
  real log_1m_c1_CQ_PMQ = log1m_inv_logit(logit_c1_CQ_PMQ);
  real log_1m_EarlyL = log1m_inv_logit(logit_EarlyL);

  log_p = log_inv_logit(logit_p);
  log_p_PMQ = log_inv_logit(logit_p_PMQ);

  log_1m_p = log1m_inv_logit(logit_p);
  log_1m_p_PMQ = log1m_inv_logit(logit_p_PMQ);
  
}

model {
  
  // ********* Prior *********
  lambda ~ gamma(Hyper_lambda_shape, Hyper_lambda_rate);
  rate_decrease ~ normal(Hyper_mean_rate_decrease,Hyper_sd_rate_decrease);

  gamma ~ gamma(Hyper_gamma_shape,Hyper_gamma_rate);
  lambda_recrud ~ gamma(Hyper_lambda_recrud_shape,Hyper_lambda_recrud_rate);

  logit_mean_p ~ normal(Hyper_logit_mean_p_mean, Hyper_logit_mean_p_sd);
  logit_sd_p ~ exponential(Hyper_logit_sd_p_lambda);
  
  logit_mean_p_PMQ ~ normal(Hyper_logit_mean_pPMQ_mean,Hyper_logit_mean_pPMQ_sd);
  logit_sd_p_PMQ ~ exponential(Hyper_logit_sd_pPMQ_lambda);

  AS_shape ~ normal(Hyper_AS_shape_mean,Hyper_AS_shape_sd);
  AS_scale ~ normal(Hyper_AS_scale_mean,Hyper_AS_scale_sd);

  CQ_shape ~ normal(Hyper_CQ_shape_mean,Hyper_CQ_shape_sd);
  CQ_scale ~ normal(Hyper_CQ_scale_mean,Hyper_CQ_scale_sd);

  // Random effects term for the reinfection after no radical cure
  // The probability of reinfection versus relapse when no radical cure
  logit_p ~ normal(logit_mean_p,logit_sd_p);
  
  // Random effects term for reinfection after radical cure
  // The probability of relapse after radical cure
  logit_p_PMQ ~ normal(logit_mean_p_PMQ,logit_sd_p_PMQ);

  // The probability of recrudescence if not reinfection
  // Same hyper prior for the three mixing weights
  logit_c1_AS ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);
  logit_c1_CQ ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);
  logit_c1_CQ_PMQ ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);

  // The probability of early relapse
  logit_EarlyL ~ normal(Early_L_logit_mean,Early_L_logit_sd);

  // ********* Likelihood *********
  // We iterate through each time to infection and add to log likelihood
  for(i in 1:Neps){
    int Ind;
    real log_probs[4];
    real reinfection_rate;
    if(Study_Period[i]==1) reinfection_rate = lambda;
    if(Study_Period[i]==2) reinfection_rate = lambda*rate_decrease;

    if(Censored[i] == 0){ // this is an observed time to new infection
      if(Drug[i] == 0){ 
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | reinfection_rate);          
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
        // Chloroquine Monotherapy
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | reinfection_rate);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lpdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lpdf(Durations[i] | lambda_recrud);

        target += log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){ 
        Ind = ID_mapped_to_PMQ_rank[ID_of_Patient[i]];
        // Chloroquine plus Primaquine
        // reInfection
        log_probs[1] = log_p_PMQ[Ind] + exponential_lpdf(Durations[i] | reinfection_rate);
        // Early reLapse
        log_probs[2] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_EarlyL + weibull_lpdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p_PMQ[Ind] + log_c1_CQ_PMQ + exponential_lpdf(Durations[i] | lambda_recrud);
        
        target += log_sum_exp(log_probs);
      }
    } 
    if(Censored[i] == 1){ // this is a right censored time to new infection
      if(Drug[i] == 0){ 
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | reinfection_rate);          
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
        // This is the Chloroquine monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | reinfection_rate);          
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lccdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lccdf(Durations[i] | lambda_recrud);
        
        target += log_sum_exp(log_probs);
      } 
      if(Drug[i] == 2){ 
        Ind = ID_mapped_to_PMQ_rank[ID_of_Patient[i]];
        // This is the Chloroquine monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p_PMQ[Ind] + exponential_lccdf(Durations[i] | reinfection_rate);          
        // Early reLapse
        log_probs[2] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_EarlyL + weibull_lccdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p_PMQ[Ind] + log_c1_CQ_PMQ + exponential_lccdf(Durations[i] | lambda_recrud);
        
        target += log_sum_exp(log_probs);
      }
    }
  }
} 

generated quantities {

  vector[Neps] log_lik;
  matrix[Neps,4] prob_labels;
  real reinfection_rates[2];
  reinfection_rates[1] = lambda;
  reinfection_rates[2] = lambda*rate_decrease;

  // This computes respective densities of each mixture component
  for(i in 1:Neps){
    int Ind;
    vector[4] prob_labels_raw;
    real reinfection_rate;

    if(Study_Period[i]==1) reinfection_rate = lambda;
    if(Study_Period[i]==2) reinfection_rate = lambda*rate_decrease;

    if(Censored[i]==0){
      if(Drug[i] == 0){ // Artesunate monotherapy
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lpdf(Durations[i] | reinfection_rate));
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
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lpdf(Durations[i] | reinfection_rate));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_EarlyL)*exp(weibull_lpdf(Durations[i] | CQ_shape, CQ_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_1m_EarlyL)*exp(exponential_lpdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p[Ind])*exp(log_c1_CQ)*exp(exponential_lpdf(Durations[i] | lambda_recrud));
      }
      if(Drug[i] == 2){ // Chloroquine + Primaquine
        Ind = ID_mapped_to_PMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p_PMQ[Ind])*exp(exponential_lpdf(Durations[i] | reinfection_rate));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p_PMQ[Ind])*exp(log_1m_c1_CQ_PMQ)*exp(log_EarlyL)*exp(weibull_lpdf(Durations[i] | CQ_shape, CQ_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p_PMQ[Ind])*exp(log_1m_c1_CQ_PMQ)*exp(log_1m_EarlyL)*exp(exponential_lpdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p_PMQ[Ind])*exp(log_c1_CQ_PMQ)*exp(exponential_lpdf(Durations[i] | lambda_recrud));
      }
    }
    if(Censored[i] == 1){ // this is a right censored time to new infection
      if(Drug[i] == 0){ // Artesunate monotherapy:
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lccdf(Durations[i] | reinfection_rate));
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
        prob_labels_raw[1] = exp(log_p[Ind])*exp(exponential_lccdf(Durations[i] | reinfection_rate));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_EarlyL)*exp(weibull_lccdf(Durations[i] | CQ_shape, CQ_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p[Ind])*exp(log_1m_c1_CQ)*exp(log_1m_EarlyL)*exp(exponential_lccdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p[Ind])*exp(log_c1_CQ)*exp(exponential_lccdf(Durations[i] | lambda_recrud));
      }
      if(Drug[i] == 2){ // Chloroquine + Primaquine
        Ind = ID_mapped_to_PMQ_rank[ID_of_Patient[i]];
        // Reinfection
        prob_labels_raw[1] = exp(log_p_PMQ[Ind])*exp(exponential_lccdf(Durations[i] | reinfection_rate));
        // Early Relapse
        prob_labels_raw[2] = exp(log_1m_p_PMQ[Ind])*exp(log_1m_c1_CQ_PMQ)*exp(log_EarlyL)*exp(weibull_lccdf(Durations[i] | CQ_shape, CQ_scale));
        // Late Relapse
        prob_labels_raw[3] = exp(log_1m_p_PMQ[Ind])*exp(log_1m_c1_CQ_PMQ)*exp(log_1m_EarlyL)*exp(exponential_lccdf(Durations[i] | gamma));
        // Recrudescence
        prob_labels_raw[4] = exp(log_1m_p_PMQ[Ind])*exp(log_c1_CQ_PMQ)*exp(exponential_lccdf(Durations[i] | lambda_recrud));

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
    real reinfection_rate;

    if(Study_Period[i]==1) reinfection_rate = lambda;
    if(Study_Period[i]==2) reinfection_rate = lambda*rate_decrease;
    if(Censored[i] == 0){ // this is an observed time to new infection
      if(Drug[i] == 0){
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | reinfection_rate);
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
        log_probs[1] = log_p[Ind] + exponential_lpdf(Durations[i] | reinfection_rate);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lpdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lpdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){
        Ind = ID_mapped_to_PMQ_rank[ID_of_Patient[i]];
        // Chloroquine + Primaquine
        // reInfection
        log_probs[1] = log_p_PMQ[Ind] + exponential_lpdf(Durations[i] | reinfection_rate);
        // Early reLapse
        log_probs[2] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_EarlyL + weibull_lpdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_1m_EarlyL + exponential_lpdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p_PMQ[Ind] + log_c1_CQ_PMQ + exponential_lpdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);

      }
    }
    if(Censored[i] == 1){ // this is a right censored time to new infection
      // This is the unobserved case (data are right censored)
      if(Drug[i] == 0){
        Ind = ID_mapped_to_noPMQ_rank[ID_of_Patient[i]];
        // This is the Artesunate monotherapy: reLapses or reInfections
        // reInfection
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | reinfection_rate);
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
        log_probs[1] = log_p[Ind] + exponential_lccdf(Durations[i] | reinfection_rate);
        // Early reLapse
        log_probs[2] = log_1m_p[Ind] + log_1m_c1_CQ + log_EarlyL + weibull_lccdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p[Ind] + log_1m_c1_CQ + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p[Ind] + log_c1_CQ + exponential_lccdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);
      }
      if(Drug[i] == 2){
        Ind = ID_mapped_to_PMQ_rank[ID_of_Patient[i]];
        // Chloroquine + Primaquine:
        // reInfection
        log_probs[1] = log_p_PMQ[Ind] + exponential_lccdf(Durations[i] | reinfection_rate);
        // Early reLapse
        log_probs[2] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_EarlyL + weibull_lccdf(Durations[i] | CQ_shape, CQ_scale);
        // Late reLapse
        log_probs[3] = log_1m_p_PMQ[Ind] + log_1m_c1_CQ_PMQ + log_1m_EarlyL + exponential_lccdf(Durations[i] | gamma);
        // recrudescence
        log_probs[4] = log_1m_p_PMQ[Ind] + log_c1_CQ_PMQ + exponential_lccdf(Durations[i] | lambda_recrud);

        log_lik[i] = log_sum_exp(log_probs);
      }
    }
  }
}

'

Timing_Model2 = stan_model(model_code = mod2)