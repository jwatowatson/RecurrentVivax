# This is model 2 but with an added seasonality component

mod3 = 
  '
data {
  // Things that are data structure and strictly data items:
  int<lower=0>                N;                        // The number of individuals
  int<lower=0>                Neps;                     // The total number of datapoints (left/right censored, and observed)
  int<lower=0>                N_noPMQ;                  // The total number of individuals who did not receive PMQ for at least one episode
  int<lower=0>                N_PMQ;                    // The total number of individuals who received PMQ for at least one episode
  real                        Durations[Neps];          // Vector of time since last infection
  real<lower=0>               WeekTime[Neps];           // Vector of Times in units of months
  int<lower=-1,upper=1>       Censored[Neps];           // -1 is left censored; 0: observed; 1: right censored
  int<lower=0,upper=2>        Drug[Neps];               // Treatment drug corresponding to each duration
  int<lower=1,upper=N>        ID_of_Episode[Neps];            // Index of individual groupings
  int<lower=0,upper=N_noPMQ>  ID_mapped_to_noPMQ_rank[N];     // Vector mapping the the individual ID (from 1 to N) to the random effects index: logit_p
  int<lower=0,upper=N_PMQ>    ID_mapped_to_PMQ_rank[N];       // Vector mapping the the individual ID (from 1 to N) to the random effects index: logit_p_PMQ
  
  // Hyper parameters that are part of the `data` section
  real<lower=0>   mu_inv_lambda;          // Mean of the prior on time to reinfection
  real<lower=0>   sigma_inv_lambda;       // SD of the prior on time to reinfection
  real<lower=0>   mu_inv_gamma;           // Mean of the prior on time to random reLapse
  real<lower=0>   sigma_inv_gamma;        // SD of the prior on time to random reLapse
  real            Early_L_logit_mean;     // mean of prior on mean of logit Early_L
  real<lower=0>   Early_L_logit_sd;       // sd of prior on mean of logit Early_L
  real            Hyper_logit_mean_p_PMQ; //
  real<lower=0>   Hyper_logit_sd_p_PMQ;   //
  real            Hyper_logit_mean_p;     // mean of prior on mean of logit p
  real<lower=0>   Hyper_logit_sd_p;       // sd of prior on mean of logit p
  real            Hyper_logit_c1_mean;
  real<lower=0>   Hyper_logit_c1_sd;
  real            beta_1_mean;            // Prior mean on amplitude of seasonality effect
  real<lower=0>   beta_1_sigma;           // Prior mean on SD of amplitude of seaosnality effect
  real            beta_0_mean;            // Prior mean time of peak reinfection season
  real<lower=0>   beta_0_sigma;           // Prior SD on time of peak reinfection season
  real<lower=0>   mu_AS_shape;      
  real<lower=0>   sigma_AS_shape;
  real<lower=0>   mu_AS_scale;
  real<lower=0>   sigma_AS_scale;
  real<lower=0>   mu_CQ_shape;
  real<lower=0>   sigma_CQ_shape;
  real<lower=0>   mu_CQ_scale;
  real<lower=0>   sigma_CQ_scale;
  real            Nweeks;
}

parameters {
  real                  logit_p[N_noPMQ];     // logit proportion of reInfections no radical cure
  real                  logit_p_PMQ[N_PMQ];   // The proportion of reInfections after radical cure
  real                  logit_c1_AS;          // proportion of recrudescences when not reinfection - AS monotherapy
  real                  logit_c1_CQ;          // proportion of recrudescences when not reinfection - CQ monotherapy
  real                  logit_c1_CQ_PMQ;      // proportion of recrudescences when not reinfection - CQ + PMQ
  real                  logit_EarlyL;         // Proportion of reLapses that are early/periodic
  real<lower=0>         inv_lambda;           // reInfection mean time
  real<lower=0>         inv_gamma;            // Late reLapse mean time
  real                  logit_mean_p;         // logit mean proportion of reInfection (no rad cure)
  real<lower=0>         logit_sd_p;           // logit sd of proportion of reInfection (no rad cure)
  real                  logit_mean_p_PMQ;     // logit mean proportion of reInfection (after rad cure)
  real<lower=0>         logit_sd_p_PMQ;       // logit sd of proportion of reInfection (after rad cure)
  real<lower=0>         Recrud_shape;         // Weibull shape parameter: recrudescence
  real<lower=0>         Recrud_scale;         // Weibull scale parameter: recrudescence
  real<lower=1>         AS_shape;             // Weibull shape parameter: Artesunate
  real<lower=0>         AS_scale;             // Weibull scale parameter: Artesunate
  real<lower=1>         CQ_shape;             // Weibull shape parameter: Chloroquine
  real<lower=0>         CQ_scale;             // Weibull scale parameter: Chloroquine
  real<lower=0,upper=6.3>      beta0;          // Seasonality centering parameter in units of weeks (rainy season about June ~ 24 weeks)
  real<lower=0>         beta1;                // Seasonality amplitude parameter (0: no seasonality)
}

transformed parameters{
  // Turn the inverse rates into rate parameters for the likelihood calculation
  real lambda = 1/inv_lambda;   
  real gamma = 1/inv_gamma;   

  // Compute the reInfection and recrudescence rates
  // We define log scale parameters from inverse logit transformation
  real log_p[N_noPMQ];
  real log_p_PMQ[N_PMQ];
  real log_1m_p[N_noPMQ];
  real log_1m_p_PMQ[N_PMQ];

  real log_c1_AS      = log_inv_logit(logit_c1_AS);
  real log_c1_CQ      = log_inv_logit(logit_c1_CQ);
  real log_c1_CQ_PMQ  = log_inv_logit(logit_c1_CQ_PMQ);
  real log_EarlyL     = log_inv_logit(logit_EarlyL);

  // The 1 minus log proportions
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
  beta0 ~ normal(beta_0_mean, beta_0_sigma);
  beta1 ~ normal(beta_1_mean, beta_1_sigma);

  inv_lambda ~ normal(mu_inv_lambda,sigma_inv_lambda);
  inv_gamma ~ normal(mu_inv_gamma,sigma_inv_gamma);

  logit_mean_p ~ normal(Hyper_logit_mean_p, Hyper_logit_sd_p);
  logit_sd_p ~ exponential(1);
  
  logit_mean_p_PMQ ~ normal(Hyper_logit_mean_p_PMQ,Hyper_logit_sd_p_PMQ);
  logit_sd_p_PMQ ~ exponential(1);

  Recrud_shape ~ normal(2, 1);
  Recrud_scale ~ normal(10, 1);

  AS_shape ~ normal(mu_AS_shape,sigma_AS_shape);
  AS_scale ~ normal(mu_AS_scale,sigma_AS_scale);

  CQ_shape ~ normal(mu_CQ_shape,sigma_CQ_shape);
  CQ_scale ~ normal(mu_CQ_scale,sigma_CQ_scale);

  // Random effects term for the reinfection after no radical cure
  // The probability of reinfection versus relapse when no radical cure
  logit_p ~ normal(logit_mean_p,logit_sd_p);

  // The probability of recrudescence if not reinfection
  // Same hyper prior for the three mixing weights
  logit_c1_AS ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);
  logit_c1_CQ ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);
  logit_c1_CQ_PMQ ~ normal(Hyper_logit_c1_mean,Hyper_logit_c1_sd);

  // The probability of early relapse
  logit_EarlyL ~ normal(Early_L_logit_mean,Early_L_logit_sd);

  // Random effects term for reinfection after radical cure
  // The probability of relapse after radical cure
  logit_p_PMQ ~ normal(logit_mean_p_PMQ,logit_sd_p_PMQ);

  // ********* Likelihood *********
  for(i in 1:Neps){

    if(Censored[i] == -1){ // This is left censored time interval (primary enrollment illness)
      // In this scenario we drop recrudescence as a possible label
      // We assume that people have not received PMQ outside of the trial 
      // The likelihood only depends on the time of enrolment
      // No Drug information 
      target += beta1*sin((2*pi()*WeekTime[i])/Nweeks + beta0) - beta1;
    }
  }
} 

generated quantities {

  vector[Neps] log_lik;
  matrix[Neps,4] prob_labels;
  real alpha_seasonal;

  // This computes respective densities of each mixture component
  for(i in 1:Neps){
    int Ind;
    vector[4] prob_labels_raw;
    
    if(Censored[i] == -1){ // Enrollment episode
      
      alpha_seasonal = exp(beta1*sin((2*pi()*WeekTime[i])/Nweeks + beta0));

      prob_labels_raw[1] = alpha_seasonal*inv_logit(logit_mean_p);
      // Early Relapse
      prob_labels_raw[2] = (1-alpha_seasonal*inv_logit(logit_mean_p))*exp(log_1m_c1_AS)*exp(log_EarlyL);
      // Late Relapse
      prob_labels_raw[3] = (1-alpha_seasonal*inv_logit(logit_mean_p))*exp(log_1m_c1_AS)*exp(log_1m_EarlyL);
      // Recrudescence
      prob_labels_raw[4] = (1-alpha_seasonal*inv_logit(logit_mean_p))*exp(log_c1_AS);
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
    real log_p_seasonal;
    real log_1m_p_seasonal;
    real log_p_PMQ_seasonal;
    real log_1m_p_PMQ_seasonal;

 
    if(Censored[i] == -1){ // primary episode no timing information
      log_lik[i] =  beta1*sin((2*pi()*WeekTime[i])/Nweeks + beta0);
    } 
  }
}

'

Timing_Model3 = stan_model(model_code = mod3)