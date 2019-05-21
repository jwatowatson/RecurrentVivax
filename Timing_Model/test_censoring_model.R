modtest = 
  '
data {
  // Things that are data structure and strictly data items:
  int<lower=0>                Neps;                       // The total number of observed durations
  real<lower=0>               Durations[Neps];            // Vector of durations
  int<lower=0,upper=1>       Censored[Neps];             // -1 is left censored; 0: observed; 1: right censored
  // Hyper parameters that are part of the `data` section
  real<lower=0>   Hyper_lambda_shape;     // Mean of the prior on time to reinfection
  real<lower=0>   Hyper_lambda_rate;      // SD of the prior on time to reinfection
}

parameters {
  real<lower=0>         lambda;             // reInfection rate
}

model {
  // ********* Prior *********
  lambda ~ gamma(Hyper_lambda_shape, Hyper_lambda_rate);
  
  // ********* Likelihood *********
  for(i in 1:Neps){

    if(Censored[i] == 0){ // this is an observed time to new infection
      target += exponential_lpdf(Durations[i] | lambda);
    } 
    if(Censored[i] == 1){ // this is a right censored time to new infection
      target += exponential_lccdf(Durations[i] | lambda);
    }
  }
} 
'


TestingModel = stan_model(model_code = modtest)
