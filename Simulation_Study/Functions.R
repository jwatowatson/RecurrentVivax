##************
generate_relapse_time = function(drug, params){
  
  type=NA
  next_event = NA
  if(runif(1) < inv.logit(params$logit_EarlyL)){ 
    # early relapse event
    if(drug %in% c('CHQ/PMQ', 'CHQ')){
      next_event = round(rweibull(1, shape = params$CQ_shape, scale = params$CQ_scale))
    } else if (drug=='AS') {
      next_event = round(rweibull(1, shape = params$AS_shape, scale = params$AS_scale))
    }
    type = 'EarlyRelapse'
  } else {
    # late relapse event
    next_event = round(rexp(1, rate = params$gamma))
    type = 'LateRelapse'
  }
  return(list(time=next_event,type=type))
}

generate_reinfection_time = function(params, Study_Period=NA){
  if(!is.na(Study_Period) & Study_Period==2) {
    lambda = params$lambda*params$rate_decrease
  } else {
    lambda = params$lambda
  }
  next_event = round(rexp(1, rate = lambda))
  return(next_event)
}

generate_reinfection_time_seasonal = function(params, week_time, 
                                              seasonal_sampling_vector, Study_Period){
  if(!is.na(Study_Period) & Study_Period==2) {
    lambda = params$lambda*params$rate_decrease
  } else {
    lambda = params$lambda
  }
  # re-adjust the seasonal sampling distribution to reflect current episode time of year
  seasonal_sampling_vector=seasonal_sampling_vector-week_time
  ind = seasonal_sampling_vector<0
  seasonal_sampling_vector[ind]=seasonal_sampling_vector[ind]+52
  # sample new event by sampling the year and then adding the weeks
  next_event = floor(rexp(1, rate = lambda)/360)*360 +
    7*sample(seasonal_sampling_vector,1)
  next_event = round(next_event)
  return(next_event)
}

generate_recrudescence_time = function(drug, params){
  next_event = round(rexp(1, rate = params$lambda_recrud))
  return(next_event)
}

simulate_dataset = function(N_PMQ, N_CQ, N_AS, 
                            FUP_time=300,
                            data_generation_function, 
                            params,
                            seasonal_sampling_vector=NA, Study_Period){
  drugs = c('CHQ/PMQ','CHQ','AS')
  Ns = c(N_PMQ, N_CQ, N_AS)
  N = sum(Ns)
  
  if(length(FUP_time) == 1){
    FUP_time = rep(FUP_time, N)
  }
  if(length(FUP_time) != N) {
    stop('FUP_time vector needs to be same length as all simulated patients')
  }
  if(length(Study_Period) != N) {
    stop('Study_Period vector needs to be same length as all simulated patients')
  }
  
  # set up data structure
  id = 1
  sim_data = data.frame(time_events=NA, censor_status=NA,
                        ID=NA, drug=NA, true_p=NA,true_event_type=NA,
                        Study_Period = NA)
  # iterate over the 3 drug types and over each patient in the drug arm
  patient_i = 1 # keep track of patient index for correct follow-up time
  for(i in 1:length(drugs)){
    drug=drugs[i]
    if(Ns[i]>0){
      for(j in 1:Ns[i]){
        dat = data_generation_function(params = params,
                                       follow_up = FUP_time[patient_i],
                                       drug = drug,
                                       ID=id, Study_Period=Study_Period[patient_i],
                                       seasonal_sampling_vector=seasonal_sampling_vector)
        id=id+1
        sim_data = rbind(sim_data,dat)
        patient_i = patient_i+1
      }
    }
    
  }
  sim_data = sim_data[-1,]
  sim_data = filter(sim_data, !(censor_status > -1 & time_events < 1))
  
  # Make data structures that will be given to the stan model
  
  N_PMQ = sum(sim_data$drug[!duplicated(sim_data$ID)] == 'CHQ/PMQ')
  # Turn drug into a numeric vector
  sim_data$numeric_drug = as.integer(revalue(sim_data$drug,
                                             c('AS'='0','CHQ'='1','CHQ/PMQ'='2')))
  
  # Create a vector that maps the nth person to 0 (always received Primaquine) 
  # or to their rank for those who didn't receive primaquine
  noPMQ_ind = which(sim_data$drug != "CHQ/PMQ")
  N_noPMQ = length(unique(sim_data$ID[noPMQ_ind])) #number of IDs without PMQ
  ID_mapped_to_noPMQ_rank = rep(0, N)
  
  index = 1
  for(id in unique(sim_data$ID)){
    ind = which(sim_data$ID==id)
    if(!2 %in% sim_data$numeric_drug[ind]){
      ID_mapped_to_noPMQ_rank[id] = index
      index = index + 1
    }
  }
  
  
  ID_mapped_to_PMQ_rank = rep(0, N)
  
  index = 1
  for(id in unique(sim_data$ID)){
    ind = which(sim_data$ID==id)
    if(2 %in% sim_data$numeric_drug[ind]){
      ID_mapped_to_PMQ_rank[id] = index
      index = index + 1
    }
  }
  # This is for the stan model
  Simdata =list(N         = as.integer(N),
                #Number of individuals
                Neps      = as.integer(nrow(sim_data)),
                #Number of durations
                N_noPMQ   = as.integer(N_noPMQ),
                # Number of individuals who do not receive PMQ
                N_PMQ     = as.integer(N_PMQ),
                # Number of individuals who do recieve PMQ
                Durations = as.double(sim_data$time_events),
                # Time to recurrence or time to censoring
                Censored  = as.integer(sim_data$censor_status),
                # If the duration is right censored or not
                Drug      = sim_data$numeric_drug,
                # drug coded as an integer
                ID_of_Patient = sim_data$ID,
                # the ID corresponding to each time interval
                ID_mapped_to_noPMQ_rank = ID_mapped_to_noPMQ_rank,
                # the index mapping PMQ individuals to their rank
                ID_mapped_to_PMQ_rank = ID_mapped_to_PMQ_rank,
                # the index mapping PMQ individuals to their rank
                Study_Period = sim_data$Study_Period
  )
  # This is to check outputs of stan model: ground simulation truth
  Simulation_truth = list(ID = sim_data$ID,
                          True_state = sim_data$true_event_type,
                          True_p = sim_data$true_p,
                          Drug = sim_data$drug)
  out = list(Simulation_truth=Simulation_truth, Simdata=Simdata)
  return(out)
}

##************ Model 1 assumptions for generating data **********
generate_patient_data_Model1 = function(params, follow_up=360, drug, ID=1,
                                        seasonal_sampling_vector=NA, Study_Period){
  time_events = c()
  censor_status = c()
  true_event_type = c()
  total_time = 0
  passed_EOF = FALSE # passed end of follow-up
  # randomly generate the propensity to relapse
  p = inv.logit(rnorm(1,mean = params$logit_mean_p,sd = params$logit_sd_p))
  if(drug=='CHQ') { 
    c1 = inv.logit(params$logit_c1_CQ) 
  } else if(drug=='AS') { 
    c1 = inv.logit(params$logit_c1_AS) 
  }
  
  R_event = NA
  while(!passed_EOF){
    if(drug == 'CHQ/PMQ'){ 
      # primaquine: 100% efficacy, reinfection is the only option
      next_event = generate_reinfection_time(params = params,Study_Period = Study_Period)
      R_event = 'Reinfection'
    } else if (drug %in% c('CHQ','AS')){
      if(runif(1) < p){
        next_event = generate_reinfection_time(params,Study_Period = Study_Period)
        R_event = 'Reinfection'
      } else {
        if(runif(1) < c1){ 
          next_event = generate_recrudescence_time(drug = drug, params = params)
          R_event = 'Recrudescence'
        } else {
          lapse = generate_relapse_time(drug = drug, params = params)
          next_event = lapse$time
          R_event = lapse$type
        }
      }
    }
    if(is.na(next_event)) stop('next_event NA this is a bug')
    total_time = total_time + next_event
    
    if(total_time <= follow_up){
      censor_status = c(censor_status, 0)
      time_events = c(time_events, next_event)
      true_event_type = c(true_event_type, R_event)
    } else {
      censor_status = c(censor_status, 1)
      passed_EOF = TRUE
      time_events = c(time_events,follow_up - (total_time-next_event))
      true_event_type = c(true_event_type, R_event)
    }
  }
  
  results = data.frame(time_events=time_events, censor_status=censor_status,
                       ID=ID, drug=drug, true_p = p,
                       true_event_type = true_event_type, Study_Period=Study_Period)
  return(results)
}

##************ Model 2 assumptions for generating data **********
generate_patient_data_Model2 = function(params, follow_up=360, drug, ID=1,
                                        seasonal_sampling_vector=NA, Study_Period=NA){
  time_events = c()
  censor_status = c()
  true_event_type = c()
  total_time = 0
  passed_EOF = FALSE # passed end of follow-up
  # randomly generate the propensity to relapse
  p = inv.logit(rnorm(1,mean = params$logit_mean_p,sd = params$logit_sd_p))
  p_PMQ = inv.logit(rnorm(1,mean = params$logit_mean_p_PMQ,sd = params$logit_sd_p_PMQ))
  
  if(drug=='CHQ' | drug == 'CHQ/PMQ') { 
    c1 = inv.logit(params$logit_c1_CQ)
  } else if(drug=='AS') { 
    c1 = inv.logit(params$logit_c1_AS)
  }
  
  while(!passed_EOF){
    if(drug == 'CHQ/PMQ'){ 
      true_p = p_PMQ
      if(runif(1) < p_PMQ){ #reinfection
        next_event = generate_reinfection_time(params,Study_Period=Study_Period)
        R_event = 'Reinfection'
      } else {
        if(runif(1) < c1){ # recrudescence
          next_event = generate_recrudescence_time(drug = drug, params = params)
          R_event = 'Recrudescence'
        } else { # relapse
          lapse = generate_relapse_time(drug = drug, params = params)
          next_event = lapse$time
          R_event = lapse$type
        }
      }
      
    } else if (drug %in% c('CHQ','AS')){
      true_p = p
      if(runif(1) < p){ #reinfection
        next_event = generate_reinfection_time(params,Study_Period=Study_Period)
        R_event = 'Reinfection'
      } else {
        if(runif(1) < c1){ #recrudescence
          next_event = generate_recrudescence_time(drug = drug, params = params)
          R_event = 'Recrudescence'
        } else { #relapse
          lapse = generate_relapse_time(drug = drug, params = params)
          next_event = lapse$time
          R_event = lapse$type
        }
      }
    }
    
    total_time = total_time + next_event
    
    if(total_time <= follow_up){
      censor_status = c(censor_status, 0)
      time_events = c(time_events, next_event)
      
    } else {
      censor_status = c(censor_status, 1)
      passed_EOF = TRUE # stop simulation
      time_events = c(time_events,follow_up - (total_time-next_event))
    }
    true_event_type = c(true_event_type, R_event)
  }
  
  results = data.frame(time_events=time_events, censor_status=censor_status,
                       ID=ID, drug=drug, true_p = true_p,
                       true_event_type = true_event_type, 
                       Study_Period=Study_Period)
  return(results)
}

##************ Model 1 assumptions for generating data with additional seasonality **********
# modify the simulation routine under the assumptions of model 1 to include seasonal reinfection
generate_patient_data_Model1_Seasonal = function(params, follow_up=360, drug, ID=1,
                                                 seasonal_sampling_vector, Study_Period){
  time_events = c()
  censor_status = c()
  true_event_type = c()
  week_time = sample(seasonal_sampling_vector,size = 1) # enrollment date
  total_time = 0
  passed_EOF = FALSE # passed end of follow-up
  # randomly generate the propensity to relapse
  p = inv.logit(rnorm(1,mean = params$logit_mean_p,sd = params$logit_sd_p))
  if(drug=='CHQ') { c1 = params$logit_c1_CQ } else if(drug=='AS') { c1 = params$logit_c1_AS}
  
  R_event = NA
  while(!passed_EOF){
    if(drug == 'CHQ/PMQ'){ 
      # primaquine: 100% efficacy, reinfection is the only option
      next_event = generate_reinfection_time_seasonal(params = params, 
                                                      week_time=week_time, 
                                                      seasonal_sampling_vector,
                                                      Study_Period=Study_Period)
      R_event = 'Reinfection'
    } else if (drug %in% c('CHQ','AS')){
      if(runif(1) < p){
        next_event = generate_reinfection_time_seasonal(params = params, 
                                                        week_time=week_time, 
                                                        seasonal_sampling_vector,
                                                        Study_Period=Study_Period)
        R_event = 'Reinfection'
      } else {
        if(runif(1) < c1){ 
          next_event = generate_recrudescence_time(drug = drug, params = params)
          R_event = 'Recrudescence'
        } else {
          lapse = generate_relapse_time(drug = drug, params = params)
          next_event = lapse$time
          R_event = lapse$type
        }
      }
    }
    week_time = (next_event - floor(next_event/360)*360)/7
    if(week_time>52) stop('week time is > 52 this is a bug')
    if(is.na(next_event)) stop('next_event is NA: this is a bug somewhere')
    total_time = total_time + next_event
    
    if(total_time <= follow_up){
      censor_status = c(censor_status, 0)
      time_events = c(time_events, next_event)
      true_event_type = c(true_event_type, R_event)
    } else {
      censor_status = c(censor_status, 1)
      passed_EOF = TRUE
      time_events = c(time_events,follow_up - (total_time-next_event))
      true_event_type = c(true_event_type, R_event)
    }
  }
  
  results = data.frame(time_events=time_events, censor_status=censor_status,
                       ID=ID, drug=drug, true_p = p,
                       true_event_type = true_event_type, Study_Period=Study_Period)
  return(results)
}

##************ Model 2 assumptions for generating data with seasonality **********
generate_patient_data_Model2_Seasonal = function(params, follow_up=360, 
                                                 drug, ID=1, Study_Period,
                                                 seasonal_sampling_vector){
  time_events = c()
  censor_status = c()
  true_event_type = c()
  week_time = sample(seasonal_sampling_vector,size = 1) # enrollment date
  total_time = 0
  passed_EOF = FALSE # passed end of follow-up
  # randomly generate the propensity to relapse
  p = inv.logit(rnorm(1,mean = params$logit_mean_p,sd = params$logit_sd_p))
  p_PMQ = inv.logit(rnorm(1,mean = params$logit_mean_p_PMQ,sd = params$logit_sd_p_PMQ))
  
  if(drug=='CHQ' | drug == 'CHQ/PMQ') { 
    c1 = params$logit_c1_CQ 
  } else if(drug=='AS') { 
    c1 = params$logit_c1_AS
  }
  
  while(!passed_EOF){
    if(drug == 'CHQ/PMQ'){ 
      true_p = p_PMQ
      if(runif(1) < p_PMQ){ #reinfection
        next_event = generate_reinfection_time_seasonal(params,week_time = week_time,
                                                        seasonal_sampling_vector = seasonal_sampling_vector,
                                                        Study_Period=Study_Period)
        R_event = 'Reinfection'
      } else {
        if(runif(1) < c1){ # recrudescence
          next_event = generate_recrudescence_time(drug = drug, params = params)
          R_event = 'Recrudescence'
        } else { # relapse
          lapse = generate_relapse_time(drug = drug, params = params)
          next_event = lapse$time
          R_event = lapse$type
        }
      }
      
    } else if (drug %in% c('CHQ','AS')){
      true_p = p
      if(runif(1) < p){ #reinfection
        next_event = generate_reinfection_time_seasonal(params,week_time = week_time,
                                                        seasonal_sampling_vector = seasonal_sampling_vector,
                                                        Study_Period=Study_Period)
        R_event = 'Reinfection'
      } else {
        if(runif(1) < c1){ #recrudescence
          next_event = generate_recrudescence_time(drug = drug, params = params)
          R_event = 'Recrudescence'
        } else { #relapse
          lapse = generate_relapse_time(drug = drug, params = params)
          next_event = lapse$time
          R_event = lapse$type
        }
      }
    }
    
    week_time = (next_event - floor(next_event/360)*360)/7
    if(week_time>52) stop('week time is > 52 this is a bug')
    
    total_time = total_time + next_event
    
    if(total_time <= follow_up){
      censor_status = c(censor_status, 0)
      time_events = c(time_events, next_event)
      
    } else {
      censor_status = c(censor_status, 1)
      passed_EOF = TRUE # stop simulation
      time_events = c(time_events,follow_up - (total_time-next_event))
    }
    true_event_type = c(true_event_type, R_event)
  }
  
  results = data.frame(time_events=time_events, censor_status=censor_status,
                       ID=ID, drug=drug, true_p = true_p,
                       true_event_type = true_event_type, Study_Period=Study_Period)
  return(results)
}


#####************* Initial value generation *************
# function that generates values from the Prior for model 1
generate_inital_values_M1 = function(Prior_params_M1, N_noPMQ, N_PMQ){
  inits = list(logit_p = rnorm(n = N_noPMQ, 
                               mean = Prior_params_M1$Hyper_logit_mean_p_mean,
                               sd = Prior_params_M1$Hyper_logit_mean_p_sd),
               logit_c1_AS = rnorm(n = 1,
                                   Prior_params_M1$Hyper_logit_c1_mean,
                                   Prior_params_M1$Hyper_logit_c1_sd),
               logit_c1_CQ = rnorm(n = 1,
                                   mean = Prior_params_M1$Hyper_logit_c1_mean,
                                   sd = Prior_params_M1$Hyper_logit_c1_sd),
               logit_EarlyL = rnorm(n = 1, 
                                    mean = Prior_params_M1$Early_L_logit_mean,
                                    sd = Prior_params_M1$Early_L_logit_sd),
               lambda = rgamma(n = 1, 
                               shape = Prior_params_M1$Hyper_lambda_shape,
                               rate = Prior_params_M1$Hyper_lambda_rate),
               gamma = rgamma(n = 1, 
                              shape = Prior_params_M1$Hyper_gamma_shape,
                              rate = Prior_params_M1$Hyper_gamma_rate),
               lambda_recrud = rgamma(n = 1, 
                                      shape = Prior_params_M1$Hyper_lambda_recrud_shape,
                                      rate = Prior_params_M1$Hyper_lambda_recrud_rate),
               logit_mean_p = rnorm(n = 1, 
                                    mean = Prior_params_M1$Hyper_logit_mean_p_mean,
                                    sd = Prior_params_M1$Hyper_logit_mean_p_sd),
               logit_sd_p = rexp(n = 1, 
                                 rate = Prior_params_M1$Hyper_logit_sd_p_lambda),
               AS_shape = rnorm(n = 1,
                                mean = Prior_params_M1$Hyper_AS_shape_mean,
                                sd = Prior_params_M1$Hyper_AS_shape_sd),
               AS_scale = rnorm(n = 1,
                                mean = Prior_params_M1$Hyper_AS_scale_mean,
                                sd = Prior_params_M1$Hyper_AS_scale_sd),
               CQ_shape = rnorm(n = 1,
                                mean = Prior_params_M1$Hyper_CQ_shape_mean,
                                sd = Prior_params_M1$Hyper_CQ_shape_sd),
               CQ_scale = rnorm(n = 1,
                                mean = Prior_params_M1$Hyper_CQ_scale_mean,
                                sd = Prior_params_M1$Hyper_CQ_scale_sd),
               rate_decrease = rtruncnorm(n = 1, a = 0, 
                                          mean = Prior_params_M1$Hyper_mean_rate_decrease, 
                                          sd = Prior_params_M1$Hyper_sd_rate_decrease))
  return(inits)
}
# function that generates values from the Prior for model 2
generate_inital_values_M2 = function(Prior_params_M2, N_noPMQ, N_PMQ){
  # the models are nested so we generate for model 1 first
  inits = generate_inital_values_M1(Prior_params_M2, N_noPMQ = N_noPMQ)
  logit_p_PMQ = list(rnorm(N_PMQ,
                           mean = Prior_params_M2$Hyper_logit_mean_pPMQ_mean,
                           sd = Prior_params_M2$Hyper_logit_mean_pPMQ_sd))
  inits = append(inits, logit_p_PMQ)
  names(inits)[length(inits)] = 'logit_p_PMQ'
  inits = c(inits,
            logit_mean_p_PMQ = rnorm(1,
                                     mean = Prior_params_M2$Hyper_logit_mean_pPMQ_mean,
                                     sd = Prior_params_M2$Hyper_logit_mean_pPMQ_sd),
            logit_sd_p_PMQ = rexp(n = 1,
                                  rate = Prior_params_M2$Hyper_logit_sd_pPMQ_lambda))
  
  return(inits)
}



###************ Plotting output of models with respect to priors *********
## Output of stan model 1
plot_output_model1 = function(thetas, Simulation_truth, Simdata_Model, Prior_params_M1){
  par(las=1)
  par(mfrow=c(2,2))
  
  # Time to reinfection
  hist(thetas$lambda, main='Reinfection rate (Period 1)', xlab='lambda', 
       freq = F, breaks = 20)
  abline(v=params_M1$lambda,col='red',lwd=3)
  abline(v=mean(thetas$lambda),col='blue',lwd=3)
  xs=quantile(x = thetas$lambda,probs = seq(0,1,by=0.01))
  lines(xs, dgamma(x = xs, 
                   Prior_params_M1$Hyper_lambda_shape,
                   Prior_params_M1$Hyper_lambda_rate),
        lwd=3,col='red')
  hist(thetas$rate_decrease, main='Decrease in reinfection rate (Period 2)', xlab='decrease proportion', 
       freq = F, breaks = 20)
  abline(v=params_M1$rate_decrease,col='red',lwd=3)
  abline(v=mean(thetas$rate_decrease),col='blue',lwd=3)
  xs=quantile(x = thetas$rate_decrease,probs = seq(0,1,by=0.01))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M1$Hyper_mean_rate_decrease,
                  Prior_params_M1$Hyper_sd_rate_decrease),
        lwd=3,col='red')
  
  # Time to late relapse
  hist(thetas$gamma, main='Mean time to late reLapse', 
       freq=F,xlab='1/gamma (days)')
  abline(v=params_M1$gamma,col='red',lwd=3)
  xs=quantile(x = thetas$gamma,probs = seq(0,1,by=0.005))
  lines(xs, dgamma(x = xs, 
                   Prior_params_M1$Hyper_gamma_shape,
                   Prior_params_M1$Hyper_gamma_rate),
        lwd=3,col='red')
  abline(v=mean(thetas$gamma),col='blue',lwd=3)
  
  # Proportion early/late relapse
  hist(thetas$logit_EarlyL, main='Logit early relapse', xlab='', freq = F)
  abline(v=params_M1$logit_EarlyL,col='red',lwd=3)
  xs=quantile(x = thetas$logit_EarlyL,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M1$Early_L_logit_mean,
                  Prior_params_M1$Early_L_logit_sd),
        lwd=3,col='red')
  abline(v=mean(thetas$logit_EarlyL),col='blue',lwd=3)
  
  par(mfrow=c(2,2))
  # recrudescence proportion
  hist(thetas$logit_c1_AS, main='Logit c1 AS', xlab='',freq=F)
  abline(v=logit(params_M1$logit_c1_AS),col='red',lwd=3)
  abline(v=mean(thetas$logit_c1_AS),col='blue',lwd=3)
  xs=quantile(x = thetas$logit_c1_AS,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M1$Hyper_logit_c1_mean,
                  Prior_params_M1$Hyper_logit_c1_sd),
        lwd=3,col='red')
  
  hist(thetas$logit_c1_CQ, main='Logit c1 CQ', xlab='',freq=F)
  lines(xs, dnorm(x = xs, 
                  Prior_params_M1$Hyper_logit_c1_mean,
                  Prior_params_M1$Hyper_logit_c1_sd),
        lwd=3,col='red')
  abline(v=logit(params_M1$logit_c1_CQ),col='red',lwd=3)
  abline(v=mean(thetas$logit_c1_CQ),col='blue',lwd=3)
  
  hist(thetas$lambda_recrud, main='Recrudescence rate', xlab='',freq=F)
  abline(v=params_M1$lambda_recrud,col='red',lwd=3)
  abline(v=mean(thetas$lambda_recrud),col='blue',lwd=3)
  xs=quantile(x = thetas$lambda_recrud,probs = seq(0,1,by=0.005))
  lines(xs, dgamma(x = xs, 
                   shape=Prior_params_M1$Hyper_lambda_recrud_shape,
                   rate=Prior_params_M1$Hyper_lambda_recrud_rate),
        lwd=3,col='red')
  
  
  par(mfrow=c(2,2))
  # shape and scale for AS early relapse
  hist(thetas$AS_shape, main='AS shape', xlab='',freq=F)
  abline(v=(params_M1$AS_shape),col='red',lwd=3)
  abline(v=mean(thetas$AS_shape),col='blue',lwd=3)
  xs=quantile(x = thetas$AS_shape,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M1$Hyper_AS_shape_mean,
                  sd=Prior_params_M1$Hyper_AS_shape_sd),
        lwd=3,col='red')
  
  hist(thetas$AS_scale, main='AS scale', xlab='',freq=F)
  abline(v=params_M1$AS_scale,col='red',lwd=3)
  abline(v=mean(thetas$AS_scale),col='blue',lwd=3)
  xs=quantile(x = thetas$AS_scale,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M1$Hyper_AS_scale_mean,
                  sd=Prior_params_M1$Hyper_AS_scale_sd),
        lwd=3,col='red')
  # shape and scale for CQ early relapse
  hist(thetas$CQ_shape, main='CQ shape', xlab='',freq=F)
  abline(v=(params_M1$CQ_shape),col='red',lwd=3)
  abline(v=mean(thetas$CQ_shape),col='blue',lwd=3)
  xs=quantile(x = thetas$CQ_shape,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M1$Hyper_CQ_shape_mean,
                  sd=Prior_params_M1$Hyper_CQ_shape_sd),
        lwd=3,col='red')
  
  hist(thetas$CQ_scale, main='CQ scale', xlab='',freq=F)
  abline(v=(params_M1$CQ_scale),col='red',lwd=3)
  abline(v=mean(thetas$CQ_scale),col='blue',lwd=3)
  xs=quantile(x = thetas$CQ_scale,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M1$Hyper_CQ_scale_mean,
                  sd=Prior_params_M1$Hyper_CQ_scale_sd),
        lwd=3,col='red')
  
  par(mfrow=c(2,2))
  # proportion reinfection
  hist(thetas$logit_mean_p, main='Logit mean p', xlab='', freq=F)
  abline(v=params_M1$logit_mean_p,col='red',lwd=3)
  abline(v=mean(thetas$logit_mean_p),col='blue',lwd=3)
  xs=quantile(x = thetas$logit_mean_p,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M1$Hyper_logit_mean_p_mean,
                  Prior_params_M1$Hyper_logit_mean_p_sd),
        lwd=3,col='red')
  
  hist(thetas$logit_sd_p, main='Logit sd p', xlab='',freq=F)
  abline(v=(params_M1$logit_sd_p),col='red',lwd=3)
  abline(v=mean(thetas$logit_sd_p),col='blue',lwd=3)
  xs=quantile(x = thetas$logit_sd_p,probs = seq(0,1,by=0.005))
  lines(xs, dexp(x = xs,Prior_params_M1$Hyper_logit_sd_p_lambda),
        lwd=3,col='red')
  # individual estimates versus true params for proportion reinfection
  p_estimates = apply(thetas$logit_p,2,mean)
  true_p_estimates= logit(Simulation_truth$True_p[!duplicated(Simulation_truth$ID) &
                                                    Simdata_Model$Drug<2])
  plot(p_estimates, true_p_estimates, xlab = 'Estimate of logit p (individual)',
       ylab = 'True logit p (individual)' )
  lines(c(-10,10),c(-10,10),lwd=2)
  plot(p_estimates, p_estimates-true_p_estimates,
       xlab = 'Estimate of logit p (individual)',
       ylab = 'Difference between estimate and true value')
  abline(h=0,lwd=2)
}

## Output of stan model 2
plot_output_model2 = function(thetas, Simulation_truth, Simdata_Model, Prior_params_M2){
  par(las=1)
  par(mfrow=c(2,2))
  
  # Time to reinfection
  hist(thetas$lambda, main='Reinfection rate', xlab='lambda', 
       freq = F, breaks = 20)
  abline(v=params_M2$lambda,col='red',lwd=3)
  abline(v=mean(thetas$lambda),col='blue',lwd=3)
  xs=quantile(x = thetas$lambda,probs = seq(0,1,by=0.01))
  lines(xs, dgamma(x = xs, 
                   Prior_params_M2$Hyper_lambda_shape,
                   Prior_params_M2$Hyper_lambda_rate),
        lwd=3,col='red')
  # Time to late relapse
  hist(thetas$gamma, main='Mean time to late reLapse', freq=F,
       xlab='1/gamma (days)')
  abline(v=params_M2$gamma,col='red',lwd=3)
  xs=quantile(x = thetas$gamma,probs = seq(0,1,by=0.005))
  lines(xs, dgamma(x = xs, 
                   Prior_params_M2$Hyper_gamma_shape,
                   Prior_params_M2$Hyper_gamma_rate),
        lwd=3,col='red')
  abline(v=mean(thetas$gamma),col='blue',lwd=3)
  
  # Proportion early/late relapse
  hist(thetas$logit_EarlyL, main='Logit early relapse', xlab='', freq = F)
  abline(v=logit(params_M2$logit_EarlyL),col='red',lwd=3)
  xs=quantile(x = thetas$logit_EarlyL,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M2$Early_L_logit_mean,
                  Prior_params_M2$Early_L_logit_sd),
        lwd=3,col='red')
  abline(v=mean(thetas$logit_EarlyL),col='blue',lwd=3)
  
  par(mfrow=c(2,2))
  # recrudescence proportion
  hist(thetas$logit_c1_AS, main='Logit c1 AS', xlab='',freq=F)
  abline(v=logit(params_M2$logit_c1_AS),col='red',lwd=3)
  abline(v=mean(thetas$logit_c1_AS),col='blue',lwd=3)
  xs=quantile(x = thetas$logit_c1_AS,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M2$Hyper_logit_c1_mean,
                  Prior_params_M2$Hyper_logit_c1_sd),
        lwd=3,col='red')
  
  hist(thetas$logit_c1_CQ, main='Logit c1 CQ', xlab='',freq=F)
  lines(xs, dnorm(x = xs, 
                  Prior_params_M2$Hyper_logit_c1_mean,
                  Prior_params_M2$Hyper_logit_c1_sd),
        lwd=3,col='red')
  abline(v=logit(params_M2$logit_c1_CQ),col='red',lwd=3)
  abline(v=mean(thetas$logit_c1_CQ),col='blue',lwd=3)
  
  hist(thetas$lambda_recrud, main='Recrudescence rate', xlab='',freq=F)
  abline(v=params_M2$lambda_recrud,col='red',lwd=3)
  abline(v=mean(thetas$lambda_recrud),col='blue',lwd=3)
  xs=quantile(x = thetas$lambda_recrud,probs = seq(0,1,by=0.005))
  lines(xs, dgamma(x = xs, 
                   shape=Prior_params_M2$Hyper_lambda_recrud_shape,
                   rate=Prior_params_M2$Hyper_lambda_recrud_rate),
        lwd=3,col='red')
  
  
  par(mfrow=c(2,2))
  # shape and scale for AS early relapse
  hist(thetas$AS_shape, main='AS shape', xlab='',freq=F)
  abline(v=(params_M2$AS_shape),col='red',lwd=3)
  abline(v=mean(thetas$AS_shape),col='blue',lwd=3)
  xs=quantile(x = thetas$AS_shape,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M2$Hyper_AS_shape_mean,
                  sd=Prior_params_M2$Hyper_AS_shape_sd),
        lwd=3,col='red')
  
  hist(thetas$AS_scale, main='AS scale', xlab='',freq=F)
  abline(v=params_M2$AS_scale,col='red',lwd=3)
  abline(v=mean(thetas$AS_scale),col='blue',lwd=3)
  xs=quantile(x = thetas$AS_scale,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M2$Hyper_AS_scale_mean,
                  sd=Prior_params_M2$Hyper_AS_scale_sd),
        lwd=3,col='red')
  # shape and scale for CQ early relapse
  hist(thetas$CQ_shape, main='CQ shape', xlab='',freq=F)
  abline(v=(params_M2$CQ_shape),col='red',lwd=3)
  abline(v=mean(thetas$CQ_shape),col='blue',lwd=3)
  xs=quantile(x = thetas$CQ_shape,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M2$Hyper_CQ_shape_mean,
                  sd=Prior_params_M2$Hyper_CQ_shape_sd),
        lwd=3,col='red')
  
  hist(thetas$CQ_scale, main='CQ scale', xlab='',freq=F)
  abline(v=(params_M2$CQ_scale),col='red',lwd=3)
  abline(v=mean(thetas$CQ_scale),col='blue',lwd=3)
  xs=quantile(x = thetas$CQ_scale,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  mean=Prior_params_M2$Hyper_CQ_scale_mean,
                  sd=Prior_params_M2$Hyper_CQ_scale_sd),
        lwd=3,col='red')
  
  par(mfrow=c(2,2))
  # proportion reinfection when no primaquine
  hist(thetas$logit_mean_p, main='Logit mean p', xlab='', freq=F)
  abline(v=params_M2$logit_mean_p,col='red',lwd=3)
  abline(v=mean(thetas$logit_mean_p),col='blue',lwd=3)
  xs=quantile(x = thetas$logit_mean_p,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M2$Hyper_logit_mean_p_mean,
                  Prior_params_M2$Hyper_logit_mean_p_sd),
        lwd=3,col='red')
  
  # proportion reinfection when given primaquine
  hist(thetas$logit_mean_p_PMQ, main='Logit mean p', xlab='', freq=F)
  abline(v=params_M2$logit_mean_p_PMQ,col='red',lwd=3)
  abline(v=mean(thetas$logit_mean_p_PMQ),col='blue',lwd=3)
  xs=quantile(x = thetas$logit_mean_p_PMQ,probs = seq(0,1,by=0.005))
  lines(xs, dnorm(x = xs, 
                  Prior_params_M2$Hyper_logit_mean_pPMQ_mean,
                  Prior_params_M2$Hyper_logit_mean_pPMQ_sd),
        lwd=3,col='red')
  
  hist(thetas$logit_sd_p, main='Logit sd p', xlab='',freq=F)
  abline(v=(params_M2$logit_sd_p),col='red',lwd=3)
  abline(v=mean(thetas$logit_sd_p),col='blue',lwd=3)
  xs=quantile(x = thetas$logit_sd_p,probs = seq(0,1,by=0.005))
  lines(xs, dexp(x = xs,Prior_params_M2$Hyper_logit_sd_p_lambda),
        lwd=3,col='red')
  # individual estimates versus true params for proportion reinfection
  p_estimates = apply(thetas$logit_p,2,mean)
  pPMQ_estimates = apply(thetas$logit_p_PMQ,2,mean)
  true_p_estimates= logit(Simulation_truth$True_p[!duplicated(Simulation_truth$ID) &
                                                    Simdata_Model$Drug<2])
  true_pPMQ_estimates= logit(Simulation_truth$True_p[!duplicated(Simulation_truth$ID) &
                                                       Simdata_Model$Drug==2])
  plot(c(p_estimates,pPMQ_estimates), 
       c(true_p_estimates, true_pPMQ_estimates),
       col = c(rep(1,length(p_estimates)),rep(2,length(pPMQ_estimates))),
       xlab = 'Estimate of logit p (individual)',
       ylab = 'True logit p (individual)' )
  lines(c(-10,10),c(-10,10),lwd=2)
  plot(c(p_estimates,pPMQ_estimates), 
       c(p_estimates,pPMQ_estimates)-c(true_p_estimates, true_pPMQ_estimates),
       col = c(rep(1,length(p_estimates)),rep(2,length(pPMQ_estimates))),
       xlab = 'Estimate of logit p (individual)',
       ylab = 'Difference between estimate and true value')
  abline(h=0,lwd=2)
}

##******* Wrapper function: simulate data; fit model; report point estimates ********
full_shebang = function(simulation_patient_function,
                        simulate_data_function,
                        stan_model,params,
                        init_values_function, 
                        N_AS, N_CQ, N_PMQ,
                        Prior_params, IT, WarmUp, Chains, thin,
                        params_interest,seasonal_sampling_vector=NA,
                        Study_Period){
  
  # simulate data
  out = simulate_dataset(N_PMQ = N_PMQ, N_CQ = N_CQ, 
                         N_AS = N_AS,FUP_time = FUP_time,
                         data_generation_function = simulation_patient_function,
                         params = params,
                         Study_Period = Study_Period,
                         seasonal_sampling_vector=seasonal_sampling_vector)
  
  Simdata_Model = out$Simdata
  Simulation_truth = out$Simulation_truth
  
  # fit stan model
  mod_Fit = sampling(stan_model,
                     data = c(Simdata_Model, Prior_params),
                     iter = IT, warmup = WarmUp,
                     chains = Chains, thin = thin, 
                     init = lapply(1:Chains, 
                                   FUN = function(x,Prior_params,N_noPMQ,N_PMQ){
                                     init_values_function(Prior_params,N_noPMQ,N_PMQ)
                                   }, 
                                   Prior_params = Prior_params, 
                                   N_noPMQ = Simdata_Model$N_noPMQ,
                                   N_PMQ = Simdata_Model$N_PMQ))
  # get theta samples
  thetas = extract(mod_Fit)
  
  summary_parameters = lapply(params_interest, FUN = function(x, thetas, params){
    theta_p_hat = mean(thetas[[x]])
    true_p = params[[x]]
    out = c(theta_p_hat=theta_p_hat,true_p=true_p)
    return(out)
  }, thetas=thetas, params = params)
  names(summary_parameters) = params_interest
  # the true failure rate in the primaquine arm
  true_failure = sum(Simulation_truth$True_state[Simulation_truth$Drug=="CHQ/PMQ"] != 'Reinfection')/sum(Simulation_truth$Drug=="CHQ/PMQ")
  
  label_probs = apply(thetas$prob_labels,c(2,3),mean)
  x=1-sum(label_probs[Simdata_Model$Drug==2,1])/sum(Simdata_Model$Drug==2)
  
  summary_stats = c(true_failure=true_failure,estimated_failure=x)
  sim_results = list(summary_stats=summary_stats,summary_parameters=summary_parameters)
  return(sim_results)
}
