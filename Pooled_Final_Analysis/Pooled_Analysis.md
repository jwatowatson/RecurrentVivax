---
title: "Pooled Analysis"
author: "Aimee Taylor and James Watson"
output:
  html_document:
    df_print: paged
    keep_md: TRUE
  html_notebook: default
  keep_md: TRUE
  pdf_document: default
---

# Preamble

Load R packages, functions and data.




Summary of the data:

```r
writeLines(sprintf('Number of individuals with at least one episode typed: %s',
                   length(unique(MS_pooled$ID))))
```

```
## Number of individuals with at least one episode typed: 164
```

```r
writeLines(sprintf('Number of episodes typed: %s',
                   length(unique(MS_pooled$Episode_Identifier))))
```

```
## Number of episodes typed: 599
```

```r
writeLines(sprintf('Number of recurrences typed: %s',
                   length(unique(MS_pooled$Episode_Identifier[MS_pooled$Episode>1]))))
```

```
## Number of recurrences typed: 435
```

```r
# We also remove MS data for which there are no recurrent data
N_episodes_typed = table(MS_pooled$ID[!duplicated(MS_pooled$Episode_Identifier)])
MS_pooled = filter(MS_pooled, ID %in% names(N_episodes_typed[N_episodes_typed>1]))

writeLines(sprintf('Number of individuals with at least two episodes typed: %s',
                   length(unique(MS_pooled$ID))))
```

```
## Number of individuals with at least two episodes typed: 159
```

```r
writeLines(sprintf('Number of recurrences typed: %s',
                   length(unique(MS_pooled$Episode_Identifier[MS_pooled$Episode>1]))))
```

```
## Number of recurrences typed: 435
```


Define the sets of microsatellite markers for the various datasets.


```r
MSs_VHX = c("PV.3.502","PV.3.27","PV.ms8","PV.1.501","PV.ms1","PV.ms5","PV.ms6")
MSs_all = c("PV.3.502","PV.3.27","PV.ms8","PV.1.501","PV.ms1","PV.ms5","PV.ms6",
            "PV.ms7","PV.ms16")
MSs_BPD = MSs_all
MSs_Main = c('PV.3.27', 'PV.3.502', 'PV.ms8') # These are typed for all episodes (the core group)
```



The approach is fully Bayesian and consists of the following:

* A prior probability vector for the recurrence state
* A likelihood based on the genetic data of being a *relapse*, a *recrudescence*, or a *reinfection* given the observed microsatellite data.

# Allele frequencies

There are a few ways of computing these. A natural first approach is to use the monoclonal data. However, some alleles are only seen in polyclonal infections, rending this approach not viable. A statistically rigorous approach would be to use a model for MS allele frequencies (e.g. Escalante 2015). At the moment we are using the empirical allele frequencies from the a specified dataset, with a Dirichlet-esque weight of 5 (5 pseudo-observations). 
Setting the weight to 0 recovers unweighted empirical allele frequencies. 


```
## Number of episodes used to compute frequencies: 159
```


## Plotting allele frequencies

These are the observed allele frequencies in the pooled data. We show 80% credible intervals (lo)

![](Pooled_Analysis_files/figure-html/AlleleFrequencies-1.png)<!-- -->


# Computing the probability of relatedness across infections

The following iterates through each individual and computes the probability of relatedness states.

## Load the time-to-event priors


```r
inds = grepl('mean_theta', colnames(Mod3_ThetaEstimates)) # Extract mean
Episode_Identifier = Mod3_ThetaEstimates$Episode_Identifier
p = data.frame(Episode_Identifier = Episode_Identifier, Mod3_ThetaEstimates[,inds],
               stringsAsFactors = F) # Reformat
colnames(p) = c('Episode_Identifier', 'C', 'L', 'I')

genetic_AND_time_data_eps = intersect(p$Episode_Identifier, MS_pooled$Episode_Identifier)
p = p[p$Episode_Identifier %in% genetic_AND_time_data_eps,]
Post_samples_matrix = Post_samples_matrix[Post_samples_matrix$Episode_Identifier %in% genetic_AND_time_data_eps,]
```


## Computation using full dataset 

We use all 9MS markers (when available).


### Full posterior computation




# Plot results

These dataframes are sorted by episode number so the columns correspond between them. We make some data.frames that store the results for ease of plotting.


```r
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
```

```
## The following `from` values were not present in `x`: AS
```

```r
thetas_9MS_Tagnostic$drug_col = mapvalues(x = thetas_9MS_Tagnostic$drug, 
                                          c('AS','CHQ','CHQ/PMQ'), mycols)
```

```
## The following `from` values were not present in `x`: AS
```


```r
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
```

## Going from time-to-event prior to posterior

Have broken it down by radical cure and no radical cure, as that is the most informative distinction here.


```r
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
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Probability of relapse, ordered from most to least likely:

```r
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
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

Probability of reinfection, ordered from most to least likely:

```r
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
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

Probability of recrudescence, ordered from most to least likely:

```r
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
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

# BPD Final Plot


```r
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
```

```
## The mean percentage of recurrences which are estimated to be relapses is 15%
```

![](Pooled_Analysis_files/figure-html/BPD_efficacy-1.png)<!-- -->


# Extra computations for VHX: too complex episodes


We remove the IDs that can be straightforwardly calculated:

```r
ind_calculated = which(MS_pooled$Episode_Identifier %in% thetas_9MS$Episode_Identifier)
IDs_calculated = unique(MS_pooled$ID[ind_calculated])
IDs_remaining = unique(MS_pooled$ID[! MS_pooled$ID %in% IDs_calculated])
```

We blow up the pooled analysis into all pairs within individuals:



Construct adjacency graphs and compute probabilities of relapse and reinfection.

```r
MS_pooled$L_or_C_state = MS_pooled$TotalEpisodes = NA
MS_pooled$L_lower = MS_pooled$L_upper = MS_pooled$L_mean = NA
MS_pooled$C_lower = MS_pooled$C_upper = MS_pooled$C_mean = NA
MS_pooled$I_lower = MS_pooled$I_upper = MS_pooled$I_mean = NA
# Arrange by complexity
# Get single rows per episode (throw away the extra MOI information)
MS_inflated = MS_inflated[!duplicated(MS_inflated$Episode_Identifier) & MS_inflated$Episode>1,]
Res$ID_True = MS_inflated$ID_True
Res$First_EpNumber = MS_inflated$First_EpNumber
Res$Second_EpNumber = MS_inflated$Second_EpNumber


## Threshold value for classification
Epsilon = 0.5

# Iterate through the ones we can calculate in one go
episodes_full_model = unique(Thetas_full_post$Episode_Identifier)
for(ep in episodes_full_model){
  ind1 = (MS_pooled$Episode_Identifier==ep)
  ind2 = (Thetas_full_post$Episode_Identifier==ep)
  MS_pooled$L_upper[ind1] = quantile(unlist(Thetas_full_post[ind2, grep('L',colnames(Thetas_full_post))]),
                                     probs=0.975)
  MS_pooled$L_lower[ind1] = quantile(unlist(Thetas_full_post[ind2, grep('L',colnames(Thetas_full_post))]),
                                     probs=0.025)
  MS_pooled$L_mean[ind1] = quantile(unlist(Thetas_full_post[ind2, grep('L',colnames(Thetas_full_post))]),
                                    probs=0.5)
  if(MS_pooled$L_upper[ind1] < Epsilon){
    MS_pooled$L_or_C_state[ind1] = 'I'
  } else if(MS_pooled$L_lower[ind1] > Epsilon){
    MS_pooled$L_or_C_state[ind1] = 'L'
  } else if(MS_pooled$L_upper[ind1] > Epsilon & MS_pooled$L_lower[ind1] < Epsilon){
    MS_pooled$L_or_C_state[ind1] = 'Uncertain'
  }
}

for(i in 1:length(IDs_remaining)){
  id = IDs_remaining[i]
  Doubles_Thetas = filter(Res, ID_True==id)
  
  for(ep in unique(Doubles_Thetas$Second_EpNumber)){
    ind1 = which(MS_pooled$ID==id & MS_pooled$Episode==ep)
    ind2 = which(Doubles_Thetas$Second_EpNumber == ep)
    
    MS_pooled$L_lower[ind1] = mean(Doubles_Thetas$L_min[ind2])
    MS_pooled$L_upper[ind1] = mean(Doubles_Thetas$L_max[ind2])
    MS_pooled$L_mean[ind1] = mean(Doubles_Thetas$L_mean[ind2])
    
    MS_pooled$C_lower[ind1] = mean(Doubles_Thetas$C_min[ind2])
    MS_pooled$C_upper[ind1] = mean(Doubles_Thetas$C_max[ind2])
    MS_pooled$C_mean[ind1] = mean(Doubles_Thetas$C_mean[ind2])

    MS_pooled$I_lower[ind1] = mean(Doubles_Thetas$I_min[ind2])
    MS_pooled$I_upper[ind1] = mean(Doubles_Thetas$I_max[ind2])
    MS_pooled$I_mean[ind1] = mean(Doubles_Thetas$I_mean[ind2])

    if(!is.na(MS_pooled$L_upper[ind1])){
      if(MS_pooled$L_upper[ind1] < MS_pooled$L_lower[ind1]){
        print(id)
      }
      if(MS_pooled$L_upper[ind1] < Epsilon){
        MS_pooled$L_or_C_state[ind1] = 'I'
      }
      if(MS_pooled$L_lower[ind1] > Epsilon){
        MS_pooled$L_or_C_state[ind1] = 'L'
      }
      if(MS_pooled$L_upper[ind1] > Epsilon & MS_pooled$L_lower[ind1] < Epsilon){
        MS_pooled$L_or_C_state[ind1] = 'Uncertain'
      }
    }
  }
  
}

MS_pooled$Drug = MS_pooled$FU = NA
for(id in MS_pooled$ID){
  ind = MS_pooled$ID==id
  MS_pooled$TotalEpisodes[ind] = max(MS_pooled$Episode[ind])
  MS_pooled$Drug[ind] = as.numeric(
    Combined_Time_Data$arm_num[Combined_Time_Data$patientid==id][1] == 'CHQ/PMQ') + 2
  MS_pooled$FU[ind] = Combined_Time_Data$FU_time[Combined_Time_Data$patientid==id][1]
}

MS_pooled$Plotting_pch_Values = as.numeric(mapvalues(MS_pooled$L_or_C_state, 
                                                     from = c('L','Uncertain','I'), to = 15:17))
MS_pooled$Plotting_col_Values = as.numeric(mapvalues(MS_pooled$L_or_C_state, 
                                                     from = c('L','Uncertain','I'), to = 1:3))
```


```r
## Time series data colored by genetic STATE: classification
mycols_states = c('red','yellow','black') # colors for states - need uncertain ones as well
mycols_drugs = brewer.pal(n=3, name = 'Set1')

MS_final = filter(MS_pooled, !is.na(L_mean))
MS_final = arrange(MS_final, Drug, desc(FU), desc(TotalEpisodes))

ids = unique(MS_final$ID)

par(las=1, bty='n', cex.axis=.3, mar=c(3,0,1,1))
plot(NA, NA, xlim = c(0,370), ylim = c(1,length(ids)),
     xaxt='n', yaxt='n')
mtext(text = 'Days from start of study', side = 1, line=2, cex=1.3)
axis(1, at = seq(0,370, by=60), cex.axis=1.5)

for(i in 1:length(ids)){
  
  id = ids[i]
  ind = which(MS_final$ID==id)
  
  # Add the follow up time line
  lines(c(0,MS_final$FU[ind[1]]), 
        c(i,i), lty=1, 
        lwd=.5, col= mycols_drugs[MS_final$Drug[ind[1]]])
  
  cols = mycols_states[MS_final$Plotting_col_Values[ind]]
  points(MS_final$timeSinceEnrolment[ind], rep(i,length(ind)), 
         pch=MS_final$Plotting_pch_Values[ind], 
         col=cols,cex=.6)
  
}
lines(x = c(0,0), y = c(0,length(ids)),lwd=3)
```

![](Pooled_Analysis_files/figure-html/CoatneyStylePLot-1.png)<!-- -->


```r
par(las=1, mfcol = c(2,1), mar=c(3,4,1,1), bty='n')
MS_final = arrange(MS_final, L_mean)
plot(1:nrow(MS_final), MS_final$L_mean,
     col = mycols_drugs[MS_final$Drug], 
     pch=20, cex=.3, yaxt='n', ylab='')
axis(2, at = c(0,.5,1))
polygon(c(1:nrow(MS_final), rev(1:nrow(MS_final))), 
        y = c(MS_final$L_lower,
              rev(MS_final$L_upper)), 
        border = NA, col = rgb(1, 0, 0,0.5))
points(1:nrow(MS_final), MS_final$L_mean, 
       col = mycols_drugs[MS_final$Drug], 
       pch=20, cex=.5)
mtext(side = 1, text = 'Recurrence index', line = 2)
mtext(side = 2, text = 'Relapse probability', line = 3,las=3)

plot(MS_final$timeSinceLastEpisode, MS_final$L_mean,
     col = mycols_drugs[MS_final$Drug], 
     pch=20, cex=1, xlab = '',yaxt='n',
     ylab = '')
axis(2, at = c(0,.5,1))
for(i in 1:nrow(MS_final)){
  
  if(abs(MS_final$L_upper[i] - MS_final$L_lower[i]) > 0.005){
    arrows(x0 = MS_final$timeSinceLastEpisode[i],
           y0 = MS_final$L_lower[i],
           x1=MS_final$timeSinceLastEpisode[i],
           y1 = MS_final$L_upper[i],
           length = 0.02,angle = 90, 
           code = 3,
           col=mycols_drugs[MS_final$Drug[i]])
  }
  i=i+1
}
mtext(side = 1, text = 'Days since last episode', line = 2)
mtext(side = 2, text = 'Relapse probability', line = 3,las=3)
```

![](Pooled_Analysis_files/figure-html/CompleteDataPlot-1.png)<!-- -->


The summaries of the final dataset:

```r
table(MS_final$Drug[!duplicated(MS_final$ID)])
```

```
## 
##  2  3 
## 80 79
```

```r
ind_CQ = MS_final$Drug==2
writeLines(sprintf('The weighted average of relapses is %s (%s-%s)',
                   round(100*sum(MS_final$L_mean[ind_CQ])/sum(ind_CQ),1),
                   round(100*sum(MS_final$L_lower[ind_CQ])/sum(ind_CQ),1),
                   round(100*sum(MS_final$L_upper[ind_CQ])/sum(ind_CQ),1)))
```

```
## The weighted average of relapses is 98.3 (92.7-100)
```

```r
writeLines(sprintf('The weighted average of recrudescences is %s (%s-%s)',
                   round(100*sum(MS_final$C_mean[ind_CQ],na.rm=T)/sum(ind_CQ),3),
                   round(100*sum(MS_final$C_lower[ind_CQ],na.rm=T)/sum(ind_CQ),3),
                   round(100*sum(MS_final$C_upper[ind_CQ],na.rm=T)/sum(ind_CQ),3)))
```

```
## The weighted average of recrudescences is 0 (0-0)
```

```r
writeLines(sprintf('The weighted average of reinfections is %s (%s-%s)',
                   round(100*sum(MS_final$I_mean[ind_CQ],na.rm=T)/sum(ind_CQ),3),
                   round(100*sum(MS_final$I_lower[ind_CQ],na.rm=T)/sum(ind_CQ),3),
                   round(100*sum(MS_final$I_upper[ind_CQ],na.rm=T)/sum(ind_CQ),3)))
```

```
## The weighted average of reinfections is 0 (0-0)
```

