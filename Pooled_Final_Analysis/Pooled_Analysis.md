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



Define the sets of microsatellite markers for the various datasets.


```r
MSs_VHX = c("PV.3.502","PV.3.27","PV.ms8","PV.1.501","PV.ms1","PV.ms5","PV.ms6")
MSs_all = c("PV.3.502","PV.3.27","PV.ms8","PV.1.501","PV.ms1","PV.ms5","PV.ms6",
            "PV.ms7","PV.ms16")
MSs_BPD = MSs_all
MSs_Main = c('PV.3.27', 'PV.3.502', 'PV.ms8') # These are typed for all episodes (the core group)
```


```r
#--------------------------------------------------------------------------
# Reformat the data s.t. there are no NA gaps in mixed infections
#--------------------------------------------------------------------------
# We also remove MS data for which there are no recurrent data
MS_data_reformated = reformat_MSdata(MSdata = MS_pooled, MSs=MSs_all)

N_episodes_typed = table(MS_pooled$ID[!duplicated(MS_pooled$Episode_Identifier)])
MS_data_reformated = filter(MS_data_reformated, ID %in% names(N_episodes_typed[N_episodes_typed>1]))
```



The approach is fully Bayesian and consists of the following:

* A prior probability vector for the recurrence state
* A likelihood based on the genetic data of being a *relapse*, a *recrudescence*, or a *reinfection* given the observed microsatellite data.

# Allele frequencies

There are a few ways of computing these. A natural first approach is to use the monoclonal data. However, some alleles are only seen in polyclonal infections, rending this approach not viable. A statistically rigorous approach would be to use a model for MS allele frequencies (e.g. Escalante 2015). At the moment we are using the empirical allele frequencies from the a specified dataset, with a Dirichlet-esque weight of 5 (5 pseudo-observations). 
Setting the weight to 0 recovers unweighted empirical allele frequencies. 


```
## Number of episodes used to compute frequencies: 164
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

genetic_AND_time_data_eps = intersect(p$Episode_Identifier, MS_data_reformated$Episode_Identifier)
p = p[p$Episode_Identifier %in% genetic_AND_time_data_eps,]
Post_samples_matrix = Post_samples_matrix[Post_samples_matrix$Episode_Identifier %in% genetic_AND_time_data_eps,]
```


## Computation using full dataset 

We use all 9MS markers (when available).


### Full posterior computation




# Plot results


```r
MS_summary = filter(MS_data_reformated, MOI_id == 1)
sum(MS_summary$Episode>1)
```

```
## [1] 435
```


These dataframes are sorted by episode number so the columns correspond between them. We make some data.frames that store the results for ease of plotting.


```r
thetas_9MS = arrange(thetas_9MS, Episode_Identifier)
thetas_9MS_Tagnostic = arrange(thetas_9MS_Tagnostic, Episode_Identifier)

Time_Estimates_1 = filter(Mod3_ThetaEstimates, 
                          Episode_Identifier %in% thetas_9MS$Episode_Identifier)
Time_Estimates_1 = arrange(Time_Estimates_1, Episode_Identifier)

thetas_9MS$drug = as.factor(Time_Estimates_1$arm_num)
# for plotting
thetas_9MS$drug_col = as.numeric(Time_Estimates_1$arm_num=='CHQ/PMQ')+1
thetas_9MS_Tagnostic$drug_col = as.numeric(Time_Estimates_1$arm_num=='CHQ/PMQ')+1

thetas_9MS_Tagnostic$drug = as.factor(Time_Estimates_1$arm_num)
```


## Going from time-to-event prior to posterior

There is some interesting correlation structure here - not quite sure what's happening exactly.
Have broken it down by radical cure and no radical cure, as that is quite a big piece of information!


```r
if(CREATE_PLOTS){
  par(mfrow=c(1,2),las=1, bty='n')
  # Time agnostic versus full posterior 
  plot(log10(thetas_9MS_Tagnostic$L), log10(thetas_9MS$L), 
       col=thetas_9MS$drug_col, main = 'Relapse',
       xlab = 'Time agnostic', ylab = 'Time included')
  lines(-10:0,-10:0)
  plot(log10(thetas_9MS_Tagnostic$I), log10(thetas_9MS$I), 
       col=thetas_9MS$drug_col, main = 'Reinfection',
       xlab = 'Time agnostic', ylab = 'Time included')
  lines(-20:0,-20:0)
  
  ##### Prior versus full posterior
  plot(log10(Time_Estimates_1$Relapse_mean_theta),
       log10(thetas_9MS$L),main = 'Relapse',
       col=thetas_9MS$drug_col,
       xlab = 'Time based prior', ylab = 'Full posterior')
  lines(-10:10,-10:10)
  plot(log10(Time_Estimates_1$ReInfection_mean_theta),
       log10(thetas_9MS$I),main = 'Reinfection',
       col=thetas_9MS$drug_col,
       xlab = 'Time based prior', ylab = 'Full posterior')
  lines(-10:10,-10:10)
}
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->![](Pooled_Analysis_files/figure-html/unnamed-chunk-10-2.png)<!-- -->

Probability of relapse, ordered from most to least likely:

```r
if(CREATE_PLOTS){
  par(las=1, bty='n')
  reLapse_ordered = sort.int(thetas_9MS$L, decreasing = TRUE, index.return = TRUE)
  plot(reLapse_ordered$x, pch=18, col = thetas_9MS$drug_col[reLapse_ordered$ix],
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
  
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
  
  reLapse_ordered_Tagn = sort.int(thetas_9MS_Tagnostic$L, decreasing = TRUE, index.return = TRUE)
  plot(reLapse_ordered_Tagn$x, pch=18, cex=.8,
       col = thetas_9MS_Tagnostic$drug_col[reLapse_ordered_Tagn$ix],
       xlab = 'Recurrence index', ylab = 'Probability of relapse state',
       main = 'Time agnostic posterior: reLapse')
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
}
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->![](Pooled_Analysis_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

Probability of reinfection, ordered from most to least likely:

```r
if(CREATE_PLOTS){
  
  par(las=1, bty='n')
  reinfection_ordered = sort.int(thetas_9MS$I, decreasing = TRUE, index.return = TRUE)
  plot(reinfection_ordered$x, pch=18, col = thetas_9MS$drug_col[reinfection_ordered$ix],
       xlab = 'Recurrence index', ylab = 'Probability of reinfection state',
       main = 'Full posterior: reInfection')
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
  
  reinfection_ordered_Tagn = sort.int(thetas_9MS_Tagnostic$I, decreasing = TRUE, index.return = TRUE)
  plot(reinfection_ordered_Tagn$x, pch=18, cex=.8,
       col = thetas_9MS_Tagnostic$drug_col[reinfection_ordered_Tagn$ix],
       xlab = 'Recurrence index', ylab = 'Probability of reinfection state',
       main = 'Time agnostic posterior: reInfection')
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
}
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->![](Pooled_Analysis_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

Probability of recrudescence, ordered from most to least likely:

```r
if(CREATE_PLOTS){
  par(las=1, bty='n')
  recrud_ordered = sort.int(thetas_9MS$C, decreasing = TRUE, index.return = TRUE)
  plot(recrud_ordered$x, pch=18, col = thetas_9MS$drug_col[recrud_ordered$ix],
       xlab = 'Recurrence index', ylab = 'Probability of recrudescence state',
       main = 'Full posterior: reCrudescence')
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
  
  recrud_ordered_Tagn = sort.int(thetas_9MS_Tagnostic$C, decreasing = TRUE, index.return = TRUE)
  plot(recrud_ordered_Tagn$x, pch=18, cex=.8,
       col = thetas_9MS_Tagnostic$drug_col[recrud_ordered_Tagn$ix],
       xlab = 'Recurrence index', ylab = 'Probability of recrudescence state',
       main = 'Time agnostic posterior: reCrudescence')
  legend('topright',col = 1:2, legend = c('No radical cure','Radical cure'),pch=18)
}
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->![](Pooled_Analysis_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

# Extra computations for VHX: too complex episodes

First we blow up the pooled analysis into all doubles

Do we need to add a check in the function for NA values?




Construct adjacency graphs and compute probabilities of relapse and reinfection.

```r
MS_pooled = MS_pooled[!duplicated(MS_pooled$Episode_Identifier),]
MS_pooled$L_or_C_state_UP = MS_pooled$L_or_C_state_LP = MS_pooled$ClusterID = MS_pooled$TotalEpisodes = NA

# Arrange by complexity
MS_inflated = arrange(MS_inflated, ID, Second_EpNumber)
# Get single rows per episode (throw away the extra MOI information)
MS_inflated = MS_inflated[!duplicated(MS_inflated$Episode_Identifier) & MS_inflated$Episode>1,]
Res$ID_True = MS_inflated$ID_True
Res$First_EpNumber = MS_inflated$First_EpNumber
Res$Second_EpNumber = MS_inflated$Second_EpNumber
MS_inflated$ID = as.factor(MS_inflated$ID)
MS_inflated$NumberClusters = NA

## Threshold value
Epsilons = c(0.1, 0.9)
ids = unique(MS_inflated$ID_True)
Graphs_lower = Graphs_upper =list()

for(i in 1:length(ids)){
  id = ids[i]
  Neps = max(MS_inflated$Second_EpNumber[MS_inflated$ID_True==id])
  Adj_Matrix = array(0, dim = c(Neps,Neps))
  diag(Adj_Matrix) = 1/2
  colnames(Adj_Matrix) = 1:Neps
  rownames(Adj_Matrix) = 1:Neps
  # Going to sum the reLapse probability and the reCrudescence probability
  Doubles_Thetas = filter(Res, ID_True==id)
  Is = Doubles_Thetas$First_EpNumber
  Js = Doubles_Thetas$Second_EpNumber
  for(k in 1:nrow(Doubles_Thetas)){
    Adj_Matrix[Is[k],Js[k]] = Doubles_Thetas$L[k]+Doubles_Thetas$C[k]
  }
  Adj_Matrix = Adj_Matrix + t(Adj_Matrix)
  
  # We're using two Epsilon values to see the uncertain ones
  Adj_Matrix_lower = Adj_Matrix_upper = Adj_Matrix
  Adj_Matrix_lower[ Adj_Matrix < Epsilons[1] ] = 0
  Adj_Matrix_lower[ Adj_Matrix >= Epsilons[1] ] = 1
  Graphs_lower[[i]] = graph_from_adjacency_matrix(Adj_Matrix_lower, 
                                                  mode = "undirected",diag = F)
  
  Adj_Matrix_upper[ Adj_Matrix < Epsilons[2] ] = 0
  Adj_Matrix_upper[ Adj_Matrix >= Epsilons[2] ] = 1
  Graphs_upper[[i]] = graph_from_adjacency_matrix(Adj_Matrix_upper, 
                                                  mode = "undirected",diag = F)
  
  
  # Add the number of clusters to the non-inflated MS dataset
  MS_pooled$TotalEpisodes[MS_pooled$ID==id] = Neps
  MS_pooled$NumberClusters[MS_pooled$ID==id] = components(Graphs_lower[[i]])$no
  # Add the membership of each cluster
  MS_pooled$ClusterID[MS_pooled$ID==id] = components(Graphs_lower[[i]])$membership
  MS_pooled$L_or_C_state_LP[MS_pooled$ID==id] = duplicated(components(Graphs_lower[[i]])$membership)
  MS_pooled$L_or_C_state_UP[MS_pooled$ID==id] = duplicated(components(Graphs_upper[[i]])$membership)
}
MS_pooled$L_or_C_state = apply(MS_pooled[, c('L_or_C_state_UP','L_or_C_state_LP')],1 , sum)
```


```r
## Time series data colored by genetic STATE
mycols = brewer.pal(5,'Set2')
par(las=1, bty='n', cex.axis=.3, mar=c(3,0,1,1))
plot(NA, NA, xlim = c(0,370), ylim = c(1,length(ids)),
     xaxt='n', yaxt='n')
#mtext(text = 'Patient ID number', side = 2, line=2,las=3, cex=1.3)
mtext(text = 'Days from start of study', side = 1, line=2, cex=1.3)
#mtext(text = paste('Threshold = ', Epsilon), side = 3, line=0)
axis(1, at = seq(0,370, by=60), cex.axis=1.5)
cc = 0
for(i in 1:length(ids)){
  
  id = ids[i]
  ind = which(MS_pooled$ID==id & MS_pooled$Episode>1)
  Neps = length(ind)
  
  drug_col = as.numeric(
    Combined_Time_Data$arm_num[Combined_Time_Data$patientid==id][1] == 'CHQ/PMQ') + 4 
  lines(c(0,Combined_Time_Data$FU_time[Combined_Time_Data$patientid==id][1]), 
        c(i,i), lty=1, 
        lwd=.5, col= mycols[drug_col])
  if(Neps > 0){
    cols = mycols[MS_pooled$L_or_C_state[ind]+1]
    points(MS_pooled$timeSinceEnrolment[ind], rep(i,Neps), pch=19, col=cols,cex=.6)
  }
}
```

![](Pooled_Analysis_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

