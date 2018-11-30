---
title: "Pooled Analysis"
author: "Aimee Taylor and James Watson"
output:
  html_document:
    df_print: paged
    keep_md: TRUE
  pdf_document: default
  keep_md: TRUE
  html_notebook: default
---

# Preamble

Load R packages, functions and data.




Summary of the data and the whole of the VHX data set versus the subset typed (in terms of number of episodes):


```
## Number of individuals with at least one episode typed: 217
```

```
## Number of episodes typed: 710
```

```
## Number of recurrences typed: 494
```

```
## 
## Overall in the dataset: breakdown by treatment group:
```

```
## 
##  AS CHQ PMQ 
##  13  90 114
```

```
## 
## Within VHX: breakdown by treatment group:
```

```
## 
##  AS CHQ PMQ 
##  13  90  34
```

```
## 
## From BPD trial there are 80 individuals with total of 167 episodes typed
```

```
## From VHX trial there are 137 individuals with total of 543 episodes typed
```



```
## Of 0 of 137 VHX individual/s selected for genotyping, 1 to 8 of their episodes were not typed
```

```
## Of 0 of 80 BPD individual/s selected for genotyping, 1 to 2 of their episodes were not typed
```

Bug in the next bit so commented out


Summary of complexity of infection based on numbers of alleles observed 

![](Pooled_Analysis_files/figure-html/COIs_VHX_BPD-1.png)<!-- -->

```
## Median COI in VHX and BPD: 1 and 1, respectively
```

```
## Number of episodes with COI >= 3: 30 of 710 (4.23 percent)
```

Define the sets of microsatellite markers for the various datasets.




# Allele frequencies

We use a multinomial-dirichlet model with subjective weight $\omega$. $\omega = 0$ recovers unweighted empirical allele frequencies. 


```
## Number of episodes used to compute frequencies: 216
```


## Plotting allele frequencies

These are the mean posterior allele frequencies (dots) and 95\% credible intervals (bars) given pooled enrollment data and $\omega=$ `D_weight_Prior`.  

![](Pooled_Analysis_files/figure-html/AlleleFrequencies-1.png)<!-- -->


# Computing the probability of relatedness across infections

The approach is Bayesian and consists of the following:

* A prior probability vector for the recurrence state from the time-to-event model
* An allele frequency estimate from the posterior distribution of allele frequencies
* A likelihood based on the genetic data of being a *relapse*, a *recrudescence*, or a *reinfection* given the observed microsatellite data.

The following iterates through each individual and computes the probability of relatedness states.


```
## Warning: package 'bindrcpp' was built under R version 3.4.4
```

```
## Number of individuals with at least two episodes typed: 212
```

```
## Number of episodes in individuals with at least two episodes: 705
```

```
## Number of recurrences typed: 493
```


## Load the time-to-event priors




## Computation using full dataset 

We use all 9MS markers (when available).



### Full posterior computation




# Plot results






## Going from time-to-event prior to posterior

Plotted by radical cure versus no radical cure, as that is the most informative distinction here.

![](Pooled_Analysis_files/figure-html/Supplementary_TimeEffect_onPosterior-1.png)<!-- -->

Probability of states, ordered from most to least likely:

![](Pooled_Analysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


# BPD Final plot



```
## The mean percentage of recurrences which are estimated to be relapses is 16%
```

![](Pooled_Analysis_files/figure-html/BPD_efficacy-1.png)<!-- -->


# Extra computations for VHX: too complex episodes


We remove the IDs that can be straightforwardly calculated:



We blow up the pooled analysis into all pairs within individuals:







Construct adjacency graphs and compute probabilities of relapse and reinfection.


```
## [1] "removing"
## [1] "removing"
```

![](Pooled_Analysis_files/figure-html/CoatneyStylePLot-1.png)<!-- -->

![](Pooled_Analysis_files/figure-html/CompleteDataPlot-1.png)<!-- -->
Individuals who appear to relapse very late (more than 300 days after last episode):


```
## The episode ids of interest are: VHX_235_3
## The episode ids of interest are: BPD_27_2
```

```
##          ID Episode Episode_Identifier Treatment MOI_id
## 60   BPD_27       1           BPD_27_1       PMQ      1
## 61   BPD_27       2           BPD_27_2       PMQ      1
## 62   BPD_27       2           BPD_27_2       PMQ      2
## 355 VHX_235       1          VHX_235_1       CHQ      1
## 356 VHX_235       1          VHX_235_1       CHQ      2
## 357 VHX_235       2          VHX_235_2       CHQ      1
## 358 VHX_235       3          VHX_235_3       CHQ      1
##     timeSinceLastEpisode timeSinceEnrolment PV.1.501 PV.3.27 PV.3.502
## 60                     0                  0        3      33        7
## 61                   308                308        3      33        7
## 62                   308                308        3      35        7
## 355                    0                  0        1       5        2
## 356                    0                  0        1       5        2
## 357                   21                 21        1       5        3
## 358                  309                330        1       5        3
##     PV.ms1 PV.ms16 PV.ms5 PV.ms6 PV.ms7 PV.ms8
## 60       4      27     24     15      5     17
## 61       4      27     24     15      5     17
## 62       4      27     24     15      5     17
## 355      3      23     13      9     10     12
## 356      3      23     13     15     10     33
## 357      4      20     13      9     10     12
## 358      4      23     11     15     10     12
```


The summaries of the final dataset:

```
## 
##   2   3 
##  99 109
```

```
## In chloroquine monotherapy individuals, the weighted average of relapses is 99.3 (96.8-99.9)
```

```
## In chloroquine monotherapy individuals, the weighted average of recrudescences is 0.3 (0.1-0.6)
```

```
## In chloroquine monotherapy individuals, the weighted average of reinfections is 0.4 (0-2.6)
```

```
## In primaquine treated individuals, the weighted average of relapses is 14.3 (12.3-16.7)
```

```
## In primaquine treated individuals, the weighted average of recrudescences is 0 (0-0.3)
```

```
## In primaquine treated individuals, the weighted average of reinfections is 85.7 (83.3-87.5)
```

