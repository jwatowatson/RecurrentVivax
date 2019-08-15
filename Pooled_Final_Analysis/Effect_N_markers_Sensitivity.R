rm(list=ls())
RUN_ANALYSIS = F
# We look at how the number of markers effects the relapse probability under null data
load('../RData/LargeFiles/Inflated_Results.bigRData')
Inflated_Results = Inflated_Results[!is.na(Inflated_Results$L),]
load('../RData/LargeFiles/APC_MSdata.bigRData')

if(RUN_ANALYSIS){
  null_ids = unique(rownames(Inflated_Results))
  Inflated_Results$N_markers = NA
  marker_cols = grep('PV', colnames(APC_MSdata))
  pb=txtProgressBar(min = 0, max = length(null_ids),style = 3)
  
  # Takes about 1 hour to run....Could be optimised
  for(ID in null_ids){
    ind = which(ID == APC_MSdata$Episode_Identifier)
    if(length(ind)>0){
      ID_APC = unique(APC_MSdata$ID[ind])
      ep_1_ind = which(APC_MSdata$ID==ID_APC & APC_MSdata$Episode==1)
      ep_2_ind = which(APC_MSdata$ID==ID_APC & APC_MSdata$Episode==2)
      ep_1_Not_NA = apply(APC_MSdata[ep_1_ind, marker_cols,drop=F],2,
                          function(x) sum(!is.na(x)) > 0)
      ep_2_Not_NA = apply(APC_MSdata[ep_2_ind, marker_cols,drop=F],2,
                          function(x) sum(!is.na(x)) > 0)
      N_markers = sum(ep_2_Not_NA & ep_1_Not_NA)
      Inflated_Results$N_markers[rownames(Inflated_Results)==ID] = N_markers
    }
    setTxtProgressBar(pb, value = which(ID==null_ids))
  }
  save(Inflated_Results, file = '../RData/LargeFiles/Inflated_Results_sensitivity_Nmarker.bigRData')
} else {
  load('../RData/LargeFiles/Inflated_Results_sensitivity_Nmarker.bigRData')
}

table(Inflated_Results$N_markers)

## Threshold value for classification
Epsilon_upper = 0.7
Epsilon_lower = 0.3
Dark2 = brewer.pal(8, 'Dark2')
transparent_pink_band = adjustcolor('blue', alpha.f = 0.2)

par(las=1, bty='n')
Inflated_Results = filter(Inflated_Results, N_markers>0)
boxplot((Inflated_Results$L) ~ as.factor(Inflated_Results$N_markers), 
        ylab='Probability of the relapse state', varwidth=T,
        xlab='Number of markers used to estimate the probability of relapse in the null data')
abline(h=(1/3), col='red',lwd=2, lty=2)
polygon(x = c(0,12,12,0),
        y = c(Epsilon_lower,
              Epsilon_lower,
              Epsilon_upper,
              Epsilon_upper),
        col = transparent_pink_band, border = NA)

table(Inflated_Results$N_markers)

writeLines(sprintf('When using 3 markers, %s%% of the comparisons for the null data will give posterior estimates greater than 1/3 (the prior value).\n When using 9 markers, %s%% of the comparisons for the null data will give posterior estimates greater than 1/3 (the prior value)',
                   round(100*sum(Inflated_Results$L>(1/3) & Inflated_Results$N_markers==3)/sum(Inflated_Results$N_markers==3)),
                   round(100*sum(Inflated_Results$L>(1/3) & Inflated_Results$N_markers==9)/sum(Inflated_Results$N_markers==9))))

writeLines(sprintf('When using 3 markers, %s%% of the comparisons for the null data will be above the upper threshold (0.7) for relapse classification.\n When using 9 markers, %s%% of the comparisons for the null data will be above the upper threshold (0.7) for relapse classification',
                   round(100*sum(Inflated_Results$L>(0.7) & Inflated_Results$N_markers==3)/sum(Inflated_Results$N_markers==3),1),
                   round(100*sum(Inflated_Results$L>(0.7) & Inflated_Results$N_markers==9)/sum(Inflated_Results$N_markers==9),1)))
