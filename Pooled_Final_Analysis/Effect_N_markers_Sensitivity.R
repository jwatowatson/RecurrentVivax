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
pdf('../Pooled_Final_Analysis/Relapse_probability_versus_Nmarkers.pdf')
par(las=1, bty='n')
boxplot(log(Inflated_Results$L) ~ as.factor(Inflated_Results$N_markers), 
        ylab='log probability of relapse state', 
        xlab='Number of markers used to estimate probability of relapse')
dev.off()
mod_test_trend = lm(log(L) ~ N_markers, data = Inflated_Results)
summary(mod_test_trend)

boxplot(log(Inflated_Results$I) ~ as.factor(Inflated_Results$N_markers), 
        ylab='log probability of reinfection state', 
        xlab='Number of markers used to estimate probability of reinfection',
        ylim=c(-2,0))

sd_NM = array(dim=9)
for(Nm in 1:9){
  ind=Inflated_Results$N_markers==Nm
  sd_NM[Nm] = sd(Inflated_Results$I[ind])
}
plot(1:9, sd_NM,xlab='Number of markers',ylab='Standard deviation of estimates')
