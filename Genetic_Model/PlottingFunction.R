#=======================================================================
# This function take as input an MS data frame and plots accordingly The order
# of the markers in the input is preserved in the plot The numbers of alleles,
# thus ranks and colours used to denote the alleles are specific to the input.
# As such, when plotting different subsets of a master data set, numbers, ranks
# and colours will change.
#
# Idea for future - set allele legend area proportional to frequency
#
# Author: Aimee Taylor 
#=======================================================================

library(RColorBrewer) # for brewer.pal()
library(fields) # for image.plot()

ColorPlot_MSdata = function(MS_data, surpress_messages = FALSE){
  
  if(!surpress_messages){
    writeLines('Please be aware: \n This function assumes all microsatellite names are preceded by PV. \n The order of markers in the input is preservered in the plot. \n Numbers of alleles, ranks and thus colours are specific to the input.')
  }
  
  # cols is a function to create ramped colours based on "Paired" brewer.pal
  cols = colorRampPalette(brewer.pal(12, "Paired"))
  
  # Extract unique episodes
  episodes_uniq = unique(MS_data$Episode_Identifier)
  no_episodes = length(episodes_uniq)           
  
  # Extract marker names (assuming all ms names are preceded by PV.)
  markers = colnames(MS_data)[grepl('PV.', colnames(MS_data))]
  no_markers = length(markers)
  Marker_data = MS_data[,markers,drop = F] # Extract genetic data
  max_ranks = apply(Marker_data, 2, function(x){length(unique(x[!is.na(x)]))}) # Extract input-specific ranks
  maxDiv = max(max_ranks) # For colour mapping
  maxMOI =  max(MS_data$MOI_id) # Maximum across entire input data 
  factorial_maxMOI = factorial(maxMOI) # For expanding to wide format 
  
  # Expand markers to represent all polyclonal infections in wide format
  markers_wide = paste(rep(markers, each = factorial_maxMOI), rep(1:factorial_maxMOI,no_markers), sep = '_')  
  no_markers_wide = length(markers_wide) # equal to no_markers * factorial_maxMOI
  
  # Reformat data s.t. columns house alternative alleles in polyclonal samples 
  Marker_data_wide_format = array(dim = c(no_episodes, no_markers_wide), 
                                  dimnames = list(episodes_uniq, markers_wide))
  for(epi in episodes_uniq){
    ind = MS_data$Episode_Identifier == epi
    Subset = MS_data[ind,] # Extract episode specific data
    for(j in 1:no_markers){
      marker_alleles = unique(Subset[,markers[j]])
      if(all(is.na(marker_alleles))){next} # If untyped skip to next ms 
      MOI_marker = length(marker_alleles[!is.na(marker_alleles)])
      To_fill_list = lapply(as.list(seq(factorial_maxMOI, 1, -factorial_maxMOI/MOI_marker)), FUN = function(x){1:x})
      for(k in 1:MOI_marker){
        To_fill = paste(rep(markers[j], each = max(To_fill_list[[k]])), To_fill_list[[k]], sep = '_')
        Marker_data_wide_format[epi, To_fill] = rep(Subset[k,markers[j]], each = max(To_fill_list[[k]]))
      }
    }
  }
  
  # Re-colour to maximise contrast per marker
  Marker_data_wide_factor = array(dim = dim(Marker_data_wide_format), 
                                  dimnames = dimnames(Marker_data_wide_format))
  for(marker in markers){
    ind = grepl(paste(marker,'_',sep = ''), colnames(Marker_data_wide_format))
    Marker_data_wide_factor[,ind] = apply(Marker_data_wide_format[, ind, drop = FALSE], 2, function(x){
      y = as.numeric(factor(x, ordered = TRUE, levels = sort(unique(Marker_data[,marker]))))
      return(y)
    })
  }
  
  # Create an index to visually group individuals
  IDs_factor <- as.factor(MS_data$ID[!duplicated(MS_data$Episode_Identifier)])
  IDs = as.numeric(IDs_factor)
  #as.numeric(do.call(rbind, strsplit(episodes_uniq, split = '_'))[,2,drop = F])
  ID_01 = c(0, rep(NA, no_episodes-1))
  for(i in 2:no_episodes){
    if(IDs[i] == IDs[i-1]){
      ID_01[i] = ID_01[i-1]
    } else { # switch
      ID_01[i] = abs(ID_01[i-1] - 1)
    }
  }
 
  # Y axis
  ID_midpoints = rep(NA, length(unique(IDs)))
  names(ID_midpoints) = unique(IDs)
  for(ID in unique(IDs)){
    inds = which(IDs == ID)
    inds_mapped_01 = (inds - 0)/(no_episodes+1 - 0)
    ID_midpoints[as.character(ID)] = median(inds_mapped_01) 
  }
  
  # Plot
  par(mfrow = c(1,1), 
      fig = c(0,0.75,0.03,0.97), # fig = c(x1, x2, y1, y2) (to leave room for the legend)
      mar = c(4,4,0,0.5)) # mar = c(bottom, left, top, right)
  
  # Add IDs to each vertical edge
  Empty_array = array(dim = dim(Marker_data_wide_factor))
  matrix_to_plot = t(rbind(-1, cbind(ID_01, Empty_array, ID_01), -1)) 
  image(matrix_to_plot, ylim=c(0,1), col = grey(c(0.25,0.75)), axes = FALSE)
  
  # Add marker data 
  COLs = seq(1, no_markers_wide, factorial_maxMOI)
  for(i in 1:length(COLs)){
    Y = Marker_data_wide_factor
    Y[,-(COLs[i]:(COLs[i]+(factorial_maxMOI)-1))] = NA
    image(t(rbind(-1, cbind(NA, Y, NA), -1)), col = c(grey(0), cols(max_ranks[i])), 
          axes = FALSE, add = TRUE)
  }
  
  # Add axes
  xaxis = seq(0, 1, length.out = no_markers_wide+2)
  xaxis_diff = (xaxis[2] - xaxis[1])
  if(maxMOI == 1){
    xaxis_at = xaxis[2:(no_markers+1)]
  } else {
    xaxis_at = xaxis[seq((1 + factorial_maxMOI/2), no_markers_wide, factorial_maxMOI)] + ifelse((factorial_maxMOI %% 2) == 0, xaxis_diff/2, 1)
  }
  axis(at = xaxis_at, side = 1, line = -0.5, labels = markers, cex.axis = 0.5, tick = F)
  axis(at = xaxis_at, side = 1, line = 0.5, labels = paste('(', max_ranks, ')', sep = ''), cex.axis = 0.5, tick = F)
  axis(side = 2, at = ID_midpoints, labels = levels(IDs_factor)[unique(IDs)], cex.axis = 0.5, las = 2, lwd.ticks = 0.25, lwd = 0)
  axis(side = 1, at = c(0, tail(xaxis, 1)), line = -0.5, labels = rep('Grouping', 2), las = 2, cex.axis = 0.5, tick = F)
  title(ylab = bquote(.(no_episodes)~'episodes (one row per episode, grouped by patient ID)'), 
        line = 3.3, cex.lab = 0.5)
  title(xlab = 'Microsatellite (number of distinct alleles)', line = 3, cex.lab = 0.5)
  
  # Add legend
  par(fig = c(0,1,0,1)) # Critical to avoid odd placement of first legend column
  legendwidth = (1-0.76)/no_markers
  for(marker in markers){
    i = which(marker == markers)
    if(i == 1){
      SMALLplot = c(0.75,0.75+legendwidth,0.03,0.97)}else{
        SMALLplot = c(0.75+legendwidth*(i-1),0.75+legendwidth*i,0.03,0.97)}
    Y = Marker_data_wide_factor
    Breaks = if(max_ranks[marker] == 1){c(0,2)}else(NULL)
    image.plot(Y[,grepl(paste(marker,'_',sep = ''), colnames(Y)),drop = FALSE], 
               col = cols(max_ranks[marker]),
               breaks = Breaks, 
               legend.only = TRUE, add = TRUE, 
               axis.args = list(at = 1:max_ranks[marker], 
                                labels = sort(unique(Marker_data[,marker])), 
                                cex.axis = 0.5, tick = FALSE,
                                line = -1.5, hadj = 0.5), 
               smallplot = SMALLplot, legend.mar = 0) 
  }
}

# # Example 
# load('../RData/GeneticModel/MS_data_PooledAnalysis.RData')
# ColorPlot_MSdata(MS_pooled)


