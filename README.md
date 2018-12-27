# RecurrentVivax

This repository has all the necessary code and data to reproduce results from **Estimating the probable cause of recurrence in Plasmodium vivax malaria: relapse, reinfection or recrudescence?** Taylor *et al*, 2018.

The full paper (preprint for the moment) can be found at: 
https://www.biorxiv.org/content/early/2018/12/25/505594

* *Genetic_Model* has the model code implementing our identity by descent approach to the estimation of relapse, recrudescence and reinfection in \textit{P. vivax} recurrences. 

* *Timing_Model* has the model code implementing our time-to-event approach to the estimation of relapse, recrudescence and reinfection in \textit{P. vivax} recurrences. 

* *Pooled_Final_Analysis* implements these models on data from two large studies on the Thailand-Myanmar border (total $n=1299$ individuals). In particular it implements the following:
    1. Full Bayesian analysis of all recurrences, where available estimate using the genetic data is based on prior probbailities from the time-to-event model
    2. Estimation of false-positive rate using genetic data alone by comparing isolates from different individuals
    3. Estimation of failure rate after supervised high-dose primaquine in the epidemiological context of the Thai-Myanmar border
    4. Prediction of failure from carboxy-primaquine (inactive, slowly eliminated metabolite) trough concentrations

* *Simulation_Study* estimates the number of microsatellites needed to reliably calculate recurrence probabilities in paired infections under certain assumptions of complexity of infection. 

All data are stored in *RData* and a microsatellite data plotting tool is given in *Plotting_MS_Data*.

Enjoy reading through the model code and model implementation! 
If you see any bugs or have any questions, the authors can be contacted at jwatowatson@gmail.com and aimee.r.taylor.85@gmail.com.
  
  