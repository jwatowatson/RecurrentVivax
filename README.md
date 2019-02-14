# RecurrentVivax

## Overview of software 

This repository has all the necessary code and data to reproduce results from **Estimating the probable cause of recurrence in Plasmodium vivax malaria: relapse, reinfection or recrudescence?** Taylor *et al*, 2018.

The full paper (preprint for the moment) can be found at: 
https://www.biorxiv.org/content/early/2018/12/25/505594

* *Genetic_Model* has the model code implementing our identity by descent approach to the estimation of relapse, recrudescence and reinfection in *P. vivax* recurrences, based on highly polymorphic microsatellite data. 

* *Timing_Model* has the model code (in *stan*) implementing our time-to-event approach to the estimation of relapse, recrudescence and reinfection in *P. vivax* recurrences, based on follow-up data post clinic level symptomatic *P. vivax* episodes. The RMarkdown script *TimingModelStan.Rmd* fits the models to pooled data from two large studies on the Thailand-Myanmar border (total *n=1299* individuals, over 1000 patient years of follow-up).

* *Pooled_Final_Analysis* provides an RMarkdown script *Pooled_Analysis.Rmd* which combines the genetic model fits and then timing model fits on the same data from the two large studies on the Thailand-Myanmar border (total *n=1299* individuals). In particular it implements the following:
    1. Full Bayesian analysis of all recurrences, and where available, estimates using the genetic data are based on prior probabilities from the time-to-event model
    2. Estimation of false-positive rate using genetic data alone by comparing isolates from different individuals
    3. Estimation of failure rate after supervised high-dose primaquine in the epidemiological context of the Thai-Myanmar border
    4. Prediction of primaquine failure from carboxy-primaquine (inactive, slowly eliminated metabolite) trough concentrations

* *Simulation_Study* estimates the number of microsatellites needed to reliably calculate recurrence probabilities in paired infections under certain assumptions of complexity of infection. This is in the RMarkdown script *SimulationStudy.Rmd*.

All data are stored in *RData* and a microsatellite data plotting tool is given in *Plotting_MS_Data*.

Enjoy reading through the model code and model implementation! 
If you see any bugs or have any questions, the authors can be contacted at jwatowatson@gmail.com and aimee.r.taylor.85@gmail.com.
  
## Installation and running guide

All software is written in R. For an installation guide please see https://cran.r-project.org/doc/manuals/R-admin.html
The R versions used to develop this code were: R version 3.4.3 (2017-11-30) and R version 3.5.0 (2018-04-23)

To run the timing model code, *rstan* is needed. For installation instructions please see: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
The rstan version we are currently using is: rstan Version 2.18.1, GitRev: 2e1f913d3ca3

Before running the analyses script please run *Install_and_load_required_packages.R*, which will load (and install if not installed already) all necessary packages. This should take a few minutes at most. Some of the non-base packages that will be loaded are as follows. A full list with version numbers can be found in session_info_txt 
* *igraph* (neat manipulation of graphs and computation of graph properties)
* *loo* (leave one out for stan models)
* *stringr* (character manipulation and pattern matching)
* *boot* (bootstrapping)
* *RColorBrewer* (get nice color palettes)
* *gdata* (some nice data manipulation functions)
* *gtools* (for calculations in log space)
* *tictoc* (timing of runs)
* *doParallel* (multicore computation) 


## Demo and fully reproducible output 

Please be aware that it takes a several days to regenerate the entire study from scratch. Some model runs are fairly computationally expensive (around 1-2 days using 6 cores on a desktop computer). We have thus made available intermediate RData, enabling results to be regenerated in a matter of minutes (see below). 

To regenerate results using intermediate RData (demo): download and put into directory RData/LargeFiles/ :

https://www.dropbox.com/sh/naslrxkyxqnvo0t/AADgAaLBta53Hc8_pX3AAzKha?dl=0

These are all the output results that are too large to be added to a github repository. The computation time to create these files by setting *RUN_MODELS___* variables to TRUE totals about 2-4 days on a standard desktop computer. 

All RMarkdown scripts can then be run in less than 10 minutes on a standard desktop computer by setting all *RUN_MODELS___* variables to FALSE (default), 
The output of the RMarkdown scripts is as follows:

* Plots of results and data as given in our paper
* Text detailing the data structures and summary statistics
* Model output when run with the *RUN_MODELS___* variables set to TRUE.

To regenerate everything from scratch, including all intermediate files, set *RUN_MODELS___* variables to TRUE in all RMarkdown scripts, then run TimingModelStan.Rmd followed by the other RMarkdown scripts in no specific order. The output is as described above. 





