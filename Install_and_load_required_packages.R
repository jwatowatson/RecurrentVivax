# This iterates through the list of necessary and packages and installs them if missing

ipak = function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

load('../RData/RPackages_List.RData')
ipak(pkgs)
sessionInfo() # See packages loaded
