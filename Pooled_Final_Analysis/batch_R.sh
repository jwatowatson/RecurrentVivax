#!/bin/bash
# James Pooled_Final_Analysis
#$ -N Vivax_Full_Model
#
#pe request for MPICH. set number of processors
#$ -pe mpich 1
#
#set to current working directory
#$ -cwd
#
#combine stderr and stdout into one file
#$ -j y
#
#Run job through bash shell
#$ -S /bin/bash
#
#Adjust MPICH procgroup for smooth shutdown
export MPICH_PROCESS_GROUP=no
#
# Print useful information
echo "Program start at: `/bin/date`."
echo "Got $NSLOTS slots."

#Run program. Use full path
/opt/openmpi/bin/mpirun -np $NSLOTS /usr/bin/R --slave -f /home/James/RecurrentVivax/Pooled_Final_Analysis/Pooled_Analysis.R

echo "Program finish at: `/bin/date`."
