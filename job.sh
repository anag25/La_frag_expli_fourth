#!/bin/bash
#SBATCH --job-name=canti_test1

## how much time is requested [ upper bound ]
#SBATCH --time=06:00:00

#specify memory needed (MB)
#SBATCH --mem=10000
 
 
##    notify me about this job
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=anag25@student.monash.edu
#SBATCH --output=testJob-%j.out
 
#SBATCH --open-mode=append
 
module load  matlab/r2016a
module list
ulimit -s 40000
  
MCR_CACHE_ROOT=$TMPDIR
export MCR_CACHE_ROOT
### on this next line is where you actually run your program
srun /run_explicit_main_monarch.sh  $MATLABROOT
   
### Here are some explanations:
### notice that when you compiled you used:
###     mcc -mv name.m
### then this produces as output the file:
###     run_name.sh
### to run this use the line above:
###     ./run_name.sh $MATLABROOT