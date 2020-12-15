#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=0:24:00
#SBATCH --partition=shas-testing
#SBATCH --ntasks=4
#SBATCH --job-name=ecosims
#SBATCH --output=ecosims.out

module purge

module load matlab/R2019b

matlab -nodesktop -nodisplay -r "ecorunner2()"