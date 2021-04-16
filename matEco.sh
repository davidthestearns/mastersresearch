#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition=shas
#SBATCH --ntasks=4
#SBATCH --job-name=ecosims
#SBATCH --output=ecosims.out

module purge

module load matlab

matlab -nodesktop -nodisplay -r "ecorunner2()"