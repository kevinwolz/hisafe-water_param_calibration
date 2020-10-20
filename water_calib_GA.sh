#!/bin/sh
#SBATCH --account=hisafe
#SBATCH --partition=defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wolzkevin@gmail.com
module purge
module load load cv-standard
module load load R/4.0.2 
cd /lustre/lecomtei/hisafe-water_param_calibration
Rscript water_calib_GA.R
