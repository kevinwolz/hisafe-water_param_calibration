#!/bin/sh
#SBATCH --account=hisafe
#SBATCH --partition=defq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wolzkevin@gmail.com
module purge
module load load cv-standard
module load load R/3.4.3 
cd /lustre/lecomtei/water_calib_GA1
Rscript water_calib_GA.R
