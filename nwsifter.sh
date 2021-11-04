#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=512M
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL
v="_*"
rm slurm*
rm $1$v 
python3 nwsifter.py $2
