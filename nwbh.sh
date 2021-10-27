#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=512M
#SBATCH -t 1:00:00
source ~/est/bin/activate
python3 nwbh.py $1 $2 $3 $4 $5
deactivate
