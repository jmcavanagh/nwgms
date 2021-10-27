#!/bin/bash
#SBATCH -n 8
#SBATCH --mem=16G
#SBATCH -t 1:00:00
source ~/est/bin/activate
python3 testspec.py $1 $2 $3 $4 $5
deactivate
