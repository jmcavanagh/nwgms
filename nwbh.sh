#!/bin/bash
#SBATCH -n 1
#SBATCH -t 1:00:00
python3 nwbh.py $1 $2 $3 $4 $5
