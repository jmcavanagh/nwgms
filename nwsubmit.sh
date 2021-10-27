#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -n 1
#SBATCH --mem=2G
nwchem $1.inp >& $1.out
