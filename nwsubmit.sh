#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 4
srun --mpi=pmix nwchem $1.inp >& $1.out
if [ -z "$2" ]
then
        rm $1.hess
        rm $1.b*
        rm $1.c
        rm $1.d*
        rm $1.zmat
        rm $1.m*
fi
