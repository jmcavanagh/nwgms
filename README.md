# nwgms
Global optimization package for use with NWChem

For use with Brown's CCV.

First, makee sure necessary packages are installed. Modify your ~/.bash_profile to include

module load nwchem/7.0.2_mvapich2-2.3.5_intel_2020.2_slurm20  
module load mpi/mvapich2-2.3.5_intel_2020.2_slurm20 intel/2020.2 cuda/11.1.1 gcc/10.2  
nwchem/6.8-openmpi  
module load python/3.7.4  


Then, install numpy with pip install numpy

Then, download the github repo. Go into the directory, for example scratch, where you want to perrform the global minimum search. Then enter the command git clone https://github.com/jmcavanagh/nwgms

Then, modify search_specs.txt. To start a global minimum search, enter the command sbatch coordinator.sh <numcalcs> <numhops> <multiplicity>. You can search several multiplicities at once.
