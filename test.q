#!/bin/bash
#
#SBATCH -J testjob
#SBATCH -A liu-compute-2020-20
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH --exclusive
#
module load Python/3.8.3-anaconda-2020.07-extras-nsc1
module load impi/.2018.1.163-eb
mpirun python3 main_mpi.py
echo "job completed"
