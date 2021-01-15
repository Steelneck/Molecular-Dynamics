#!/bin/bash
#
#SBATCH -J testjob
#SBATCH -A LIU-compute-2020-20
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --exclusive
#
module load impi/.2018.1.163
module loard Python/3.6.4-nsc2-intel-2018a-eb

time mpirun python3 MPI_run.py

echo "job completed"
