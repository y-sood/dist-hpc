#!/bin/bash
#SBATCH -A uppmax2025-2-247
#SBATCH -N 2
#SBATCH -n 32
#SBATCH --ntasks-per-node=16
#SBATCH -p node
#SBATCH --time=00:30:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --job-name=scaling_tests

mpirun -np 32 ./main 3200000
mpirun -np 16 ./main 1600000
mpirun -np 8 ./main 800000
mpirun -np 4 ./main 400000
mpirun -np 2 ./main 200000
mpirun -np 1 ./main 100000

mpirun -np 32 ./main 500000
mpirun -np 16 ./main 500000
mpirun -np 8 ./main 500000
mpirun -np 4 ./main 500000
mpirun -np 2 ./main 500000
mpirun -np 1 ./main 500000
