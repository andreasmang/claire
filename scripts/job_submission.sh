#!/bin/bash

#### sbatch parameters
#SBATCH -J clairegpu
#SBATCH -o stdout.o
#SBATCH -e stderr.e
#SBATCH -p v100
#SBATCH -N 8
#SBATCH -n 32
#SBATCH -t 01:00:00
#SBATCH --mail-user=naveen@ices.utexas.edu
#SBATCH --mail-type=all

bin=/home/04716/naveen15/claire-dev/bingpu/test
for p in 1 2 4 8 16 32; do
    mpirun -np $p -pami_noib -gpu $bin -interp -nx 256 > /home/04716/naveen15/claire-dev/scripts/interp_nx-256_p-${p}.log
done

