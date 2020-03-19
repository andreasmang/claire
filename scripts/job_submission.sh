#!/bin/bash

#### sbatch parameters
#SBATCH -J clairegpu
#SBATCH -o stdout.o
#SBATCH -e stderr.e
#SBATCH -p v100
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH --mail-user=naveen@ices.utexas.edu
#SBATCH --mail-type=all

nx=512x512x256
case=random
bin=/home/04716/naveen15/claire-dev/bingpu/test
for p in 1 2 4 8 16 32; do
    mpirun -np $p -gpu -pami_noib $bin -interp -nx ${nx} > /home/04716/naveen15/claire-dev/scripts/multi-node-random/interp_${case}_nx-${nx}_p-${p}.log
done

