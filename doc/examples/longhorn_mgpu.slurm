#!/bin/bash

#SBATCH -J claire
#SBATCH -p v100
#SBATCH --exclusive
#SBATCH --gpus-per-task=1
#SBATCH -t 00:05:00
#SBATCH -o stdout.o
#SBATCH -e stderr.e
#SBATCH -N 1 # number of nodes
#SBATCH -n 2 # number of GPUs per node

BINDIR=${HOME}/claire/bin
module load xl spectrum_mpi
source ${HOME}/claire/deps/env_source.sh

CLAIREOPT= # add your options here

export MY_SPECTRUM_OPTIONS="--gpu --aff on"
HOSTS=`srun hostname | sort -n |  sed -r "s/^(.*?)\.longhorn\.tacc\.utexas\.edu$/\1/" | tr "\n" "," | sed "s/,$//"`
echo $HOSTS
mpirun --gpu --aff on --host $HOSTS -np 2 $BINDIR/claire $CLAIREOPT
