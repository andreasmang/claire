#!/bin/bash

#### sbatch parameters
#SBATCH -J clairegpu
#SBATCH -o cuda_fwdsolve_benchmark.o
#SBATCH -e cuda_fwdsolve_benchmark.e
#SBATCH -p vis
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:02:00
#SBATCH --mail-user=naveen@ices.utexas.edu
#SBATCH --mail-type=all


source ${HOME}/claire/external/libs/environment_vars.sh 

#### define paths
BDIR=${HOME}/claire/bin
RDIR=${WORK}/claire_results
LDIR=benchmark_forward_solver
RDIR=$RDIR/$LDIR
WORK=${WORK}


#### submitt job
#ibrun tacc_affinity $BDIR/claire -x $RDIR/ -mr ${WORK}/nirep_updated/na01_128.nii.gz -mt ${WORK}/nirep_updated/na02_128.nii.gz -nx 128x128x128 -pdesolver sl -regnorm h2s -beta 1.000000e-02 -opttol 1.000000e-01 -nt 4 -verbosity 2 -x $RDIR/ -velocity -monitordefgrad -detdefgrad
#ibrun tacc_affinity /work/04716/naveen15/stampede2/claire/bin/clairetools -computeravensmap -labels 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33 -ifile ${WORK}/nirep_updated/na01_seg_128.nii.gz -v1 ${RDIR}/velocity-field-x1.nii.gz -v2 ${RDIR}/velocity-field-x2.nii.gz -v3 ${RDIR}/velocity-field-x3.nii.gz -xfile ${RDIR}/mr_to_mt_ravens_map/mr_to_mt_ravens_transport.nii.gz -r2t

#ibrun tacc_affinity /work/04716/naveen15/stampede2/claire/bin/clairetools -computeravensmap -labels 1 -ifile ${RDIR}/temp.nii.gz -v1 ${RDIR}/velocity-field-x1.nii.gz -v2 ${RDIR}/velocity-field-x2.nii.gz -v3 ${RDIR}/velocity-field-x3.nii.gz -xfile ${RDIR}/temp_ravens_transport.nii.gz

#ibrun tacc_affinity $BDIR/claire -mrc 2 ${RDIR}/0010Y02_GM.nii.gz ${RDIR}/0010Y02_WM.nii.gz -mtc 2 ${RDIR}/0010Y01_GM.nii.gz ${RDIR}/0010Y01_WM.nii.gz -regnorm h2s -beta 1.000000e-02 -opttol 1.000000e-01 -verbosity 2 -x $RDIR/atlas- -velocity -monitordefgrad -detdefgrad
#ibrun tacc_affinity ${BDIR}/claire -mrc 2 ${WORK}/nirep_updated/na02.nii.gz ${WORK}/nirep_updated/na04.nii.gz -mtc 2 ${WORK}/nirep_updated/na01.nii.gz $WORK/nirep_updated/na03.nii.gz -x ${RDIR}/temp- -regnorm h2s -beta 1.000000e-02 -opttol 1.000000e-01 -nt 4 -verbosity 2 -velocity -monitordefgrad -detdefgrad
#ibrun tacc_affinity $BDIR/claire -x $RDIR/ -mrc 2 ${RDIR}/\{0010Y02_GM.nii.gz,0010Y02_WM.nii.gz\} -mtc 2 ${RDIR}/\{0010Y01_GM.nii.gz,0010Y01_WM.nii.gz\} -nx 240x240x155 -pdesolver sl -regnorm h2s -beta 1.000000e-02 -opttol 1.000000e-01 -nt 4 -verbosity 2 -x $RDIR/atlas -velocity -monitordefgrad -detdefgrad

ibrun tacc_affinity ${BDIR}/benchmark -terror -np 4x4 -x ${RDIR}/ -nx 64

