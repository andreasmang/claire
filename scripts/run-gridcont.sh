#!/bin/bash
set -x

# reference image and template image (filename without extension)
mT=na02
mR=na01

betav=1E-2
betaw=1E-4

maxit=1
krylovmaxit=1

# directory where to find refernce and template image
DDIR=/scratch/data/medical_images/brain/nirep/netcdf

# file extension
ext=nc
np=16

defpara='-pdesolver sl -regnorm h1s -opttol 1.000000e-01 -nt 4 -xresult -usenc -logworkload -verbosity 2 -ric -betaw 1.000000e-04'

# directory for binary
BDIR=~/code/develop/cold/bin

# output directory
GDIR=/scratch/andreas/results/compare-regmethods/current-study/nirep/cold/fastsolve_2/h1s-ric/nirep-grid-continuation


LDIR=${mT}-to-${mR}-pdesolver-sl-nt-4-maxit-${maxit}-regnorm-h1s-betav-${betav}-grid-continuation-gtol-0.1

source ../external/libs/environment_vars.sh

RDIR=${GDIR}/${LDIR}
if [[ ! -e ${RDIR} ]]; then
    mkdir -p ${RDIR}
elif [[ ! -d ${RDIR} ]]; then
    echo "${RDIR} already exists but is not a directory" 1>&2
fi

#### submitt job

if [[ ! -e ${RDIR}/velocity-field-x1.${ext} ]]; then
	mpirun -n ${np} ${BDIR}/runcoldreg	-x ${RDIR}/ -mr ${DDIR}/${mR}-nx-256x300x256.${ext} -mt ${DDIR}/${mT}-nx-256x300x256.${ext} \
								-nx 256x300x256 ${defpara} -betav ${betav}  -betaw ${betaw} -maxit ${maxit} -krylovmaxit ${krylovmaxit} > ${RDIR}/log.txt
fi

nt=4

if [ ! -f ${RDIR}/${mT}_seg.${ext} ]; then
	ln -s ${DDIR}/${mT}_seg-nx-256x300x256.${ext} ${RDIR}/${mT}_seg.${ext}
fi

if [ ! -f ${RDIR}/${mT}.${ext} ]; then
	ln -s ${DDIR}/${mT}-nx-256x300x256.${ext} ${RDIR}/${mT}.${ext}
fi

if [ ! -f ${RDIR}/${mT}_seg-transported.${ext} ]; then
	mpirun -n ${np} ${BDIR}/regtools -ifile ${RDIR}/${mT}_seg.${ext} -ivecx1 ${RDIR}/velocity-field-x1.${ext} -ivecx2 ${RDIR}/velocity-field-x2.${ext} -ivecx3 ${RDIR}/velocity-field-x3.${ext} -tlabelmap -nt ${nt} -verbosity 2 -usenc -nx 256x300x256
fi

if [ ! -f ${RDIR}/${mT}-transported.nc ]; then
	mpirun -n ${np} ${BDIR}/regtools -ifile ${RDIR}/${mT}.${ext} -ivecx1 ${RDIR}/velocity-field-x1.${ext} -ivecx2 ${RDIR}/velocity-field-x2.${ext} -ivecx3 ${RDIR}/velocity-field-x3.${ext} -tscafield -nt ${nt} -verbosity 2 -usenc -nx 256x300x256
fi

if [ ! -f ${RDIR}/sl-nt=${nt}_det-deformation-grad.${ext} ]; then
	mpirun -n ${np} ${BDIR}/regtools -i ${RDIR}/ -x ${RDIR}/sl-nt=${nt}_ -detdefgrad -nt ${nt} -verbosity 2 -usenc -nx 256x300x256
fi

if [ ! -f ${RDIR}/sl-nt=${nt}_deformation-map-x1.${ext} ]; then
	mpirun -n ${np} ${BDIR}/regtools -i ${RDIR}/ -x ${RDIR}/sl-nt=${nt}_ -defmap -nt ${nt} -verbosity 2 -usenc -nx 256x300x256
fi


if [ ! -f ${RDIR}/sl-nt=${nt}_inverse-det-deformation-grad.${ext} ]; then
	mpirun -n ${np} ${BDIR}/regtools -i ${RDIR}/ -x ${RDIR}/sl-nt=${nt}_ -invdetdefgrad -nt ${nt} -verbosity 2 -usenc -nx 256x300x256
fi



