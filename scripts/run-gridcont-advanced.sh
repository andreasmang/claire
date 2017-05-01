#!/bin/bash
set -x

# reference image and template image (filename without extension)
mT=na02
mR=na01

# directory where to find refernce and template image
DDIR=/scratch/data/medical_images/brain/nirep/netcdf

# file extension
ext=nc

defpara='-pdesolver sl -regnorm h1s -precond 2level -pcsolver fpcg -pcsolvermaxit 10 -opttol 1.000000e-01 -nt 4 -xresult -usenc -verbosity 2 -ric -betaw 1.000000e-04'

# number of time steps

# drid levels
gridsizel1=64x75x64
gridsizel2=128x150x128
gridsizel3=256x300x256

npl1=10
npl2=40
npl3=40

# directory for binary
BDIR=~/code/develop/cold/bin

# output directory
GDIR=/scratch/andreas/results/compare-regmethods/current-study/nirep/cold/fastsolve_2/h1s-ric/nirep-grid-continuation


LDIR=${mT}-to-${mR}

# regularization weights
betavl1=5.000000e-3
betavl2=5.000000e-3
betavl3=5.000000e-3


# iterations per level
maxitl1=10
maxitl2=5
maxitl3=2

# max number of krylov iterations per level
krylovmaxitl1=3
krylovmaxitl2=3
krylovmaxitl3=3


# construct output filename
LDIRL1=${LDIR}-pdesolver-maxit-${maxitl1}-krylovmaxit-${krylovmaxitl1}-pc-2level-precond-regnorm-h1s-betav-${betavl1}-grid-continuation-gtol-0.1-${gridsizel1}
LDIRL2=${LDIR}-pdesolver-maxit-${maxitl2}-krylovmaxit-${krylovmaxitl2}-pc-2level-precond-regnorm-h1s-betav-${betavl2}-grid-continuation-gtol-0.1-${gridsizel2}
LDIRL3=${LDIR}-pdesolver-maxit-${maxitl3}-krylovmaxit-${krylovmaxitl3}-pc-2level-precond-regnorm-h1s-betav-${betavl3}-grid-continuation-gtol-0.1-${gridsizel3}


source ../external/libs/environment_vars.sh

########################################
##### level one
########################################
RDIRL1=${GDIR}/${LDIRL1}
if [[ ! -e ${RDIRL1} ]]; then
    mkdir -p ${RDIRL1}
elif [[ ! -d ${RDIRL1} ]]; then
    echo "${RDIRL1} already exists but is not a directory" 1>&2
fi

#### submitt job

if [[ ! -e ${RDIRL1}/velocity-field-x1.${ext} ]]; then
	mpirun -n ${npl1} ${BDIR}/runcoldreg	-x ${RDIRL1}/ -mr ${DDIR}/${mR}-nx-${gridsizel1}.${ext} -mt ${DDIR}/${mT}-nx-${gridsizel1}.${ext} \
								-nx ${gridsizel1} ${defpara} -betav ${betavl1} -maxit ${maxitl1} -krylovmaxit ${krylovmaxitl1} > ${RDIRL1}/log.txt
fi

if [[ ! -e ${RDIRL1}/resampled_velocity-field-x1.${ext} ]]; then
	mpirun -n ${npl1} ${BDIR}/regtools -resample	-ivecx1 ${RDIRL1}/velocity-field-x1.${ext} \
										-ivecx2 ${RDIRL1}/velocity-field-x2.${ext} \
										-ivecx3 ${RDIRL1}/velocity-field-x3.${ext} -rscale 2 -verbosity 2 -nx $gridsizel1
fi




########################################
##### level two
########################################
RDIRL2=${GDIR}/${LDIRL2}
if [[ ! -e ${RDIRL2} ]]; then
    mkdir -p ${RDIRL2}
elif [[ ! -d ${RDIRL2} ]]; then
    echo "${RDIRL2} already exists but is not a directory" 1>&2
fi

if [[ ! -e ${RDIRL2}/initialguess ]]; then
    mkdir -p ${RDIRL2}/initialguess
elif [[ ! -d ${RDIRL2}/initialguess ]]; then
    echo "${RDIRL2}/initialguess already exists but is not a directory" 1>&2
fi


if [[ ! -e ${RDIRL2}/initialguess/velocity-field-x1.${ext} ]]; then
	ln -s ${RDIRL1}/resampled_velocity-field-x1.${ext} ${RDIRL2}/initialguess/velocity-field-x1.${ext}
	ln -s ${RDIRL1}/resampled_velocity-field-x2.${ext} ${RDIRL2}/initialguess/velocity-field-x2.${ext}
	ln -s ${RDIRL1}/resampled_velocity-field-x3.${ext} ${RDIRL2}/initialguess/velocity-field-x3.${ext}
fi

if [[ ! -e ${RDIRL2}/velocity-field-x1.${ext} ]]; then
	mpirun -n ${npl2} ${BDIR}/runcoldreg -x ${RDIRL2}/ -mr ${DDIR}/${mR}-nx-${gridsizel2}.${ext} -mt ${DDIR}/${mT}-nx-${gridsizel2}.${ext} \
										-nx ${gridsizel2} ${defpara} -betav ${betavl2} -maxit ${maxitl2} -krylovmaxit ${krylovmaxitl2}    \
										-vx1 ${RDIRL2}/initialguess/velocity-field-x1.${ext} \
										-vx2 ${RDIRL2}/initialguess/velocity-field-x2.${ext} \
									 -vx3 ${RDIRL2}/initialguess/velocity-field-x3.${ext} > ${RDIRL2}/log.txt
fi

if [[ ! -e ${RDIRL2}/resampled_velocity-field-x1.${ext} ]]; then
	mpirun -n ${npl2} ${BDIR}/regtools -resample -ivecx1 ${RDIRL2}/velocity-field-x1.${ext} -ivecx2 ${RDIRL2}/velocity-field-x2.${ext} -ivecx3 ${RDIRL2}/velocity-field-x3.${ext} -rscale 2 -verbosity 2 -nx $gridsizel2
fi



########################################
##### level three
########################################
RDIRL3=${GDIR}/${LDIRL3}
if [[ ! -e ${RDIRL3} ]]; then
    mkdir -p ${RDIRL3}
elif [[ ! -d ${RDIRL3} ]]; then
    echo "${RDIRL3} already exists but is not a directory" 1>&2
fi

if [[ ! -e ${RDIRL3}/initialguess ]]; then
    mkdir -p ${RDIRL3}/initialguess
elif [[ ! -d ${RDIRL3}/initialguess ]]; then
    echo "${RDIRL3}/initialguess already exists but is not a directory" 1>&2
fi


if [[ ! -e ${RDIRL3}/initialguess/velocity-field-x1.${ext} ]]; then
	ln -s ${RDIRL2}/resampled_velocity-field-x1.${ext} ${RDIRL3}/initialguess/velocity-field-x1.${ext}
	ln -s ${RDIRL2}/resampled_velocity-field-x2.${ext} ${RDIRL3}/initialguess/velocity-field-x2.${ext}
	ln -s ${RDIRL2}/resampled_velocity-field-x3.${ext} ${RDIRL3}/initialguess/velocity-field-x3.${ext}
fi

if [[ ! -e ${RDIRL3}/velocity-field-x1.nc ]]; then
	mpirun -n ${npl3} ${BDIR}/runcoldreg -x ${RDIRL3}/ -mr ${DDIR}/${mR}-nx-${gridsizel3}.${ext} -mt ${DDIR}/${mT}-nx-${gridsizel3}.${ext} \
										-nx ${gridsizel3} ${defpara} -betav ${betavl3} -maxit ${maxitl3} -krylovmaxit ${krylovmaxitl3} \
										-vx1 ${RDIRL3}/initialguess/velocity-field-x1.${ext} \
										-vx2 ${RDIRL3}/initialguess/velocity-field-x2.${ext} \
										-vx3 ${RDIRL3}/initialguess/velocity-field-x3.${ext} > ${RDIRL3}/log.txt
fi


nt=4

if [ ! -f ${RDIRL3}/${mT}_seg.${ext} ]; then
	ln -s ${DDIR}/${mT}_seg-nx-${gridsizel3}.${ext} ${RDIRL3}/${mT}_seg.${ext}
fi

if [ ! -f ${RDIRL3}/${mT}.${ext} ]; then
	ln -s ${DDIR}/${mT}-nx-${gridsizel3}.${ext} ${RDIRL3}/${mT}.${ext}
fi

if [ ! -f ${RDIRL3}/${mT}_seg-transported.${ext} ]; then
	mpirun -np ${npl3} ${BDIR}/regtools -ifile ${RDIRL3}/${mT}_seg.${ext} -ivecx1 ${RDIRL3}/velocity-field-x1.${ext} -ivecx2 ${RDIRL3}/velocity-field-x2.${ext} -ivecx3 ${RDIRL3}/velocity-field-x3.${ext} -tlabelmap -nt ${nt} -verbosity 2 -usenc -nx ${gridsizel3}
fi

if [ ! -f ${RDIRL3}/${mT}-transported.nc ]; then
	mpirun -np ${npl3} ${BDIR}/regtools -ifile ${RDIRL3}/${mT}.${ext} -ivecx1 ${RDIRL3}/velocity-field-x1.${ext} -ivecx2 ${RDIRL3}/velocity-field-x2.${ext} -ivecx3 ${RDIRL3}/velocity-field-x3.${ext} -tscafield -nt ${nt} -verbosity 2 -usenc -nx ${gridsizel3}
fi

if [ ! -f ${RDIRL3}/sl-nt=${nt}_det-deformation-grad.${ext} ]; then
	mpirun -np ${npl3} ${BDIR}/regtools -i ${RDIRL3}/ -x ${RDIRL3}/sl-nt=${nt}_ -detdefgrad -nt ${nt} -verbosity 2 -usenc -nx ${gridsizel3}
fi

if [ ! -f ${RDIRL3}/sl-nt=${nt}_deformation-map-x1.${ext} ]; then
	mpirun -np ${npl3} ${BDIR}/regtools -i ${RDIRL3}/ -x ${RDIRL3}/sl-nt=${nt}_ -defmap -nt ${nt} -verbosity 2 -usenc -nx ${gridsizel3}
fi


if [ ! -f ${RDIRL3}/sl-nt=${nt}_inverse-det-deformation-grad.${ext} ]; then
	mpirun -np ${npl3} ${BDIR}/regtools -i ${RDIRL3}/ -x ${RDIRL3}/sl-nt=${nt}_ -invdetdefgrad -nt ${nt} -verbosity 2 -usenc -nx ${gridsizel3}
fi



