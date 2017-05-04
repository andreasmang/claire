#!/bin/bash
set -x

# reference image and template image (filename without extension)
mT=0013Y02
mR=0014Y02-to-0013Y02

betav=1E-4
betaw=1E-4

maxit=20
krylovmaxit=10

# directory where to find refernce and template image
DDIR=/scratch/data/medical_images/brain/698_templates/updated

# file extension
ext=nc
np=40

#defpara='-pdesolver sl -regnorm h1s -opttol 1.000000e-01 -nt 4 -xresult -usenc -verbosity 2 -ric -betaw 1.000000e-04'
defpara='-pdesolver sl -regnorm h2s -opttol 1.000000e-01 -nt 4 -xresult -usenc -debug -disablesmoothing'

# directory for binary
BDIR=~/code/develop/cold/bin

# output directory
GDIR=/scratch/andreas/results/compare-regmethods/current-study/multi-component


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
	mpirun -n ${np} ${BDIR}/runcoldreg	-mrc 4	${DDIR}/${mR}-wm-128x128x128-sigma-2.${ext}  \
												${DDIR}/${mR}-gm-128x128x128-sigma-2.${ext}  \
												${DDIR}/${mR}-csf-128x128x128-sigma-2.${ext} \
												${DDIR}/${mR}-glm-128x128x128-sigma-2.${ext} \
										-mtc 4	${DDIR}/${mT}-wm-128x128x128-sigma-2.${ext}  \
												${DDIR}/${mT}-gm-128x128x128-sigma-2.${ext}  \
												${DDIR}/${mT}-csf-128x128x128-sigma-2.${ext} \
												${DDIR}/${mT}-glm-128x128x128-sigma-2.${ext} \
												-x ${RDIR}/ -nx 128x128x128 ${defpara} -betav ${betav} -maxit ${maxit} -krylovmaxit ${krylovmaxit} -xdetdefgrad -xdeftemplate -xmr -xmt # > ${RDIR}/log.txt
fi



