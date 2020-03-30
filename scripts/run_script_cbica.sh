#!/bin/bash

# Script to 
# 1. Warp subject to template space
# Pre-requisites: affine registration and affine transform file is needed as input


# fixed arguments
CLAIRE_DIR=${HOME}/claire
CLAIRE_SCRIPT_DIR=$CLAIRE_DIR/scripts
CLAIRE_BDIR=$CLAIRE_DIR/bin


myline="----------------------------------------------------------------------------------"

for i in "$@"
do
case $i in
    --subject=*)
    SUBJECT="${i#*=}"    
	echo $SUBJECT
    shift # past argument=value
    ;;
    --subject_label=*)
    SUBJECT_LAB="${i#*=}"
    shift # past argument=value
    ;;
    --template=*)
    TEMPLATE="${i#*=}"
    shift # past argument=value
    ;;
    --N=*)
    N="${i#*=}"
    shift # past argument=value
    ;;
    --affine=*)
    AFFINE="${i#*=}"
    shift # past argument=value
    ;;
    --x=*)
    OUTPUT_DIR="${i#*=}"
    shift # past argument=value
    ;;
    --mode=*)
    MODE="${i#*=}"
    shift # past argument=value
    ;;
    --i=*)
    INPUT_DIR="${i#*=}"
    shift # past argument=value
    ;;
    --labels=*)
    LABELS="${i#*=}"
    shift # past argument=value
    ;;
    --distance_metric=*)
    DISTANCE="${i#*=}"
    shift # past argument=value
    ;;
    --parameter_cont=*)
    PARAM_CONT="${i#*=}"
    shift # past argument=value
    ;;
    --parameter=*)
    PARAM="${i#*=}"
    shift # past argument=value
    ;;
    --scale_factor=*)
    SCALE_FACTOR="${i#*=}"
    shift # past argument=value
    ;;
    --reg_res=*)
    REG_RES="${i#*=}"
    shift # past argument=value
    ;;
    --help)
    echo "script to run claire to warp subject image to template image" 
    echo ${myline}    
    echo " options for this script are"
    echo ${myline}
    echo "     --help                           print this message"
    echo ${myline}
    echo "     --subject=*.nii.gz               absolute path to subject image before affine registration(moving image)"
    echo "     --subject_label=*.nii.gz         absolute path to subject segmentation before affine registration(moving image segmentation)"
    echo "     --labels=                        comma separated label ids in the subject(moving) image segmentation. e.g. 10,50,150,250 (needed)"
    echo "     --template=*.nii.gz              absolute path to the template image(fixed image)"
    echo "     --N=lxmxn                        size of template image. eg. 182x218x218 (needed)"
    echo "     --reg_res=lxmxn                  resolution in which to run claire. This should typically be higher or the same resolution"
    echo "                                      as the input image. If the input image resolution (--N) is smaller than 256x256x256, consider giving"
    echo "                                      --reg_res=256x256x256 as the input, otherwise if the input image resolution is higher than"
    echo "                                      256x256x256, then try giving the input image dimensions as the input to --reg_res. Upsampling"
    echo "                                      in claire is usually required to avoid unwanted aliasing effects"
    echo "     --affine=*.txt                   absolute path to the Affine.txt file genereated from ANTS affine registration of subject(moving) to"
    echo "                                      template(fixed)"
    echo "     --x=/path/to/output              path to the output directory"
    echo "     --mode=<1,2>                     run mode. 1: warp subject image to template image (default)"
    echo "                                                2: apply claire velocities to a new image"
    echo "                                                   (need to provide --i with this mode)"
    echo "     --i=/path/to/claire/results      path to directory where claire velocity fields are stored."
    echo "     --distance_metric=<sl2,ncc>      distance metric for claire. options: 1. sl2 - sum of squared differences"
    echo "                                                                           2. ncc - normalized cross correlation (default)"    
    echo "     --parameter_cont=<0,1>           whether to perform parameter continuation for diffeomorphic registration (suggested). default: 1"
    echo "     --parameter=<float>              parameter(regularization) value for diffeomorhic registration if --parameter_cont=0 is given"
    echo "                                      default: 1e-2"
    echo "     --scale_factor=                  scaling factor for RAVENS (i.e Jacobian Determinant in Template space). default: 1"
    echo ${myline}
    echo ""
    exit;
    shift # past argument=value
    ;;
    *)
    # unknown option
    ;;
esac
shift
done


###########################################
if [ -z "$SUBJECT" ]; then
    echo "--subject is needed"
    exit
fi
SUBJECT_FNAME="${SUBJECT##*/}"

###########################################
if [ -z "$SUBJECT_LAB" ]; then
    echo "--subject_label is needed"
    exit  
fi
SUBJECT_LAB_FNAME="${SUBJECT_LAB##*/}"

###########################################
if [ -z "$TEMPLATE" ]; then
    echo "--template is needed"
    exit
fi
TEMPLATE_FNAME="${TEMPLATE##*/}"

###########################################
if [ -z "$AFFINE" ]; then
    echo "--affine is needed"
    exit
fi

###########################################
if [ -z "$LABELS" ]; then
    echo "--labels is needed"
    exit
fi

###########################################
if [ -z "$OUTPUT_DIR" ]; then
    echo "--x not provided, seetting to current directory"
    OUTPUT_DIR=$PWD
fi

###########################################
if [ -z "$MODE" ]; then
    echo "--mode not provided, setting mode to 1"
    MODE=1
fi

###########################################
if [ "$MODE" == "2" ]; then 
    if [ -z "$INPUT_DIR" ]; then
        echo "with with mode 2, --i is needed"
        exit
    fi
fi

###########################################
if [ -z "$DISTANCE" ]; then
    echo "--distance_metric not provided, setting to ncc (normalized cross correlation"
    DISTANCE=ncc
fi

###########################################
if [ -z "$PARAM_CONT" ]; then
    echo "--parameter_cont not provided, setting to 1(binary search)"
    PARAM_CONT=binary
fi

###########################################
if [ -z "$PARAM" ]; then
    echo "--parameter not provided, setting to 1e-2"
    PARAM=1e-2
fi

###########################################
if [ -z "$SCALE_FACTOR" ]; then
    echo "--scale_factor not provided, setting to 1"
    SCALE_FACTOR=1
fi

###########################################
if [ -z "$REG_RES" ]; then
    echo "--reg_res needs to be provided"
    exit;
fi

###########################################
if [ -z "$N" ]; then
    echo "--N needs to be provided"
    exit;
fi

echo "subject $SUBJECT"
echo "subject_label $SUBJECT_LAB"
echo "template $TEMPLATE"
echo "affine $AFFINE"
echo "mode $MODE"
echo "labels $LABELS"
echo "output $OUTPUT_DIR"

mkdir -p $OUTPUT_DIR

if [ $MODE == 1 ]; then


# apply AFFINE transform to subject T1
antsApplyTransforms \
    -d 3 \
    -i ${SUBJECT} \
    -r ${TEMPLATE} \
    -o ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz \
    -n Linear \
    -t ${AFFINE}

# apply AFFINE transform to subject labels
antsApplyTransforms \
    -d 3 \
    -i ${SUBJECT_LAB} \
    -r ${TEMPLATE} \
    -o ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_affine_warped.nii.gz \
    -n GenericLabel[Linear] \
    -t ${AFFINE}

# make the input image periodic
python ${CLAIRE_SCRIPT_DIR}/utils.py -mode 1 -input_image ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz -output_image ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz -nz 5

python ${CLAIRE_SCRIPT_DIR}/utils.py -mode 1 -input_image ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_affine_warped.nii.gz -output_image ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_affine_warped.nii.gz -nz 5


# get jacobian determinant of affine transformation
# first compose affine transform
ComposeMultiTransform \
    3 \
    ${OUTPUT_DIR}/affine.nii.gz \
    -R ${TEMPLATE} \
    ${AFFINE}
# now compute affine jacobian determinant
CreateJacobianDeterminantImage \
    3 \
    ${OUTPUT_DIR}/affine.nii.gz \
    ${OUTPUT_DIR}/affine_detdefgrad.nii.gz

if [ ! $REG_RES == ${N} ]; then
    # resample the input images (template and affine warped image to required resolution)
    ResampleImage \
    3 \
    ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz \
    ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped_upsampled.nii.gz \
    ${REG_RES} \
    1,0 \
    0 \
    7

    ResampleImage \
    3 \
    $TEMPLATE \
    $OUTPUT_DIR/${TEMPLATE_FNAME%.nii.gz}_upsampled.nii.gz \
    ${REG_RES} \
    1,0 \
    0 \
    7
else
    ln -sf ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz  ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped_upsampled.nii.gz 
    ln -sf $TEMPLATE  $OUTPUT_DIR/${TEMPLATE_FNAME%.nii.gz}_upsampled.nii.gz
fi


# run claire on the AFFINE warped T1 image
if [ ! -f ${OUTPUT_DIR}/velocity-field-x1.nii.gz ]; then
    if [ "$PARAM_CONT" == binary ]; then
    	mpirun --prefix $OPENMPI -np $NSLOTS \
        $CLAIRE_BDIR/claire \
        -x $OUTPUT_DIR/ \
        -mt ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped_upsampled.nii.gz \
        -mr $OUTPUT_DIR/${TEMPLATE_FNAME%.nii.gz}_upsampled.nii.gz \
        -regnorm h1s-div \
        -train binary \
        -jbound 0.2 \
        -format nifti \
        -velocity \
        -verbosity 2 \
        -distance ${DISTANCE} \
        > $OUTPUT_DIR/solver_log.txt
    else
        mpirun --prefix $OPENMPI -np $NSLOTS \
        $CLAIRE_BDIR/claire \
        -x $OUTPUT_DIR/ \
        -mt ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped_upsampled.nii.gz \
        -mr $OUTPUT_DIR/${TEMPLATE_FNAME%.nii.gz}_upsampled.nii.gz \
        -regnorm h1s-div \
        -beta ${PARAM} \
        -format nifti \
        -velocity \
        -verbosity 2 \
        -distance ${DISTANCE} \
        > $OUTPUT_DIR/solver_log.txt
fi
fi


# filter the velocity field to get rid of aliasing (remove the 5 highest frequencies)
python $CLAIRE_SCRIPT_DIR/utils.py -mode 2 -input_image $OUTPUT_DIR/ -output_image $OUTPUT_DIR/ 


# if resampling was done then downsample the velocity field
if [ ! $REG_RES == ${N} ]; then
    ResampleImage \
    3 \
    $OUTPUT_DIR/filtered_velocity-field-x1.nii.gz \
    $OUTPUT_DIR/downsampled_velocity-field-x1.nii.gz \
    $N \
    1,0 \
    4 \
    7

    ResampleImage \
    3 \
    $OUTPUT_DIR/filtered_velocity-field-x2.nii.gz \
    $OUTPUT_DIR/downsampled_velocity-field-x2.nii.gz \
    $N \
    1,0 \
    4 \
    7

    ResampleImage \
    3 \
    $OUTPUT_DIR/filtered_velocity-field-x3.nii.gz \
    $OUTPUT_DIR/downsampled_velocity-field-x3.nii.gz \
    $N \
    1,0 \
    4 \
    7
else
    # simple sylink to the original velocity if no upsampling was done
    ln -sf $OUTPUT_DIR/filtered_velocity-field-x1.nii.gz $OUTPUT_DIR/downsampled_velocity-field-x1.nii.gz
    ln -sf $OUTPUT_DIR/filtered_velocity-field-x2.nii.gz $OUTPUT_DIR/downsampled_velocity-field-x2.nii.gz
    ln -sf $OUTPUT_DIR/filtered_velocity-field-x3.nii.gz $OUTPUT_DIR/downsampled_velocity-field-x3.nii.gz
fi

# run apply the claire warp field on the affine registered subject image
if [ -f ${OUTPUT_DIR}/downsampled_velocity-field-x1.nii.gz ]; then
	mpirun --prefix $OPENMPI -np $NSLOTS \
    $CLAIRE_BDIR/clairetools \
    -deformimage \
    -ifile ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz \
    -xfile ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_rTemplate.nii.gz \
    -v1 $OUTPUT_DIR/downsampled_velocity-field-x1.nii.gz \
    -v2 $OUTPUT_DIR/downsampled_velocity-field-x2.nii.gz \
    -v3 $OUTPUT_DIR/downsampled_velocity-field-x3.nii.gz \
    > $OUTPUT_DIR/transport_log.txt

	mpirun --prefix $OPENMPI -np $NSLOTS \
    $CLAIRE_BDIR/clairetools \
    -tlabelmap \
    -labels ${LABELS} \
    -ifile ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_affine_warped.nii.gz \
    -xfile ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate.nii.gz \
    -v1 $OUTPUT_DIR/downsampled_velocity-field-x1.nii.gz \
    -v2 $OUTPUT_DIR/downsampled_velocity-field-x2.nii.gz \
    -v3 $OUTPUT_DIR/downsampled_velocity-field-x3.nii.gz \
    >> $OUTPUT_DIR/transport_log.txt

	mpirun --prefix $OPENMPI -np $NSLOTS \
    $CLAIRE_BDIR/clairetools \
    -invdetdefgrad \
    -x ${OUTPUT_DIR}/ \
    -v1 $OUTPUT_DIR/downsampled_velocity-field-x1.nii.gz \
    -v2 $OUTPUT_DIR/downsampled_velocity-field-x2.nii.gz \
    -v3 $OUTPUT_DIR/downsampled_velocity-field-x3.nii.gz \
    >> $OUTPUT_DIR/transport_log.txt
fi

# compute jacobian determinant of affine+CLAIRE by multiplying the 2 jacobians
if [ -f ${OUTPUT_DIR}/det-deformation-grad.nii.gz ]; then
	3dcalc \
    -prefix ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_JacDet.nii.gz \
    -a ${OUTPUT_DIR}/det-deformation-grad.nii.gz \
    -b ${OUTPUT_DIR}/affine_detdefgrad.nii.gz \
    -expr "a*b" \
    -verbose \
    -nscale \
    -float;
fi

for l in $(echo ${LABELS} | tr "," " "); do
    #echo "computing ravens for ${l}"
    # extract the label image
    ThresholdImage \
    3 \
    ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate.nii.gz \
    ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate_${l}.nii.gz \
    ${l} ${l} \
    1 0

    # compute ravens map by masking the jacobian determinant with the label image
    3dcalc \
    -prefix ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_RAVENS_${l}.nii.gz \
    -a ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_JacDet.nii.gz \
    -b ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate_${l}.nii.gz \
    -expr "a*b*$SCALE_FACTOR" \
    -verbose \
    -nscale \
    -float;
done

fi


if [ $MODE == 2 ]; then

    # if it is indeed supplied then use it
    # apply AFFINE transform to subject T1
    antsApplyTransforms \
        -d 3 \
        -i ${SUBJECT} \
        -r ${TEMPLATE} \
        -o ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz \
        -n Linear \
        -t ${AFFINE}

    # apply AFFINE transform to subject labels
    antsApplyTransforms \
        -d 3 \
        -i ${SUBJECT_LAB} \
        -r ${TEMPLATE} \
        -o ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_affine_warped.nii.gz \
        -n GenericLabel[Linear] \
        -t ${AFFINE}

    # get jacobian determinant of affine transformation
    # first compose affine transform
    ComposeMultiTransform \
        3 \
        ${OUTPUT_DIR}/affine.nii.gz \
        -R ${TEMPLATE} \
        ${AFFINE}

    # now compute affine jacobian determinant
    CreateJacobianDeterminantImage \
         3 \
         ${OUTPUT_DIR}/affine.nii.gz \
         ${OUTPUT_DIR}/affine_detdefgrad.nii.gz

    # run apply the claire warp field on the affine registered subject image
    # use the velocity field from the input directory of claire results
    if [ -f ${INPUT_DIR}/downsampled_velocity-field-x1.nii.gz ]; then
        mpirun --prefix $OPENMPI -np $NSLOTS \
            $CLAIRE_BDIR/clairetools \
            -deformimage \
            -ifile ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_affine_warped.nii.gz \
            -xfile ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_rTemplate.nii.gz \
            -v1 $INPUT_DIR/downsampled_velocity-field-x1.nii.gz \
            -v2 $INPUT_DIR/downsampled_velocity-field-x2.nii.gz \
            -v3 $INPUT_DIR/downsampled_velocity-field-x3.nii.gz \
            > $OUTPUT_DIR/transport_log.txt

        mpirun --prefix $OPENMPI -np $NSLOTS \
            $CLAIRE_BDIR/clairetools \
            -tlabelmap \
            -labels ${LABELS} \
            -ifile ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_affine_warped.nii.gz \
            -xfile ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate.nii.gz \
            -v1 $INPUT_DIR/downsampled_velocity-field-x1.nii.gz \
            -v2 $INPUT_DIR/downsampled_velocity-field-x2.nii.gz \
            -v3 $INPUT_DIR/downsampled_velocity-field-x3.nii.gz \
            >> $OUTPUT_DIR/transport_log.txt

        mpirun --prefix $OPENMPI -np $NSLOTS \
            $CLAIRE_BDIR/clairetools \
            -invdetdefgrad \
            -x ${OUTPUT_DIR}/ \
            -v1 $INPUT_DIR/downsampled_velocity-field-x1.nii.gz \
            -v2 $INPUT_DIR/downsampled_velocity-field-x2.nii.gz \
            -v3 $INPUT_DIR/downsampled_velocity-field-x3.nii.gz \
            >> $OUTPUT_DIR/transport_log.txt
    else
	echo "velocity field does not exist in the directory of the claire results. Compute it first using mode 1. Exiting"
	exit 1
    fi


    # compute jacobian determinant of affine+CLAIRE by multiplying the 2 jacobians
    if [ -f ${OUTPUT_DIR}/det-deformation-grad.nii.gz ]; then
        3dcbinary        -prefix ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_JacDet.nii.gz \
        -a ${OUTPUT_DIR}/det-deformation-grad.nii.gz \
        -b ${OUTPUT_DIR}/affine_detdefgrad.nii.gz \
        -expr "a*b" \
        -verbose \
        -nscale \
        -float;
    fi

    for l in $(echo ${LABELS} | tr "," " "); do
        # extract binarybel image
        ThresholdImage \
        3 \
        ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate.nii.gz \
        ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate_${l}.nii.gz \
        ${l} ${l} \
        1 0

        # compute ravens map by masking the jacobian determinant with the label image
        3dcalc \
        -prefix ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_RAVENS_${l}.nii.gz \
        -a ${OUTPUT_DIR}/${SUBJECT_FNAME%.nii.gz}_JacDet.nii.gz \
        -b ${OUTPUT_DIR}/${SUBJECT_LAB_FNAME%.nii.gz}_rTemplate_${l}.nii.gz \
        -expr "a*b*$SCALE_FACTOR" \
        -verbose \
        -nscale \
        -float;
    done
fi
