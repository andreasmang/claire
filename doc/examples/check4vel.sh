#!/bin/bash

v1=velocity-field-x1.nii.gz
v2=velocity-field-x2.nii.gz
v3=velocity-field-x3.nii.gz

# run the registration first if velocities do not exist
if [ -f ${v1} ] && [ -f ${v2} ] && [ -f ${v3} ]; then
   echo "velocities found"
else
	# execute the registration to get some velocities
	./runclaire06.sh
fi
