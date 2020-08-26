import numpy as np
import os
import nibabel as nib
import argparse
from copy import deepcopy


def writeNII(img, filename, affine=None, ref_image=None):
    '''
    function to write a nifti image, creates a new nifti object
    '''
    if ref_image is not None:
        data = nib.Nifti1Image(img, affine=ref_image.affine, header=ref_image.header);
        data.header['datatype'] = 64
        #data.header['glmax'] = np.amax(img.get_fdata().flatten())
        #data.header['glmin'] = np.amin(img.get_fdata().flatten())
    elif affine is not None:
        data = nib.Nifti1Image(img, affine=affine);
    else:
        data = nib.Nifti1Image(img, np.eye(4))

    nib.save(data, filename);

def MakeImagePeriodic(args):
    input_image = nib.load(args.input_image)
    img = input_image.get_fdata()
    img[:,:,0:args.nz] = 0;
    writeNII(img, args.output_image, ref_image=input_image);

def FilterSingleImage(img, nf):
    img_k = np.fft.fftn(img)

    nf1 = nf - 1;
    N = img.shape
    NM = []
    for i,n in enumerate(list(N)):
        if n%2 == 0:
            NM.append(n//2)
        else:
            NM.append((n-1)//2)

    img_k[NM[0]-nf1-1:NM[0]+2+nf1, :, :] = 0
    img_k[:, NM[1]-nf1-1:NM[1]+2+nf1, :] = 0
    img_k[:, :, NM[2]-nf1-1:NM[2]+2+nf1] = 0

    filtered_img = np.fft.ifftn(img_k)

    return np.real(filtered_img)


def FilterImage(args):
    # this is the path to the claire output folder
    input_path = args.input_image;
    # this is the path to the claire output folder as well
    output_path = args.output_image;
    v1 = nib.load(os.path.join(input_path, 'velocity-field-x1.nii.gz'))
    ref = v1
    v1 = v1.get_fdata()
    v2 = nib.load(os.path.join(input_path, 'velocity-field-x2.nii.gz')).get_fdata()
    v3 = nib.load(os.path.join(input_path, 'velocity-field-x3.nii.gz')).get_fdata()

    nf = 5
    filtered_v1 = FilterSingleImage(v1, nf)
    filtered_v2 = FilterSingleImage(v2, nf)
    filtered_v3 = FilterSingleImage(v3, nf)

    writeNII(filtered_v1, os.path.join(input_path, 'filtered_velocity-field-x1.nii.gz'), ref_image=ref)
    writeNII(filtered_v2, os.path.join(input_path, 'filtered_velocity-field-x2.nii.gz'), ref_image=ref)
    writeNII(filtered_v3, os.path.join(input_path, 'filtered_velocity-field-x3.nii.gz'), ref_image=ref)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='process input images')
    parser.add_argument ('-mode', type=int, help = 'run modes. 1 = zero out the first nz slices in axial direction', required=True);
    parser.add_argument ('-input_image', type=str, help='path to input image', required=True);
    parser.add_argument ('-output_image', type=str, help='path to output image', required=True);
    parser.add_argument ('-nz', type=int, help='zero out the first nz slices in axial direction. Give this argument with -mode 1');
    args = parser.parse_args();


    if args.mode == 1:
        MakeImagePeriodic(args)

    if args.mode == 2:
        FilterImage(args)
