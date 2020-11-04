#!/usr/bin/python3

import nibabel as nib
import nibabel.processing as nipr
import numpy as np
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input nii.gz", required=True)
parser.add_argument("-i2", help="input2 nii.gz")
parser.add_argument("-o", help="input nii.gz")
parser.add_argument("-s", help="smoothing size XxYxZ")
parser.add_argument("-n", help="rescale size XxYxZ")
parser.add_argument("-l", help="pointwise lambda v, pos, stat:")
parser.add_argument("-l2", help="pointwise lambda v1, v2, pos, stat1, stat2:")

xyz_re = re.compile("(?P<x>\d+)x(?P<y>\d+)x(?P<z>\d+)")

def parse_xyz(arg):
  m = xyz_re.match(arg)
  if m:
    return int(m.group('x')), int(m.group('y')), int(m.group('z'))
  else:
    return None

args = parser.parse_args()

nifti = nib.load(args.i)

nx, ny, nz = nifti.header['dim'][1:4]
min_val = np.min(nifti.dataobj)
max_val = np.max(nifti.dataobj)

print("Input Data %ix%ix%i [%f - %f]"%(nx,ny,nz,min_val, max_val))

if args.n:
  nx_new, ny_new, nz_new = parse_xyz(args.n)
  print("scale to %ix%ix%i"%(nx_new,ny_new,nz_new))
  p_x, p_y, p_z = nifti.header['pixdim'][1:4]
  nifti = nipr.resample_to_output(nifti,[p_x*nx/nx_new,p_y*ny/ny_new,p_z*nz/nz_new],order=3)
  nx = nx_new
  ny = ny_new
  nz = nz_new
  
if args.l:
  data = np.array(nifti.dataobj)
  kernel = eval("lambda v, pos, stat : " + args.l)
  min_val = np.min(nifti.dataobj)
  max_val = np.max(nifti.dataobj)
  for x in range(nx):
    for y in range(ny):
      for z in range(nz):
        data[x,y,z] = kernel(data[x,y,z],(x,y,z),(nx,ny,nz,min_val,max_val))
  nifti = nib.Nifti1Image(data, nifti.affine, nifti.header)

if args.s:
  data = np.array(nifti.dataobj)
  lx, ly, lz = parse_xyz(args.s)

  for x in range(lx):
    for y in range(ny):
      for z in range(nz):
        data[x,y,z] *= float(x)/lx
        data[-1-x,y,z] *= float(x)/lx
        
  for x in range(nx):
    for y in range(ly):
      for z in range(nz):
        data[x,y,z] *= float(y)/ly
        data[x,ny-1-y,z] *= float(y)/ly
        
  for x in range(nx):
    for y in range(ny):
      for z in range(lz):
        data[x,y,z] *= float(z)/lz
        data[x,y,nz-1-z] *= float(z)/lz

  nifti = nib.Nifti1Image(data, nifti.affine, nifti.header)
  
if args.l2:
  nifti2 = nib.load(args.i2)
  data = np.array(nifti.dataobj)
  data2 = np.array(nifti2.dataobj)
  kernel = eval("lambda v1, v2, pos, stat1, stat2 : " + args.l2)
  min_val1 = np.min(nifti.dataobj)
  min_val2 = np.min(nifti2.dataobj)
  max_val1 = np.max(nifti.dataobj)
  max_val2 = np.max(nifti2.dataobj)
  for x in range(nx):
    for y in range(ny):
      for z in range(nz):
        data[x,y,z] = kernel(data[x,y,z],data2[x,y,z],(x,y,z),(nx,ny,nz,min_val2,max_val2),(nx,ny,nz,min_val2,max_val2))
  nifti = nib.Nifti1Image(data, nifti.affine, nifti.header)


if args.o:
  nx, ny, nz = nifti.header['dim'][1:4]
  min_val = np.min(nifti.dataobj)
  max_val = np.max(nifti.dataobj)
  print("Output Data %ix%ix%i [%f - %f]"%(nx,ny,nz,min_val, max_val))
  nib.save(nifti, args.o)
