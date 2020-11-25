from mpi4py import MPI
import pyclaire
#import nibabel as nib



obj = pyclaire.claire()

config = "-synthetic 4 -betacont 5e-4 -regnorm h1s-div -precond invreg -diffpde finite -iporder 1 -nx 256x256x256 -verbosity 2"

obj.setParameters(config.split(" "))

obj.runRegistration()
