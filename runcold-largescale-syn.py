# this is a python script, that creates a bash script to control sbatch and submits it
# (maverick: max: 32 nodes; 20 cores per node)
import time
from subprocess import call
import os

wch=0; # max wall clock hours
wcm=40; # max wall clock minutes

regprob = "syn";

cluster = "stampede2";
#cluster = "lonestar";
#cluster = "maverick";

if cluster == "stampede2":
    ncpn = 12; # stampede2 max 68 cores per node
elif cluster == "lonestar":
    ncpn = 12; # lonestar: 4104 cores
else:
    ncpn = 20; # maverick

mpitasks = 32;
nnodes = 1 + int(mpitasks/ncpn);
reqmpitasks = nnodes*ncpn;
ompthreads = 3;

betav = 1E-2;
betaw = 1E-4;

scale = 1;
greltol = 1E-2;
regnorm = "h1s"; # regularization norm
opttype = "gnls"; # optimization type
precond = 'invreg';
pdesolver = "sl";
nx1 = int(scale*256);
nx2 = int(scale*300);
nx3 = int(scale*256);
nt = int(4);

maxit = 5;
krylovmaxit = 10;

warmstart = False;
reducebetav = False;
gridcont = False;
learnbeta = False;

regmodel = 'ric';

gridstr = regprob + "-nx-" + str(nx1) + "x" + str(nx2) + "x" + str(nx3);
basedir = os.environ['WORK'] + "/cold/results/" + gridstr;

bindir  = os.environ['WORK'] + "/cold/bin/"
datadir = os.environ['WORK'] + "/data/nirep/"
libdir  = os.environ['WORK'] + "/cold/external/libs/"


outdir = basedir;
lxdir = "smooth";
lxdir = lxdir + "-nodes-" + str(nnodes);
lxdir = lxdir + "-mpitasks-" + str(mpitasks);
lxdir = lxdir + "-threads-" + str(ompthreads);
lxdir = lxdir + "-pdesolver-" + pdesolver;
lxdir = lxdir + "-nt-" + str(nt);
lxdir = lxdir + "-regnorm-" + regnorm;

if reducebetav:
    lxdir = lxdir + "-betavcont";
    lxdir = lxdir + "-betav-" + str(betavtarget);
elif learnbeta:
    lxdir = lxdir + "-betav-estimated";
elif gridcont:
    lxdir = lxdir + "-grid-continuation";
else:
    lxdir = lxdir + "-betav-" + str(betav);

if regmodel == "ric":
    lxdir=lxdir + "-betaw-" + str(betaw);

lxdir=lxdir + "-gtol-" + str(greltol);

if regmodel == "ric":
    lxdir=lxdir + "-ric";

# make output directory
if not os.path.exists(outdir + "/" + lxdir):
    os.makedirs(outdir + "/" + lxdir);
else:
    print "run already performed"
    quit(); # if already computed, go home

# construct bash file name
bashjobfile = "job-submission.sh";
bashfilename = outdir + "/" + lxdir + "/" + bashjobfile;

# if output tab line exists, delete it
if os.path.isfile(bashfilename):
    os.remove(bashfilename);


# create bash file
bashfile = open(bashfilename,'w');

# heading
bashfile.write("#!/bin/bash\n\n");
bashfile.write("\n");
bashfile.write("#### sbatch parameters\n");
bashfile.write("#SBATCH -J coldreg\n");
bashfile.write("#SBATCH -o coldreg.o%j\n");
bashfile.write("#SBATCH -e coldreg.e%j\n");
bashfile.write("#SBATCH -p normal\n");
bashfile.write("#SBATCH -N " + str(nnodes) + "\n");
bashfile.write("#SBATCH -n " + str(reqmpitasks) + "\n");

if cluster == "lonestar":
    if nnodes > 171:
        bashfile.write("#SBATCH -p large\n"); # large or normal
    else:
        bashfile.write("#SBATCH -p normal\n"); # large or normal

bashfile.write("#SBATCH -t " + str(wch) + ":" + str(wcm) + ":00\n");
bashfile.write("#SBATCH --mail-user=naveen@ices.utexas.edu\n");
bashfile.write("#SBATCH --mail-type=all\n");
bashfile.write("#SBATCH --mail-type=fail\n");
bashfile.write("\n");
bashfile.write("\n");
bashfile.write("\n#### define paths\n");
bashfile.write("export OMP_NUM_THREADS="+str(ompthreads)+ "\n" );
bashfile.write("BDIR=" + bindir + "\n");
bashfile.write("RDIR=" + outdir + "\n");
bashfile.write("LDIR=" + lxdir + "\n");
bashfile.write("RDIR=$RDIR/$LDIR\n")
bashfile.write("\n");
bashfile.write("\n");


# construct command line
commandline="ibrun -o 0 -n " + str(mpitasks) + " tacc_affinity $BDIR/runcoldreg";


commandline = commandline + " -x $RDIR/";
commandline = commandline + " -nx " + str(nx1) + "x" + str(nx2) + "x" + str(nx3);
commandline = commandline + " -pdesolver " + pdesolver;
commandline = commandline + " -regnorm " + regnorm;
commandline = commandline + " -maxit " + str(maxit);
commandline = commandline + " -betav " + str(betav);
commandline = commandline + " -usenc";
commandline = commandline + " -krylovmaxit " + str(krylovmaxit);
if regmodel == "ric":
    commandline = commandline + " -betaw " + '{0:.6e}'.format(betaw);
    commandline = commandline + " -ric";

commandline = commandline + " -opttol " + '{0:.6e}'.format(greltol);
commandline = commandline + " -nt " + str(nt);
commandline = commandline + " -logworkload ";

commandline = commandline + "\n";
commandline = commandline;

bashfile.write("#### submitt job\n");
bashfile.write(commandline);

bashfile.write("\n");

# write out done
bashfile.close();

# change to directory
os.chdir(outdir + "/" + lxdir);

# submitt the job
call(["sbatch",bashjobfile]);

time.sleep(1);


