# CLAIRE: Quick Installation Guide

Here, we only provide a minimal installation guide. We provide scripts to download and install the dependencies using generic settings that have worked on our systems. If this brief installation guide does not work for you, please consult the [detailed installation guide](README-INSTALL.md) found in [README-INSTALL.md](README-INSTALL.md).

## Contents
1. [One Shot](#oneshot)
2. [Step by Step](#stepbystep)

## One Shot <a name="oneshot"></a>


```bash
cd deps
./get_libs.sh
./build_libs.sh --build
cd ..
make -j
./bin/claire -synthetic 0
```

For GPU build on TACC's Longhorn system (need to be on the [gpu](https://github.com/andreasmang/claire/tree/gpu) branch of CLAIRE)

```bash
cd deps
./get_libs.sh
./build_libs.sh --build --gpu=V100 --enableCUDA --POWER9
cd ..
make -f makefile_p9 -j
mpirun -np 1 ./bindev/claire -synthetic 0
```



## Step by Step  <a name="stepbystep"></a>


### Step 1) Installing Dependencies

To install the dependencies go to the top level directory of CLAIRE in your command window and execute the following commands within your command window:

```bash
cd deps
./get_libs.sh
./build_libs.sh --build
```

The first bash script downloads the libraries you need to install. The second bash script compiles these libraries using default settings. To add these dependencies to your environment type the following into your command line and press return:

```bash
source environment_vars.sh
```

This step needs to be done every time you log out of your computer and would like to recompile CLAIRE. As an alternative, you can add the content of `environment_vars.sh` to your `.bashrc` or `bash_profile`.

### Step 2) Compiling CLAIRE

Assuming that your in the top level directory of CLAIRE, all you need to do is to type

```bash
make -j
```

### Step 3) Executing CLAIRE

If you would like to verify if CLAIRE has been installed correctly run the following command in your command window:

```bash
./bin/claire -synthetic 0
```


### Troubleshooting

* If you have trouble **installing** CLAIRE consult the [detailed installation guide](README-INSTALL.md) found in [README-INSTALL.md](README-INSTALL.md).
* If you have trouble executing CLAIRE take a look at the [examples](README-RUNME.md) described in [README-RUNME.md](README-RUNME.md).
