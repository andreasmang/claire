# Installing and Running COLDREG


## Before Compiling

Make sure that the standard **MPI wrappers** for `mpicc` and `mpicxx` are available on your system (either by loading the appropriate modules and/or by setting up the appropriate `PATH` and `LD_LIBRARY_PATH` definitions). Add the following to your `~/.bashrc`:

```bash
export PATH=/path/to/mpicxx:/path/to/mpicc:${PATH}
export LD_LIBRARY_PATH=/path/to/mpi/lib:${LD_LIBRARY_PATH}
```


## Build COLDREG

Before you can build COLDREG you need to 

* make sure that you have installed the **external dependencies** (visit [doc/README-EXTLIBS.md](README-EXTLIBS.md) to learn more)

* check the [makefile](makefile) before building the code:
	* if you use an **intel compiler** (`icc`) set the `USEINTEL` flag to `yes`
	* if you use a **GNU compiler** (`gcc`) set the `USEINTEL` flag to `no`
	* if you use a **IntelMPI** (`impi`) set the `USEINTELMPI` flag to `yes` (if not, set it to `no`)
* make sure all paths needed in the makefile are available on your system (to check, you can do `env` in your bash); to add the pathes necessary to find and link against the library you can `source libs/environment_vars.sh` or add the content of `libs/environment_vars.sh` to your `~/.bashrc`

To build the code using the `make` system do:

```bash
make -j
```


## Run COLDREG

If everything compiled correctly, you can run a test example by doing:

```bash
./bin/runcoldreg
```

If you get an error message that indicates that the PETSc library could not be found, you probably forgot to

```bash
source libs/environment_vars.sh
```

To get a general idea on how to run the binary do: 

```bash
./bin/runcoldreg -help
```

For more advanced options do:

```bash
./bin/runcoldreg -advanced
```
