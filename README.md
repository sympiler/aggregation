# LBC
Load-balance Level Coarsening

## Install

### Linux
First export METIS and MKL library that are only required for the triangular example:
```bash
export MKLROOT <path to MKL>
export METISROOT <path to METIS>
```

Then create a build directory and use camke to build the LBC code:

```bash
mkdir build
cd build
cmake ..
make
```

### Mac
TODO: Setting the C and CXX compilers


