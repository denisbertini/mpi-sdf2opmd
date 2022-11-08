# mpi-sdf2opmd

Converts SDF to openPMD compliant format via openpmd-api interface.

## Overview

**mpi-sdf2opmd** converts [Epoch](https://github.com/Warwick-Plasma/epoch)  `SDF` data format to
[openPMD standard schema](https://github.com/openPMD/openPMD-standard).

It  uses internally the functionality of the [openPMD-api](https://openpmd-api.readthedocs.io/en/0.14.5/) library which provides
a reference API for the [openPMD schema](https://github.com/openPMD/openPMD-standard). 
Furthermore the  openPMD-api library implementing  various I/O backends i.e

- [JSON](https://de.wikipedia.org/wiki/JavaScript_Object_Notation) (
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [ADIOS1](https://github.com/ornladios/ADIOS)
- [ADIOS2](https://adios2.readthedocs.io/en/latest/)

the converter supports also automatically  all these formats.



## Usage

To convert a 1D/2D or 3D SDF file: 
```
mpirun -n N sdf2opmd_1d/2d/3d -f <sdf_file> -m <fields> -s <particles> -o <output_format>

```
where

- `<sdf_file>` is the full name of the input SDF file including full path

- `<fields>` is the list of relevant fields to be converted. For example to convert all electrical fields components use
the options `-m ex:ey:ez` 

- `<particles>` is the list of relevant particles to be converted in the format <species1:species2:species3 ...>.
For example  to convert the electrons and proton species use the option `-s electron:proton`

- `<output_format>` defines the output format. Use the option `-o hdf5` to convert in HSDF5 format or `-o adios` to convert in ADIOS2
format.

