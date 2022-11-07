# mpi-sdf2opmd
SDF to openPMD compliant format via openpmd-api interface.

## Overview
`mpi-sdf2opmd` converts [Epoch](https://github.com/Warwick-Plasma/epoch)  `SDF` data format to
[openPMD standard schema](https://github.com/openPMD/openPMD-standard).
It  uses internally the functionality of the [openPMD-api](https://openpmd-api.readthedocs.io/en/0.14.5/) library which provides
a reference API for the [openPMD schema](https://github.com/openPMD/openPMD-standard). 
Furthermore the  openPMD-api library implementing  
various I/O backends i.e

- [JSON](https://de.wikipedia.org/wiki/JavaScript_Object_Notation)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [ADIOS1](https://github.com/ornladios/ADIOS)
- [ADIOS2](https://adios2.readthedocs.io/en/latest/)

the converter will also automatically supports all this formats.



