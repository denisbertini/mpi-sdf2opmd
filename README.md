# mpi-sdf2opmd

Converts SDF to openPMD compliant format via openpmd-api interface.

## Overview

**mpi-sdf2opmd** converts [Epoch](https://github.com/Warwick-Plasma/epoch)  `SDF` data format to
[openPMD standard schema](https://github.com/openPMD/openPMD-standard).

It  uses internally the functionality of the [openPMD-api](https://openpmd-api.readthedocs.io/en/0.14.5/) library which provides
a reference API for the [openPMD schema](https://github.com/openPMD/openPMD-standard). 
Furthermore the  openPMD-api library implementing  various I/O backends i.e

- [JSON](https://de.wikipedia.org/wiki/JavaScript_Object_Notation) (only for serial workflow)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [ADIOS1](https://github.com/ornladios/ADIOS)
- [ADIOS2](https://adios2.readthedocs.io/en/latest/)

the converter supports also automatically  all these formats.

**mpi-sdf2opmd** uses parallel I/O for both reading and writing, making convertion of large datasets efficient on 
HPC systems.

## Usage

To convert a 1D/2D or 3D SDF file: 
```
mpirun -n N sdf2opmd_1d/2d/3d -p <dir> -f <sdf_file> -m <fields> -d <derived_fields> -s <particles> -o <output_format>
```
where

- `<dir>` is the directory where the Epoch SDF output files are stored. if only the `-p` option is used all the SDF files from the input directory `<dir>` will be converted. 

- `<sdf_files>` is a list of  input SDF file to be converted from the input directory <dir>. Example: `-f "0012.sdf:0014.sdf" option will tell the converter to only convert files `0012.sdf` and `0014.sdf` from the input directory.

- `<fields>` is the list of relevant fields to be converted. For example to convert all electrical fields components use
the options `-m ex:ey:ez` 


- `<derived_fields>` is the list of relevant derived data fields to be converted. For example to convert charge density and number density calculated for all species use
the options `-d charge_density:number_density` 



- `<particles>` is the list of relevant particles to be converted in the format <species1:species2:species3 ...>.
For example  to convert the electrons and proton species use the option `-s electron:proton`

- `<output_format>` defines the output format. Use the option `-o hdf5` to convert in HSDF5 format or `-o adios` to convert in ADIOS2
format.

To obtain the SDF input file layout:
```
mpirun -n N sdf2opmd_1d/2d/3d -i all -p <dir> -f <sdf_files>
```
This will not perform any conversion but gives the meta-data contents of the input SDF file.
Example output:

```
sdf2opmd_1d: dry run:  reading SDF blocks info: all from file: /lustre/rz/dbertini/sdf2opmd/sim/epoch1d/data/0100.sdf
 Input file contains:           29  blocks
 Block info:            1  : Run_info                                                         : run_info                         :            7
 Block info:            2  : CPUs/Original rank                                               : cpu_rank                         :           20
 Block info:            3  : Wall-time                                                        : elapsed_time                     :            5
 Block info:            4  : Electric Field/Ex                                                : ex                               :            3
 Block info:            5  : Electric Field/Ey                                                : ey                               :            3
 Block info:            6  : Electric Field/Ez                                                : ez                               :            3
 Block info:            7  : Magnetic Field/Bx                                                : bx                               :            3
 Block info:            8  : Magnetic Field/By                                                : by                               :            3
 Block info:            9  : Magnetic Field/Bz                                                : bz                               :            3
 Block info:           10  : Current/Jx                                                       : jx                               :            3
 Block info:           11  : CPU split/electron_r                                             : cpu/electron_r                   :           20
 Block info:           12  : CPU split/electron_l                                             : cpu/electron_l                   :           20
 Block info:           13  : Particles/Px/electron_r                                          : px/electron_r                    :            4
 Block info:           14  : Particles/Px/electron_l                                          : px/electron_l                    :            4
 Block info:           15  : Grid/Particles/electron_r                                        : grid/electron_r                  :            2
 Block info:           16  : Grid/Particles/electron_l                                        : grid/electron_l                  :            2
 Block info:           17  : Derived/Average_Particle_Energy                                  : ekbar                            :            3
 Block info:           18  : Derived/Charge_Density                                           : charge_density                   :            3
 Block info:           19  : Derived/Number_Density                                           : number_density                   :            3
 Block info:           20  : Derived/Number_Density/electron_r                                : number_density/electron_r        :            3
 Block info:           21  : Derived/Number_Density/electron_l                                : number_density/electron_l        :            3
 Block info:           22  : Derived/Temperature                                              : temperature                      :            3
 Block info:           23  : Derived/Temperature/electron_r                                   : temperature/electron_r           :            3
 Block info:           24  : Derived/Temperature/electron_l                                   : temperature/electron_l           :            3
 Block info:           25  : Grid/Grid                                                        : grid                             :            1
 Block info:           26  : Grid/x_px/electron_r                                             : grid/x_px/electron_r             :            1
 Block info:           27  : dist_fn/x_px/electron_r                                          : x_px/electron_r                  :            3
 Block info:           28  : Grid/x_px/electron_l                                             : grid/x_px/electron_l             :            1
 Block info:           29  : dist_fn/x_px/electron_l                                          : x_px/electron_l                  :            3

```
The third column corresponds to the EPOCH `block_id` metadata (Fields, Derived and Particles) that can be given as input to the converter.

**Limitation:** distribution function output blocks i.e  `dist_fn` are not supported for conversion. 


## Working with pp-containers
The converter has been tested with the standard [pp-containers](https://git.gsi.de/d.bertini/pp-containers).

### Compilation

- first launch the container on baremetal submit node `virgo2.hpc.gsi.de`, for example:

```
export CONT=/lustre/rz/dbertini/containers/prod/rlx8_ompi_ucx.sif
singularity exec $CONT bash -l
```

- execute setup script from `mpi-sdf2opmd` top directory 

```
. ./setup_mpi -c cont
```

- compile on a separate build directory
```
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<path_to_install_dir>/sdf2opmd <path_to_converter_source>/mpi-sdf2opmd
make
make install
```

### Running 

To launch a converter job from the baremetal submit node you will need typically the following scripts

- `submit.sh`

```
sbatch --tasks 2 --ntasks-per-core 1 --cpus-per-task 1 --no-requeue --job-name r_mpi --mem-per-cpu 4000 --mail-type ALL --mail-user d.bertini@gsi.de --partition high_mem --time 0-08:00:00 -D ./data -o %j.out.log -e %j.err.log  ./run-file.sh
```


- `run-file.sh`

```
#!/bin/bash
export OMP_NUM_THREADS=1
export APPTAINER_DISABLE_CACHE=1

# Standard RLX8 container
export CONT=/lustre/rz/dbertini/containers/prod/rlx8_ompi_ucx.sif

export APPTAINER_BINDPATH=/lustre/rz/dbertini/,/cvmfs
export OMPI_MCA_io=romio321

# Select files to convert
files="0012.sdf:0014.sdf"

srun --export=ALL  singularity exec $CONT ./convert.sh $files

```

- `convert.sh`

```
#!/bin/bash
export SIMDIR=/lustre/rz/dbertini/sdf2opmd
export PATH=$SIMDIR/bin:$PATH
export LD_LIBRARY_PATH=$SIMDIR/lib:$LD_LIBRARY_PATH

# Define SDF files from args
sdf_files=$1

# Define the derived data that should be included in the conversion
export derived_data=charge_density:number_density

# Convert all SDF files from input directory  to HDF5 format 
$SIMDIR/bin/sdf2opmd_2d -p $SIMDIR/sim/epoch2d/data  -m ex:ey:ez -d $derived_data -s electron_r:electron_l -o hdf5

# Convert SDF to ADIOS2 format including selection
$SIMDIR/bin/sdf2opmd_2d -p $SIMDIR/sim/epoch2d/data -f $sdf_files -m ex:ey:ez -d $derived_data -s electron_r:electron_l -o adios

```

### Converted File Outputs

The converted output file is written in the `sbatch` output directory. In the above example, the converted files will be written in the `data` sub-directory   `sbatch -D ./data ...`  with the extention corresponding to either
- `HDF5: .h5` 
- `ADIOS: .bp`


Example:
```
.
├── convert_mult_file.sh
├── convert.sh
├── data
│   ├── 0000.h5
│   ├── 0001.h5
│   ├── 0002.h5
│   ├── 0003.h5
│   ├── 0004.h5
│   ├── 0005.h5
│   ├── 0006.h5
│   ├── 0007.h5
│   ├── 0008.h5
│   ├── 0009.h5
│   ├── 0100.h5
│   ├── 11860289.err.log
│   ├── 11860289.out.log
│   └── convert.sh -> ../convert.sh
├── run-file.sh
├── setup_mpi.sh
└── submit.sh

```


### Compression via Particle Merging

The converter can on demand perform data compression. For that purpose, the  particle merging algorithm as described in
[M. Vranic et al. (2005)](https://www.sciencedirect.com/science/article/abs/pii/S0010465515000405?via%3Dihub) is used.

The algorithm implementation is similar to the one used in the 
[Smilei pic code](https://smileipic.github.io/Smilei/Understand/particle_merging.html) and uses a cartesian momentum
discretization.

To activate the compression, just add `-z cartesian` as additional option to the converter (2D/3D only):

```
# Convert SDF to ADIOS2 format including selection
$SIMDIR/bin/sdf2opmd_2d -p $SIMDIR/sim/epoch2d/data -f $sdf_files -m ex:ey:ez -d $derived_data -s electron_r:electron_l -z cartesian -o adios

```
