#module load comm_libs/openmpi/nvhpc-5.0.5
export VERNO=25.5
module load nvhpc
module load misc/nvhpc-build/${VERNO}/netcdf-c/
module load misc/nvhpc-build/${VERNO}/netcdf-fortran/
module load hdf5
export PYTHON=python3
