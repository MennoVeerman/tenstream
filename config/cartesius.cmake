# Cartesius supercomputer at SARA Amsterdam
# modules loaded:

# module unload compilerwrappers
# module load hdf5 netcdf cmake python git
# module load blas mkl

# Petsc installed with:
# module unload compilerwrappers
# base_opt="--with-fortran --with-fortran-interfaces --with-shared-libraries=1 --with-blas-lapack-dir=$SURFSARA_MKL_ROOT"
# intel_compilers="--with-cc=$(which mpiicc) --with-fc=$(which mpiifort) --with-cxx=$(which mpiicpc)"
# export PETSC_ARCH='prod_icc'
# ./configure $intel_compilers $base_opt $batch_compile_opts \
# '--with-debugging=0' COPTFLAGS='-mkl' FOPTFLAGS='-mkl' && make

#set(NETCDF_DIR      "$ENV{SURFSARA_NETCDF_ROOT}")
#set(NETCDF_DIR_F90  "$ENV{SURFSARA_NETCDF_ROOT}")
set(NETCDF_DIR      "$ENV{EBROOTNETCDF}")
set(NETCDF_DIR_F90  "$ENV{EBROOTNETCDFMINFORTRAN}")
#set(BLA_PREFER_PKGCONFIG "$ENV{MKLROOT}")

set(USER_C_FLAGS "-nowarn -std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend-source -g -mkl -lsvml")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xCORE-AVX2 -fp-model source -fno-omit-frame-pointer")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")


#set(USER_C_FLAGS "-nowarn -std=c99")
#set(USER_Fortran_FLAGS "-cpp -traceback -extend-source -g -mkl ")
#set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xCORE-AVX2 -fp-model source -fno-omit-frame-pointer")
#set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")
#set(USER_C_FLAGS "-nowarn  -ftz -O3 -xCORE-AVX2 -std=c99 -I/hpc/sw/petsc-3.6.3-intel-impi5-par/avx2/include  -L/sw/arch/RedHatEnterpriseServer7/EB_production/2019/software/ifort/2018.3.222-GCC-7.3.0-2.30/lib/intel64_lin -lsvml")
#set(USER_Fortran_FLAGS "-cpp -ftz -extend-source -mkl -nostandard-realloc-lhs -L/sw/arch/RedHatEnterpriseServer7/EB_production/2019/software/ifort/2018.3.222-GCC-7.3.0-2.30/lib/intel64_lin -lsvml")
#set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xCORE-AVX2 -fp-model fast=2 -fno-omit-frame-pointer")#-I/hpc/sw/petsc-3.6.3-intel-impi5-par/avx2/include")
#set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -traceback -check nopointers -check noarg_temp_created -I/hpc/sw/petsc-3.6.3-intel-impi5-par/avx2/include -L/sw/arch/RedHatEnterpriseServer7/EB_production/2019/software/ifort/2018")

set(CMAKE_C_COMPILER   "mpiicc")
set(CMAKE_Fortran_COMPILER   "mpiifort")
set(Fortran_COMPILER_WRAPPER "mpiifort")
