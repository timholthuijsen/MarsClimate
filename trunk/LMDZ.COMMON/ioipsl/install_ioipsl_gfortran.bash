#!/bin/bash
# script to download and install the latest version of IOIPSL
# using gfortran
# You'll probably have to change path to NetCDF library
# below to adapt this script to your computer.

#0. Preliminary stuff 
if (( $# == 0 ))
then
  # default behavior: get latest version of IOIPSL
  rev="HEAD"
else
  # but otherwise the first argument of the script can be the version to use
  if (( $# == 1 ))
  then
    rev=$1
  else
    echo "Error, invalid script arguments"
    echo "Usage:"
    echo "$0 rev"
    echo " where optional rev is the IOIPSL revision number"
    exit
  fi
fi

# Where is the NetCDF library root located? Hopefully nf-config can tell us
# but you might need an appropriate "module load netcdf***" beforehand
NETCDF_HOME=/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/netcdf

# cleanup possible previous attempt:
\rm -rf ../../IOIPSL modipsl 

# 1. Get IOIPSL 
# move up at same level as LMDZ.COMMON , etc.
cd ../..
svn co --revision $rev http://forge.ipsl.jussieu.fr/igcmg/svn/IOIPSL/trunk IOIPSL

# 2. Set correct settings: make some arch.* files
# Ideally these arch files should be the same as the GCM's
# arch.env file (add any suitable module load here)
cd IOIPSL/arch
echo "" > arch-gfortran.env
echo "export NETCDF_HOME=$NETCDF_HOME" >> arch-gfortran.env
# arch.fcm file 
echo '%COMPILER            gfortran' > arch-gfortran.fcm
echo '%LINK                gfortran' >> arch-gfortran.fcm
echo '%AR                  ar' >> arch-gfortran.fcm
echo '%MAKE                make' >> arch-gfortran.fcm
echo '%FPP_FLAGS           -P -traditional' >> arch-gfortran.fcm
echo '%BASE_FFLAGS         -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fno-align-commons -fcray-pointer' >> arch-gfortran.fcm
echo '%PROD_FFLAGS         -O3' >> arch-gfortran.fcm
echo '%DEV_FFLAGS          -O -Wall -fbounds-check' >> arch-gfortran.fcm
echo '%DEBUG_FFLAGS        -ffpe-trap=invalid,zero,overflow -fbounds-check -g3 -O0 -fstack-protector-all -finit-real=snan -fbacktrace' >> arch-gfortran.fcm
echo '%MPI_FFLAGS          ' >> arch-gfortran.fcm
echo '%OMP_FFLAGS          ' >> arch-gfortran.fcm
echo '%BASE_LD             ' >> arch-gfortran.fcm
echo '%MPI_LD              ' >> arch-gfortran.fcm
echo '%OMP_LD              ' >> arch-gfortran.fcm
# arch.path file
echo "NETCDF_INCDIR=\"-I$NETCDF_HOME/include\"" > arch-gfortran.path
echo "NETCDF_LIBDIR=\"-L$NETCDF_HOME/lib\"" >> arch-gfortran.path
echo 'NETCDF_LIB="-lnetcdff"' >> arch-gfortran.path
echo '' >> arch-gfortran.path
echo 'HDF5_INCDIR=""' >> arch-gfortran.path
echo 'HDF5_LIBDIR=""' >> arch-gfortran.path
echo 'HDF5_LIB=""' >> arch-gfortran.path
echo '' >> arch-gfortran.path
echo 'MPI_INCDIR=""' >> arch-gfortran.path
echo 'MPI_LIBDIR=""' >> arch-gfortran.path
echo 'MPI_LIB=""' >> arch-gfortran.path

## 3. build ioipsl:
cd ..
./makeioipsl_fcm -arch gfortran -job 8 > makeioipsl.out 2>&1

## 4. Check if the library was indeed built:
whereami=`pwd -P`
if [[ -f lib/libioipsl.a ]] 
  then
  echo "OK: ioipsl library is in ${whereami}/lib"
else
  echo "Something went wrong... check messages in ${whereami}/makeioipsl.out"
  exit
fi

## 5. Comply with old setup and make appropriate links
cd ../LMDZ.COMMON/ioipsl
mkdir modipsl 
cd modipsl
# lib + module files
mkdir lib
cd lib
ln -s ../../../../IOIPSL/lib/libioipsl.a .
ln -s ../../../../IOIPSL/inc/* .
cd ..
# rebuild utility
mkdir bin
cd bin
ln -s ../../../../IOIPSL/bin/* .

