#!/bin/bash
# script to download and install the latest version of IOIPSL
# using ifort
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
NETCDF_HOME=$(nf-config --prefix)

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
echo "" > arch-ifort.env
echo "export NETCDF_HOME=$NETCDF_HOME" >> arch-ifort.env
# arch.fcm file 
echo '%COMPILER            ifort' > arch-ifort.fcm
echo '%LINK                ifort' >> arch-ifort.fcm
echo '%AR                  ar' >> arch-ifort.fcm
echo '%MAKE                make' >> arch-ifort.fcm
echo '%FPP_FLAGS           -P -traditional' >> arch-ifort.fcm
echo '%BASE_FFLAGS         -real-size 64 -ip -mkl=parallel -fp-model precise -align common' >> arch-ifort.fcm
echo '%PROD_FFLAGS         -O3' >> arch-ifort.fcm
echo '%DEV_FFLAGS          -O2' >> arch-ifort.fcm
echo '%DEBUG_FFLAGS        -fpe0 -g -no-ftz -check -traceback -ftrapuv -fp-stack-check -O0' >> arch-ifort.fcm
echo '%MPI_FFLAGS          ' >> arch-ifort.fcm
echo '%OMP_FFLAGS          ' >> arch-ifort.fcm
echo '%BASE_LD             ' >> arch-ifort.fcm
echo '%MPI_LD              ' >> arch-ifort.fcm
echo '%OMP_LD              ' >> arch-ifort.fcm
# arch.path file
echo "NETCDF_INCDIR=\"-I$NETCDF_HOME/include\"" > arch-ifort.path
echo "NETCDF_LIBDIR=\"-L$NETCDF_HOME/lib\"" >> arch-ifort.path
echo 'NETCDF_LIB="-lnetcdff"' >> arch-ifort.path
echo '' >> arch-ifort.path
echo 'HDF5_INCDIR=""' >> arch-ifort.path
echo 'HDF5_LIBDIR=""' >> arch-ifort.path
echo 'HDF5_LIB=""' >> arch-ifort.path
echo '' >> arch-ifort.path
echo 'MPI_INCDIR=""' >> arch-ifort.path
echo 'MPI_LIBDIR=""' >> arch-ifort.path
echo 'MPI_LIB=""' >> arch-ifort.path

## 3. build ioipsl:
cd ..
./makeioipsl_fcm -arch ifort -job 8 > makeioipsl.out 2>&1

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

