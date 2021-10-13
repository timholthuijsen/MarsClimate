#!/bin/bash
# script to download and install the latest version of IOIPSL on Occigen
#

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

# cleanup possible previous attempt:
\rm -rf ../../IOIPSL modipsl 

# 1. Get IOIPSL 
# move up at same level as LMDZ.COMMON , etc.
cd ../..
svn co --revision $rev http://forge.ipsl.jussieu.fr/igcmg/svn/IOIPSL/trunk IOIPSL

# 2. Set correct settings: copy over arch.* files
cd IOIPSL
cp -f ../LMDZ.COMMON/arch/arch-X64_OCCIGEN.env arch
cp -f ../LMDZ.COMMON/arch/arch-X64_OCCIGEN.path arch
cp -f ../LMDZ.COMMON/arch/arch-X64_OCCIGEN.fcm arch

## 3. build ioipsl:
# but first make a small correction to makeioipsl_fcm
sed -i -e s:'$HDF5_LIBDIR $HDF5_LIB':'':1 makeioipsl_fcm
./makeioipsl_fcm -arch X64_OCCIGEN -job 8 > makeioipsl.out 2>&1

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
