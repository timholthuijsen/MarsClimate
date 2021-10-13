#!/bin/bash
# script to download and install the latest version of IOIPSL on Jean Zay
#

#0. Preliminary stuff 
module load subversion

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
cp -f ../LMDZ.COMMON/arch/arch-X64_JEANZAY-pgi.env arch
cp -f ../LMDZ.COMMON/arch/arch-X64_JEANZAY-pgi.path arch
cp -f ../LMDZ.COMMON/arch/arch-X64_JEANZAY-pgi.fcm arch

## 3. build ioipsl:
./makeioipsl_fcm -arch X64_JEANZAY-pgi -job 8 > makeioipsl.out 2>&1

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

