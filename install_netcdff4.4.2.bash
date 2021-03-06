#!/bin/bash 
#########################################################################
# Script to automatically install the NetCDF 4.4.2 Fortran Library
#########################################################################
# Defaults
#########################################################################
install_dir=$(pwd)/netcdf-fortran-4.4.2
compiler_suite=gnu
#########################################################################
#  Options 
#########################################################################
while (($# > 0))
   do
   case $1 in
     "-h") cat <<........fin
    $0 [ -prefix path ]       where (path) to install
                              (default: $install_dir)
       [ -compiler gnu | intel ] compiler suite (Fortran, C, C++) to use
                              (default: $compiler_suite)
........fin
     exit ;;
     "-prefix") install_dir=$2 ; shift ; shift ;;
     "-compiler") compiler_suite=$2 ; shift ; shift ;;
     *) echo "Error, bad argument $1" ; $0 -h ; exit
   esac
done

# Install directory (get full path)
mkdir -p $install_dir
install_dir=$(cd $install_dir ; pwd -P )
cd $install_dir

rm -rf netcdf-fortran-4.4.2*
wget -nv http://www.lmd.jussieu.fr/~lmdz/Distrib/netcdf-fortran-4.4.2.tar.gz
tar xzf netcdf-fortran-4.4.2.tar.gz ; cd netcdf-fortran-4.4.2 

if [[ ${compiler_suite} == "gnu" ]] ; then
  f_compiler="gfortran"
  c_compiler="gcc"
  cxx_compiler="g++"
elif [[ ${compiler_suite} == "intel" ]] ; then
  f_compiler="ifort"
  c_compiler="icc"
  cxx_compiler="icpc"
else
  echo "unknown compiler family $compiler_suite"
  echo "might as well stop here"
  exit
fi

export FC=$f_compiler
export F90=$f_compiler
export CC=$c_compiler
export CXX=$cxx_compiler
if [[ ${f_compiler} == "gfortran" ]] ; then
  export FFLAGS=" -O2 -fPIC"
  export FCFLAGS="-O2 -ffree-form -fPIC"
  export CPPFLAGS=""
  export CFLAGS="-O2 -fPIC"
  export CXXFLAGS="-O2 -fPIC"
elif [[ ${f_compiler} == "ifort" ]] ; then
  export CPP="icc -E"
  export FFLAGS="-O2 -ip -fpic"
  export FCFLAGS="-O2 -ip -fpic"
  export CPPFLAGS=""
  export CFLAGS="-O2 -ip -fpic"
  export CXXFLAGS="-O2 -ip -fpic"
else
  echo "unknown compiler $f_compiler"
  echo "might as well stop here"
  exit
fi

./configure --enable-shared --prefix=$install_dir > configure.log 2>&1 

make > make.log 2>&1

make test > make_test.log 2>&1

make install > make_install.log 2>&1

if [[ -f $install_dir/bin/nc-config ]] ; then
  echo "successfully installed the netcdf library in $install_dir"
fi
