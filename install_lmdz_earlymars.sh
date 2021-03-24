#!/bin/bash

# stop as soon as a comand returns an error:
set -e

## Default values for various script options
installdir="trunk"
svn_planeto="HEAD"
arch="local"
compiler_suite="gnu"

machine=`hostname`
if [[ ${machine:7:10} == "occigen" ]] ; then
  arch=X64_OCCIGEN
fi
if [[ ${machine:0:5} == "ada33" ]] ; then
  arch=X64_ADA
fi

# 0.1 Check that prerequisits are there
for logiciel in svn wget tar gzip make ; do
if [ "`which $logiciel`" = "" ] ; then
echo You must first install $logiciel on your system
exit
fi
done

if [ "`which fcm`" = "" ] ; then
echo "You must first install fcm and add it to your PATH, e.g.:"
echo "cd ; svn co http://forge.ipsl.jussieu.fr/fcm/svn/PATCHED/FCM_V1.2"
echo "and add FCM_V1.2/bin to your PATH environment variable,"
echo "e.g. add in your .bashrc the following line:"
echo 'export PATH=$PATH:$HOME/FCM_V1.2/bin'
exit
fi

# 0.2 Get user options from command line
while (($# > 0))
   do
   case $1 in
    "-h") cat <<........end
    Script usage:
    $0 [option1 value1 option2 value2 ...]
    where available options are:
    -installdir name of subdirectory to install in (default: $installdir)
    -planeto_rev revision number of planeto repository (default: $svn_planeto)
    -compiler compiler family (default: $compiler_suite ; gnu | intel )
........end
     exit ;;
    "-installdir") installdir=$2 ; shift ; shift ;;
    "-planeto_rev") svn_planeto==$2 ; shift ; shift ;;
    "-compiler") compiler_suite=$2 ; shift ; shift ;;
    *) $0 -h ; exit
   esac
done

# check compiler is available
if [[ ${compiler_suite} == "gnu" ]] ; then
  f_compiler="gfortran"
  c_compiler="gcc"
elif [[ ${compiler_suite} == "intel" ]] ; then
  f_compiler="ifort"
  c_compiler="icc"
else
  echo "unknown compiler family $compiler_suite"
  echo "might as well stop here"
  exit
fi
for logiciel in $f_compiler $c_compiler ; do
if [ "`which $logiciel`" = "" ] ; then
echo You must first install $logiciel on your system
exit
fi
done

rootdir=`pwd -P`

# 1. Check out "Planeto" repository
cd $rootdir
svn co --revision $svn_planeto http://svn.lmd.jussieu.fr/Planeto/trunk --depth empty  $installdir
cd $installdir
svn update LMDZ.COMMON LMDZ.GENERIC DOC

# 1.2 Download bench test case
cd $rootdir
wget -nv http://www.lmd.jussieu.fr/~lmdz/planets/generic/bench_earlymars_32x32x15_b32x36.tar.gz

# 1.3. Check out and build the NetCDF library
cd $rootdir
wget -nv http://www.lmd.jussieu.fr/~lmdz/pub/import/install_netcdff4.4.2.bash
chmod u=rwx install_netcdff4.4.2.bash
netcdflog=`pwd`/netcdf.log
./install_netcdff4.4.2.bash -compiler $compiler_suite > $netcdflog 2>&1
export NETCDF=$rootdir/netcdf-fortran-4.4.2
export NETCDFINCLUDE=$NETCDF/include
export NETCDFDIR=$NETCDF/lib

# 2.1 Compile IOIPSL
cd $rootdir/$installdir/LMDZ.COMMON/ioipsl
if  [[ ${compiler_suite} == "gnu" ]] ; then
sed -i -e s:'NETCDF_HOME=$(nf-config --prefix)':"NETCDF_HOME=$NETCDF":1 install_ioipsl_gfortran.bash
./install_ioipsl_gfortran.bash > install_ioipsl_gfortran.out 2>&1
else
sed -i -e s:'NETCDF_HOME=$(nf-config --prefix)':"NETCDF_HOME=$NETCDF":1 install_ioipsl_ifort.bash
./install_ioipsl_ifort.bash > install_ioipsl_ifort.out 2>&1
fi

# 2.2 Compile LMDZ
cd $rootdir/$installdir/LMDZ.COMMON
# 2.2.1 Make arch files
cd arch
echo 'ROOT=$PWD' > arch-local.path
echo "" >> arch-local.path
echo 'NETCDF_LIBDIR="-L'$NETCDF'/lib"' >> arch-local.path
echo 'NETCDF_LIB="-lnetcdf -lnetcdff"' >> arch-local.path
echo 'NETCDF_INCDIR="-I'$NETCDF'/include"' >> arch-local.path
echo "" >> arch-local.path
echo 'IOIPSL_INCDIR="-I$ROOT/ioipsl/modipsl/lib"' >> arch-local.path
echo 'IOIPSL_LIBDIR="-L$ROOT/ioipsl/modipsl/lib"' >> arch-local.path
echo 'IOIPSL_LIB="-lioipsl"' >> arch-local.path
echo "" >> arch-local.path
echo 'XIOS_INCDIR="-I$ROOT/../XIOS/inc"' >> arch-local.path
echo 'XIOS_LIBDIR="-L$ROOT/../XIOS/lib"' >> arch-local.path
echo 'XIOS_LIB="-lxios -lstdc++"' >> arch-local.path

if  [[ ${compiler_suite} == "gnu" ]] ; then
cp arch-gfortran.fcm arch-local.fcm
else
cp arch-linux-ifort.fcm arch-local.fcm
fi

echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'${NETCDF}'/lib' > arch-local.env

# 2.2.2 Compile GCM
cd $rootdir/$installdir/LMDZ.COMMON
./makelmdz_fcm -arch $arch -p std -s 2 -d 32x32x15 -b 32x36 -j 8 gcm > makelmdz_fcm.out 2>&1

if [ -f bin/gcm_32x32x15_phystd_seq.e ] ; then
  echo "success! Executable is in $rootdir/$installdir/LMDZ.COMMON/bin"
else
  echo "something went wrong"
  echo "check file $rootdir/$installdir/LMDZ.COMMON/makelmdz_fcm.out"
fi

# 3. Run Bench
cd $rootdir
tar xvzf bench_earlymars_32x32x15_b32x36.tar.gz
cp $installdir/LMDZ.COMMON/bin/gcm_32x32x15_phystd_seq.e bench_earlymars_32x32x15_b32x36/
cd bench_earlymars_32x32x15_b32x36/

echo 'source '$rootdir'/'$installdir'/LMDZ.COMMON/arch.env' > run_gcm.sh
echo './gcm_32x32x15_phystd_seq.e > gcm.out 2>&1' >> run_gcm.sh
chmod u=rwx run_gcm.sh
./run_gcm.sh
echo "Bench testcase run in $rootdir/bench_earlymars_32x32x15_b32x36/ using run_gcm.sh"
echo "Outputs are in gcm.out which ends with:"
tail gcm.out
