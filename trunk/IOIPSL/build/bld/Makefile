# Automatic Makefile

FCM_BLD_CFG = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/IOIPSL
FCM_BLD_TARGETS = libioipsl.a flio_rbld.exe rebuild

export PERL5LIB := /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/IOIPSL/tools/FCM_V1.2/lib

export FCM_ROOTDIR = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/IOIPSL/build
export FCM_BINDIR = $(FCM_ROOTDIR)/bin
export FCM_BLDDIR = $(FCM_ROOTDIR)/bld
export FCM_CACHEDIR = $(FCM_ROOTDIR)/.cache
export FCM_CFGDIR = $(FCM_ROOTDIR)/cfg
export FCM_DONEDIR = $(FCM_ROOTDIR)/done
export FCM_ETCDIR = $(FCM_ROOTDIR)/etc
export FCM_FLAGSDIR = $(FCM_ROOTDIR)/flags
export FCM_INCDIR = $(FCM_ROOTDIR)/inc
export FCM_LIBDIR = $(FCM_ROOTDIR)/lib
export FCM_OBJDIR = $(FCM_ROOTDIR)/obj
export FCM_PPSRCDIR = $(FCM_ROOTDIR)/ppsrc
export FCM_SRCDIR = $(FCM_ROOTDIR)/src
export FCM_TMPDIR = $(FCM_ROOTDIR)/tmp

export FCM_ROOTPATH = $(FCM_ROOTDIR)
export FCM_BINPATH = $(FCM_BINDIR)
export FCM_BLDPATH = $(FCM_BLDDIR)
export FCM_CACHEPATH = $(FCM_CACHEDIR)
export FCM_CFGPATH = $(FCM_CFGDIR)
export FCM_DONEPATH = $(FCM_DONEDIR)
export FCM_ETCPATH = $(FCM_ETCDIR)
export FCM_FLAGSPATH = $(FCM_FLAGSDIR)
export FCM_INCPATH = $(FCM_INCDIR)
export FCM_LIBPATH = $(FCM_LIBDIR)
export FCM_OBJPATH = $(FCM_OBJDIR)
export FCM_PPSRCPATH = $(FCM_PPSRCDIR)
export FCM_SRCPATH = $(FCM_SRCDIR)
export FCM_TMPPATH = $(FCM_TMPDIR)

export AR = ar
export ARFLAGS = rs
export CC = cc
export CC_COMPILE = -c
export CC_DEFINE = -D
export CC_INCLUDE = -I
export CC_OUTPUT = -o
export CFLAGS = 
export CPPKEYS = 
export FC = gfortran
export FC_COMPILE = -c
export FC_DEFINE = -D
export FC_INCLUDE = -I
export FC_OUTPUT = -o
export FFLAGS = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fno-align-commons -fcray-pointer -O3 -I/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/netcdf/include
export FPPKEYS = CPP_NONE %FPP_DEF
export LD = gfortran
export LDFLAGS =  -L/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/netcdf/lib -lnetcdff
export LD_LIBLINK = -l
export LD_LIBSEARCH = -L
export LD_OUTPUT = -o
export SHELL = /bin/bash

export FCM_VERBOSE ?= 1

vpath % $(FCM_BINPATH)
vpath %.bash $(FCM_BINPATH)
vpath %.csh $(FCM_BINPATH)
vpath %.ksh $(FCM_BINPATH)
vpath %.pl $(FCM_BINPATH)
vpath %.pm $(FCM_BINPATH)
vpath %.pro $(FCM_BINPATH)
vpath %.py $(FCM_BINPATH)
vpath %.sh $(FCM_BINPATH)
vpath %.tcl $(FCM_BINPATH)
vpath %.exe $(FCM_BINPATH)
vpath %.mk $(FCM_BLDPATH)
vpath %.done $(FCM_DONEPATH)
vpath %.etc $(FCM_DONEPATH)
vpath %.idone $(FCM_DONEPATH)
vpath %.flags $(FCM_FLAGSPATH)
vpath %.h $(FCM_INCPATH)
vpath %.h90 $(FCM_INCPATH)
vpath %.inc $(FCM_INCPATH)
vpath %.interface $(FCM_INCPATH)
vpath %.a $(FCM_LIBPATH)
vpath %.o $(FCM_OBJPATH)
vpath %.dummy $(FCM_DONEDIR)

.PHONY : all clean

all : $(FCM_BLD_TARGETS)

clean : 
	rm -rf $(FCM_BINDIR) $(FCM_LIBDIR) $(FCM_OBJDIR) $(FCM_DONEDIR)
	rm -f lib__*.a *.o

$(FCM_DONEDIR)/FCM_CP.dummy :
	touch $@

OBJECTS__src = calendar.o defprec.o errioipsl.o flincom.o fliocom.o getincom.o histcom.o ioipsl.o mathelp.o restcom.o stringop.o

libsrc.a : $(OBJECTS__src)
	fcm_internal archive $@ $(^F)

OBJECTS__rebuild =

librebuild.a : $(OBJECTS__rebuild)
	fcm_internal archive $@ $(^F)

export OBJECTS = $(OBJECTS__rebuild) $(OBJECTS__src)

libioipsl.a : $(OBJECTS)
	fcm_internal archive $@ $(^F)

CC.flags :
	touch $(FCM_FLAGSDIR)/$@

CFLAGS.flags : CC.flags
	touch $(FCM_FLAGSDIR)/$@

CFLAGS__rebuild.flags CFLAGS__src.flags : CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__rebuild.flags CPPKEYS__src.flags : CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FC.flags :
	touch $(FCM_FLAGSDIR)/$@

FFLAGS.flags : FC.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__rebuild.flags FFLAGS__src.flags : FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__rebuild.flags FPPKEYS__src.flags : FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LD.flags :
	touch $(FCM_FLAGSDIR)/$@

LD__rebuild.flags LD__src.flags : LD.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS.flags :
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__rebuild.flags LDFLAGS__src.flags : LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

include $(FCM_BLDDIR)/rebuild.mk
include $(FCM_BLDDIR)/src.mk

# EOF
