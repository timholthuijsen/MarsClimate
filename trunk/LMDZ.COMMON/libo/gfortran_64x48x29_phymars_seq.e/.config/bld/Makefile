# Automatic Makefile

FCM_BLD_CFG = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON
FCM_BLD_TARGETS = liblmdz.a gcm_64x48x29_phymars_seq.e

export PERL5LIB := /mnt/d/OneDrive/MarsClimate/MetJacky/FCM_V1.2/lib

export FCM_ROOTDIR = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libo/gfortran_64x48x29_phymars_seq.e/.config
export FCM_BINDIR = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/bin
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
export CC = gcc
export CC_COMPILE = -c
export CC_DEFINE = -D
export CC_INCLUDE = -I
export CC_OUTPUT = -o
export CFLAGS = -O3 -I/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/ioipsl/modipsl/lib -I/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libf/aeronomars -I/include
export CPPKEYS = CPP_PHYS CPP_IOIPSL LAPACK
export FC = gfortran
export FC_COMPILE = -c
export FC_DEFINE = -D
export FC_INCLUDE = -I
export FC_OUTPUT = -o
export FFLAGS = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fno-align-commons -O3  -I/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/ioipsl/modipsl/lib -I/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libf/aeronomars -I/include
export FFLAGS__muphy = -c -I/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/ioipsl/modipsl/lib -I/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libf/aeronomars -I/include -O3 
export FPPKEYS = CPP_PHYS CPP_IOIPSL NC_DOUBLE
export LD = gfortran
export LDFLAGS =   -L/mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/ioipsl/modipsl/lib -lioipsl -L/lib -lnetcdf -lnetcdff
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
vpath %_64x48x29_phymars_seq.e $(FCM_BINPATH)
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

OBJECTS__strataer =

libstrataer.a : $(OBJECTS__strataer)
	fcm_internal archive $@ $(^F)

OBJECTS__sisvat =

libsisvat.a : $(OBJECTS__sisvat)
	fcm_internal archive $@ $(^F)

OBJECTS__rrtm =

librrtm.a : $(OBJECTS__rrtm)
	fcm_internal archive $@ $(^F)

OBJECTS__phys = aerave.o aeropacity_mod.o aeroptproperties.o albedocaps.o blackl.o blendrad.o calcstormfract_mod.o call_dayperi.o calldrag_noro_mod.o callradite_mod.o callsedim_mod.o calltherm_interface.o check_fields_mod.o co2cloud_mod.o co2condens_mod.o co2condens_mod4micro.o co2sat.o co2snow.o comcstfi_h.o comgeomfi_h.o comm_wrf.o compute_dtau_mod.o comsaison_h.o comsoil_h.o comtherm_h.o conc_mod.o conf_phys.o convadj.o cvmgp.o cvmgt.o datafile_mod.o def_var.o density_co2_ice_mod.o dimphy.o dimradmars_mod.o drag_noro_mod.o dust_param_mod.o dust_rad_adjust_mod.o dust_scaling_mod.o dustdevil.o dustlift.o eofdump_mod.o flusv.o geticecover.o getslopes.o growthrate.o gwprofil_mod.o gwstress_mod.o hdo_surfex_mod.o improvedclouds_mod.o improvedco2clouds_mod.o iniorbit.o inistats.o initracer.o iniwrite.o iniwritesoil.o interp_line.o iostart.o lwb.o lwdiff.o lwflux.o lwi.o lwmain_mod.o lwtt.o lwu.o lwxb.o lwxd.o lwxn.o massflowrateco2.o mkstats.o mucorr.o mufract.o newsedim_mod.o nir_leedat.o nirco2abs.o planckdp.o mzesc110.o nlte_setup.o nlte_tcool.o nltecool.o nlthermeq.o nonoro_gwd_ran_mod.o nuclea.o nucleaco2.o orbite.o orodrag_mod.o orosetup.o param_slope.o pbl_parameters.o phyetat0_mod.o phyredem.o phys_state_var_init_mod.o physiq_mod.o planete_h.o planetwide_mod.o read_dust_scenario.o rocketduststorm_mod.o scopyi.o sig.o simpleclouds.o slope_mod.o soil.o soil_settings.o soil_tifeedback.o solang.o solarlong.o statto_mod.o suaer.o sugwd.o surfdat_h.o surfini.o swmain_mod.o swr_fouquart.o swr_toon.o swrayleigh.o tabfi.o tcondco2.o tcondwater.o thermcell_dqup.o thermcell_main_mars.o time_phylmdz_mod.o topmons_mod.o tracer_mod.o turb_mod.o updaterad.o updatereffrad_mod.o vdif_cd.o vdif_kc.o vdifc_mod.o vlz_fi.o watercloud_mod.o watersat_mod.o writediagfi.o writediagmicrofi.o writediagsoil.o writeg1d.o wstats.o xios_output_mod.o yamada4.o yomlw_h.o

libphys.a : $(OBJECTS__phys)
	fcm_internal archive $@ $(^F)

OBJECTS__phy_common = abort_physic.o geometry_mod.o init_print_control_mod.o ioipsl_getin_p_mod.o mod_grid_phy_lmdz.o mod_phys_lmdz_mpi_data.o mod_phys_lmdz_mpi_transfert.o mod_phys_lmdz_omp_data.o mod_phys_lmdz_omp_transfert.o mod_phys_lmdz_para.o mod_phys_lmdz_transfert_para.o physics_distribution_mod.o print_control_mod.o regular_lonlat_mod.o vertical_layers_mod.o write_field_phy.o xios_writefield_mod.o

libphy_common.a : $(OBJECTS__phy_common)
	fcm_internal archive $@ $(^F)

OBJECTS__muphy =

libmuphy.a : $(OBJECTS__muphy)
	fcm_internal archive $@ $(^F)

OBJECTS__misc = arth_m.o assert_eq_m.o assert_m.o cbrt.o chfev.o coefpoly_m.o scopy.o fdump.o formcoord.o handle_err_m.o i1mach.o interpolation.o ioipsl_errioipsl.o ioipsl_getincom.o ioipsl_stringop.o ismax.o ismin.o j4save.o juldate.o minmax.o minmax2.o netcdf95.o new_unit_m.o nf95_def_var_m.o nf95_get_att_m.o nf95_get_var_m.o nf95_gw_var_m.o nf95_put_att_m.o nf95_put_var_m.o nrtype.o pchdf.o pchfe.o pchfe_95_m.o pchsp.o pchsp_95_m.o pres2lev_mod.o q_sat.o ran1.o regr1_conserv_m.o regr1_lint_m.o regr3_lint_m.o simple.o slopes_m.o sort.o sort_mod.o vampir.o write_field.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o xgetua.o

libmisc.a : $(OBJECTS__misc)
	fcm_internal archive $@ $(^F)

OBJECTS__grid =

libgrid.a : $(OBJECTS__grid)
	fcm_internal archive $@ $(^F)

OBJECTS__filtrez = acc.o eigen.o eigen_sort.o filtreg.o filtreg_mod.o inifgn.o jacobi.o mkl_dft_type.o mkl_dfti.o mod_fft.o mod_fft_fftw.o mod_fft_mathkeisan.o mod_fft_mkl.o mod_fft_wrapper.o mod_filtre_fft.o mod_filtre_fft_loc.o timer_filtre.o

libfiltrez.a : $(OBJECTS__filtrez)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn_phys_sub = avg_horiz_mod.o caldyn0.o callphysiq_mod.o datareadnc.o defrun_new.o grid_noro1.o ini_archive.o iniphysiq_mod.o interp_vert.o lect_start_archive.o mvc_horiz_mod.o readhead_nc.o scal_wind.o wind_scal.o write_archive.o

libdyn_phys_sub.a : $(OBJECTS__dyn_phys_sub)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn_phys = calfis.o calfis_p.o gr_dyn_fi.o gr_dyn_fi_p.o gr_fi_dyn.o gr_fi_dyn_p.o grid_atob_m.o inigeomphy_mod.o mod_interface_dyn_phys.o

libdyn_phys.a : $(OBJECTS__dyn_phys)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn3d_common = adaptdt.o advn.o advx.o advxp.o advy.o advyp.o advz.o advzp.o bernoui.o com_io_dyn_mod.o comconst_mod.o comuforc_h.o comvert_mod.o conf_dat_m.o conf_planete.o control_mod.o convflu.o convmas.o coordij.o covcont.o covnat.o cpdet_mod.o defrun.o diagedyn.o disvert.o disvert_noterre.o diverg.o diverg_gam.o divergf.o divergst.o divgrad.o divgrad2.o dump2d.o dynetat0.o dynredem0.o dynredem_mod.o ener_mod.o enercin.o exner_hyb_m.o exner_milieu_m.o extrapol.o flumass.o fxhyp_m.o fxy.o fxysinus.o fyhyp_m.o geopot.o gr_ecrit_fi.o gr_int_dyn.o gr_u_scal.o gr_v_scal.o grad.o gradiv.o gradiv2.o grilles_gcm_netcdf_sub.o heavyside.o infotrac.o iniacademic.o iniconst.o inidissip.o inigeom.o inigrads.o initdynav.o initfluxsto.o inithist.o initial0.o inter_barxy_m.o interpost.o interpre.o invert_lat.o invert_zoom_x_m.o iso_verif_nonan_nostop.o laplacien.o laplacien_gam.o laplacien_rot.o laplacien_rotgam.o limx.o limy.o limz.o massbar.o massbarxy.o massdair.o misc_mod.o nxgrad.o nxgrad_gam.o nxgradst.o nxgraro2.o nxgrarot.o pbar.o pentes_ini.o ppm3d.o prather.o pression.o principal_cshift_m.o psextbar.o relaxetat0.o rotat.o rotat_nfil.o rotatf.o rotatst.o serre_mod.o sortvarc.o sortvarc0.o test_period.o tidal_forces.o tourpot.o traceurpole.o ugeostr.o vitvert.o writedynav.o writehist.o

libdyn3d_common.a : $(OBJECTS__dyn3d_common)
	fcm_internal archive $@ $(^F)

OBJECTS__dyn = abort_gcm.o addfi.o advect.o advtrac.o bilan_dyn.o caladvtrac.o caldyn.o check_isotopes_seq.o conf_gcm.o dissip.o dteta1.o dudv1.o dudv2.o fluxstokenc.o friction.o getparam.o groupe.o groupeun.o guide_mod.o iniinterp_horiz.o integrd.o interp_horiz.o leapfrog.o leapfrog_nogcm.o logic_mod.o mod_const_mpi.o moyzon_mod.o planetary_operations.o qminimum.o sponge_mod.o sw_case_williamson91_6.o temps_mod.o tetaleveli1j.o tetaleveli1j1.o top_bound.o vlsplt.o vlspltqs.o wrgrads.o

libdyn.a : $(OBJECTS__dyn)
	fcm_internal archive $@ $(^F)

OBJECTS__dust =

libdust.a : $(OBJECTS__dust)
	fcm_internal archive $@ $(^F)

OBJECTS__cosp =

libcosp.a : $(OBJECTS__cosp)
	fcm_internal archive $@ $(^F)

OBJECTS__cloud =

libcloud.a : $(OBJECTS__cloud)
	fcm_internal archive $@ $(^F)

OBJECTS__chem =

libchem.a : $(OBJECTS__chem)
	fcm_internal archive $@ $(^F)

OBJECTS__aerono = calchim_mod.o chemthermos.o chemthermos_readini.o concentrations.o conduction.o deposition.o dtridgl.o euvheat.o hrtherm.o inichim_newstart.o intrplf.o inv.o iono_h.o jthermcalc.o jthermcalc_e107.o moldiff.o moldiff_red.o moldiffcoeff.o moldiffcoeff_red.o molvis.o param_read.o param_read_e107.o param_v4_h.o paramfoto_compact.o perosat.o photochemistry.o photolysis.o photolysis_mod.o photolysis_online.o read_phototable.o surfacearea.o thermosphere.o types_asis.o

libaerono.a : $(OBJECTS__aerono)
	fcm_internal archive $@ $(^F)

export OBJECTS = $(OBJECTS__aerono) $(OBJECTS__chem) $(OBJECTS__cloud) $(OBJECTS__cosp) $(OBJECTS__dust) $(OBJECTS__dyn) $(OBJECTS__dyn3d_common) $(OBJECTS__dyn_phys) $(OBJECTS__dyn_phys_sub) $(OBJECTS__filtrez) $(OBJECTS__grid) $(OBJECTS__misc) $(OBJECTS__muphy) $(OBJECTS__phy_common) $(OBJECTS__phys) $(OBJECTS__rrtm) $(OBJECTS__sisvat) $(OBJECTS__strataer)

liblmdz.a : $(OBJECTS)
	fcm_internal archive $@ $(^F)

CC.flags :
	touch $(FCM_FLAGSDIR)/$@

CFLAGS.flags : CC.flags
	touch $(FCM_FLAGSDIR)/$@

CFLAGS__aerono.flags CFLAGS__chem.flags CFLAGS__cloud.flags CFLAGS__cosp.flags CFLAGS__dust.flags CFLAGS__dyn.flags CFLAGS__dyn3d_common.flags CFLAGS__dyn_phys.flags CFLAGS__dyn_phys_sub.flags CFLAGS__filtrez.flags CFLAGS__grid.flags CFLAGS__misc.flags CFLAGS__muphy.flags CFLAGS__phy_common.flags CFLAGS__phys.flags CFLAGS__rrtm.flags CFLAGS__sisvat.flags CFLAGS__strataer.flags : CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__aerono.flags CPPKEYS__chem.flags CPPKEYS__cloud.flags CPPKEYS__cosp.flags CPPKEYS__dust.flags CPPKEYS__dyn.flags CPPKEYS__dyn3d_common.flags CPPKEYS__dyn_phys.flags CPPKEYS__dyn_phys_sub.flags CPPKEYS__filtrez.flags CPPKEYS__grid.flags CPPKEYS__misc.flags CPPKEYS__muphy.flags CPPKEYS__phy_common.flags CPPKEYS__phys.flags CPPKEYS__rrtm.flags CPPKEYS__sisvat.flags CPPKEYS__strataer.flags : CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FC.flags :
	touch $(FCM_FLAGSDIR)/$@

FFLAGS.flags : FC.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__aerono.flags FFLAGS__chem.flags FFLAGS__cloud.flags FFLAGS__cosp.flags FFLAGS__dust.flags FFLAGS__dyn.flags FFLAGS__dyn3d_common.flags FFLAGS__dyn_phys.flags FFLAGS__dyn_phys_sub.flags FFLAGS__filtrez.flags FFLAGS__grid.flags FFLAGS__misc.flags FFLAGS__muphy.flags FFLAGS__phy_common.flags FFLAGS__phys.flags FFLAGS__rrtm.flags FFLAGS__sisvat.flags FFLAGS__strataer.flags : FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS.flags :
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__aerono.flags FPPKEYS__chem.flags FPPKEYS__cloud.flags FPPKEYS__cosp.flags FPPKEYS__dust.flags FPPKEYS__dyn.flags FPPKEYS__dyn3d_common.flags FPPKEYS__dyn_phys.flags FPPKEYS__dyn_phys_sub.flags FPPKEYS__filtrez.flags FPPKEYS__grid.flags FPPKEYS__misc.flags FPPKEYS__muphy.flags FPPKEYS__phy_common.flags FPPKEYS__phys.flags FPPKEYS__rrtm.flags FPPKEYS__sisvat.flags FPPKEYS__strataer.flags : FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LD.flags :
	touch $(FCM_FLAGSDIR)/$@

LD__aerono.flags LD__chem.flags LD__cloud.flags LD__cosp.flags LD__dust.flags LD__dyn.flags LD__dyn3d_common.flags LD__dyn_phys.flags LD__dyn_phys_sub.flags LD__filtrez.flags LD__grid.flags LD__misc.flags LD__muphy.flags LD__phy_common.flags LD__phys.flags LD__rrtm.flags LD__sisvat.flags LD__strataer.flags : LD.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS.flags :
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__aerono.flags LDFLAGS__chem.flags LDFLAGS__cloud.flags LDFLAGS__cosp.flags LDFLAGS__dust.flags LDFLAGS__dyn.flags LDFLAGS__dyn3d_common.flags LDFLAGS__dyn_phys.flags LDFLAGS__dyn_phys_sub.flags LDFLAGS__filtrez.flags LDFLAGS__grid.flags LDFLAGS__misc.flags LDFLAGS__muphy.flags LDFLAGS__phy_common.flags LDFLAGS__phys.flags LDFLAGS__rrtm.flags LDFLAGS__sisvat.flags LDFLAGS__strataer.flags : LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

include $(FCM_BLDDIR)/aerono.mk
include $(FCM_BLDDIR)/dyn.mk
include $(FCM_BLDDIR)/dyn3d_common.mk
include $(FCM_BLDDIR)/dyn_phys.mk
include $(FCM_BLDDIR)/dyn_phys_sub.mk
include $(FCM_BLDDIR)/filtrez.mk
include $(FCM_BLDDIR)/grid.mk
include $(FCM_BLDDIR)/misc.mk
include $(FCM_BLDDIR)/phy_common.mk
include $(FCM_BLDDIR)/phys.mk

# EOF
