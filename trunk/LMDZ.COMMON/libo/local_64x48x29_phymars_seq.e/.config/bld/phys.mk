# Automatic Make rule for phys

PPSRCDIR0__phys = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libo/local_64x48x29_phymars_seq.e/.config/ppsrc/phys

SRCDIR0__phys = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libf/phymars

FFLAGS__phys__aerave.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aerave.done: \
          aerave.o
	touch $(FCM_DONEDIR)/$@

aerave.o: \
          $(PPSRCDIR0__phys)/aerave.f \
          FFLAGS__phys__aerave.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__aeropacity_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aeropacity_mod.done: \
          aeropacity_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          comgeomfi_h.done \
          density_co2_ice_mod.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          dust_scaling_mod.done \
          geometry_mod.done \
          ioipsl_getin_p_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

aeropacity_mod.o: \
          $(PPSRCDIR0__phys)/aeropacity_mod.f \
          FFLAGS__phys__aeropacity_mod.flags \
          callkeys.h \
          comcstfi_h.o \
          comgeomfi_h.o \
          density_co2_ice_mod.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          dust_scaling_mod.o \
          geometry_mod.o \
          ioipsl_getin_p_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__aeroptproperties.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

aeroptproperties.done: \
          aeroptproperties.o \
          callkeys.h.idone \
          dimradmars_mod.done
	touch $(FCM_DONEDIR)/$@

aeroptproperties.o: \
          $(PPSRCDIR0__phys)/aeroptproperties.f \
          FFLAGS__phys__aeroptproperties.flags \
          callkeys.h \
          dimradmars_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__albedocaps.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

albedocaps.done: \
          albedocaps.o \
          datafile_mod.done \
          geometry_mod.done \
          ioipsl_getin_p_mod.done \
          surfdat_h.done
	touch $(FCM_DONEDIR)/$@

albedocaps.o: \
          $(PPSRCDIR0__phys)/albedocaps.f90 \
          FFLAGS__phys__albedocaps.flags \
          datafile_mod.o \
          geometry_mod.o \
          ioipsl_getin_p_mod.o \
          surfdat_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__blackl.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

blackl.done: \
          blackl.o
	touch $(FCM_DONEDIR)/$@

blackl.o: \
          $(PPSRCDIR0__phys)/blackl.f \
          FFLAGS__phys__blackl.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__blendrad.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

blendrad.done: \
          blendrad.o \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

blendrad.o: \
          $(PPSRCDIR0__phys)/blendrad.f \
          FFLAGS__phys__blendrad.flags \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__calcstormfract_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

calcstormfract_mod.done: \
          calcstormfract_mod.o \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

calcstormfract_mod.o: \
          $(PPSRCDIR0__phys)/calcstormfract_mod.f90 \
          FFLAGS__phys__calcstormfract_mod.flags \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__calldrag_noro_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

calldrag_noro_mod.done: \
          calldrag_noro_mod.o \
          dimradmars_mod.done \
          drag_noro_mod.done \
          surfdat_h.done
	touch $(FCM_DONEDIR)/$@

calldrag_noro_mod.o: \
          $(PPSRCDIR0__phys)/calldrag_noro_mod.f \
          FFLAGS__phys__calldrag_noro_mod.flags \
          dimradmars_mod.o \
          drag_noro_mod.o \
          surfdat_h.o
	fcm_internal compile:F phys $< $@

callkeys.h: \
          $(SRCDIR0__phys)/callkeys.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

callkeys.h.idone: \
          $(SRCDIR0__phys)/callkeys.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__callradite_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

callradite_mod.done: \
          callradite_mod.o \
          aeropacity_mod.done \
          callkeys.h.idone \
          comcstfi_h.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          lwmain_mod.done \
          swmain_mod.done \
          time_phylmdz_mod.done \
          updatereffrad_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

callradite_mod.o: \
          $(PPSRCDIR0__phys)/callradite_mod.f \
          FFLAGS__phys__callradite_mod.flags \
          aeropacity_mod.o \
          callkeys.h \
          comcstfi_h.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          lwmain_mod.o \
          swmain_mod.o \
          time_phylmdz_mod.o \
          updatereffrad_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__callsedim_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

callsedim_mod.done: \
          callsedim_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          ioipsl_getin_p_mod.done \
          newsedim_mod.done \
          tracer_mod.done \
          updaterad.done
	touch $(FCM_DONEDIR)/$@

callsedim_mod.o: \
          $(PPSRCDIR0__phys)/callsedim_mod.f \
          FFLAGS__phys__callsedim_mod.flags \
          callkeys.h \
          comcstfi_h.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          ioipsl_getin_p_mod.o \
          newsedim_mod.o \
          tracer_mod.o \
          updaterad.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__calltherm_interface.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

calltherm_interface.done: \
          calltherm_interface.o \
          comcstfi_h.done \
          comtherm_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

calltherm_interface.o: \
          $(PPSRCDIR0__phys)/calltherm_interface.f90 \
          FFLAGS__phys__calltherm_interface.flags \
          comcstfi_h.o \
          comtherm_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__call_dayperi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

call_dayperi.done: \
          call_dayperi.o
	touch $(FCM_DONEDIR)/$@

call_dayperi.o: \
          $(PPSRCDIR0__phys)/call_dayperi.f \
          FFLAGS__phys__call_dayperi.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__check_fields.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

check_fields_mod.done: \
          check_fields_mod.o \
          dimphy.done
	touch $(FCM_DONEDIR)/$@

check_fields_mod.o: \
          $(PPSRCDIR0__phys)/check_fields.f90 \
          FFLAGS__phys__check_fields.flags \
          dimphy.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__co2cloud.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

co2cloud_mod.done: \
          co2cloud_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          conc_mod.done \
          datafile_mod.done \
          dimradmars_mod.done \
          improvedco2clouds_mod.done \
          ioipsl_getincom.done \
          microphys.h.idone \
          newsedim_mod.done \
          tracer_mod.done \
          updaterad.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

co2cloud_mod.o: \
          $(PPSRCDIR0__phys)/co2cloud.f90 \
          FFLAGS__phys__co2cloud.flags \
          callkeys.h \
          comcstfi_h.o \
          conc_mod.o \
          datafile_mod.o \
          dimradmars_mod.o \
          improvedco2clouds_mod.o \
          ioipsl_getincom.o \
          microphys.h \
          newsedim_mod.o \
          tracer_mod.o \
          updaterad.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__co2condens_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

co2condens_mod.done: \
          co2condens_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          dust_param_mod.done \
          geometry_mod.done \
          planete_h.done \
          surfdat_h.done \
          tracer_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

co2condens_mod.o: \
          $(PPSRCDIR0__phys)/co2condens_mod.f \
          FFLAGS__phys__co2condens_mod.flags \
          callkeys.h \
          comcstfi_h.o \
          dust_param_mod.o \
          geometry_mod.o \
          planete_h.o \
          surfdat_h.o \
          tracer_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__co2condens_mod4micro.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

co2condens_mod4micro.done: \
          co2condens_mod4micro.o \
          callkeys.h.idone \
          comcstfi_h.done \
          geometry_mod.done \
          planete_h.done \
          surfdat_h.done \
          tracer_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

co2condens_mod4micro.o: \
          $(PPSRCDIR0__phys)/co2condens_mod4micro.f90 \
          FFLAGS__phys__co2condens_mod4micro.flags \
          callkeys.h \
          comcstfi_h.o \
          geometry_mod.o \
          planete_h.o \
          surfdat_h.o \
          tracer_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__co2sat.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

co2sat.done: \
          co2sat.o
	touch $(FCM_DONEDIR)/$@

co2sat.o: \
          $(PPSRCDIR0__phys)/co2sat.f \
          FFLAGS__phys__co2sat.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__co2snow.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

co2snow.done: \
          co2snow.o \
          geometry_mod.done \
          surfdat_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

co2snow.o: \
          $(PPSRCDIR0__phys)/co2snow.f \
          FFLAGS__phys__co2snow.flags \
          geometry_mod.o \
          surfdat_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comcstfi_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comcstfi_h.done: \
          comcstfi_h.o
	touch $(FCM_DONEDIR)/$@

comcstfi_h.o: \
          $(PPSRCDIR0__phys)/comcstfi_h.f90 \
          FFLAGS__phys__comcstfi_h.flags
	fcm_internal compile:F phys $< $@

comg1d.h: \
          $(SRCDIR0__phys)/comg1d.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

comg1d.h.idone: \
          $(SRCDIR0__phys)/comg1d.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__comgeomfi_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comgeomfi_h.done: \
          comgeomfi_h.o
	touch $(FCM_DONEDIR)/$@

comgeomfi_h.o: \
          $(PPSRCDIR0__phys)/comgeomfi_h.f90 \
          FFLAGS__phys__comgeomfi_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comm_wrf.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comm_wrf.done: \
          comm_wrf.o
	touch $(FCM_DONEDIR)/$@

comm_wrf.o: \
          $(PPSRCDIR0__phys)/comm_wrf.f90 \
          FFLAGS__phys__comm_wrf.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__compute_dtau_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

compute_dtau_mod.done: \
          compute_dtau_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          geometry_mod.done \
          time_phylmdz_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

compute_dtau_mod.o: \
          $(PPSRCDIR0__phys)/compute_dtau_mod.f90 \
          FFLAGS__phys__compute_dtau_mod.flags \
          callkeys.h \
          comcstfi_h.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          geometry_mod.o \
          time_phylmdz_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comsaison_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comsaison_h.done: \
          comsaison_h.o
	touch $(FCM_DONEDIR)/$@

comsaison_h.o: \
          $(PPSRCDIR0__phys)/comsaison_h.f90 \
          FFLAGS__phys__comsaison_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comsoil_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comsoil_h.done: \
          comsoil_h.o
	touch $(FCM_DONEDIR)/$@

comsoil_h.o: \
          $(PPSRCDIR0__phys)/comsoil_h.f90 \
          FFLAGS__phys__comsoil_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__comtherm_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

comtherm_h.done: \
          comtherm_h.o
	touch $(FCM_DONEDIR)/$@

comtherm_h.o: \
          $(PPSRCDIR0__phys)/comtherm_h.f90 \
          FFLAGS__phys__comtherm_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__conc_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

conc_mod.done: \
          conc_mod.o
	touch $(FCM_DONEDIR)/$@

conc_mod.o: \
          $(PPSRCDIR0__phys)/conc_mod.f90 \
          FFLAGS__phys__conc_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__conf_phys.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

conf_phys.done: \
          conf_phys.o \
          aeropacity_mod.done \
          calchim_mod.done \
          callkeys.h.idone \
          co2condens_mod.done \
          datafile_mod.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          ioipsl_getin_p_mod.done \
          microphys.h.idone \
          surfdat_h.done \
          time_phylmdz_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

conf_phys.o: \
          $(PPSRCDIR0__phys)/conf_phys.f \
          FFLAGS__phys__conf_phys.flags \
          aeropacity_mod.o \
          calchim_mod.o \
          callkeys.h \
          co2condens_mod.o \
          datafile_mod.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          ioipsl_getin_p_mod.o \
          microphys.h \
          surfdat_h.o \
          time_phylmdz_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__convadj.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

convadj.done: \
          convadj.o \
          callkeys.h.idone \
          comcstfi_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

convadj.o: \
          $(PPSRCDIR0__phys)/convadj.f \
          FFLAGS__phys__convadj.flags \
          callkeys.h \
          comcstfi_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__cvmgp.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

cvmgp.done: \
          cvmgp.o
	touch $(FCM_DONEDIR)/$@

cvmgp.o: \
          $(PPSRCDIR0__phys)/cvmgp.f \
          FFLAGS__phys__cvmgp.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__cvmgt.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

cvmgt.done: \
          cvmgt.o
	touch $(FCM_DONEDIR)/$@

cvmgt.o: \
          $(PPSRCDIR0__phys)/cvmgt.f \
          FFLAGS__phys__cvmgt.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__datafile_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

datafile_mod.done: \
          datafile_mod.o
	touch $(FCM_DONEDIR)/$@

datafile_mod.o: \
          $(PPSRCDIR0__phys)/datafile_mod.f90 \
          FFLAGS__phys__datafile_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__def_var.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

def_var.done: \
          def_var.o
	touch $(FCM_DONEDIR)/$@

def_var.o: \
          $(PPSRCDIR0__phys)/def_var.f90 \
          FFLAGS__phys__def_var.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__density_co2_ice.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

density_co2_ice_mod.done: \
          density_co2_ice_mod.o
	touch $(FCM_DONEDIR)/$@

density_co2_ice_mod.o: \
          $(PPSRCDIR0__phys)/density_co2_ice.f90 \
          FFLAGS__phys__density_co2_ice.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dimphy.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dimphy.done: \
          dimphy.o
	touch $(FCM_DONEDIR)/$@

dimphy.o: \
          $(PPSRCDIR0__phys)/dimphy.f90 \
          FFLAGS__phys__dimphy.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dimradmars_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dimradmars_mod.done: \
          dimradmars_mod.o
	touch $(FCM_DONEDIR)/$@

dimradmars_mod.o: \
          $(PPSRCDIR0__phys)/dimradmars_mod.f90 \
          FFLAGS__phys__dimradmars_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__drag_noro_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

drag_noro_mod.done: \
          drag_noro_mod.o \
          comcstfi_h.done \
          dimradmars_mod.done \
          orodrag_mod.done
	touch $(FCM_DONEDIR)/$@

drag_noro_mod.o: \
          $(PPSRCDIR0__phys)/drag_noro_mod.f \
          FFLAGS__phys__drag_noro_mod.flags \
          comcstfi_h.o \
          dimradmars_mod.o \
          orodrag_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dustdevil.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dustdevil.done: \
          dustdevil.o \
          comcstfi_h.done \
          surfdat_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

dustdevil.o: \
          $(PPSRCDIR0__phys)/dustdevil.f \
          FFLAGS__phys__dustdevil.flags \
          comcstfi_h.o \
          surfdat_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dustlift.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dustlift.done: \
          dustlift.o \
          comcstfi_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

dustlift.o: \
          $(PPSRCDIR0__phys)/dustlift.f \
          FFLAGS__phys__dustlift.flags \
          comcstfi_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dust_param_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dust_param_mod.done: \
          dust_param_mod.o
	touch $(FCM_DONEDIR)/$@

dust_param_mod.o: \
          $(PPSRCDIR0__phys)/dust_param_mod.f90 \
          FFLAGS__phys__dust_param_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dust_rad_adjust_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dust_rad_adjust_mod.done: \
          dust_rad_adjust_mod.o \
          dust_param_mod.done \
          geometry_mod.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

dust_rad_adjust_mod.o: \
          $(PPSRCDIR0__phys)/dust_rad_adjust_mod.f90 \
          FFLAGS__phys__dust_rad_adjust_mod.flags \
          dust_param_mod.o \
          geometry_mod.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__dust_scaling_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

dust_scaling_mod.done: \
          dust_scaling_mod.o \
          dimradmars_mod.done \
          dust_param_mod.done \
          dust_rad_adjust_mod.done
	touch $(FCM_DONEDIR)/$@

dust_scaling_mod.o: \
          $(PPSRCDIR0__phys)/dust_scaling_mod.f90 \
          FFLAGS__phys__dust_scaling_mod.flags \
          dimradmars_mod.o \
          dust_param_mod.o \
          dust_rad_adjust_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__eofdump_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

eofdump_mod.done: \
          eofdump_mod.o \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          nrtype.done \
          time_phylmdz_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

eofdump_mod.o: \
          $(PPSRCDIR0__phys)/eofdump_mod.f90 \
          FFLAGS__phys__eofdump_mod.flags \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          nrtype.o \
          time_phylmdz_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__flusv.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

flusv.done: \
          flusv.o \
          dimradmars_mod.done
	touch $(FCM_DONEDIR)/$@

flusv.o: \
          $(PPSRCDIR0__phys)/flusv.f \
          FFLAGS__phys__flusv.flags \
          dimradmars_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__geticecover.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

geticecover.done: \
          geticecover.o
	touch $(FCM_DONEDIR)/$@

geticecover.o: \
          $(PPSRCDIR0__phys)/geticecover.f90 \
          FFLAGS__phys__geticecover.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__getslopes.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

getslopes.done: \
          getslopes.o \
          comcstfi_h.done \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          slope_mod.done
	touch $(FCM_DONEDIR)/$@

getslopes.o: \
          $(PPSRCDIR0__phys)/getslopes.f90 \
          FFLAGS__phys__getslopes.flags \
          comcstfi_h.o \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          slope_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__growthrate.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

growthrate.done: \
          growthrate.o \
          comcstfi_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

growthrate.o: \
          $(PPSRCDIR0__phys)/growthrate.f \
          FFLAGS__phys__growthrate.flags \
          comcstfi_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__gwprofil_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

gwprofil_mod.done: \
          gwprofil_mod.o \
          dimradmars_mod.done
	touch $(FCM_DONEDIR)/$@

gwprofil_mod.o: \
          $(PPSRCDIR0__phys)/gwprofil_mod.f \
          FFLAGS__phys__gwprofil_mod.flags \
          dimradmars_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__gwstress_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

gwstress_mod.done: \
          gwstress_mod.o \
          dimradmars_mod.done \
          yoegwd.h.idone
	touch $(FCM_DONEDIR)/$@

gwstress_mod.o: \
          $(PPSRCDIR0__phys)/gwstress_mod.f \
          FFLAGS__phys__gwstress_mod.flags \
          dimradmars_mod.o \
          yoegwd.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__hdo_surfex_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

hdo_surfex_mod.done: \
          hdo_surfex_mod.o \
          callkeys.h.idone \
          geometry_mod.done \
          surfdat_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

hdo_surfex_mod.o: \
          $(PPSRCDIR0__phys)/hdo_surfex_mod.f \
          FFLAGS__phys__hdo_surfex_mod.flags \
          callkeys.h \
          geometry_mod.o \
          surfdat_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__improvedclouds_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

improvedclouds_mod.done: \
          improvedclouds_mod.o \
          comcstfi_h.done \
          conc_mod.done \
          tracer_mod.done \
          updaterad.done \
          watersat_mod.done
	touch $(FCM_DONEDIR)/$@

improvedclouds_mod.o: \
          $(PPSRCDIR0__phys)/improvedclouds_mod.f \
          FFLAGS__phys__improvedclouds_mod.flags \
          comcstfi_h.o \
          conc_mod.o \
          tracer_mod.o \
          updaterad.o \
          watersat_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__improvedco2clouds_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

improvedco2clouds_mod.done: \
          improvedco2clouds_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          conc_mod.done \
          datafile_mod.done \
          density_co2_ice_mod.done \
          microphys.h.idone \
          tracer_mod.done \
          updaterad.done
	touch $(FCM_DONEDIR)/$@

improvedco2clouds_mod.o: \
          $(PPSRCDIR0__phys)/improvedco2clouds_mod.f90 \
          FFLAGS__phys__improvedco2clouds_mod.flags \
          callkeys.h \
          comcstfi_h.o \
          conc_mod.o \
          datafile_mod.o \
          density_co2_ice_mod.o \
          microphys.h \
          tracer_mod.o \
          updaterad.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniorbit.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniorbit.done: \
          iniorbit.o \
          comcstfi_h.done \
          planete_h.done
	touch $(FCM_DONEDIR)/$@

iniorbit.o: \
          $(PPSRCDIR0__phys)/iniorbit.f \
          FFLAGS__phys__iniorbit.flags \
          comcstfi_h.o \
          planete_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__inistats.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

inistats.done: \
          inistats.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          nrtype.done \
          regular_lonlat_mod.done \
          statto_mod.done \
          time_phylmdz_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

inistats.o: \
          $(PPSRCDIR0__phys)/inistats.f \
          FFLAGS__phys__inistats.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          nrtype.o \
          regular_lonlat_mod.o \
          statto_mod.o \
          time_phylmdz_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__initracer.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

initracer.done: \
          initracer.o \
          callkeys.h.idone \
          comcstfi_h.done \
          dust_param_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

initracer.o: \
          $(PPSRCDIR0__phys)/initracer.f \
          FFLAGS__phys__initracer.flags \
          callkeys.h \
          comcstfi_h.o \
          dust_param_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniwrite.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniwrite.done: \
          iniwrite.o \
          comcstfi_h.done \
          comsoil_h.done \
          mod_grid_phy_lmdz.done \
          regular_lonlat_mod.done \
          time_phylmdz_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

iniwrite.o: \
          $(PPSRCDIR0__phys)/iniwrite.f \
          FFLAGS__phys__iniwrite.flags \
          comcstfi_h.o \
          comsoil_h.o \
          mod_grid_phy_lmdz.o \
          regular_lonlat_mod.o \
          time_phylmdz_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iniwritesoil.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iniwritesoil.done: \
          iniwritesoil.o \
          comcstfi_h.done \
          comsoil_h.done \
          mod_grid_phy_lmdz.done \
          regular_lonlat_mod.done
	touch $(FCM_DONEDIR)/$@

iniwritesoil.o: \
          $(PPSRCDIR0__phys)/iniwritesoil.f90 \
          FFLAGS__phys__iniwritesoil.flags \
          comcstfi_h.o \
          comsoil_h.o \
          mod_grid_phy_lmdz.o \
          regular_lonlat_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__interp_line.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

interp_line.done: \
          interp_line.o
	touch $(FCM_DONEDIR)/$@

interp_line.o: \
          $(PPSRCDIR0__phys)/interp_line.f \
          FFLAGS__phys__interp_line.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__iostart.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

iostart.done: \
          iostart.o \
          comsoil_h.done \
          dimphy.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

iostart.o: \
          $(PPSRCDIR0__phys)/iostart.f90 \
          FFLAGS__phys__iostart.flags \
          comsoil_h.o \
          dimphy.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwb.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwb.done: \
          lwb.o \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwb.o: \
          $(PPSRCDIR0__phys)/lwb.f \
          FFLAGS__phys__lwb.flags \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwdiff.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwdiff.done: \
          lwdiff.o \
          comcstfi_h.done \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwdiff.o: \
          $(PPSRCDIR0__phys)/lwdiff.f \
          FFLAGS__phys__lwdiff.flags \
          comcstfi_h.o \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwflux.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwflux.done: \
          lwflux.o \
          callkeys.h.idone \
          comg1d.h.idone \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwflux.o: \
          $(PPSRCDIR0__phys)/lwflux.f \
          FFLAGS__phys__lwflux.flags \
          callkeys.h \
          comg1d.h \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwi.done: \
          lwi.o \
          comcstfi_h.done \
          dimradmars_mod.done \
          time_phylmdz_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwi.o: \
          $(PPSRCDIR0__phys)/lwi.f \
          FFLAGS__phys__lwi.flags \
          comcstfi_h.o \
          dimradmars_mod.o \
          time_phylmdz_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwmain_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwmain_mod.done: \
          lwmain_mod.o \
          callkeys.h.idone \
          comg1d.h.idone \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwmain_mod.o: \
          $(PPSRCDIR0__phys)/lwmain_mod.f \
          FFLAGS__phys__lwmain_mod.flags \
          callkeys.h \
          comg1d.h \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwtt.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwtt.done: \
          lwtt.o \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwtt.o: \
          $(PPSRCDIR0__phys)/lwtt.f \
          FFLAGS__phys__lwtt.flags \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwu.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwu.done: \
          lwu.o \
          callkeys.h.idone \
          comcstfi_h.done \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwu.o: \
          $(PPSRCDIR0__phys)/lwu.f \
          FFLAGS__phys__lwu.flags \
          callkeys.h \
          comcstfi_h.o \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwxb.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwxb.done: \
          lwxb.o \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwxb.o: \
          $(PPSRCDIR0__phys)/lwxb.f \
          FFLAGS__phys__lwxb.flags \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwxd.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwxd.done: \
          lwxd.o \
          callkeys.h.idone \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwxd.o: \
          $(PPSRCDIR0__phys)/lwxd.f \
          FFLAGS__phys__lwxd.flags \
          callkeys.h \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__lwxn.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

lwxn.done: \
          lwxn.o \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

lwxn.o: \
          $(PPSRCDIR0__phys)/lwxn.f \
          FFLAGS__phys__lwxn.flags \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__massflowrateco2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

massflowrateco2.done: \
          massflowrateco2.o \
          comcstfi_h.done \
          microphys.h.idone \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

massflowrateco2.o: \
          $(PPSRCDIR0__phys)/massflowrateco2.f90 \
          FFLAGS__phys__massflowrateco2.flags \
          comcstfi_h.o \
          microphys.h \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

microphys.h: \
          $(SRCDIR0__phys)/microphys.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

microphys.h.idone: \
          $(SRCDIR0__phys)/microphys.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__mkstat.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mkstats.done: \
          mkstats.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          statto_mod.done
	touch $(FCM_DONEDIR)/$@

mkstats.o: \
          $(PPSRCDIR0__phys)/mkstat.f90 \
          FFLAGS__phys__mkstat.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          statto_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__mucorr.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mucorr.done: \
          mucorr.o
	touch $(FCM_DONEDIR)/$@

mucorr.o: \
          $(PPSRCDIR0__phys)/mucorr.f \
          FFLAGS__phys__mucorr.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__mufract.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mufract.done: \
          mufract.o
	touch $(FCM_DONEDIR)/$@

mufract.o: \
          $(PPSRCDIR0__phys)/mufract.f \
          FFLAGS__phys__mufract.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__newsedim_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

newsedim_mod.done: \
          newsedim_mod.o \
          comcstfi_h.done
	touch $(FCM_DONEDIR)/$@

newsedim_mod.o: \
          $(PPSRCDIR0__phys)/newsedim_mod.f \
          FFLAGS__phys__newsedim_mod.flags \
          comcstfi_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nirco2abs.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nirco2abs.done: \
          nirco2abs.o \
          callkeys.h.idone \
          comcstfi_h.done \
          comgeomfi_h.done \
          nirdata.h.idone \
          time_phylmdz_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

nirco2abs.o: \
          $(PPSRCDIR0__phys)/nirco2abs.f \
          FFLAGS__phys__nirco2abs.flags \
          callkeys.h \
          comcstfi_h.o \
          comgeomfi_h.o \
          nirdata.h \
          time_phylmdz_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

nirdata.h: \
          $(SRCDIR0__phys)/nirdata.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

nirdata.h.idone: \
          $(SRCDIR0__phys)/nirdata.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__nir_leedat.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nir_leedat.done: \
          nir_leedat.o \
          datafile_mod.done \
          nirdata.h.idone
	touch $(FCM_DONEDIR)/$@

nir_leedat.o: \
          $(PPSRCDIR0__phys)/nir_leedat.f \
          FFLAGS__phys__nir_leedat.flags \
          datafile_mod.o \
          nirdata.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nltecool.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nltecool.done: \
          nltecool.o \
          conc_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

nltecool.o: \
          $(PPSRCDIR0__phys)/nltecool.f \
          FFLAGS__phys__nltecool.flags \
          conc_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

nltedata.h: \
          $(SRCDIR0__phys)/nltedata.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

nltedata.h.idone: \
          $(SRCDIR0__phys)/nltedata.h
	touch $(FCM_DONEDIR)/$@

nlteparams.h: \
          $(SRCDIR0__phys)/nlteparams.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

nlteparams.h.idone: \
          $(SRCDIR0__phys)/nlteparams.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__nlte_aux.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

planckdp.done: \
          planckdp.o \
          nlte_commons.h.idone \
          nlte_paramdef.h.idone
	touch $(FCM_DONEDIR)/$@

planckdp.o: \
          $(PPSRCDIR0__phys)/nlte_aux.f \
          FFLAGS__phys__nlte_aux.flags \
          nlte_commons.h \
          nlte_paramdef.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nlte_calc.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

mzesc110.done: \
          mzesc110.o \
          nlte_commons.h.idone \
          nlte_paramdef.h.idone
	touch $(FCM_DONEDIR)/$@

mzesc110.o: \
          $(PPSRCDIR0__phys)/nlte_calc.f \
          FFLAGS__phys__nlte_calc.flags \
          nlte_commons.h \
          nlte_paramdef.h
	fcm_internal compile:F phys $< $@

nlte_commons.h: \
          $(SRCDIR0__phys)/nlte_commons.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

nlte_commons.h.idone: \
          $(SRCDIR0__phys)/nlte_commons.h
	touch $(FCM_DONEDIR)/$@

nlte_paramdef.h: \
          $(SRCDIR0__phys)/nlte_paramdef.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

nlte_paramdef.h.idone: \
          $(SRCDIR0__phys)/nlte_paramdef.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__nlte_setup.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nlte_setup.done: \
          nlte_setup.o \
          datafile_mod.done \
          nlte_commons.h.idone \
          nlte_paramdef.h.idone
	touch $(FCM_DONEDIR)/$@

nlte_setup.o: \
          $(PPSRCDIR0__phys)/nlte_setup.f \
          FFLAGS__phys__nlte_setup.flags \
          datafile_mod.o \
          nlte_commons.h \
          nlte_paramdef.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nlte_tcool.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nlte_tcool.done: \
          nlte_tcool.o \
          chimiedata.h.idone \
          conc_mod.done \
          nlte_commons.h.idone \
          nlte_paramdef.h.idone \
          param_v4_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

nlte_tcool.o: \
          $(PPSRCDIR0__phys)/nlte_tcool.f \
          FFLAGS__phys__nlte_tcool.flags \
          chimiedata.h \
          conc_mod.o \
          nlte_commons.h \
          nlte_paramdef.h \
          param_v4_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nlthermeq.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nlthermeq.done: \
          nlthermeq.o \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

nlthermeq.o: \
          $(PPSRCDIR0__phys)/nlthermeq.f \
          FFLAGS__phys__nlthermeq.flags \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nonoro_gwd_ran_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nonoro_gwd_ran_mod.done: \
          nonoro_gwd_ran_mod.o \
          assert_m.done \
          callkeys.h.idone \
          comcstfi_h.done \
          geometry_mod.done \
          ioipsl_getin_p_mod.done \
          vertical_layers_mod.done \
          yoegwd.h.idone
	touch $(FCM_DONEDIR)/$@

nonoro_gwd_ran_mod.o: \
          $(PPSRCDIR0__phys)/nonoro_gwd_ran_mod.f90 \
          FFLAGS__phys__nonoro_gwd_ran_mod.flags \
          assert_m.o \
          callkeys.h \
          comcstfi_h.o \
          geometry_mod.o \
          ioipsl_getin_p_mod.o \
          vertical_layers_mod.o \
          yoegwd.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nuclea.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nuclea.done: \
          nuclea.o \
          callkeys.h.idone \
          comcstfi_h.done
	touch $(FCM_DONEDIR)/$@

nuclea.o: \
          $(PPSRCDIR0__phys)/nuclea.f \
          FFLAGS__phys__nuclea.flags \
          callkeys.h \
          comcstfi_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__nucleaco2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

nucleaco2.done: \
          nucleaco2.o \
          callkeys.h.idone \
          comcstfi_h.done \
          microphys.h.idone
	touch $(FCM_DONEDIR)/$@

nucleaco2.o: \
          $(PPSRCDIR0__phys)/nucleaco2.f \
          FFLAGS__phys__nucleaco2.flags \
          callkeys.h \
          comcstfi_h.o \
          microphys.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__orbite.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

orbite.done: \
          orbite.o \
          comcstfi_h.done \
          planete_h.done
	touch $(FCM_DONEDIR)/$@

orbite.o: \
          $(PPSRCDIR0__phys)/orbite.f \
          FFLAGS__phys__orbite.flags \
          comcstfi_h.o \
          planete_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__orodrag_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

orodrag_mod.done: \
          orodrag_mod.o \
          comcstfi_h.done \
          dimradmars_mod.done \
          gwprofil_mod.done \
          gwstress_mod.done \
          yoegwd.h.idone
	touch $(FCM_DONEDIR)/$@

orodrag_mod.o: \
          $(PPSRCDIR0__phys)/orodrag_mod.f \
          FFLAGS__phys__orodrag_mod.flags \
          comcstfi_h.o \
          dimradmars_mod.o \
          gwprofil_mod.o \
          gwstress_mod.o \
          yoegwd.h
	fcm_internal compile:F phys $< $@

FFLAGS__phys__orosetup.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

orosetup.done: \
          orosetup.o \
          comcstfi_h.done \
          dimradmars_mod.done
	touch $(FCM_DONEDIR)/$@

orosetup.o: \
          $(PPSRCDIR0__phys)/orosetup.f \
          FFLAGS__phys__orosetup.flags \
          comcstfi_h.o \
          dimradmars_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__param_slope.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

param_slope.done: \
          param_slope.o
	touch $(FCM_DONEDIR)/$@

param_slope.o: \
          $(PPSRCDIR0__phys)/param_slope.f90 \
          FFLAGS__phys__param_slope.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__pbl_parameters.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

pbl_parameters.done: \
          pbl_parameters.o \
          comcstfi_h.done
	touch $(FCM_DONEDIR)/$@

pbl_parameters.o: \
          $(PPSRCDIR0__phys)/pbl_parameters.f \
          FFLAGS__phys__pbl_parameters.flags \
          comcstfi_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__phyetat0_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

phyetat0_mod.done: \
          phyetat0_mod.o \
          callkeys.h.idone \
          compute_dtau_mod.done \
          dust_param_mod.done \
          dust_rad_adjust_mod.done \
          ioipsl_getin_p_mod.done \
          iostart.done \
          nonoro_gwd_ran_mod.done \
          surfdat_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

phyetat0_mod.o: \
          $(PPSRCDIR0__phys)/phyetat0_mod.f90 \
          FFLAGS__phys__phyetat0_mod.flags \
          callkeys.h \
          compute_dtau_mod.o \
          dust_param_mod.o \
          dust_rad_adjust_mod.o \
          ioipsl_getin_p_mod.o \
          iostart.o \
          nonoro_gwd_ran_mod.o \
          surfdat_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__phyredem.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

phyredem.done: \
          phyredem.o \
          callkeys.h.idone \
          comcstfi_h.done \
          compute_dtau_mod.done \
          comsoil_h.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          dust_rad_adjust_mod.done \
          geometry_mod.done \
          iostart.done \
          mod_grid_phy_lmdz.done \
          nonoro_gwd_ran_mod.done \
          planete_h.done \
          surfdat_h.done \
          time_phylmdz_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

phyredem.o: \
          $(PPSRCDIR0__phys)/phyredem.f90 \
          FFLAGS__phys__phyredem.flags \
          callkeys.h \
          comcstfi_h.o \
          compute_dtau_mod.o \
          comsoil_h.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          dust_rad_adjust_mod.o \
          geometry_mod.o \
          iostart.o \
          mod_grid_phy_lmdz.o \
          nonoro_gwd_ran_mod.o \
          planete_h.o \
          surfdat_h.o \
          time_phylmdz_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__physiq_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

physiq_mod.done: \
          physiq_mod.o \
          calchim_mod.done \
          calcstormfract_mod.done \
          calldrag_noro_mod.done \
          callkeys.h.idone \
          callradite_mod.done \
          callsedim_mod.done \
          check_fields_mod.done \
          co2cloud_mod.done \
          co2condens_mod.done \
          co2condens_mod4micro.done \
          comcstfi_h.done \
          comg1d.h.idone \
          comgeomfi_h.done \
          compute_dtau_mod.done \
          comsaison_h.done \
          comsoil_h.done \
          conc_mod.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          eofdump_mod.done \
          geometry_mod.done \
          ioipsl_getin_p_mod.done \
          iono_h.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_omp_data.done \
          nlteparams.h.idone \
          nonoro_gwd_ran_mod.done \
          param_v4_h.done \
          phyetat0_mod.done \
          phyredem.done \
          planete_h.done \
          planetwide_mod.done \
          rocketduststorm_mod.done \
          slope_mod.done \
          surfdat_h.done \
          time_phylmdz_mod.done \
          topmons_mod.done \
          tracer_mod.done \
          turb_mod.done \
          vdifc_mod.done \
          vertical_layers_mod.done \
          watercloud_mod.done \
          watersat_mod.done
	touch $(FCM_DONEDIR)/$@

physiq_mod.o: \
          $(PPSRCDIR0__phys)/physiq_mod.f \
          FFLAGS__phys__physiq_mod.flags \
          calchim_mod.o \
          calcstormfract_mod.o \
          calldrag_noro_mod.o \
          callkeys.h \
          callradite_mod.o \
          callsedim_mod.o \
          check_fields_mod.o \
          co2cloud_mod.o \
          co2condens_mod.o \
          co2condens_mod4micro.o \
          comcstfi_h.o \
          comg1d.h \
          comgeomfi_h.o \
          compute_dtau_mod.o \
          comsaison_h.o \
          comsoil_h.o \
          conc_mod.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          eofdump_mod.o \
          geometry_mod.o \
          ioipsl_getin_p_mod.o \
          iono_h.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_omp_data.o \
          nlteparams.h \
          nonoro_gwd_ran_mod.o \
          param_v4_h.o \
          phyetat0_mod.o \
          phyredem.o \
          planete_h.o \
          planetwide_mod.o \
          rocketduststorm_mod.o \
          slope_mod.o \
          surfdat_h.o \
          time_phylmdz_mod.o \
          topmons_mod.o \
          tracer_mod.o \
          turb_mod.o \
          vdifc_mod.o \
          vertical_layers_mod.o \
          watercloud_mod.o \
          watersat_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__phys_state_var_init_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

phys_state_var_init_mod.done: \
          phys_state_var_init_mod.o \
          calchim_mod.done \
          co2cloud_mod.done \
          comcstfi_h.done \
          comgeomfi_h.done \
          compute_dtau_mod.done \
          comsaison_h.done \
          comsoil_h.done \
          conc_mod.done \
          dimradmars_mod.done \
          dust_param_mod.done \
          dust_rad_adjust_mod.done \
          init_print_control_mod.done \
          nonoro_gwd_ran_mod.done \
          rocketduststorm_mod.done \
          slope_mod.done \
          surfdat_h.done \
          time_phylmdz_mod.done \
          topmons_mod.done \
          tracer_mod.done \
          turb_mod.done \
          watercloud_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

phys_state_var_init_mod.o: \
          $(PPSRCDIR0__phys)/phys_state_var_init_mod.f90 \
          FFLAGS__phys__phys_state_var_init_mod.flags \
          calchim_mod.o \
          co2cloud_mod.o \
          comcstfi_h.o \
          comgeomfi_h.o \
          compute_dtau_mod.o \
          comsaison_h.o \
          comsoil_h.o \
          conc_mod.o \
          dimradmars_mod.o \
          dust_param_mod.o \
          dust_rad_adjust_mod.o \
          init_print_control_mod.o \
          nonoro_gwd_ran_mod.o \
          rocketduststorm_mod.o \
          slope_mod.o \
          surfdat_h.o \
          time_phylmdz_mod.o \
          topmons_mod.o \
          tracer_mod.o \
          turb_mod.o \
          watercloud_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__planete_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

planete_h.done: \
          planete_h.o
	touch $(FCM_DONEDIR)/$@

planete_h.o: \
          $(PPSRCDIR0__phys)/planete_h.f90 \
          FFLAGS__phys__planete_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__planetwide_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

planetwide_mod.done: \
          planetwide_mod.o \
          dimphy.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done
	touch $(FCM_DONEDIR)/$@

planetwide_mod.o: \
          $(PPSRCDIR0__phys)/planetwide_mod.f90 \
          FFLAGS__phys__planetwide_mod.flags \
          dimphy.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__read_dust_scenario.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

read_dust_scenario.done: \
          read_dust_scenario.o \
          callkeys.h.idone \
          datafile_mod.done \
          dust_param_mod.done \
          geometry_mod.done \
          planete_h.done
	touch $(FCM_DONEDIR)/$@

read_dust_scenario.o: \
          $(PPSRCDIR0__phys)/read_dust_scenario.f90 \
          FFLAGS__phys__read_dust_scenario.flags \
          callkeys.h \
          datafile_mod.o \
          dust_param_mod.o \
          geometry_mod.o \
          planete_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__rocketduststorm_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

rocketduststorm_mod.done: \
          rocketduststorm_mod.o \
          callkeys.h.idone \
          callradite_mod.done \
          comcstfi_h.done \
          comsaison_h.done \
          dimradmars_mod.done \
          surfdat_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

rocketduststorm_mod.o: \
          $(PPSRCDIR0__phys)/rocketduststorm_mod.f90 \
          FFLAGS__phys__rocketduststorm_mod.flags \
          callkeys.h \
          callradite_mod.o \
          comcstfi_h.o \
          comsaison_h.o \
          dimradmars_mod.o \
          surfdat_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__scopyi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

scopyi.done: \
          scopyi.o
	touch $(FCM_DONEDIR)/$@

scopyi.o: \
          $(PPSRCDIR0__phys)/scopyi.f \
          FFLAGS__phys__scopyi.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__sig.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

sig.done: \
          sig.o
	touch $(FCM_DONEDIR)/$@

sig.o: \
          $(PPSRCDIR0__phys)/sig.f \
          FFLAGS__phys__sig.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__simpleclouds.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

simpleclouds.done: \
          simpleclouds.o \
          comcstfi_h.done \
          dimradmars_mod.done \
          tracer_mod.done \
          updaterad.done \
          watersat_mod.done
	touch $(FCM_DONEDIR)/$@

simpleclouds.o: \
          $(PPSRCDIR0__phys)/simpleclouds.f \
          FFLAGS__phys__simpleclouds.flags \
          comcstfi_h.o \
          dimradmars_mod.o \
          tracer_mod.o \
          updaterad.o \
          watersat_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__slope_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

slope_mod.done: \
          slope_mod.o
	touch $(FCM_DONEDIR)/$@

slope_mod.o: \
          $(PPSRCDIR0__phys)/slope_mod.f90 \
          FFLAGS__phys__slope_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__soil.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

soil.done: \
          soil.o \
          comsoil_h.done \
          surfdat_h.done
	touch $(FCM_DONEDIR)/$@

soil.o: \
          $(PPSRCDIR0__phys)/soil.f \
          FFLAGS__phys__soil.flags \
          comsoil_h.o \
          surfdat_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__soil_settings.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

soil_settings.done: \
          soil_settings.o \
          comsoil_h.done \
          iostart.done
	touch $(FCM_DONEDIR)/$@

soil_settings.o: \
          $(PPSRCDIR0__phys)/soil_settings.f \
          FFLAGS__phys__soil_settings.flags \
          comsoil_h.o \
          iostart.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__soil_tifeedback.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

soil_tifeedback.done: \
          soil_tifeedback.o \
          comsoil_h.done \
          surfdat_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

soil_tifeedback.o: \
          $(PPSRCDIR0__phys)/soil_tifeedback.f \
          FFLAGS__phys__soil_tifeedback.flags \
          comsoil_h.o \
          surfdat_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__solang.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

solang.done: \
          solang.o
	touch $(FCM_DONEDIR)/$@

solang.o: \
          $(PPSRCDIR0__phys)/solang.f \
          FFLAGS__phys__solang.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__solarlong.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

solarlong.done: \
          solarlong.o \
          comcstfi_h.done \
          planete_h.done
	touch $(FCM_DONEDIR)/$@

solarlong.o: \
          $(PPSRCDIR0__phys)/solarlong.f \
          FFLAGS__phys__solarlong.flags \
          comcstfi_h.o \
          planete_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__statto_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

statto_mod.done: \
          statto_mod.o
	touch $(FCM_DONEDIR)/$@

statto_mod.o: \
          $(PPSRCDIR0__phys)/statto_mod.f90 \
          FFLAGS__phys__statto_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__suaer.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

suaer.done: \
          suaer.o \
          callkeys.h.idone \
          datafile_mod.done \
          dimradmars_mod.done
	touch $(FCM_DONEDIR)/$@

suaer.o: \
          $(PPSRCDIR0__phys)/suaer.f90 \
          FFLAGS__phys__suaer.flags \
          callkeys.h \
          datafile_mod.o \
          dimradmars_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__sugwd.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

sugwd.done: \
          sugwd.o
	touch $(FCM_DONEDIR)/$@

sugwd.o: \
          $(PPSRCDIR0__phys)/sugwd.f \
          FFLAGS__phys__sugwd.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__surfdat_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

surfdat_h.done: \
          surfdat_h.o
	touch $(FCM_DONEDIR)/$@

surfdat_h.o: \
          $(PPSRCDIR0__phys)/surfdat_h.f90 \
          FFLAGS__phys__surfdat_h.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__surfini.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

surfini.done: \
          surfini.o \
          callkeys.h.idone \
          comcstfi_h.done \
          datafile_mod.done \
          geometry_mod.done \
          ioipsl_getin_p_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          surfdat_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

surfini.o: \
          $(PPSRCDIR0__phys)/surfini.f \
          FFLAGS__phys__surfini.flags \
          callkeys.h \
          comcstfi_h.o \
          datafile_mod.o \
          geometry_mod.o \
          ioipsl_getin_p_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          surfdat_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__swmain_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

swmain_mod.done: \
          swmain_mod.o \
          callkeys.h.idone \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

swmain_mod.o: \
          $(PPSRCDIR0__phys)/swmain_mod.f \
          FFLAGS__phys__swmain_mod.flags \
          callkeys.h \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__swrayleigh.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

swrayleigh.done: \
          swrayleigh.o \
          comcstfi_h.done
	touch $(FCM_DONEDIR)/$@

swrayleigh.o: \
          $(PPSRCDIR0__phys)/swrayleigh.f \
          FFLAGS__phys__swrayleigh.flags \
          comcstfi_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__swr_fouquart.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

swr_fouquart.done: \
          swr_fouquart.o \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

swr_fouquart.o: \
          $(PPSRCDIR0__phys)/swr_fouquart.f \
          FFLAGS__phys__swr_fouquart.flags \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__swr_toon.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

swr_toon.done: \
          swr_toon.o \
          dimradmars_mod.done \
          yomlw_h.done
	touch $(FCM_DONEDIR)/$@

swr_toon.o: \
          $(PPSRCDIR0__phys)/swr_toon.f \
          FFLAGS__phys__swr_toon.flags \
          dimradmars_mod.o \
          yomlw_h.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__tabfi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

tabfi.done: \
          tabfi.o \
          comcstfi_h.done \
          comsoil_h.done \
          dimradmars_mod.done \
          ioipsl_getin_p_mod.done \
          iostart.done \
          mod_phys_lmdz_para.done \
          planete_h.done \
          surfdat_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

tabfi.o: \
          $(PPSRCDIR0__phys)/tabfi.f \
          FFLAGS__phys__tabfi.flags \
          comcstfi_h.o \
          comsoil_h.o \
          dimradmars_mod.o \
          ioipsl_getin_p_mod.o \
          iostart.o \
          mod_phys_lmdz_para.o \
          planete_h.o \
          surfdat_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__tcondco2.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

tcondco2.done: \
          tcondco2.o \
          comcstfi_h.done \
          conc_mod.done
	touch $(FCM_DONEDIR)/$@

tcondco2.o: \
          $(PPSRCDIR0__phys)/tcondco2.f90 \
          FFLAGS__phys__tcondco2.flags \
          comcstfi_h.o \
          conc_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__tcondwater.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

tcondwater.done: \
          tcondwater.o
	touch $(FCM_DONEDIR)/$@

tcondwater.o: \
          $(PPSRCDIR0__phys)/tcondwater.f90 \
          FFLAGS__phys__tcondwater.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_dqup.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_dqup.done: \
          thermcell_dqup.o
	touch $(FCM_DONEDIR)/$@

thermcell_dqup.o: \
          $(PPSRCDIR0__phys)/thermcell_dqup.f90 \
          FFLAGS__phys__thermcell_dqup.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__thermcell_main_mars.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

thermcell_main_mars.done: \
          thermcell_main_mars.o \
          comcstfi_h.done \
          comtherm_h.done \
          planetwide_mod.done
	touch $(FCM_DONEDIR)/$@

thermcell_main_mars.o: \
          $(PPSRCDIR0__phys)/thermcell_main_mars.f90 \
          FFLAGS__phys__thermcell_main_mars.flags \
          comcstfi_h.o \
          comtherm_h.o \
          planetwide_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__time_phylmdz_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

time_phylmdz_mod.done: \
          time_phylmdz_mod.o
	touch $(FCM_DONEDIR)/$@

time_phylmdz_mod.o: \
          $(PPSRCDIR0__phys)/time_phylmdz_mod.f90 \
          FFLAGS__phys__time_phylmdz_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__topmons_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

topmons_mod.done: \
          topmons_mod.o \
          callradite_mod.done \
          comcstfi_h.done \
          comsaison_h.done \
          dimradmars_mod.done \
          surfdat_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

topmons_mod.o: \
          $(PPSRCDIR0__phys)/topmons_mod.f90 \
          FFLAGS__phys__topmons_mod.flags \
          callradite_mod.o \
          comcstfi_h.o \
          comsaison_h.o \
          dimradmars_mod.o \
          surfdat_h.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__tracer_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

tracer_mod.done: \
          tracer_mod.o
	touch $(FCM_DONEDIR)/$@

tracer_mod.o: \
          $(PPSRCDIR0__phys)/tracer_mod.f90 \
          FFLAGS__phys__tracer_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__turb_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

turb_mod.done: \
          turb_mod.o
	touch $(FCM_DONEDIR)/$@

turb_mod.o: \
          $(PPSRCDIR0__phys)/turb_mod.f90 \
          FFLAGS__phys__turb_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__updaterad.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

updaterad.done: \
          updaterad.o \
          comcstfi_h.done \
          density_co2_ice_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

updaterad.o: \
          $(PPSRCDIR0__phys)/updaterad.f90 \
          FFLAGS__phys__updaterad.flags \
          comcstfi_h.o \
          density_co2_ice_mod.o \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__updatereffrad_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

updatereffrad_mod.done: \
          updatereffrad_mod.o \
          callkeys.h.idone \
          dimradmars_mod.done \
          dust_param_mod.done \
          tracer_mod.done \
          updaterad.done
	touch $(FCM_DONEDIR)/$@

updatereffrad_mod.o: \
          $(PPSRCDIR0__phys)/updatereffrad_mod.f \
          FFLAGS__phys__updatereffrad_mod.flags \
          callkeys.h \
          dimradmars_mod.o \
          dust_param_mod.o \
          tracer_mod.o \
          updaterad.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vdifc_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vdifc_mod.done: \
          vdifc_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          compute_dtau_mod.done \
          dust_param_mod.done \
          hdo_surfex_mod.done \
          microphys.h.idone \
          surfdat_h.done \
          tracer_mod.done \
          turb_mod.done \
          watersat_mod.done
	touch $(FCM_DONEDIR)/$@

vdifc_mod.o: \
          $(PPSRCDIR0__phys)/vdifc_mod.f \
          FFLAGS__phys__vdifc_mod.flags \
          callkeys.h \
          comcstfi_h.o \
          compute_dtau_mod.o \
          dust_param_mod.o \
          hdo_surfex_mod.o \
          microphys.h \
          surfdat_h.o \
          tracer_mod.o \
          turb_mod.o \
          watersat_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vdif_cd.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vdif_cd.done: \
          vdif_cd.o \
          comcstfi_h.done \
          turb_mod.done
	touch $(FCM_DONEDIR)/$@

vdif_cd.o: \
          $(PPSRCDIR0__phys)/vdif_cd.f \
          FFLAGS__phys__vdif_cd.flags \
          comcstfi_h.o \
          turb_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vdif_kc.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vdif_kc.done: \
          vdif_kc.o \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

vdif_kc.o: \
          $(PPSRCDIR0__phys)/vdif_kc.f \
          FFLAGS__phys__vdif_kc.flags \
          tracer_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__vlz_fi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

vlz_fi.done: \
          vlz_fi.o
	touch $(FCM_DONEDIR)/$@

vlz_fi.o: \
          $(PPSRCDIR0__phys)/vlz_fi.f \
          FFLAGS__phys__vlz_fi.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__watercloud_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

watercloud_mod.done: \
          watercloud_mod.o \
          callkeys.h.idone \
          dimradmars_mod.done \
          improvedclouds_mod.done \
          ioipsl_getin_p_mod.done \
          tracer_mod.done \
          updaterad.done \
          watersat_mod.done
	touch $(FCM_DONEDIR)/$@

watercloud_mod.o: \
          $(PPSRCDIR0__phys)/watercloud_mod.f \
          FFLAGS__phys__watercloud_mod.flags \
          callkeys.h \
          dimradmars_mod.o \
          improvedclouds_mod.o \
          ioipsl_getin_p_mod.o \
          tracer_mod.o \
          updaterad.o \
          watersat_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__watersat_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

watersat_mod.done: \
          watersat_mod.o
	touch $(FCM_DONEDIR)/$@

watersat_mod.o: \
          $(PPSRCDIR0__phys)/watersat_mod.f \
          FFLAGS__phys__watersat_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writediagfi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writediagfi.done: \
          writediagfi.o \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          surfdat_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writediagfi.o: \
          $(PPSRCDIR0__phys)/writediagfi.f \
          FFLAGS__phys__writediagfi.flags \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          surfdat_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writediagmicrofi.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writediagmicrofi.done: \
          writediagmicrofi.o \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          surfdat_h.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writediagmicrofi.o: \
          $(PPSRCDIR0__phys)/writediagmicrofi.f \
          FFLAGS__phys__writediagmicrofi.flags \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          surfdat_h.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writediagsoil.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writediagsoil.done: \
          writediagsoil.o \
          comsoil_h.done \
          geometry_mod.done \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writediagsoil.o: \
          $(PPSRCDIR0__phys)/writediagsoil.f90 \
          FFLAGS__phys__writediagsoil.flags \
          comsoil_h.o \
          geometry_mod.o \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__writeg1d.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

writeg1d.done: \
          writeg1d.o \
          time_phylmdz_mod.done
	touch $(FCM_DONEDIR)/$@

writeg1d.o: \
          $(PPSRCDIR0__phys)/writeg1d.f \
          FFLAGS__phys__writeg1d.flags \
          time_phylmdz_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__wstats.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

wstats.done: \
          wstats.o \
          mod_grid_phy_lmdz.done \
          mod_phys_lmdz_para.done \
          statto_mod.done
	touch $(FCM_DONEDIR)/$@

wstats.o: \
          $(PPSRCDIR0__phys)/wstats.f90 \
          FFLAGS__phys__wstats.flags \
          mod_grid_phy_lmdz.o \
          mod_phys_lmdz_para.o \
          statto_mod.o
	fcm_internal compile:F phys $< $@

FFLAGS__phys__xios_output_mod.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

xios_output_mod.done: \
          xios_output_mod.o
	touch $(FCM_DONEDIR)/$@

xios_output_mod.o: \
          $(PPSRCDIR0__phys)/xios_output_mod.f90 \
          FFLAGS__phys__xios_output_mod.flags
	fcm_internal compile:F phys $< $@

FFLAGS__phys__yamada4.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

yamada4.done: \
          yamada4.o \
          callkeys.h.idone \
          tracer_mod.done \
          turb_mod.done
	touch $(FCM_DONEDIR)/$@

yamada4.o: \
          $(PPSRCDIR0__phys)/yamada4.f \
          FFLAGS__phys__yamada4.flags \
          callkeys.h \
          tracer_mod.o \
          turb_mod.o
	fcm_internal compile:F phys $< $@

yoegwd.h: \
          $(SRCDIR0__phys)/yoegwd.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

yoegwd.h.idone: \
          $(SRCDIR0__phys)/yoegwd.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__phys__yomlw_h.flags: \
          FFLAGS__phys.flags
	touch $(FCM_FLAGSDIR)/$@

yomlw_h.done: \
          yomlw_h.o \
          dimradmars_mod.done
	touch $(FCM_DONEDIR)/$@

yomlw_h.o: \
          $(PPSRCDIR0__phys)/yomlw_h.f90 \
          FFLAGS__phys__yomlw_h.flags \
          dimradmars_mod.o
	fcm_internal compile:F phys $< $@

