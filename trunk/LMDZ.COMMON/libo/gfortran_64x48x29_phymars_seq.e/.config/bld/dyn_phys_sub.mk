# Automatic Make rule for dyn_phys_sub

PPSRCDIR0__dyn_phys_sub = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libo/gfortran_64x48x29_phymars_seq.e/.config/ppsrc/dyn_phys_sub

SRCDIR0__dyn_phys_sub = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libf/dynphy_lonlat/phymars

FFLAGS__dyn_phys_sub__avg_horiz_mod.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

avg_horiz_mod.done: \
          avg_horiz_mod.o
	touch $(FCM_DONEDIR)/$@

avg_horiz_mod.o: \
          $(PPSRCDIR0__dyn_phys_sub)/avg_horiz_mod.f \
          FFLAGS__dyn_phys_sub__avg_horiz_mod.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__caldyn0.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

caldyn0.done: \
          caldyn0.o \
          comvert_mod.done
	touch $(FCM_DONEDIR)/$@

caldyn0.o: \
          $(PPSRCDIR0__dyn_phys_sub)/caldyn0.f \
          FFLAGS__dyn_phys_sub__caldyn0.flags \
          comvert_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__callphysiq_mod.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

callphysiq_mod.done: \
          callphysiq_mod.o \
          control_mod.done \
          mod_grid_phy_lmdz.done \
          physiq_mod.done
	touch $(FCM_DONEDIR)/$@

callphysiq_mod.o: \
          $(PPSRCDIR0__dyn_phys_sub)/callphysiq_mod.f90 \
          FFLAGS__dyn_phys_sub__callphysiq_mod.flags \
          control_mod.o \
          mod_grid_phy_lmdz.o \
          physiq_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__datareadnc.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

datareadnc.done: \
          datareadnc.o \
          avg_horiz_mod.done \
          comconst_mod.done \
          comgeom.h.idone \
          datafile_mod.done \
          dimensions.h.idone \
          ioipsl_getincom.done \
          mvc_horiz_mod.done \
          paramet.h.idone
	touch $(FCM_DONEDIR)/$@

datareadnc.o: \
          $(PPSRCDIR0__dyn_phys_sub)/datareadnc.f \
          FFLAGS__dyn_phys_sub__datareadnc.flags \
          avg_horiz_mod.o \
          comconst_mod.o \
          comgeom.h \
          datafile_mod.o \
          dimensions.h \
          ioipsl_getincom.o \
          mvc_horiz_mod.o \
          paramet.h
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__defrun_new.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

defrun_new.done: \
          defrun_new.o \
          control_mod.done \
          ioipsl_getincom.done \
          logic_mod.done \
          serre_mod.done \
          sponge_mod.done
	touch $(FCM_DONEDIR)/$@

defrun_new.o: \
          $(PPSRCDIR0__dyn_phys_sub)/defrun_new.f \
          FFLAGS__dyn_phys_sub__defrun_new.flags \
          control_mod.o \
          ioipsl_getincom.o \
          logic_mod.o \
          serre_mod.o \
          sponge_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__grid_noro1.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

grid_noro1.done: \
          grid_noro1.o \
          comconst_mod.done
	touch $(FCM_DONEDIR)/$@

grid_noro1.o: \
          $(PPSRCDIR0__dyn_phys_sub)/grid_noro1.f \
          FFLAGS__dyn_phys_sub__grid_noro1.flags \
          comconst_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__iniphysiq_mod.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

iniphysiq_mod.done: \
          iniphysiq_mod.o \
          comgeomfi_h.done \
          dimphy.done \
          geometry_mod.done \
          infotrac.done \
          inigeomphy_mod.done \
          iniprint.h.idone \
          mod_phys_lmdz_para.done \
          phys_state_var_init_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

iniphysiq_mod.o: \
          $(PPSRCDIR0__dyn_phys_sub)/iniphysiq_mod.f90 \
          FFLAGS__dyn_phys_sub__iniphysiq_mod.flags \
          comgeomfi_h.o \
          dimphy.o \
          geometry_mod.o \
          infotrac.o \
          inigeomphy_mod.o \
          iniprint.h \
          mod_phys_lmdz_para.o \
          phys_state_var_init_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__ini_archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

ini_archive.done: \
          ini_archive.o \
          comconst_mod.done \
          comsoil_h.done \
          comvert_mod.done \
          ener_mod.done \
          logic_mod.done \
          serre_mod.done
	touch $(FCM_DONEDIR)/$@

ini_archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/ini_archive.f \
          FFLAGS__dyn_phys_sub__ini_archive.flags \
          comconst_mod.o \
          comsoil_h.o \
          comvert_mod.o \
          ener_mod.o \
          logic_mod.o \
          serre_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__interp_vert.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

interp_vert.done: \
          interp_vert.o
	touch $(FCM_DONEDIR)/$@

interp_vert.o: \
          $(PPSRCDIR0__dyn_phys_sub)/interp_vert.f \
          FFLAGS__dyn_phys_sub__interp_vert.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__lect_start_archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

lect_start_archive.done: \
          lect_start_archive.o \
          comconst_mod.done \
          comsoil_h.done \
          comvert_mod.done \
          infotrac.done \
          planete_h.done
	touch $(FCM_DONEDIR)/$@

lect_start_archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/lect_start_archive.f \
          FFLAGS__dyn_phys_sub__lect_start_archive.flags \
          comconst_mod.o \
          comsoil_h.o \
          comvert_mod.o \
          infotrac.o \
          planete_h.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__mvc_horiz_mod.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

mvc_horiz_mod.done: \
          mvc_horiz_mod.o
	touch $(FCM_DONEDIR)/$@

mvc_horiz_mod.o: \
          $(PPSRCDIR0__dyn_phys_sub)/mvc_horiz_mod.f \
          FFLAGS__dyn_phys_sub__mvc_horiz_mod.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__newstart.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__dyn_phys_sub__newstart.flags: \
          LDFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LD__dyn_phys_sub__newstart.flags: \
          LD__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

newstart.o: \
          $(PPSRCDIR0__dyn_phys_sub)/newstart.f \
          FFLAGS__dyn_phys_sub__newstart.flags \
          clesph0.h \
          co2cloud_mod.o \
          comconst_mod.o \
          comdissnew.h \
          comgeom2.h \
          comsoil_h.o \
          comvert_mod.o \
          control_mod.o \
          datafile_mod.o \
          dimensions.h \
          dimradmars_mod.o \
          dust_param_mod.o \
          ener_mod.o \
          exner_hyb_m.o \
          filtreg_mod.o \
          geometry_mod.o \
          infotrac.o \
          iniphysiq_mod.o \
          ioipsl_getincom.o \
          iostart.o \
          mod_const_mpi.o \
          mod_phys_lmdz_para.o \
          paramet.h \
          phyetat0_mod.o \
          phyredem.o \
          serre_mod.o \
          surfdat_h.o \
          temps_mod.o \
          tracer_mod.o \
          turb_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

newstart_64x48x29_phymars_seq.e: \
          newstart.o \
          LD__dyn_phys_sub__newstart.flags \
          LDFLAGS__dyn_phys_sub__newstart.flags \
          $(OBJECTS) \
          clesph0.h.idone \
          co2cloud_mod.done \
          comconst_mod.done \
          comdissnew.h.idone \
          comgeom2.h.idone \
          comsoil_h.done \
          comvert_mod.done \
          control_mod.done \
          datafile_mod.done \
          dimensions.h.idone \
          dimradmars_mod.done \
          dust_param_mod.done \
          ener_mod.done \
          exner_hyb_m.done \
          filtreg_mod.done \
          geometry_mod.done \
          infotrac.done \
          iniphysiq_mod.done \
          ioipsl_getincom.done \
          iostart.done \
          mod_const_mpi.done \
          mod_phys_lmdz_para.done \
          paramet.h.idone \
          phyetat0_mod.done \
          phyredem.done \
          serre_mod.done \
          surfdat_h.done \
          temps_mod.done \
          tracer_mod.done \
          turb_mod.done
	fcm_internal load dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__readhead_NC.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

readhead_nc.done: \
          readhead_nc.o \
          comconst_mod.done \
          comvert_mod.done \
          ener_mod.done \
          temps_mod.done
	touch $(FCM_DONEDIR)/$@

readhead_nc.o: \
          $(PPSRCDIR0__dyn_phys_sub)/readhead_NC.f \
          FFLAGS__dyn_phys_sub__readhead_NC.flags \
          comconst_mod.o \
          comvert_mod.o \
          ener_mod.o \
          temps_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__scal_wind.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

scal_wind.done: \
          scal_wind.o
	touch $(FCM_DONEDIR)/$@

scal_wind.o: \
          $(PPSRCDIR0__dyn_phys_sub)/scal_wind.f \
          FFLAGS__dyn_phys_sub__scal_wind.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__start2archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__dyn_phys_sub__start2archive.flags: \
          LDFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LD__dyn_phys_sub__start2archive.flags: \
          LD__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

start2archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/start2archive.f \
          FFLAGS__dyn_phys_sub__start2archive.flags \
          comconst_mod.o \
          comdissip.h \
          comgeom.h \
          comsoil_h.o \
          comvert_mod.o \
          control_mod.o \
          dimensions.h \
          exner_hyb_m.o \
          filtreg_mod.o \
          infotrac.o \
          iniphysiq_mod.o \
          mod_const_mpi.o \
          paramet.h \
          phyetat0_mod.o \
          surfdat_h.o \
          temps_mod.o
	fcm_internal compile:F dyn_phys_sub $< $@

start2archive_64x48x29_phymars_seq.e: \
          start2archive.o \
          LD__dyn_phys_sub__start2archive.flags \
          LDFLAGS__dyn_phys_sub__start2archive.flags \
          $(OBJECTS) \
          comconst_mod.done \
          comdissip.h.idone \
          comgeom.h.idone \
          comsoil_h.done \
          comvert_mod.done \
          control_mod.done \
          dimensions.h.idone \
          exner_hyb_m.done \
          filtreg_mod.done \
          infotrac.done \
          iniphysiq_mod.done \
          mod_const_mpi.done \
          paramet.h.idone \
          phyetat0_mod.done \
          surfdat_h.done \
          temps_mod.done
	fcm_internal load dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__wind_scal.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

wind_scal.done: \
          wind_scal.o
	touch $(FCM_DONEDIR)/$@

wind_scal.o: \
          $(PPSRCDIR0__dyn_phys_sub)/wind_scal.f \
          FFLAGS__dyn_phys_sub__wind_scal.flags
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__write_archive.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

write_archive.done: \
          write_archive.o \
          comsoil_h.done
	touch $(FCM_DONEDIR)/$@

write_archive.o: \
          $(PPSRCDIR0__dyn_phys_sub)/write_archive.f \
          FFLAGS__dyn_phys_sub__write_archive.flags \
          comsoil_h.o
	fcm_internal compile:F dyn_phys_sub $< $@

FFLAGS__dyn_phys_sub__xvik.flags: \
          FFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__dyn_phys_sub__xvik.flags: \
          LDFLAGS__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

LD__dyn_phys_sub__xvik.flags: \
          LD__dyn_phys_sub.flags
	touch $(FCM_FLAGSDIR)/$@

xvik.o: \
          $(PPSRCDIR0__dyn_phys_sub)/xvik.f \
          FFLAGS__dyn_phys_sub__xvik.flags \
          comconst_mod.o \
          comdissip.h \
          comgeom2.h \
          dimensions.h \
          filtreg_mod.o \
          paramet.h
	fcm_internal compile:F dyn_phys_sub $< $@

xvik_64x48x29_phymars_seq.e: \
          xvik.o \
          LD__dyn_phys_sub__xvik.flags \
          LDFLAGS__dyn_phys_sub__xvik.flags \
          $(OBJECTS) \
          comconst_mod.done \
          comdissip.h.idone \
          comgeom2.h.idone \
          dimensions.h.idone \
          filtreg_mod.done \
          paramet.h.idone
	fcm_internal load dyn_phys_sub $< $@

