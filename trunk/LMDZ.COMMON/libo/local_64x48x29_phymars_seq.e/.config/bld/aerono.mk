# Automatic Make rule for aerono

PPSRCDIR0__aerono = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libo/local_64x48x29_phymars_seq.e/.config/ppsrc/aerono

SRCDIR0__aerono = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/LMDZ.COMMON/libf/aeronomars

FFLAGS__aerono__calchim_mod.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

calchim_mod.done: \
          calchim_mod.o \
          callkeys.h.idone \
          comcstfi_h.done \
          conc_mod.done \
          iono_h.done \
          photolysis_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

calchim_mod.o: \
          $(PPSRCDIR0__aerono)/calchim_mod.f90 \
          FFLAGS__aerono__calchim_mod.flags \
          callkeys.h \
          comcstfi_h.o \
          conc_mod.o \
          iono_h.o \
          photolysis_mod.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__chemthermos.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

chemthermos.done: \
          chemthermos.o \
          comcstfi_h.done \
          param_v4_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

chemthermos.o: \
          $(PPSRCDIR0__aerono)/chemthermos.f90 \
          FFLAGS__aerono__chemthermos.flags \
          comcstfi_h.o \
          param_v4_h.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__chemthermos_readini.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

chemthermos_readini.done: \
          chemthermos_readini.o \
          param_v4_h.done
	touch $(FCM_DONEDIR)/$@

chemthermos_readini.o: \
          $(PPSRCDIR0__aerono)/chemthermos_readini.f \
          FFLAGS__aerono__chemthermos_readini.flags \
          param_v4_h.o
	fcm_internal compile:F aerono $< $@

chimiedata.h: \
          $(SRCDIR0__aerono)/chimiedata.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

chimiedata.h.idone: \
          $(SRCDIR0__aerono)/chimiedata.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__aerono__concentrations.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

concentrations.done: \
          concentrations.o \
          conc_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

concentrations.o: \
          $(PPSRCDIR0__aerono)/concentrations.f \
          FFLAGS__aerono__concentrations.flags \
          conc_mod.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__conduction.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

conduction.done: \
          conduction.o \
          conc_mod.done
	touch $(FCM_DONEDIR)/$@

conduction.o: \
          $(PPSRCDIR0__aerono)/conduction.f \
          FFLAGS__aerono__conduction.flags \
          conc_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__deposition.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

deposition.done: \
          deposition.o \
          conc_mod.done \
          surfdat_h.done
	touch $(FCM_DONEDIR)/$@

deposition.o: \
          $(PPSRCDIR0__aerono)/deposition.f \
          FFLAGS__aerono__deposition.flags \
          conc_mod.o \
          surfdat_h.o
	fcm_internal compile:F aerono $< $@

diffusion.h: \
          $(SRCDIR0__aerono)/diffusion.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

diffusion.h.idone: \
          $(SRCDIR0__aerono)/diffusion.h
	touch $(FCM_DONEDIR)/$@

FFLAGS__aerono__dtridgl.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

dtridgl.done: \
          dtridgl.o
	touch $(FCM_DONEDIR)/$@

dtridgl.o: \
          $(PPSRCDIR0__aerono)/dtridgl.f \
          FFLAGS__aerono__dtridgl.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__euvheat.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

euvheat.done: \
          euvheat.o \
          conc_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

euvheat.o: \
          $(PPSRCDIR0__aerono)/euvheat.f90 \
          FFLAGS__aerono__euvheat.flags \
          conc_mod.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__hrtherm.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

hrtherm.done: \
          hrtherm.o \
          callkeys.h.idone \
          param_v4_h.done
	touch $(FCM_DONEDIR)/$@

hrtherm.o: \
          $(PPSRCDIR0__aerono)/hrtherm.f \
          FFLAGS__aerono__hrtherm.flags \
          callkeys.h \
          param_v4_h.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__inichim_newstart.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

inichim_newstart.done: \
          inichim_newstart.o \
          callkeys.h.idone \
          datafile_mod.done \
          dust_param_mod.done \
          mod_grid_phy_lmdz.done \
          tracer_mod.done \
          vertical_layers_mod.done
	touch $(FCM_DONEDIR)/$@

inichim_newstart.o: \
          $(PPSRCDIR0__aerono)/inichim_newstart.f90 \
          FFLAGS__aerono__inichim_newstart.flags \
          callkeys.h \
          datafile_mod.o \
          dust_param_mod.o \
          mod_grid_phy_lmdz.o \
          tracer_mod.o \
          vertical_layers_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__intrplf.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

intrplf.done: \
          intrplf.o
	touch $(FCM_DONEDIR)/$@

intrplf.o: \
          $(PPSRCDIR0__aerono)/intrplf.f \
          FFLAGS__aerono__intrplf.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__inv.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

inv.done: \
          inv.o
	touch $(FCM_DONEDIR)/$@

inv.o: \
          $(PPSRCDIR0__aerono)/inv.f \
          FFLAGS__aerono__inv.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__iono_h.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

iono_h.done: \
          iono_h.o
	touch $(FCM_DONEDIR)/$@

iono_h.o: \
          $(PPSRCDIR0__aerono)/iono_h.f90 \
          FFLAGS__aerono__iono_h.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__jthermcalc.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

jthermcalc.done: \
          jthermcalc.o \
          callkeys.h.idone \
          comsaison_h.done \
          param_v4_h.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

jthermcalc.o: \
          $(PPSRCDIR0__aerono)/jthermcalc.f \
          FFLAGS__aerono__jthermcalc.flags \
          callkeys.h \
          comsaison_h.o \
          param_v4_h.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__jthermcalc_e107.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

jthermcalc_e107.done: \
          jthermcalc_e107.o \
          callkeys.h.idone \
          comsaison_h.done \
          param_v4_h.done
	touch $(FCM_DONEDIR)/$@

jthermcalc_e107.o: \
          $(PPSRCDIR0__aerono)/jthermcalc_e107.f \
          FFLAGS__aerono__jthermcalc_e107.flags \
          callkeys.h \
          comsaison_h.o \
          param_v4_h.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__moldiff.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

moldiff.done: \
          moldiff.o \
          comcstfi_h.done \
          conc_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

moldiff.o: \
          $(PPSRCDIR0__aerono)/moldiff.f \
          FFLAGS__aerono__moldiff.flags \
          comcstfi_h.o \
          conc_mod.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__moldiffcoeff.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

moldiffcoeff.done: \
          moldiffcoeff.o \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

moldiffcoeff.o: \
          $(PPSRCDIR0__aerono)/moldiffcoeff.f \
          FFLAGS__aerono__moldiffcoeff.flags \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__moldiffcoeff_red.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

moldiffcoeff_red.done: \
          moldiffcoeff_red.o \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

moldiffcoeff_red.o: \
          $(PPSRCDIR0__aerono)/moldiffcoeff_red.f \
          FFLAGS__aerono__moldiffcoeff_red.flags \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__moldiff_red.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

moldiff_red.done: \
          moldiff_red.o \
          diffusion.h.idone \
          geometry_mod.done \
          planetwide_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

moldiff_red.o: \
          $(PPSRCDIR0__aerono)/moldiff_red.f90 \
          FFLAGS__aerono__moldiff_red.flags \
          diffusion.h \
          geometry_mod.o \
          planetwide_mod.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__molvis.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

molvis.done: \
          molvis.o \
          conc_mod.done
	touch $(FCM_DONEDIR)/$@

molvis.o: \
          $(PPSRCDIR0__aerono)/molvis.f \
          FFLAGS__aerono__molvis.flags \
          conc_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__paramfoto_compact.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

paramfoto_compact.done: \
          paramfoto_compact.o \
          callkeys.h.idone \
          iono_h.done \
          param_v4_h.done
	touch $(FCM_DONEDIR)/$@

paramfoto_compact.o: \
          $(PPSRCDIR0__aerono)/paramfoto_compact.f \
          FFLAGS__aerono__paramfoto_compact.flags \
          callkeys.h \
          iono_h.o \
          param_v4_h.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__param_read.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

param_read.done: \
          param_read.o \
          datafile_mod.done \
          param_v4_h.done
	touch $(FCM_DONEDIR)/$@

param_read.o: \
          $(PPSRCDIR0__aerono)/param_read.f \
          FFLAGS__aerono__param_read.flags \
          datafile_mod.o \
          param_v4_h.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__param_read_e107.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

param_read_e107.done: \
          param_read_e107.o \
          callkeys.h.idone \
          datafile_mod.done \
          param_v4_h.done
	touch $(FCM_DONEDIR)/$@

param_read_e107.o: \
          $(PPSRCDIR0__aerono)/param_read_e107.f \
          FFLAGS__aerono__param_read_e107.flags \
          callkeys.h \
          datafile_mod.o \
          param_v4_h.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__param_v4_h.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

param_v4_h.done: \
          param_v4_h.o
	touch $(FCM_DONEDIR)/$@

param_v4_h.o: \
          $(PPSRCDIR0__aerono)/param_v4_h.f90 \
          FFLAGS__aerono__param_v4_h.flags
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__perosat.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

perosat.done: \
          perosat.o \
          comcstfi_h.done \
          conc_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

perosat.o: \
          $(PPSRCDIR0__aerono)/perosat.f \
          FFLAGS__aerono__perosat.flags \
          comcstfi_h.o \
          conc_mod.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__photochemistry.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

photochemistry.done: \
          photochemistry.o \
          callkeys.h.idone \
          comcstfi_h.done \
          param_v4_h.done \
          photolysis_mod.done \
          tracer_mod.done \
          types_asis.done
	touch $(FCM_DONEDIR)/$@

photochemistry.o: \
          $(PPSRCDIR0__aerono)/photochemistry.f90 \
          FFLAGS__aerono__photochemistry.flags \
          callkeys.h \
          comcstfi_h.o \
          param_v4_h.o \
          photolysis_mod.o \
          tracer_mod.o \
          types_asis.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__photolysis.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

photolysis.done: \
          photolysis.o \
          comcstfi_h.done
	touch $(FCM_DONEDIR)/$@

photolysis.o: \
          $(PPSRCDIR0__aerono)/photolysis.f90 \
          FFLAGS__aerono__photolysis.flags \
          comcstfi_h.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__photolysis_mod.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

photolysis_mod.done: \
          photolysis_mod.o \
          datafile_mod.done
	touch $(FCM_DONEDIR)/$@

photolysis_mod.o: \
          $(PPSRCDIR0__aerono)/photolysis_mod.f90 \
          FFLAGS__aerono__photolysis_mod.flags \
          datafile_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__photolysis_online.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

photolysis_online.done: \
          photolysis_online.o \
          photolysis_mod.done
	touch $(FCM_DONEDIR)/$@

photolysis_online.o: \
          $(PPSRCDIR0__aerono)/photolysis_online.f \
          FFLAGS__aerono__photolysis_online.flags \
          photolysis_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__read_phototable.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

read_phototable.done: \
          read_phototable.o \
          chimiedata.h.idone \
          datafile_mod.done \
          ioipsl_getin_p_mod.done
	touch $(FCM_DONEDIR)/$@

read_phototable.o: \
          $(PPSRCDIR0__aerono)/read_phototable.f90 \
          FFLAGS__aerono__read_phototable.flags \
          chimiedata.h \
          datafile_mod.o \
          ioipsl_getin_p_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__surfacearea.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

surfacearea.done: \
          surfacearea.o \
          callkeys.h.idone \
          comcstfi_h.done \
          conc_mod.done \
          tracer_mod.done
	touch $(FCM_DONEDIR)/$@

surfacearea.o: \
          $(PPSRCDIR0__aerono)/surfacearea.f \
          FFLAGS__aerono__surfacearea.flags \
          callkeys.h \
          comcstfi_h.o \
          conc_mod.o \
          tracer_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__thermosphere.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

thermosphere.done: \
          thermosphere.o \
          comcstfi_h.done \
          conc_mod.done
	touch $(FCM_DONEDIR)/$@

thermosphere.o: \
          $(PPSRCDIR0__aerono)/thermosphere.f \
          FFLAGS__aerono__thermosphere.flags \
          comcstfi_h.o \
          conc_mod.o
	fcm_internal compile:F aerono $< $@

FFLAGS__aerono__types_asis.flags: \
          FFLAGS__aerono.flags
	touch $(FCM_FLAGSDIR)/$@

types_asis.done: \
          types_asis.o
	touch $(FCM_DONEDIR)/$@

types_asis.o: \
          $(PPSRCDIR0__aerono)/types_asis.f90 \
          FFLAGS__aerono__types_asis.flags
	fcm_internal compile:F aerono $< $@

