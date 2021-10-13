# Automatic Make rule for rebuild

SRCDIR0__rebuild = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/IOIPSL/rebuild

rebuild.etc : \
          $(SRCDIR0__rebuild)/AA_make \
          $(SRCDIR0__rebuild)/AA_make.ldef \
          $(FCM_DONEDIR)/FCM_CP.dummy
	cp $^ $(FCM_ETCDIR)
	touch $(FCM_DONEDIR)/$@

FFLAGS__rebuild__flio_rbld.flags: \
          FFLAGS__rebuild.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__rebuild__flio_rbld.flags: \
          LDFLAGS__rebuild.flags
	touch $(FCM_FLAGSDIR)/$@

LD__rebuild__flio_rbld.flags: \
          LD__rebuild.flags
	touch $(FCM_FLAGSDIR)/$@

flio_rbld.exe: \
          flio_rbld.o \
          LD__rebuild__flio_rbld.flags \
          LDFLAGS__rebuild__flio_rbld.flags \
          $(OBJECTS) \
          ioipsl.done \
          defprec.done
	fcm_internal load rebuild $< $@

flio_rbld.o: \
          $(SRCDIR0__rebuild)/flio_rbld.f90 \
          FFLAGS__rebuild__flio_rbld.flags \
          ioipsl.o \
          defprec.o
	fcm_internal compile:F rebuild $< $@

rebuild: \
          $(SRCDIR0__rebuild)/rebuild
	cp $< $(FCM_BINDIR)
	chmod u+w $(FCM_BINDIR)/$@

