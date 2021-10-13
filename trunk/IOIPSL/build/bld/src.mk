# Automatic Make rule for src

SRCDIR0__src = /mnt/d/OneDrive/MarsClimate/MetJacky/NextAttempt/trunk/IOIPSL/src

src.etc : \
          $(SRCDIR0__src)/AA_make \
          $(SRCDIR0__src)/AA_make.ldef \
          $(SRCDIR0__src)/def.prec \
          $(FCM_DONEDIR)/FCM_CP.dummy
	cp $^ $(FCM_ETCDIR)
	touch $(FCM_DONEDIR)/$@

FFLAGS__src__calendar.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

calendar.done: \
          calendar.o \
          errioipsl.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

calendar.o: \
          $(SRCDIR0__src)/calendar.f90 \
          FFLAGS__src__calendar.flags \
          errioipsl.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__defprec.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

defprec.done: \
          defprec.o
	touch $(FCM_DONEDIR)/$@

defprec.o: \
          $(SRCDIR0__src)/defprec.f90 \
          FFLAGS__src__defprec.flags
	fcm_internal compile:F src $< $@

FFLAGS__src__errioipsl.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

errioipsl.done: \
          errioipsl.o
	touch $(FCM_DONEDIR)/$@

errioipsl.o: \
          $(SRCDIR0__src)/errioipsl.f90 \
          FFLAGS__src__errioipsl.flags
	fcm_internal compile:F src $< $@

FFLAGS__src__flincom.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

flincom.done: \
          flincom.o \
          calendar.done \
          errioipsl.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

flincom.o: \
          $(SRCDIR0__src)/flincom.f90 \
          FFLAGS__src__flincom.flags \
          calendar.o \
          errioipsl.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__fliocom.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

fliocom.done: \
          fliocom.o \
          calendar.done \
          defprec.done \
          errioipsl.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

fliocom.o: \
          $(SRCDIR0__src)/fliocom.f90 \
          FFLAGS__src__fliocom.flags \
          calendar.o \
          defprec.o \
          errioipsl.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__getincom.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

getincom.done: \
          getincom.o \
          errioipsl.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

getincom.o: \
          $(SRCDIR0__src)/getincom.f90 \
          FFLAGS__src__getincom.flags \
          errioipsl.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__histcom.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

histcom.done: \
          histcom.o \
          calendar.done \
          errioipsl.done \
          fliocom.done \
          mathelp.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

histcom.o: \
          $(SRCDIR0__src)/histcom.f90 \
          FFLAGS__src__histcom.flags \
          calendar.o \
          errioipsl.o \
          fliocom.o \
          mathelp.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__ioipsl.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

ioipsl.done: \
          ioipsl.o \
          calendar.done \
          errioipsl.done \
          flincom.done \
          fliocom.done \
          getincom.done \
          histcom.done \
          mathelp.done \
          restcom.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

ioipsl.o: \
          $(SRCDIR0__src)/ioipsl.f90 \
          FFLAGS__src__ioipsl.flags \
          calendar.o \
          errioipsl.o \
          flincom.o \
          fliocom.o \
          getincom.o \
          histcom.o \
          mathelp.o \
          restcom.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__mathelp.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

mathelp.done: \
          mathelp.o \
          errioipsl.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

mathelp.o: \
          $(SRCDIR0__src)/mathelp.f90 \
          FFLAGS__src__mathelp.flags \
          errioipsl.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__restcom.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

restcom.done: \
          restcom.o \
          calendar.done \
          defprec.done \
          errioipsl.done \
          fliocom.done \
          mathelp.done \
          stringop.done
	touch $(FCM_DONEDIR)/$@

restcom.o: \
          $(SRCDIR0__src)/restcom.f90 \
          FFLAGS__src__restcom.flags \
          calendar.o \
          defprec.o \
          errioipsl.o \
          fliocom.o \
          mathelp.o \
          stringop.o
	fcm_internal compile:F src $< $@

FFLAGS__src__stringop.flags: \
          FFLAGS__src.flags
	touch $(FCM_FLAGSDIR)/$@

stringop.done: \
          stringop.o
	touch $(FCM_DONEDIR)/$@

stringop.o: \
          $(SRCDIR0__src)/stringop.f90 \
          FFLAGS__src__stringop.flags
	fcm_internal compile:F src $< $@

