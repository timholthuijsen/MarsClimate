      SUBROUTINE conf_phys(ngrid,nlayer,nq)
 
!=======================================================================
!
!   purpose:
!   -------
!
!   Initialisation for the physical parametrisations of the LMD 
!   martian atmospheric general circulation modele.
!
!   author: Frederic Hourdin 15 / 10 /93
!   -------
!   modified: Sebastien Lebonnois 11/06/2003 (new callphys.def)
!             Ehouarn Millour (oct. 2008) tracers are now identified
!              by their names and may not be contiguously
!              stored in the q(:,:,:,:) array
!             E.M. (june 2009) use getin routine to load parameters
!             adapted to the mesoscale use - Aymeric Spiga - 01/2007-07/2011
!             separated inifis into conf_phys and phys_state_var_init (A. Spiga)
!
!
!   arguments:
!   ----------
!
!   input:
!   ------
!
!    nq                    Number of tracers
!
!=======================================================================
!
!-----------------------------------------------------------------------
!   declarations:
!   -------------
      USE ioipsl_getin_p_mod, ONLY : getin_p
      use tracer_mod, only : nuice_sed, ccn_factor, nuiceco2_sed,
     &                       nuice_ref,nuiceco2_ref
      use surfdat_h, only: albedo_h2o_cap,albedo_h2o_frost,
     &                     frost_albedo_threshold, inert_h2o_ice
      use time_phylmdz_mod, only: ecritphy,day_step,iphysiq,ecritstart,
     &                            daysec,dtphys
      use dimradmars_mod, only: naerkind, name_iaer,
     &                      ini_scatterers,tauvis
      use datafile_mod, only: datadir
      use calchim_mod, only: ichemistry
      use co2condens_mod, only: scavco2cond
      use dust_param_mod, only: dustbin, doubleq, submicron, active,
     &                          lifting, freedust, callddevil,
     &                          dustscaling_mode
      use aeropacity_mod, only: iddist, topdustref
      IMPLICIT NONE
      include "callkeys.h"
      include "microphys.h"

      INTEGER,INTENT(IN) :: ngrid,nlayer,nq
      INTEGER ierr,j
      character(len=20),parameter :: modname="conf_phys"
 
      CHARACTER ch1*12
#ifndef MESOSCALE
      ! read in some parameters from "run.def" for physics,
      ! or shared between dynamics and physics.
      ecritphy=240 ! default value
      call getin_p("ecritphy",ecritphy) ! frequency of outputs in physics,
                                      ! in dynamical steps
      day_step=960 ! default value
      call getin_p("day_step",day_step) ! number of dynamical steps per day
      iphysiq=20 ! default value
      call getin_p("iphysiq",iphysiq) ! call physics every iphysiq dyn step
      ecritstart=0 ! default value
      call getin_p("ecritstart",ecritstart) ! write a restart every ecristart steps
#endif

! --------------------------------------------------------------
!  Reading the "callphys.def" file controlling some key options
! --------------------------------------------------------------
     
      ! check that 'callphys.def' file is around
      OPEN(99,file='callphys.def',status='old',form='formatted'
     &     ,iostat=ierr)
      CLOSE(99)
      
      IF(ierr.EQ.0) THEN
         PRINT*
         PRINT*
         PRINT*,'--------------------------------------------'
         PRINT*,' conf_phys: Parameters for the physics (callphys.def)'
         PRINT*,'--------------------------------------------'

         write(*,*) "Directory where external input files are:"
         ! default path is set in datafile_mod
         call getin_p("datadir",datadir) 
         write(*,*) " datadir = ",trim(datadir)

         write(*,*) "Initialize physics with startfi.nc file ?"
         startphy_file=.true.
         call getin_p("startphy_file",startphy_file)
         write(*,*) "startphy_file", startphy_file
         
         write(*,*) "Run with or without tracer transport ?"
         tracer=.false. ! default value
         call getin_p("tracer",tracer)
         write(*,*) " tracer = ",tracer

         write(*,*) "Diurnal cycle ?"
         write(*,*) "(if diurnal=False, diurnal averaged solar heating)"
         diurnal=.true. ! default value
         call getin_p("diurnal",diurnal)
         write(*,*) " diurnal = ",diurnal

         write(*,*) "Seasonal cycle ?"
         write(*,*) "(if season=False, Ls stays constant, to value ",
     &   "set in 'start'"
         season=.true. ! default value
         call getin_p("season",season)
         write(*,*) " season = ",season

         write(*,*) "Write some extra output to the screen ?"
         lwrite=.false. ! default value
         call getin_p("lwrite",lwrite)
         write(*,*) " lwrite = ",lwrite

         write(*,*) "Save statistics in file stats.nc ?"
#ifdef MESOSCALE
         callstats=.false. ! default value
#else
         callstats=.true. ! default value
#endif
         call getin_p("callstats",callstats)
         write(*,*) " callstats = ",callstats

         write(*,*) "Save EOF profiles in file 'profiles' for ",
     &              "Climate Database?"
         calleofdump=.false. ! default value
         call getin_p("calleofdump",calleofdump)
         write(*,*) " calleofdump = ",calleofdump

         write(*,*) "Dust scenario: 1=constant dust (read from startfi",
     &   " or set as tauvis); 2=Viking scenario; =3 MGS scenario,",
     &   "=6 cold (low dust) scenario; =7 warm (high dust) scenario ",
     &   "=24,25 ... 30 :Mars Year 24, ... or 30 from TES assimilation"
         iaervar=3 ! default value
         call getin_p("iaervar",iaervar)
         write(*,*) " iaervar = ",iaervar

         write(*,*) "Reference (visible) dust opacity at 610 Pa ",
     &   "(matters only if iaervar=1)"
         ! NB: default value of tauvis is set/read in startfi.nc file
         call getin_p("tauvis",tauvis)
         write(*,*) " tauvis = ",tauvis

         write(*,*) "Dust vertical distribution:"
         write(*,*) "(=1 top set by topdustref parameter;",
     & " =2 Viking scenario; =3 MGS scenario)"
         iddist=3 ! default value
         call getin_p("iddist",iddist)
         write(*,*) " iddist = ",iddist

         write(*,*) "Dust top altitude (km). (Matters only if iddist=1)"
         topdustref= 90.0 ! default value
         call getin_p("topdustref",topdustref)
         write(*,*) " topdustref = ",topdustref

         write(*,*) "Prescribed surface thermal flux (H/(rho*cp),K m/s)"
         tke_heat_flux=0. ! default value
         call getin_p("tke_heat_flux",tke_heat_flux)
         write(*,*) " tke_heat_flux = ",tke_heat_flux
         write(*,*) " 0 means the usual schemes are computing"

         write(*,*) "call radiative transfer ?"
         callrad=.true. ! default value
         call getin_p("callrad",callrad)
         write(*,*) " callrad = ",callrad

         write(*,*) "call slope insolation scheme ?",
     &              "(matters only if callrad=T)"
#ifdef MESOSCALE
         callslope=.true. ! default value
#else
         callslope=.false. ! default value (not supported yet)
#endif
         call getin_p("callslope",callslope)
         write(*,*) " callslope = ",callslope

         write(*,*) "call NLTE radiative schemes ?",
     &              "(matters only if callrad=T)"
         callnlte=.false. ! default value
         call getin_p("callnlte",callnlte)
         write(*,*) " callnlte = ",callnlte
         
         nltemodel=0    !default value
         write(*,*) "NLTE model?"
         write(*,*) "0 -> old model, static O"
         write(*,*) "1 -> old model, dynamic O"
         write(*,*) "2 -> new model"
         write(*,*) "(matters only if callnlte=T)"
         call getin_p("nltemodel",nltemodel)
         write(*,*) " nltemodel = ",nltemodel

         write(*,*) "call CO2 NIR absorption ?",
     &              "(matters only if callrad=T)"
         callnirco2=.false. ! default value
         call getin_p("callnirco2",callnirco2)
         write(*,*) " callnirco2 = ",callnirco2

         write(*,*) "New NIR NLTE correction ?",
     $              "0-> old model (no correction)",
     $              "1-> new correction",
     $              "(matters only if callnirco2=T)"
#ifdef MESOSCALE
         nircorr=0      !default value. this is OK below 60 km.
#else
         nircorr=0      !default value
#endif
         call getin_p("nircorr",nircorr)
         write(*,*) " nircorr = ",nircorr

         write(*,*) "call turbulent vertical diffusion ?"
         calldifv=.true. ! default value
         call getin_p("calldifv",calldifv)
         write(*,*) " calldifv = ",calldifv

         write(*,*) "call thermals ?"
         calltherm=.false. ! default value
         call getin_p("calltherm",calltherm)
         write(*,*) " calltherm = ",calltherm

         write(*,*) "call convective adjustment ?"
         calladj=.true. ! default value
         call getin_p("calladj",calladj)
         write(*,*) " calladj = ",calladj

         if (calltherm .and. calladj) then
          print*,'!!! PLEASE NOTE !!!'
          print*,'convective adjustment is on'
          print*,'but since thermal plume model is on'
          print*,'convadj is only activated above the PBL'
         endif
        
         write(*,*) "used latest version of yamada scheme?"
         callyamada4=.true. ! default value
         call getin_p("callyamada4",callyamada4)
         write(*,*) " callyamada4 = ",callyamada4

         if (calltherm .and. .not.callyamada4) then
          print*,'!!!! WARNING WARNING WARNING !!!!'
          print*,'if calltherm=T we strongly advise that '
          print*,'you set the flag callyamada4 to T '
          print*,'!!!! WARNING WARNING WARNING !!!!'
         endif
 
         write(*,*) "call Richardson-based surface layer ?"
         callrichsl=.false. ! default value
         call getin_p("callrichsl",callrichsl)
         write(*,*) " callrichsl = ",callrichsl

         if (calltherm .and. .not.callrichsl) then
          print*,'WARNING WARNING WARNING'
          print*,'if calltherm=T we strongly advise that '
          print*,'you use the new surface layer scheme '
          print*,'by setting callrichsl=T '
         endif

         if (calladj .and. callrichsl .and. (.not. calltherm)) then
          print*,'You should not be calling the convective adjustment
     & scheme with the Richardson surface-layer and without the thermals
     &. This approach is not
     & physically consistent and can lead to unrealistic friction
     & values.'
          print*,'If you want to use the Ri. surface-layer, either
     & activate thermals OR de-activate the convective adjustment.'
          call abort_physic(modname,
     &     "Richardson layer must be used with thermals",1)
         endif

         write(*,*) "call CO2 condensation ?"
         callcond=.true. ! default value
         call getin_p("callcond",callcond)
         write(*,*) " callcond = ",callcond

         write(*,*)"call thermal conduction in the soil ?"
         callsoil=.true. ! default value
         call getin_p("callsoil",callsoil)
         write(*,*) " callsoil = ",callsoil
         

         write(*,*)"call Lott's gravity wave/subgrid topography ",
     &             "scheme ?"
         calllott=.true. ! default value
         call getin_p("calllott",calllott)
         write(*,*)" calllott = ",calllott

         write(*,*)"call Lott's non-oro GWs parameterisation ",
     &             "scheme ?"
         calllott_nonoro=.false. ! default value
         call getin_p("calllott_nonoro",calllott_nonoro)
         write(*,*)" calllott_nonoro = ",calllott_nonoro

! rocket dust storm injection scheme
         write(*,*)"call rocket dust storm parametrization"
         rdstorm=.false. ! default value
         call getin_p("rdstorm",rdstorm)
         write(*,*)" rdstorm = ",rdstorm
! rocket dust storm detrainment coefficient        
        coeff_detrainment=0.05 ! default value
        call getin_p("coeff_detrainment",coeff_detrainment)
        write(*,*)" coeff_detrainment = ",coeff_detrainment

! entrainment by slope wind scheme
         write(*,*)"call slope wind lifting parametrization"
         slpwind=.false. ! default value
         call getin_p("slpwind",slpwind)
         write(*,*)" slpwind = ",slpwind

! latent heat release from ground water ice sublimation/condensation
         write(*,*)"latent heat release during sublimation", 
     &              " /condensation of ground water ice"
         latentheat_surfwater=.true. ! default value
         call getin_p("latentheat_surfwater",latentheat_surfwater)
         write(*,*)" latentheat_surfwater = ",latentheat_surfwater

         write(*,*)"rad.transfer is computed every iradia",
     &             " physical timestep"
         iradia=1 ! default value
         call getin_p("iradia",iradia)
         write(*,*)" iradia = ",iradia
         

         write(*,*)"Output of the exchange coefficient mattrix ?",
     &             "(for diagnostics only)"
         callg2d=.false. ! default value
         call getin_p("callg2d",callg2d)
         write(*,*)" callg2d = ",callg2d

         write(*,*)"Rayleigh scattering : (should be .false. for now)"
         rayleigh=.false.
         call getin_p("rayleigh",rayleigh)
         write(*,*)" rayleigh = ",rayleigh


! TRACERS:

! dustbin
         write(*,*)"Transported dust ? (if >0, use 'dustbin' dust bins)"
         dustbin=0 ! default value
         call getin_p("dustbin",dustbin)
         write(*,*)" dustbin = ",dustbin
! active
         write(*,*)"Radiatively active dust ? (matters if dustbin>0)"
         active=.false. ! default value
         call getin_p("active",active)
         write(*,*)" active = ",active

! Test of incompatibility:
! if active is used, then dustbin should be > 0

         if (active.and.(dustbin.lt.1)) then
           print*,'if active is used, then dustbin should > 0'
           call abort_physic(modname,
     &          "active option requires dustbin < 0",1)
         endif
! doubleq
         write(*,*)"use mass and number mixing ratios to predict",
     &             " dust size ?"
         doubleq=.false. ! default value
         call getin_p("doubleq",doubleq)
         write(*,*)" doubleq = ",doubleq
! submicron
         submicron=.false. ! default value
         call getin_p("submicron",submicron)
         write(*,*)" submicron = ",submicron

! Test of incompatibility:
! if doubleq is used, then dustbin should be 2

         if (doubleq.and.(dustbin.ne.2)) then
           print*,'if doubleq is used, then dustbin should be 2'
           call abort_physic(modname,
     &          "doubleq option requires dustbin = 2",1)
         endif
         if (doubleq.and.submicron.and.(nq.LT.3)) then
           print*,'If doubleq is used with a submicron tracer,'
           print*,' then the number of tracers has to be'
           print*,' larger than 3.'
           call abort_physic(modname,
     &          "submicron option requires dustbin > 2",1)
         endif

! lifting
         write(*,*)"dust lifted by GCM surface winds ?"
         lifting=.false. ! default value
         call getin_p("lifting",lifting)
         write(*,*)" lifting = ",lifting

! Test of incompatibility:
! if lifting is used, then dustbin should be > 0

         if (lifting.and.(dustbin.lt.1)) then
           print*,'if lifting is used, then dustbin should > 0'
           call abort_physic(modname,
     &          "lifting option requires dustbin > 0",1)
         endif

! dust injection scheme
        dustinjection=0 ! default: no injection scheme
        call getin_p("dustinjection",dustinjection)
        write(*,*)" dustinjection = ",dustinjection
! dust injection scheme coefficient        
        coeff_injection=1. ! default value
        call getin_p("coeff_injection",coeff_injection)
        write(*,*)" coeff_in,jection = ",coeff_injection
! timing for dust injection        
        ti_injection=10. ! default value
        tf_injection=12. ! default value
        call getin_p("ti_injection",ti_injection)
        write(*,*)" ti_injection = ",ti_injection
        call getin_p("tf_injection",tf_injection)
        write(*,*)" tf_injection = ",tf_injection

! free evolving dust
! freedust=true just says that there is no lifting and no dust opacity scaling.
         write(*,*)"dust lifted by GCM surface winds ?"
         freedust=.false. ! default value
         call getin_p("freedust",freedust)
         write(*,*)" freedust = ",freedust
         if (freedust.and..not.doubleq) then
           print*,'freedust should be used with doubleq !'
           call abort_physic(modname,
     &          "freedust option requires doubleq",1)
         endif

! dust rescaling mode (if any)
         if (freedust) then
           dustscaling_mode=0
         else
           dustscaling_mode=1 ! GCMv5.3 style
         endif
         call getin_p("dustscaling_mode",dustscaling_mode)
         write(*,*) "dustscaling_mode=",dustscaling_mode

#ifndef MESOSCALE
         ! this test is valid in GCM case
         ! ... not in mesoscale case, for we want to activate mesoscale lifting
         if (freedust.and.dustinjection.eq.0)then
           if(lifting) then
             print*,'if freedust is used and dustinjection = 0, 
     &      then lifting should not be used'
             call abort_physic(modname,
     &          "freedust option with dustinjection = 0"//
     &          " requires lifting to be false",1)
           endif
         endif
#endif
         if (dustinjection.eq.1)then
           if(.not.lifting) then
             print*,"if dustinjection=1, then lifting should be true"
             call abort_physic(modname,
     &          "dustinjection=1 requires lifting",1)
           endif
           if(.not.freedust) then
             print*,"if dustinjection=1, then freedust should be true"
             call abort_physic(modname,
     &          "dustinjection=1 requires freedust",1)
           endif
         endif
! rocket dust storm and entrainment by slope wind
! Test of incompatibility:
! if rdstorm or slpwind is used, then doubleq should be true
         if ((rdstorm.or.slpwind).and..not.doubleq) then
           print*,'if rdstorm or slpwind is used, then doubleq 
     &            should be used !'
           call abort_physic(modname,
     &          "rdstorm or slpwind requires doubleq",1)
         endif
         if ((rdstorm.or.slpwind).and..not.active) then
           print*,'if rdstorm or slpwind is used, then active 
     &            should be used !'
           call abort_physic(modname,
     &          "rdstorm or slpwind requires activ",1)
         endif
         if (rdstorm.and..not.lifting) then
           print*,'if rdstorm is used, then lifting 
     &            should be used !'
           call abort_physic(modname,
     &          "rdstorm requires lifting",1)
         endif
         if ((rdstorm.or.slpwind).and..not.freedust) then
           print*,'if rdstorm or slpwind is used, then freedust 
     &            should be used !'
           call abort_physic(modname,
     &          "rdstorm or slpwind requires freedust",1)
         endif
         if (rdstorm.and.(dustinjection.eq.0)) then
           print*,'if rdstorm is used, then dustinjection
     &            should be used !'
           call abort_physic(modname,
     &          "rdstorm requires dustinjection",1)
         endif
! Dust IR opacity
         write(*,*)" Wavelength for infrared opacity of dust ?"
         write(*,*)" Choices are:"
         write(*,*)" tes  --- > 9.3 microns  [default]"
         write(*,*)" mcs  --- > 21.6 microns"
         !
         ! WARNING WARNING WARNING WARNING WARNING WARNING
         !
         ! BEFORE ADDING A NEW VALUE, BE SURE THAT THE
         ! CORRESPONDING WAVELENGTH IS IN THE LOOKUP TABLE,
         ! OR AT LEAST NO TO FAR, TO AVOID FALLACIOUS INTERPOLATIONS.
         !
         dustiropacity="tes" !default value - is expected to shift to mcs one day
         call getin_p("dustiropacity",dustiropacity)
         write(*,*)" dustiropacity = ",trim(dustiropacity)
         select case (trim(dustiropacity))
           case ("tes")
             dustrefir = 9.3E-6
           case ("mcs")
             dustrefir = 21.6E-6
           case default
              write(*,*) trim(dustiropacity),
     &                  " is not a valid option for dustiropacity"
             call abort_physic(modname,
     &          "invalid dustiropacity option value",1)
         end select

! callddevil
         write(*,*)" dust lifted by dust devils ?"
         callddevil=.false. !default value
         call getin_p("callddevil",callddevil)
         write(*,*)" callddevil = ",callddevil

! Test of incompatibility:
! if dustdevil is used, then dustbin should be > 0

         if (callddevil.and.(dustbin.lt.1)) then
           print*,'if dustdevil is used, then dustbin should > 0'
           call abort_physic(modname,
     &          "callddevil requires dustbin > 0",1)
         endif
! sedimentation
         write(*,*) "Gravitationnal sedimentation ?"
         sedimentation=.true. ! default value
         call getin_p("sedimentation",sedimentation)
         write(*,*) " sedimentation = ",sedimentation
! activice
         write(*,*) "Radiatively active transported atmospheric ",
     &              "water ice ?"
         activice=.false. ! default value
         call getin_p("activice",activice)
         write(*,*) " activice = ",activice
! water
         write(*,*) "Compute water cycle ?"
         water=.false. ! default value
         call getin_p("water",water)
         write(*,*) " water = ",water
! hdo
         write(*,*) "Compute hdo cycle ?"
         hdo=.false. ! default value
         call getin_p("hdo",hdo)
         write(*,*) " hdo = ",hdo

         write(*,*) "Use fractionation for hdo?"
         hdofrac=.true. ! default value
         call getin_p("hdofrac",hdofrac)
         write(*,*) " hdofrac = ",hdofrac

! Activeco2ice
         write(*,*) "Radiatively active transported atmospheric ",
     &              "Co2 ice ?"
         activeco2ice=.false. ! default value
         call getin_p("activeco2ice",activeco2ice)
         write(*,*) " activeco2ice = ",activeco2ice

! sub-grid cloud fraction: fixed
         write(*,*) "Fixed cloud fraction?"
         CLFfixval=1.0 ! default value
         call getin_p("CLFfixval",CLFfixval)
         write(*,*) "CLFfixval=",CLFfixval
! sub-grid cloud fraction: varying
         write(*,*) "Use partial nebulosity?"
         CLFvarying=.false. ! default value
         call getin_p("CLFvarying",CLFvarying)
         write(*,*)"CLFvarying=",CLFvarying

!CO2 clouds scheme?
         write(*,*) "Compute CO2 clouds (implies microphysical scheme)?"
         co2clouds=.false. ! default value
         call getin_p("co2clouds",co2clouds)
         write(*,*) " co2clouds = ",co2clouds
!Can water ice particles serve as CCN for CO2clouds
         write(*,*) "Use water ice as CO2 clouds CCN ?"
         co2useh2o=.false. ! default value
         call getin_p("co2useh2o",co2useh2o)
         write(*,*) " co2useh2o = ",co2useh2o
!Do we allow a supply of meteoritic paricles to serve as CO2 ice CCN?
         write(*,*) "Supply meteoritic particle for CO2 clouds ?"
         meteo_flux=.false. !Default value
         call getin_p("meteo_flux",meteo_flux)
         write(*,*)  " meteo_flux = ",meteo_flux
!Do we allow a sub-grid temperature distribution for the CO2 microphysics
         write(*,*) "sub-grid temperature distribution for CO2 clouds?"
         CLFvaryingCO2=.false. !Default value
         call getin_p("CLFvaryingCO2",CLFvaryingCO2)
         write(*,*)  " CLFvaryingCO2 = ",CLFvaryingCO2
!Amplitude of the sub-grid temperature distribution for the CO2 microphysics
         write(*,*) "sub-grid temperature amplitude for CO2 clouds?"
         spantCO2=0 !Default value
         call getin_p("spantCO2",spantCO2)
         write(*,*)  " spantCO2 = ",spantCO2
!Do you want to filter the sub-grid T distribution by a Saturation index?
         write(*,*) "filter sub-grid temperature by Saturation index?"
         satindexco2=.true.
         call getin_p("satindexco2",satindexco2)
         write(*,*)  " satindexco2 = ",satindexco2


! thermal inertia feedback
         write(*,*) "Activate the thermal inertia feedback ?"
         tifeedback=.false. ! default value
         call getin_p("tifeedback",tifeedback)
         write(*,*) " tifeedback = ",tifeedback

! Test of incompatibility:

         if (tifeedback.and..not.water) then
           print*,'if tifeedback is used,'
           print*,'water should be used too'
           call abort_physic(modname,
     &          "tifeedback requires water",1)
         endif

         if (tifeedback.and..not.callsoil) then
           print*,'if tifeedback is used,'
           print*,'callsoil should be used too'
           call abort_physic(modname,
     &          "tifeedback requires callsoil",1)
         endif

         if (activice.and..not.water) then
           print*,'if activice is used, water should be used too'
           call abort_physic(modname,
     &          "activeice requires water",1)
         endif

         if (water.and..not.tracer) then
           print*,'if water is used, tracer should be used too'
           call abort_physic(modname,
     &          "water requires tracer",1)
         endif

         if (hdo.and..not.water) then
           print*,'if hdo is used, water should be used too'
           call abort_physic(modname,
     &          "hd2 requires tracer",1)
         endif

         
         if (activeco2ice.and..not.co2clouds) then
          print*,'if activeco2ice is used, co2clouds should be used too'
          call abort_physic(modname,
     &          "activeco2ice requires co2clouds",1)
         endif

! water ice clouds effective variance distribution for sedimentaion       
        write(*,*) "Sed effective variance for water ice clouds ?"
        nuice_sed=0.45 
        call getin_p("nuice_sed",nuice_sed)
        write(*,*) "water_param nueff Sedimentation:", nuice_sed
              
        write(*,*) "Sed effective variance for CO2 clouds ?"
        nuiceco2_sed=0.45 
        call getin_p("nuiceco2_sed",nuiceco2_sed)
        write(*,*) "CO2 nueff Sedimentation:", nuiceco2_sed
  
        write(*,*) "REF effective variance for CO2 clouds ?"
        nuiceco2_ref=0.45 
        call getin_p("nuiceco2_ref",nuiceco2_ref)
        write(*,*) "CO2 nueff Sedimentation:", nuiceco2_ref

        write(*,*) "REF effective variance for water clouds ?"
        nuice_ref=0.45 
        call getin_p("nuice_ref",nuice_ref)
        write(*,*) "CO2 nueff Sedimentation:", nuice_ref


! ccn factor if no scavenging         
        write(*,*) "water param CCN reduc. factor ?"
        ccn_factor = 4.5
        call getin_p("ccn_factor",ccn_factor)
        write(*,*)" ccn_factor = ",ccn_factor
        write(*,*)"Careful: only used when microphys=F, otherwise"
        write(*,*)"the contact parameter is used instead;"

       ! microphys
        write(*,*)"Microphysical scheme for water-ice clouds?"
        microphys=.false.       ! default value
        call getin_p("microphys",microphys)
        write(*,*)" microphys = ",microphys

      ! supersat
        write(*,*)"Allow super-saturation of water vapor?"
        supersat=.true.         ! default value
        call getin_p("supersat",supersat)
        write(*,*)"supersat = ",supersat

! microphysical parameter contact       
        write(*,*) "water contact parameter ?"
        mteta  = 0.95 ! default value
        temp_dependant_m  = .false. ! default value
        call getin_p("temp_dependant_m",temp_dependant_m)
        if (temp_dependant_m) then
           print*,'You have chosen a temperature-dependant water'
           print*,'contact parameter ! From Maattanen et al. 2014'
        else if (.not.temp_dependant_m) then
           print*,'Water contact parameter is constant'
           call getin_p("mteta",mteta)
           write(*,*) "mteta = ", mteta
        endif
        
! scavenging
        write(*,*)"Dust scavenging by H2O/CO2 snowfall ?"
        scavenging=.false.      ! default value
        call getin_p("scavenging",scavenging)
        write(*,*)" scavenging = ",scavenging
         

! Test of incompatibility:
! if scavenging is used, then dustbin should be > 0

        if ((microphys.and..not.doubleq).or.
     &       (microphys.and..not.water)) then
           print*,'if microphys is used, then doubleq,'
           print*,'and water must be used!'
           call abort_physic(modname,
     &          "microphys requires water and doubleq",1)
        endif
        if (microphys.and..not.scavenging) then
           print*,''
           print*,'----------------WARNING-----------------'
           print*,'microphys is used without scavenging !!!'
           print*,'----------------WARNING-----------------'
           print*,''
        endif
        
        if ((scavenging.and..not.microphys).or.
     &       (scavenging.and.(dustbin.lt.1)))then
           print*,'if scavenging is used, then microphys'
           print*,'must be used!'
           call abort_physic(modname,
     &          "scavenging requires microphys",1)
        endif

! Instantaneous scavenging by CO2
! -> expected to be replaced by scavenging with microphysics (flag scavenging) one day
        write(*,*)"Dust scavenging by instantaneous CO2 snowfall ?"
        scavco2cond=.false.      ! default value
        call getin_p("scavco2cond",scavco2cond)
        write(*,*)" scavco2cond = ",scavco2cond
! Test of incompatibility:
! if scavco2cond is used, then dustbin should be > 0
        if (scavco2cond.and.(dustbin.lt.1))then
           print*,'if scavco2cond is used, then dustbin should be > 0'
           call abort_physic(modname,
     &          "scavco2cond requires dustbin > 0",1)
        endif
! if co2clouds is used, then there is no need for scavco2cond
        if (co2clouds.and.scavco2cond) then
           print*,''
           print*,'----------------WARNING-----------------'
           print*,'     microphys scavenging is used so    '
	   print*,'        no need for scavco2cond !!!     '
           print*,'----------------WARNING-----------------'
           print*,''
	   call abort_physic(modname,
     &          "incompatible co2cloud and scavco2cond options",1)
        endif
	
! Test of incompatibility:

         write(*,*) "Permanent water caps at poles ?",
     &               " .true. is RECOMMENDED"
         write(*,*) "(with .true., North cap is a source of water ",
     &   "and South pole is a cold trap)"
         caps=.true. ! default value
         call getin_p("caps",caps)
         write(*,*) " caps = ",caps

! JN : now separated between albedo_h2o_cap and 
! albedo_h2o_frost. Retrocompatible with old
! callphys.def with albedo_h2o_ice
         write(*,*) "water ice albedo ? Old settings use ",
     &              "albedo_h2o_ice, new settings use ", 
     &              "albedo_h2o_cap and albedo_h2o_frost " 
         albedo_h2o_cap=0.35
         albedo_h2o_frost=0.35
         call getin_p("albedo_h2o_ice",albedo_h2o_cap)
         albedo_h2o_frost=albedo_h2o_cap
         call getin_p("albedo_h2o_cap",albedo_h2o_cap)
         write(*,*) " albedo_h2o_cap = ",albedo_h2o_cap
         call getin_p("albedo_h2o_frost",albedo_h2o_frost)
         write(*,*) " albedo_h2o_frost = ",albedo_h2o_frost

! Northern polar cap albedo (JN 2021)
         write(*,*)"Watercaptag albedo is unchanged by water frost",
     &              " deposition (default is false)"
         cap_albedo=.false. ! default value
         call getin_p("cap_albedo",cap_albedo)
         write(*,*)"cap_albedo = ",cap_albedo

! inert_h2o_ice
         write(*,*) "water ice thermal inertia ?"
         inert_h2o_ice=2400 ! (J.m^-2.K^-1.s^-1/2)
         call getin_p("inert_h2o_ice",inert_h2o_ice)
         write(*,*) " inert_h2o_ice = ",inert_h2o_ice
! frost_albedo_threshold
         write(*,*) "frost thickness threshold for albedo ?"
         frost_albedo_threshold=0.005 ! 5.4 mic (i.e 0.005 kg.m-2)
         call getin_p("frost_albedo_threshold",
     &    frost_albedo_threshold)
         write(*,*) " frost_albedo_threshold = ",
     &            frost_albedo_threshold

! call Titus crocus line -- DEFAULT IS NONE
         write(*,*) "Titus crocus line ?"
         tituscap=.false.  ! default value
         call getin_p("tituscap",tituscap)
         write(*,*) "tituscap",tituscap
                     
! Chemistry:
         write(*,*) "photochemistry: include chemical species"
         photochem=.false. ! default value
         call getin_p("photochem",photochem)
         write(*,*) " photochem = ",photochem
         
         write(*,*) "Compute chemistry (if photochem is .true.)",
     &   "every ichemistry physics step (default: ichemistry=1)"
         ichemistry=1
         call getin_p("ichemistry",ichemistry)
         write(*,*) " ichemistry = ",ichemistry


! SCATTERERS
         write(*,*) "how many scatterers?"
         naerkind=1 ! default value
         call getin_p("naerkind",naerkind)
         write(*,*)" naerkind = ",naerkind

! Test of incompatibility
c        Logical tests for radiatively active water-ice clouds:
         IF ( (activice.AND.(.NOT.water)).OR.
     &        (activice.AND.(naerkind.LT.2)) ) THEN
           WRITE(*,*) 'If activice is TRUE, water has to be set'
           WRITE(*,*) 'to TRUE, and "naerkind" must be at least'
           WRITE(*,*) 'equal to 2.'
           call abort_physic(modname,
     &          "radiatively active dust and water"//
     &          " require naerkind > 1",1)
         ENDIF

!------------------------------------------
!------------------------------------------
! once naerkind is known allocate arrays
! -- we do it here and not in phys_var_init
! -- because we need to know naerkind first
         CALL ini_scatterers(ngrid,nlayer)
!------------------------------------------
!------------------------------------------


c        Please name the different scatterers here ----------------
         name_iaer(1) = "dust_conrath"   !! default choice is good old Conrath profile
         IF (doubleq.AND.active) name_iaer(1) = "dust_doubleq" !! two-moment scheme

         if (nq.gt.1) then
           ! trick to avoid problems compiling with 1 tracer
           ! and picky compilers who know name_iaer(2) is out of bounds
           j=2
           IF (rdstorm.AND..NOT.activice.AND..NOT.slpwind) then
             name_iaer(j) = "stormdust_doubleq" !! storm dust two-moment scheme
             j = j+1
           END IF

           IF (rdstorm.AND.water.AND.activice.AND..NOT.slpwind) then
             name_iaer(j) = "stormdust_doubleq" 
             j = j+1
           END IF

           IF (slpwind.AND..NOT.activice.AND..NOT.rdstorm) then
             name_iaer(j) = "topdust_doubleq" !! storm dust two-moment scheme
             j = j+1
           END IF
 
           IF (slpwind.AND.water.AND.activice.AND..NOT.rdstorm) then
             name_iaer(j) =  "topdust_doubleq"
             j = j+1
           END IF

           IF (rdstorm.AND.slpwind.AND..NOT.activice) THEN 
             name_iaer(j) = "stormdust_doubleq"
             name_iaer(j+1) = "topdust_doubleq"
             j = j+2
           ENDIF

           IF (rdstorm.AND.slpwind.AND.water.AND.activice) THEN 
             name_iaer(j) = "stormdust_doubleq"
             name_iaer(j+1) = "topdust_doubleq"
             j = j+2
           ENDIF

           IF (water.AND.activice) then
            name_iaer(j) = "h2o_ice"      !! radiatively-active clouds
            j = j+1
           END IF

           IF (co2clouds.AND.activeco2ice) then
             name_iaer(j) = "co2_ice" !! radiatively-active co2 clouds
             j = j+1
           ENDIF

           IF (submicron.AND.active) then
             name_iaer(j) = "dust_submicron" !! JBM experimental stuff
             j = j+1
           ENDIF
         endif ! of if (nq.gt.1)
c        ----------------------------------------------------------

! THERMOSPHERE

         write(*,*) "call thermosphere ?"
         callthermos=.false. ! default value
         call getin_p("callthermos",callthermos)
         write(*,*) " callthermos = ",callthermos
         

         write(*,*) " water included without cycle ",
     &              "(only if water=.false.)"
         thermoswater=.false. ! default value
         call getin_p("thermoswater",thermoswater)
         write(*,*) " thermoswater = ",thermoswater

         write(*,*) "call thermal conduction ?",
     &    " (only if callthermos=.true.)"
         callconduct=.false. ! default value
         call getin_p("callconduct",callconduct)
         write(*,*) " callconduct = ",callconduct

         write(*,*) "call EUV heating ?",
     &   " (only if callthermos=.true.)"
         calleuv=.false.  ! default value
         call getin_p("calleuv",calleuv)
         write(*,*) " calleuv = ",calleuv

         write(*,*) "call molecular viscosity ?",
     &   " (only if callthermos=.true.)"
         callmolvis=.false. ! default value
         call getin_p("callmolvis",callmolvis)
         write(*,*) " callmolvis = ",callmolvis

         write(*,*) "call molecular diffusion ?",
     &   " (only if callthermos=.true.)"
         callmoldiff=.false. ! default value
         call getin_p("callmoldiff",callmoldiff)
         write(*,*) " callmoldiff = ",callmoldiff
         

         write(*,*) "call thermospheric photochemistry ?",
     &   " (only if callthermos=.true.)"
         thermochem=.false. ! default value
         call getin_p("thermochem",thermochem)
         write(*,*) " thermochem = ",thermochem

         write(*,*) "Method to include solar variability"
         write(*,*) "0-> fixed value of E10.7 (fixed_euv_value); ",
     &          "1-> daily evolution of E10.7 (for given solvaryear)"
         solvarmod=1
         call getin_p("solvarmod",solvarmod)
         write(*,*) " solvarmod = ",solvarmod

         write(*,*) "Fixed euv (for solvarmod==0) 10.7 value?"
         write(*,*) " (min=80 , ave=140, max=320)"
         fixed_euv_value=140 ! default value
         call getin_p("fixed_euv_value",fixed_euv_value)
         write(*,*) " fixed_euv_value = ",fixed_euv_value
         
         write(*,*) "Solar variability as observed for MY: "
         write(*,*) "Only if solvarmod=1"
         solvaryear=24
         call getin_p("solvaryear",solvaryear)
         write(*,*) " solvaryear = ",solvaryear

         write(*,*) "UV heating efficiency:",
     &   "measured values between 0.19 and 0.23 (Fox et al. 1996)",
     &   "lower values may be used to compensate low 15 um cooling"
         euveff=0.21 !default value
         call getin_p("euveff",euveff)
         write(*,*) " euveff = ", euveff


         if (.not.callthermos) then
           if (thermoswater) then
             print*,'if thermoswater is set, callthermos must be true'
             call abort_physic(modname,
     &          "thermoswater requires callthermos",1)
           endif          
           if (callconduct) then
             print*,'if callconduct is set, callthermos must be true'
             call abort_physic(modname,
     &          "callconduct requires callthermos",1)
           endif        
           if (calleuv) then
             print*,'if calleuv is set, callthermos must be true'
             call abort_physic(modname,
     &          "calleuv requires callthermos",1)
           endif         
           if (callmolvis) then
             print*,'if callmolvis is set, callthermos must be true'
             call abort_physic(modname,
     &          "callmolvis requires callthermos",1)
           endif        
           if (callmoldiff) then
             print*,'if callmoldiff is set, callthermos must be true'
             call abort_physic(modname,
     &          "callmoldiff requires callthermos",1)
           endif          
           if (thermochem) then
             print*,'if thermochem is set, callthermos must be true'
             call abort_physic(modname,
     &          "thermochem requires callthermos",1)
           endif          
        endif

! Test of incompatibility:
! if photochem is used, then water should be used too

         if (photochem.and..not.water) then
           print*,'if photochem is used, water should be used too'
           call abort_physic(modname,
     &          "photochem requires water",1)
         endif

! if callthermos is used, then thermoswater should be used too 
! (if water not used already)

         if (callthermos .and. .not.water) then
           if (callthermos .and. .not.thermoswater) then
             print*,'if callthermos is used, water or thermoswater 
     &               should be used too'
             call abort_physic(modname,
     &          "callthermos requires water or thermoswater",1)
           endif
         endif

         PRINT*,'--------------------------------------------'
         PRINT*
         PRINT*
      ELSE
         write(*,*)
         write(*,*) 'Cannot read file callphys.def. Is it here ?'
         call abort_physic(modname,
     &          "missing callphys.def file",1)
      ENDIF

8000  FORMAT(t5,a12,l8)
8001  FORMAT(t5,a12,i8)

      PRINT*
      PRINT*,'conf_phys: daysec',daysec
      PRINT*
      PRINT*,'conf_phys: The radiative transfer is computed:'
      PRINT*,'           each ',iradia,' physical time-step'
      PRINT*,'        or each ',iradia*dtphys,' seconds'
      PRINT*
! --------------------------------------------------------------
!  Managing the Longwave radiative transfer
! --------------------------------------------------------------

!     In most cases, the run just use the following values :
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      callemis=.true.     
!     ilwd=10*int(daysec/dtphys) ! bug before 22/10/01       
      ilwd=1
      ilwn=1 !2
      ilwb=1 !2
      linear=.true.        
      ncouche=3
      alphan=0.4
      semi=0

!     BUT people working hard on the LW may want to read them in 'radia.def' 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      OPEN(99,file='radia.def',status='old',form='formatted'
     .     ,iostat=ierr)
      IF(ierr.EQ.0) THEN
         write(*,*) 'conf_phys: Reading radia.def !!!'
         READ(99,fmt='(a)') ch1
         READ(99,*) callemis
         WRITE(*,8000) ch1,callemis

         READ(99,fmt='(a)') ch1
         READ(99,*) iradia
         WRITE(*,8001) ch1,iradia

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwd
         WRITE(*,8001) ch1,ilwd

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwn
         WRITE(*,8001) ch1,ilwn

         READ(99,fmt='(a)') ch1
         READ(99,*) linear
         WRITE(*,8000) ch1,linear

         READ(99,fmt='(a)') ch1
         READ(99,*) ncouche
         WRITE(*,8001) ch1,ncouche

         READ(99,fmt='(a)') ch1
         READ(99,*) alphan
         WRITE(*,*) ch1,alphan

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwb
         WRITE(*,8001) ch1,ilwb


         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callg2d
         WRITE(*,8000) ch1,callg2d

         READ(99,fmt='(a)') ch1
         READ(99,*) semi
         WRITE(*,*) ch1,semi
      end if
      CLOSE(99)

      END
