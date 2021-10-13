










      subroutine thermosphere(ngrid,nlayer,nq,
     &     pplev,pplay,dist_sol,
     $     mu0,ptimestep,ptime,zday,tsurf,zzlev,zzlay,
     &     pt,pq,pu,pv,pdt,pdq,
     $     zdteuv,zdtconduc,zdumolvis,zdvmolvis,zdqmoldiff,
     $     PhiEscH,PhiEscH2,PhiEscD)

      use conc_mod, only: rnew, cpnew
      USE comcstfi_h, only: r, cpp
      implicit none

!
! For Fortran 77/Fortran 90 compliance always use line continuation
! symbols '&' in columns 73 and 6
!
! NB: to keep commons aligned, it is better to split them in groups
!     of given types (logical, integer, real, ...)

      COMMON/callkeys_l/callrad,calldifv,calladj,callcond,callsoil      &
     &   ,season,diurnal,lwrite,calllott,callstats,calleofdump          &
     &   ,callnirco2,callnlte,callthermos,callconduct,calleuv           &
     &   ,callmolvis,callmoldiff,thermochem,thermoswater,callemis       &
     &   ,callg2d,linear,rayleigh,tracer                                &
     &   ,scavenging,sedimentation                                      &
     &   ,activice,water,tifeedback,microphys,supersat,caps,photochem   &
     &   ,calltherm,callrichsl,callslope,tituscap,callyamada4,co2clouds &
     &   ,co2useh2o,meteo_flux,activeco2ice,CLFvaryingCO2,spantCO2      &
     &   ,CLFvarying,satindexco2,rdstorm,slpwind,calllott_nonoro        &
     &   ,latentheat_surfwater,gwd_convective_source,startphy_file      &
     &   ,hdo,hdofrac,cap_albedo,temp_dependant_m
     
      COMMON/callkeys_i/iradia,iaervar,ilwd,ilwb,ilwn,ncouche           &
     &   ,nltemodel,nircorr,solvarmod,solvaryear,dustinjection
     
      COMMON/callkeys_r/semi,alphan,euveff,                             &
     &   tke_heat_flux,dustrefir,fixed_euv_value,CLFfixval,             &
     &   coeff_injection,ti_injection,tf_injection,coeff_detrainment
     
      LOGICAL callrad,calldifv,calladj,callcond,callsoil,               &
     &   season,diurnal,lwrite,calllott,calllott_nonoro                 &
     &   ,callstats,calleofdump                                         &
     &   ,callnirco2,callnlte,callthermos,callconduct,                  &
     &    calleuv,callmolvis,callmoldiff,thermochem,thermoswater        &
     &   ,calltherm,callrichsl,callslope,tituscap,callyamada4

      COMMON/aeroutput/dustiropacity

      logical startphy_file

      logical callemis
      logical callg2d
      logical linear
      logical gwd_convective_source

      real semi
      real alphan
      real fixed_euv_value
      real euveff
      real tke_heat_flux
      real coeff_injection ! dust injection scheme coefficient
      real ti_injection ! local time of beginning injection
      real tf_injection ! local time of end injection
      real coeff_detrainment ! rocket dust detrainment coefficient
      real CLFfixval

      integer iaervar
      integer iradia
      integer ilwd
      integer ilwb
      integer ilwn
      integer ncouche
      integer solvarmod   ! model for solar EUV variation
      integer solvaryear  ! mars year for realisticly varying solar EUV 
      integer dustinjection ! dust injection scheme number 

      logical rayleigh
      logical tracer
      logical scavenging
      logical rdstorm ! rocket dust storm parametrization
      logical slpwind ! entrainment by slope wind parametrization
      logical latentheat_surfwater ! latent heat release from ground water ice sublimation/condensation
      logical cap_albedo ! polar cap albedo remains unchanged by water frost deposition
      logical temp_dependant_m ! temperature-dependant water contact parameter
      logical sedimentation
      logical activice,tifeedback,supersat,caps
      logical co2clouds,co2useh2o,meteo_flux,CLFvaryingCO2,satindexco2
      logical activeco2ice
      integer spantCO2
      logical CLFvarying
      logical water
      logical hdo
      logical hdofrac
      logical microphys
      logical photochem
      integer nltemodel
      integer nircorr

      character(len=100) dustiropacity
      real               dustrefir 
  
      integer swrtype ! type of short wave (solar wavelength) radiative
      ! transfer to use 1: Fouquart 2: Toon.
      parameter (swrtype=2)
!      parameter (swrtype=2)

      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of advected tracers
      REAL,INTENT(in) :: pplay(ngrid,nlayer)
      REAL,INTENT(in) :: pplev(ngrid,nlayer+1)
      REAL,INTENT(in) :: zzlay(ngrid,nlayer)
      REAL,INTENT(in) :: zzlev(ngrid,nlayer+1)
      REAL,INTENT(in) :: pt(ngrid,nlayer)
      REAL,INTENT(in) :: zday
      REAL,INTENT(in) :: dist_sol
      REAL,INTENT(in) :: mu0(ngrid)
      REAL,INTENT(in) :: pq(ngrid,nlayer,nq)
      REAL,INTENT(in) :: ptimestep
      REAL,INTENT(in) :: ptime
      REAL,INTENT(in) :: tsurf(ngrid)
      REAL,INTENT(in) :: pu(ngrid,nlayer),pv(ngrid,nlayer)
      REAL,INTENT(in) :: pdt(ngrid,nlayer),pdq(ngrid,nlayer,nq)

      REAL,INTENT(out) :: zdteuv(ngrid,nlayer)
      REAL,INTENT(out) :: zdtconduc(ngrid,nlayer)
      REAL,INTENT(out) :: zdumolvis(ngrid,nlayer)
      REAL,INTENT(out) :: zdvmolvis(ngrid,nlayer)
      REAL,INTENT(out) :: zdqmoldiff(ngrid,nlayer,nq)
      REAL*8,INTENT(out) :: PhiEscH,PhiEscH2,PhiEscD

      INTEGER :: l,ig
      logical,save :: firstcall=.true.

      if (firstcall) then
        if (.not. tracer) then
          do l=1,nlayer
            do ig=1,ngrid
              rnew(ig,l)=r
              cpnew(ig,l)=cpp
            enddo
          enddo
        endif
        firstcall= .false.
      endif

      ! initialize tendencies to zero in all cases
      ! (tendencies are added later on, even if parametrization is not called)
      zdteuv(1:ngrid,1:nlayer)=0
      zdtconduc(1:ngrid,1:nlayer)=0
      zdumolvis(1:ngrid,1:nlayer)=0
      zdvmolvis(1:ngrid,1:nlayer)=0
      zdqmoldiff(1:ngrid,1:nlayer,1:nq)=0
      
      if (calleuv) then
        call euvheat(ngrid,nlayer,nq,pt,pdt,pplev,pplay,zzlay,
     $               mu0,ptimestep,ptime,zday,pq,pdq,zdteuv)
      endif

      if (callconduct) THEN
        call conduction(ngrid,nlayer,ptimestep,pplay,pplev,pt,zdteuv,
     $                   tsurf,zzlev,zzlay,zdtconduc)
      endif

      if (callmolvis) THEN
        call molvis(ngrid,nlayer,ptimestep,pplay,pplev,pt,
     &                zdteuv,zdtconduc,pu,
     $                   tsurf,zzlev,zzlay,zdumolvis)
        call molvis(ngrid,nlayer,ptimestep,pplay,pplev,pt,
     &                zdteuv,zdtconduc,pv,
     $                   tsurf,zzlev,zzlay,zdvmolvis)
      endif

      if (callmoldiff) THEN
        call moldiff_red(ngrid,nlayer,nq,
     &                   pplay,pplev,pt,pdt,pq,pdq,ptimestep,
     &                   zzlay,zdteuv,zdtconduc,zdqmoldiff,
     &                   PhiEscH,PhiEscH2,PhiEscD)
      endif

      end


