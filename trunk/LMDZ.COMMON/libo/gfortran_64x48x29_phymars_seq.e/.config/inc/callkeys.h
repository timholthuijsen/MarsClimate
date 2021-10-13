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
