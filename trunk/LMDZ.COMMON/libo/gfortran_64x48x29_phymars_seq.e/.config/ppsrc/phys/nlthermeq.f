










      subroutine nlthermeq(ngrid, nlayer, pplev, pplay)
c
c  Compute the number of layers nlaylte (stored in module yomlw_h)
c  over which local thermodynamic equilibrium
c  radiation scheme should be run to be sure of covering at least to a
c  height greater than (pressure lower than) p=pminte, set in nlteparams.h.
c  The maximum layer needed is found for the worst possible case.
c  Stephen Lewis 6/2000
c  Modified Y. Wanherdrick/ F. Forget 09/2000
      use yomlw_h, only: nlaylte
      implicit none
c-----------------------------------------------------------------------
c   INCLUDE 'nlteparams.h'
c
c   Parameters which govern the transition from LTE to NLTE radiation
c   tendencies.
c-----------------------------------------------------------------------

      real ptrans           ! central pressure for transition (Pa)
      parameter (ptrans = 0.1)
      real zw               ! half-width for transition (scale heights)
      parameter (zw = 0.5)
      real pminte           ! pressure up to which LTE is calculated (Pa)
      parameter (pminte = 0.4*ptrans)
c     almost one scale height above transition in worst case is very safe
      real zwi
      parameter (zwi = 2./zw)

c-----------------------------------------------------------------------
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

c
c     Input:
      integer ngrid, nlayer
      real pplev(ngrid, nlayer+1)
      real pplay(ngrid, nlayer)
c
c     Local:
      integer igpmax, ismax
      logical firstcall
      data firstcall /.true./
      save firstcall, igpmax
c
      if(firstcall) then
c     Find the location of maximum surface pressure.
c     Location won't vary much so only do it at the start;
c     with no topography location would vary, but this is only
c     needed for an estimate so any point would do in that case.
!!    AS: can be problem w MESOSCALE nesting (ignored for the moment)
         igpmax = ismax(ngrid, pplev, 1)
         write(*, 10) ptrans
         write(*, 20) zw
         write(*, 30) pminte
         firstcall = .false.
      endif
c
      IF(callnlte) then
c       Find first layer above pminte at this location
        do nlaylte = nlayer, 1, -1
           if (pplay(igpmax, nlaylte).gt.pminte)  go to 100
        enddo
      ELSE
        nlaylte=nlayer
      END IF
  100 write(*,*) 'LTE rad. calculations up to layer ',  nlaylte
c
      return
c
   10 format(' nlthermeq: transition to NLTE centred at ',f6.2,'Pa')
   20 format('               half-width (scale heights) ',f6.2)
   30 format('          suggested LTE coverage at least ',f6.2,'Pa')
      end
