










c=======================================================================
c     Program for simulate the impact of the CO2 snow fall on
c     the surface infrared emission (emissivity)
c-----------------------------------------------------------------------
c     Algorithme:
c     1 - Initialization
c     2 - Compute the surface emissivity
c-----------------------------------------------------------------------
c     Author: F. Forget 
c-----------------------------------------------------------------------
c     Reference paper:
c     Francois Forget and James B. Pollack
c
c     "Thermal infrared observations of the condensing Martian polar
c     caps: CO2 ice temperatures and radiative budget"
c
c     JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 101, NO. E7, PAGES 16,
c     865-16,879, JULY 25, 1996
c=======================================================================

      SUBROUTINE co2snow(ngrid, nlayer, ptimestep, emisref, condsub,
     &                   pplev, pcondicea, pcondices, pfallice,
     &                   pemisurf)
      
      use surfdat_h, only: iceradius, dtemisice
      use geometry_mod, only: latitude ! grid point latitudes (rad)
      use time_phylmdz_mod, only: daysec

      IMPLICIT NONE

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

c=======================================================================
c     Variables
c=======================================================================
c     Inputs:
c     -------
      integer, intent(in) ::
     &   ngrid, ! number of atmospheric columns
     &   nlayer ! number of atmospheric layers
                        
      real, intent(in) ::
     &   ptimestep,               ! timestep of the physics (s)
     &   emisref(ngrid),          ! grd or ice emissivity without snow 
     &   pplev(ngrid,nlayer+1),   ! interlayer pressure (Pa)
     &   pcondicea(ngrid,nlayer), ! CO2 condensation rate in layers
                                  !   (kg/m2/s)
     &   pcondices(ngrid),        ! CO2 condensation rate on the surface
                                  !   (kg/m2/s)        
     &   pfallice(ngrid,nlayer+1) ! falling CO2 ice (kg/m2/s)

      logical, intent(in) ::
     &   condsub(ngrid) ! true if there is CO2 condensation or
                        ! sublimation in the column

c     Output:
c     -------
      real, intent(out) ::
     &   pemisurf(ngrid) ! surface emissivity

c     local:
c     ------
      integer ::
     &   l,
     &   ig,
     &   icap

      real ::
     &   sumdaer,
     &   zdemisurf

      real, parameter ::
     &   alpha = 0.45

c     saved:
c     ------
      real, save ::
     &   Kscat(2) ! Kscat: coefficient for decreasing the surface 
c                 !        emissivity
c                 ! 
c                 ! Kscat = (0.001/3.) * alpha / iceradius
c                 !          with 0.3 < alpha < 0.6
c                 !
c                 ! alpha set to 0.45 (coeff from emis = f (tau)) and
c                 ! iceradius the mean radius of the scaterring
c                 ! particles (200.e-6 < iceradius < 10.e-6)
c                 !
c                 ! (2) = N and S hemispheres

      logical, save ::
     &   firstcall = .true.
c=======================================================================
c BEGIN
c=======================================================================
c 1 - Initialization
c=======================================================================
      if (firstcall) then
        Kscat(1) = (0.001/3.) * alpha / iceradius(1)
        Kscat(2) = (0.001/3.) * alpha / iceradius(2)
                   
        firstcall = .false.
      end if
c=======================================================================
c 2 - Compute the surface emissivity
c=======================================================================
      do ig = 1, ngrid
        if (condsub(ig)) then

          if (latitude(ig).lt.0.) then
              icap = 2 ! Southern hemisphere
          else
              icap = 1 ! Northern Hemisphere
          end if

          ! compute zdemisurf using an integrated form for numerical
          ! safety instead
          zdemisurf =
     &     (emisref(ig) - pemisurf(ig)) / (dtemisice(icap) * daysec)
     &     +
     &     ( emisref(ig) * 
     &              ( (pemisurf(ig)/emisref(ig))**(-3) +
     &                3.*Kscat(icap)*pfallice(ig,1)*ptimestep )**(-1/3.)
     &       - pemisurf(ig)
     &      ) / ptimestep
          pemisurf(ig) = pemisurf(ig) + zdemisurf * ptimestep 
          if (pemisurf(ig).lt.0.1) then
            write(*,*)'ds co2snow: emis < 0.1 !!!'
            write(*,*)'ig =' , ig
            write(*,*)'pemisurf(ig)', pemisurf(ig)
            write(*,*)'zdemisurf*ptimestep', zdemisurf*ptimestep
          end if
        else ! if condsub(ig) is false
           pemisurf(ig) = emisref(ig)
        end if
      end do ! ngrid
c=======================================================================
c END
c=======================================================================
      return
      end 
