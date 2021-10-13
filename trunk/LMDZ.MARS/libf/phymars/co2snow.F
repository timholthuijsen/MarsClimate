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

#include "callkeys.h"

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
