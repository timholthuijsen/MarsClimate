










      MODULE hdo_surfex_mod

      IMPLICIT NONE

      CONTAINS

      subroutine hdo_surfex(ngrid,nlay,nq,ptimestep,
     &                      zt,zq,pqsurf,
     &                      old_h2o_vap,pdqsdif,dwatercap_dif,
     &                      hdoflux)

      use tracer_mod, only: igcm_h2o_vap, igcm_h2o_ice,
     &                      igcm_hdo_vap, igcm_hdo_ice,
     &                      qperemin
      use surfdat_h, only: watercaptag
      use geometry_mod, only: longitude_deg,latitude_deg

      implicit none
c------------------------------------------------------------------
c               Routine to compute the fluxes between air and surface
c               for HDO, based of the fluxes for H2O
c           L. Rossi.; M. Vals 2019
c------------------------------------------------------------------
      include "callkeys.h"

c------------------------------------------------------------------
c     Arguments:
c     ---------
c     Inputs:
      INTEGER, INTENT(IN) :: ngrid,nlay
      INTEGER, INTENT(IN) :: nq                 ! nombre de traceurs

      REAL, INTENT(IN) :: ptimestep             ! pas de temps physique (s)
      REAL, INTENT(IN) :: zt(ngrid,nlay)       ! local value of temperature
      REAL, INTENT(IN) :: zq(ngrid,nlay,nq)    ! local value of tracers
      REAL, INTENT(IN) :: pqsurf(ngrid,nq)
      REAL, INTENT(IN) :: old_h2o_vap(ngrid)     ! traceur d'eau avant
                                           !traitement de l'eau (kg/kg)
      REAL, INTENT(IN) :: dwatercap_dif(ngrid)  ! trend related to permanent ice
      REAL, INTENT(INOUT) :: pdqsdif(ngrid,nq)    ! tendance towards surface 
                                 !   (kg/kg.s-1)

c     Output:
      REAL, INTENT(OUT) :: hdoflux(ngrid)       ! value of vapour flux of HDO

c------------------------------------------------------------------
c     Local variables:

      REAL alpha_c(ngrid)  ! fractionation factor
      REAL extrasublim ! sublimation in excess of surface ice
      REAL tmpratio(ngrid)   ! D/H ratio in flux to surf
      REAL h2oflux(ngrid)       ! value of vapour flux of H2O
                                      ! same sign as pdqsdif

      INTEGER ig,l

      REAL DoH_vap(ngrid)

c-----------------------------------------------------------------------
c    Calculation of the fluxes for HDO
        
        alpha_c(1:ngrid)=0.

        DO ig=1,ngrid
              
            h2oflux(ig) = pdqsdif(ig,igcm_h2o_ice) + 
     &          dwatercap_dif(ig)

            !! IF Sublimation
            if (h2oflux(ig).le.0.) then

               if (pqsurf(ig,igcm_h2o_ice).gt.qperemin) then
                pdqsdif(ig,igcm_hdo_ice) =
     &            pdqsdif(ig,igcm_h2o_ice)*
     &             (pqsurf(ig,igcm_hdo_ice)/
     &             pqsurf(ig,igcm_h2o_ice) )
               else
                pdqsdif(ig,igcm_hdo_ice) = 0.
               endif

                pdqsdif(ig,igcm_hdo_ice)=
     &             max(pdqsdif(ig,igcm_hdo_ice),
     &            -pqsurf(ig,igcm_hdo_ice)/ptimestep)

                hdoflux(ig) = pdqsdif(ig,igcm_hdo_ice)

             if(watercaptag(ig)) then

              !if we sublimate more than qsurf
              if ((-h2oflux(ig)*ptimestep)
     &           .gt.pqsurf(ig,igcm_h2o_ice)) then

C               dwatercap_dif is how much we sublimate in excess of
C               pqsurf for H2O                        
C               hdoflux(ig) is the flux of HDO from atm. to surf.
c               The D/H of the old ice is supposed to be 5 SMOW
c               We need D/H of the flux to be 5, so we need
c               dwatercap_dif* 5 * 2 * 155.76e-6 (=1 SMOW)
                    hdoflux(ig)= hdoflux(ig)
     &                   +(dwatercap_dif(ig)*(2.*155.76e-6)*5.)
                endif
             endif ! watercap

            else ! condensation

               if (hdofrac) then !do we use fractionation?
c               alpha_c(ig) = exp(16288./zt(ig,1)**2.-9.34e-2)
                alpha_c = exp(13525./zt(ig,1)**2.-5.59e-2) !Lamb
               else
                alpha_c(ig) = 1.
               endif

               if (old_h2o_vap(ig).gt.qperemin) then
                         pdqsdif(ig,igcm_hdo_ice)=
     &                      alpha_c(ig)*pdqsdif(ig,igcm_h2o_ice)*
     &                      (zq(ig,1,igcm_hdo_vap)/
     &                          old_h2o_vap(ig))
                else
                 pdqsdif(ig,igcm_hdo_ice)= 0.
                endif

               if (hdofrac) then !do we use fractionation?
                pdqsdif(ig,igcm_hdo_ice)= 
     &      min(  pdqsdif(ig,igcm_hdo_ice), 
     &           (zq(ig,1,igcm_hdo_vap)/ptimestep) )
               endif

                hdoflux(ig)=pdqsdif(ig,igcm_hdo_ice)

            endif !sublimation

        ENDDO ! of DO ig=1,ngrid

c           CALL WRITEDIAGFI(ngrid,'extrasublim',
c    &                       'extrasublimation',
c    &                       ' ',2,tmpratio)
c           CALL WRITEDIAGFI(ngrid,'alpha_c_s',
c    &                       'alpha_c_s',
c    &                       ' ',2,alpha_c)

       return
      end subroutine hdo_surfex
c------------------------------------------------------------------

      end module hdo_surfex_mod
