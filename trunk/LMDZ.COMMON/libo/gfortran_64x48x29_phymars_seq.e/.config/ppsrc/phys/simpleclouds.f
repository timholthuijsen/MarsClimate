










      subroutine simpleclouds(ngrid,nlay,ptimestep,
     &             pplay,pzlay,pt,pdt,
     &             pq,pdq,pdqcloud,pdtcloud,
     &             nq,tau,rice)
      USE updaterad
      USE watersat_mod, ONLY: watersat
      use tracer_mod, only: igcm_h2o_vap, igcm_h2o_ice,
     &                      igcm_hdo_vap, igcm_hdo_ice,
     &                      qperemin
      USE comcstfi_h
      use dimradmars_mod, only: naerkind

      implicit none
c------------------------------------------------------------------
c  This routine is used to form clouds when a parcel of the GCM is
c    saturated. It is a simplified scheme, and there is almost no
c    microphysics involved. When the air is saturated, water-ice
c    clouds form on a fraction of the dust particles, specified by
c    the constant called "ccn_factor". There is no supersaturation,
c    and no nucleation rates computed. A more accurate scheme can
c    be found in the routine called "improvedclouds.F".

c  Modif de zq si saturation dans l'atmosphere
c  si zq(ig,l)> zqsat(ig,l) ->    zq(ig,l)=zqsat(ig,l)
c  Le test est effectue de bas en haut. L'eau condensee
c    (si saturation) est remise dans la couche en dessous.
c  L'eau condensee dans la couche du bas est deposee a la surface

c  Authors: Franck Montmessin (water ice scheme)
c           Francois Forget (changed nuclei density & outputs)
c           Ehouarn Millour (sept.2008, tracers are now handled
c                                   by name and not fixed index)
c           J.-B. Madeleine (developed a single routine called
c                            simpleclouds.F, and corrected calculations
c                            of the typical CCN profile, Oct. 2011)
c------------------------------------------------------------------
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

c------------------------------------------------------------------
c     Arguments:
c     ---------
c     Inputs:
      INTEGER ngrid,nlay
      integer nq                 ! nombre de traceurs
      REAL ptimestep             ! pas de temps physique (s)
      REAL pplay(ngrid,nlay)     ! pression au milieu des couches (Pa)
      REAL pzlay(ngrid,nlay)     ! altitude at the middle of the layers
      REAL pt(ngrid,nlay)        ! temperature at the middle of the
                                 !   layers (K)
      REAL pdt(ngrid,nlay)       ! tendance temperature des autres
                                 !   param.
      real pq(ngrid,nlay,nq)     ! traceur (kg/kg)
      real pdq(ngrid,nlay,nq)    ! tendance avant condensation
                                 !   (kg/kg.s-1)
      REAL tau(ngrid,naerkind)   ! Column dust optical depth at each point

c     Output:
      REAL rice(ngrid,nlay)      ! Ice mass mean radius (m)
                                 ! (r_c in montmessin_2004)
      real pdqcloud(ngrid,nlay,nq) ! tendance de la condensation
                                   !   H2O(kg/kg.s-1)
      REAL pdtcloud(ngrid,nlay)    ! tendance temperature due
                                   !   a la chaleur latente

c------------------------------------------------------------------
c     Local variables:

      REAL rhocloud(ngrid,nlay)  ! Cloud density (kg.m-3)

      INTEGER ig,l

      REAL zq(ngrid,nlay,nq)    ! local value of tracers
      REAL zq0(ngrid,nlay,nq)   ! local initial value of tracers
      REAL zt(ngrid,nlay)       ! local value of temperature
      REAL zqsat(ngrid,nlay)    ! saturation
      REAL*8 dzq                      ! masse de glace echangee (kg/kg)
      REAL lw                         !Latent heat of sublimation (J.kg-1) 
      REAL,PARAMETER :: To=273.15     ! reference temperature, T=273.15 K
      real rdusttyp(ngrid,nlay) ! Typical dust geom. mean radius (m)
      REAL ccntyp(ngrid,nlay)
                                      ! Typical dust number density (#/kg)
      REAL alpha_c(ngrid,nlay) !!MARGAUX: alpha_c as a spatial variable

c     CCN reduction factor
c      REAL, PARAMETER :: ccn_factor = 4.5  !! comme TESTS_JB // 1. avant
      REAL DoH_vap(ngrid,nlay)
      REAL DoH_ice(ngrid,nlay)
c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

c    On "update" la valeur de q(nq) (water vapor) et temperature.
c    On effectue qqes calculs preliminaires sur les couches : 

      do l=1,nlay
        do ig=1,ngrid
          zq(ig,l,igcm_h2o_vap)=
     &      pq(ig,l,igcm_h2o_vap)+pdq(ig,l,igcm_h2o_vap)*ptimestep
          zq(ig,l,igcm_h2o_vap)=max(zq(ig,l,igcm_h2o_vap),1.E-30) ! FF 12/2004 
          zq0(ig,l,igcm_h2o_vap)=zq(ig,l,igcm_h2o_vap)
          zt(ig,l)=pt(ig,l)+ pdt(ig,l)*ptimestep

          zq(ig,l,igcm_h2o_ice)=
     &      pq(ig,l,igcm_h2o_ice)+pdq(ig,l,igcm_h2o_ice)*ptimestep
          zq(ig,l,igcm_h2o_ice)=max(zq(ig,l,igcm_h2o_ice),0.) ! FF 12/2004 
          zq0(ig,l,igcm_h2o_ice)=zq(ig,l,igcm_h2o_ice)

          if (hdo) then
          zq(ig,l,igcm_hdo_vap)=
     &      pq(ig,l,igcm_hdo_vap)+pdq(ig,l,igcm_hdo_vap)*ptimestep
          zq(ig,l,igcm_hdo_vap)=max(zq(ig,l,igcm_hdo_vap),1e-30) ! FF 12/2004 
          zq0(ig,l,igcm_hdo_vap)=zq(ig,l,igcm_hdo_vap)

          zq(ig,l,igcm_hdo_ice)=
     &      pq(ig,l,igcm_hdo_ice)+pdq(ig,l,igcm_hdo_ice)*ptimestep
          zq(ig,l,igcm_hdo_ice)=max(zq(ig,l,igcm_hdo_ice),1e-30) ! FF 12/2004 
          zq0(ig,l,igcm_hdo_ice)=zq(ig,l,igcm_hdo_ice)

          endif !hdo
        enddo
      enddo

      pdqcloud(1:ngrid,1:nlay,1:nq)=0
      pdtcloud(1:ngrid,1:nlay)=0

      alpha_c(1:ngrid,1:nlay)=0.

c     ----------------------------------------------
c
c     Rapport de melange a saturation dans la couche l : -------
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call watersat(ngrid*nlay,zt,pplay,zqsat)

c     taux de condensation (kg/kg/s-1) dans les differentes couches
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do l=1,nlay
        do ig=1,ngrid

          if (zq(ig,l,igcm_h2o_vap).ge.zqsat(ig,l))then  !  Condensation
            dzq=zq(ig,l,igcm_h2o_vap)-zqsat(ig,l)

          elseif(zq(ig,l,igcm_h2o_vap).lt.zqsat(ig,l))then  ! Sublimation
            dzq=-min(zqsat(ig,l)-zq(ig,l,igcm_h2o_vap),
     &               zq(ig,l,igcm_h2o_ice))
          endif
            
c         Water Mass change
c         ~~~~~~~~~~~~~~~~~
          zq(ig,l,igcm_h2o_ice)=zq(ig,l,igcm_h2o_ice)+dzq
          zq(ig,l,igcm_h2o_vap)=zq(ig,l,igcm_h2o_vap)-dzq
        
        enddo ! of do ig=1,ngrid
      enddo ! of do l=1,nlay


c     Tendance finale
c     ~~~~~~~~~~~~~~~
      do l=1, nlay
        do ig=1,ngrid
          pdqcloud(ig,l,igcm_h2o_vap)=(zq(ig,l,igcm_h2o_vap)
     &                            -zq0(ig,l,igcm_h2o_vap))/ptimestep
          pdqcloud(ig,l,igcm_h2o_ice) =
     &      (zq(ig,l,igcm_h2o_ice) - zq0(ig,l,igcm_h2o_ice))/ptimestep

          lw=(2834.3-0.28*(zt(ig,l)-To)-0.004*(zt(ig,l)-To)**2)*1.e+3
          pdtcloud(ig,l)=-pdqcloud(ig,l,igcm_h2o_vap)*lw/cpp

          if (hdo) then
            if (pdqcloud(ig,l,igcm_h2o_ice).gt.0.0) then !condens

                if (hdofrac) then ! do we use fractionation?
c               alpha_c(ig,l) = exp(16288./zt(ig,l)**2.-9.34d-2)
                alpha_c(ig,l) = exp(13525./zt(ig,l)**2.-5.59d-2)  !Lamb
                else
                alpha_c(ig,l) = 1.d0
                endif
               
                if (zq0(ig,l,igcm_h2o_vap).gt.qperemin) then
                  pdqcloud(ig,l,igcm_hdo_ice)= 
     &              pdqcloud(ig,l,igcm_h2o_ice)*alpha_c(ig,l)*
     &         ( zq0(ig,l,igcm_hdo_vap)
     &                 /zq0(ig,l,igcm_h2o_vap) )
                else
                   pdqcloud(ig,l,igcm_hdo_ice)= 0.0
                endif

                pdqcloud(ig,l,igcm_hdo_ice) =
     &            min(pdqcloud(ig,l,igcm_hdo_ice),
     &              zq0(ig,l,igcm_hdo_vap)/ptimestep)

                pdqcloud(ig,l,igcm_hdo_vap)= 
     &               -pdqcloud(ig,l,igcm_hdo_ice)       

          else  ! sublimation

             if (zq0(ig,l,igcm_h2o_ice).gt.qperemin) then 
                pdqcloud(ig,l,igcm_hdo_ice)= 
     &               pdqcloud(ig,l,igcm_h2o_ice)*
     &      ( zq0(ig,l,igcm_hdo_ice)
     &              /zq0(ig,l,igcm_h2o_ice) )
             else
                pdqcloud(ig,l,igcm_hdo_ice)= 0.
             endif

              pdqcloud(ig,l,igcm_hdo_ice) =
     &          max(pdqcloud(ig,l,igcm_hdo_ice),
     &            -zq0(ig,l,igcm_hdo_ice)/ptimestep)

              pdqcloud(ig,l,igcm_hdo_vap)= 
     &             -pdqcloud(ig,l,igcm_hdo_ice)        

            endif ! condensation/sublimation

          endif ! hdo

        enddo ! of do ig=1,ngrid
      enddo ! of do l=1,nlay

c     ice crystal radius
      do l=1, nlay
        do ig=1,ngrid
          call updaterice_typ(zq(ig,l,igcm_h2o_ice),
     &       tau(ig,1),pzlay(ig,l),rice(ig,l))
        end do
      end do

c     if (hdo) then
c           CALL WRITEDIAGFI(ngrid,'alpha_c',
c    &                       'alpha_c',
c    &                       ' ',3,alpha_c) 
c     endif !hdo
c------------------------------------------------------------------
      return
      end
