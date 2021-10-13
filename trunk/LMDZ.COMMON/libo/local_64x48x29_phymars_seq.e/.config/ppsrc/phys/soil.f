










      subroutine soil(ngrid,nsoil,firstcall,
     &          therm_i,
     &          timestep,tsurf,tsoil,
     &          capcal,fluxgrd)

      use comsoil_h, only: layer, mlayer, volcapa,
     &                     mthermdiff, thermdiff, coefq,
     &                     coefd, alph, beta, mu
      use surfdat_h, only: watercaptag, inert_h2o_ice

      implicit none

!-----------------------------------------------------------------------
!  Author: Ehouarn Millour
!
!  Purpose: Compute soil temperature using an implict 1st order scheme
!  
!  Note: depths of layers and mid-layers, soil thermal inertia and 
!        heat capacity are commons in comsoil_h
!-----------------------------------------------------------------------

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

c-----------------------------------------------------------------------
!  arguments
!  ---------
!  inputs:
      integer ngrid	! number of (horizontal) grid-points 
      integer nsoil	! number of soil layers 
      logical firstcall ! identifier for initialization call 
      real therm_i(ngrid,nsoil) ! thermal inertia
      real timestep	    ! time step
      real tsurf(ngrid)   ! surface temperature
! outputs:
      real tsoil(ngrid,nsoil) ! soil (mid-layer) temperature
      real capcal(ngrid) ! surface specific heat
      real fluxgrd(ngrid) ! surface diffusive heat flux

! local variables:
      integer ig,ik

! 0. Initialisations and preprocessing step
      if (firstcall.or.tifeedback) then
      ! note: firstcall is set to .true. or .false. by the caller
      !       and not changed by soil.F 
! 0.1 Build mthermdiff(:), the mid-layer thermal diffusivities
      do ig=1,ngrid
        if (watercaptag(ig)) then
          do ik=0,nsoil-1
! If we have permanent ice, we use the water ice thermal inertia from ground to last layer.
              mthermdiff(ig,ik)=inert_h2o_ice*inert_h2o_ice/volcapa
          enddo
        else
          do ik=0,nsoil-1
           mthermdiff(ig,ik)=therm_i(ig,ik+1)*therm_i(ig,ik+1)/volcapa
          enddo
        endif
      enddo


! 0.2 Build thermdiff(:), the "interlayer" thermal diffusivities
      do ig=1,ngrid
        do ik=1,nsoil-1
      thermdiff(ig,ik)=((layer(ik)-mlayer(ik-1))*mthermdiff(ig,ik)
     &                +(mlayer(ik)-layer(ik))*mthermdiff(ig,ik-1))
     &                    /(mlayer(ik)-mlayer(ik-1))
!	write(*,*),'soil: ik: ',ik,' thermdiff:',thermdiff(ig,ik)
	enddo
      enddo

! 0.3 Build coefficients mu, q_{k+1/2}, d_k, alpha_k and capcal
      ! mu
      mu=mlayer(0)/(mlayer(1)-mlayer(0))

      ! q_{1/2}
      coefq(0)=volcapa*layer(1)/timestep
	! q_{k+1/2}
        do ik=1,nsoil-1
          coefq(ik)=volcapa*(layer(ik+1)-layer(ik))
     &                 /timestep
	enddo

      do ig=1,ngrid
	! d_k
	do ik=1,nsoil-1
	  coefd(ig,ik)=thermdiff(ig,ik)/(mlayer(ik)-mlayer(ik-1))
	enddo
	
	! alph_{N-1}
	alph(ig,nsoil-1)=coefd(ig,nsoil-1)/
     &                  (coefq(nsoil-1)+coefd(ig,nsoil-1))
        ! alph_k
        do ik=nsoil-2,1,-1
	  alph(ig,ik)=coefd(ig,ik)/(coefq(ik)+coefd(ig,ik+1)*
     &                              (1.-alph(ig,ik+1))+coefd(ig,ik))
	enddo

        ! capcal
! Cstar
        capcal(ig)=volcapa*layer(1)+
     &              (thermdiff(ig,1)/(mlayer(1)-mlayer(0)))*
     &              (timestep*(1.-alph(ig,1)))
! Cs
        capcal(ig)=capcal(ig)/(1.+mu*(1.0-alph(ig,1))*
     &                         thermdiff(ig,1)/mthermdiff(ig,0))
!      write(*,*)'soil: ig=',ig,' capcal(ig)=',capcal(ig)
      enddo ! of do ig=1,ngrid
            
      endif ! of if (firstcall.or.tifeedback)

!  1. Compute soil temperatures
      IF (.not.firstcall) THEN
! First layer:
      do ig=1,ngrid
        tsoil(ig,1)=(tsurf(ig)+mu*beta(ig,1)*
     &                         thermdiff(ig,1)/mthermdiff(ig,0))/
     &              (1.+mu*(1.0-alph(ig,1))*
     &               thermdiff(ig,1)/mthermdiff(ig,0))
      enddo
! Other layers:
      do ik=1,nsoil-1
        do ig=1,ngrid
	  tsoil(ig,ik+1)=alph(ig,ik)*tsoil(ig,ik)+beta(ig,ik)
	enddo
      enddo
      
      ENDIF! of if (.not.firstcall)

!  2. Compute beta coefficients (preprocessing for next time step)
! Bottom layer, beta_{N-1}
      do ig=1,ngrid
        beta(ig,nsoil-1)=coefq(nsoil-1)*tsoil(ig,nsoil)
     &                   /(coefq(nsoil-1)+coefd(ig,nsoil-1))
      enddo
! Other layers
      do ik=nsoil-2,1,-1
        do ig=1,ngrid
	  beta(ig,ik)=(coefq(ik)*tsoil(ig,ik+1)+
     &                 coefd(ig,ik+1)*beta(ig,ik+1))/
     &                 (coefq(ik)+coefd(ig,ik+1)*(1.0-alph(ig,ik+1))
     &                  +coefd(ig,ik))
	enddo
      enddo

!  3. Compute surface diffusive flux & calorific capacity
      do ig=1,ngrid
! Cstar
!        capcal(ig)=volcapa(ig,1)*layer(ig,1)+
!     &              (thermdiff(ig,1)/(mlayer(ig,1)-mlayer(ig,0)))*
!     &              (timestep*(1.-alph(ig,1)))
! Fstar
        fluxgrd(ig)=(thermdiff(ig,1)/(mlayer(1)-mlayer(0)))*
     &              (beta(ig,1)+(alph(ig,1)-1.0)*tsoil(ig,1))

!        mu=mlayer(ig,0)/(mlayer(ig,1)-mlayer(ig,0))
!        capcal(ig)=capcal(ig)/(1.+mu*(1.0-alph(ig,1))*
!     &                         thermdiff(ig,1)/mthermdiff(ig,0))
! Fs
        fluxgrd(ig)=fluxgrd(ig)+(capcal(ig)/timestep)*
     &              (tsoil(ig,1)*(1.+mu*(1.0-alph(ig,1))*
     &                         thermdiff(ig,1)/mthermdiff(ig,0))
     &               -tsurf(ig)-mu*beta(ig,1)*
     &                          thermdiff(ig,1)/mthermdiff(ig,0))
      enddo

      end

