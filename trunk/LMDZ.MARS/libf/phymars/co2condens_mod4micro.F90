!======================================================================================================================!
! Module: CO2 condensation for the CO2 cloud microphysics =============================================================!
!----------------------------------------------------------------------------------------------------------------------!
! Authors: Christophe Mathé, Anni Määttänen
! Date: 16/04/2020
!----------------------------------------------------------------------------------------------------------------------!
! Contains subroutines:
!     - co2condens4micro: condensation/sublimation of CO2 ice on the ground and compute pressure change resulting
!
!     - vl1d: Van-Leer scheme
!======================================================================================================================!
module co2condens_mod4micro

implicit none

contains
!======================================================================================================================!
! SUBROUTINE: co2condens4micro ========================================================================================!
!======================================================================================================================!
! Subject:
!---------
!   Condensation/sublimation of CO2 ice on the ground and compute pressure change resulting
!----------------------------------------------------------------------------------------------------------------------!
! Comments:
!----------
!   Adapted from co2condens_mod.F
!----------------------------------------------------------------------------------------------------------------------!
! Paper:
!-------
!  Forget et al. (2008), "Non condensable gas enrichment and depletion in the Martian polar regions."
!----------------------------------------------------------------------------------------------------------------------!
! Algorithm:
!-----------
!   1. Initialization
!   2. Firstcall
!   3. Compute CO2 Volume mixing ratio
!   4. Set zcondicea, zfallice from co2clouds condensation rate and set zt
!   5. Main co2condens
!     5.1. Forecast of ground temperature ztsrf and frost temperature ztcondsol
!     5.2. Check if we have condensation/sublimation on the ground
!     5.3. Compute zfallheat
!     5.4. Compute direct condensation/sublimation of CO2 ice
!       5.4.a. If there is not enough CO2 tracer in 1st layer to condense
!       5.4.b. If the entire CO2 ice layer sublimes (including what has just condensed in the atmosphere)
!     5.5. Changing CO2 ice amount and pressure
!     5.6. Surface albedo and emissivity of the surface below the snow (emisref)
!       5.6.a. Check that amont of CO2 ice is not problematic
!       5.6.b. Set albedo and emissivity of the surface
!       5.6.c. Set pemisurf to emissiv when there is bare surface (needed for co2snow)
!     5.7. Correction to account for redistribution between sigma or hybrid layers when changing surface pressure (and
!         warming/cooling of the CO2 currently changing phase).
!       5.7.a. Mass fluxes through the sigma levels (kg.m-2.s-1)  (>0 when up)
!       5.7.b. Mass of each layer at the end of timestep
!       5.7.c. Corresponding fluxes for T, U, V and Q (averaging operator for transport)
!         5.7.c.i. Value transfert at the surface interface when condensation/sublimation
!         5.7.c.ii. Van Leer scheme
!         5.7.c.iii Compute tendencies on T, U, V, Q
!   6. CO2 snow / clouds scheme
!   7. Extra special case for surface temperature tendency pdtsrfc
!----------------------------------------------------------------------------------------------------------------------!
  subroutine co2condens4micro(ngrid, nlayer, nq, ptimestep, pcapcal, pplay, pplev, ptsrf, pt, pphi, pdt, pdu, pdv, &
                              pdtsrf, pu, pv, pq, pdq, piceco2, psolaralb, pemisurf, pdtc, pdtsrfc, pdpsrf, pduc, &
                              pdvc, pdqc, fluxsurf_sw, zls, zdqssed_co2, pcondicea_co2microp)

  use tracer_mod, only: noms, igcm_co2, igcm_co2_ice, igcm_h2o_vap, igcm_h2o_ice, igcm_dust_number, igcm_dust_mass, &
                        igcm_ccnco2_number, igcm_ccnco2_mass, igcm_ccn_mass

  use surfdat_h, only: emissiv, phisfi

  use geometry_mod, only: latitude, longitude_deg, latitude_deg

  use planete_h, only: obliquit

  use comcstfi_h, only: cpp, g, r, pi

#ifndef MESOSCALE
  use vertical_layers_mod, only: ap, bp
#endif

  implicit none

  include "callkeys.h"
!----------------------------------------------------------------------------------------------------------------------!
! VARIABLE DECLARATION
!----------------------------------------------------------------------------------------------------------------------!
! Input arguments:
!-----------------
  integer, intent(in) :: &
     nq,    &! number of tracers
     ngrid, &! number of atmospheric columns
     nlayer  ! number of vertical layers

  real, intent(in) :: &
     ptimestep,                       &! physics timestep (s)
     pcapcal(ngrid),                  &! surface specific heat
     pplay(ngrid,nlayer),             &! mid-layer pressure (Pa)
     pplev(ngrid,nlayer+1),           &! inter-layer pressure (Pa)
     ptsrf(ngrid),                    &! surface temperature (K)
     pphi(ngrid,nlayer),              &! geopotential (m2.s-2)
     pt(ngrid,nlayer),                &! atmospheric temperature (K)
     pu(ngrid,nlayer),                &! zonal wind (m/s)
     pv(ngrid,nlayer),                &! meridional wind (m/s)
     pdt(ngrid,nlayer),               &! tendency on temperature from previous physical processes (K/s)
     pdu(ngrid,nlayer),               &! tendency on zonal wind from previous physical processes(m/s2)
     pdv(ngrid,nlayer),               &! tendency on meridional wind from previous physical processes (m/s2)
     pdtsrf(ngrid),                   &! tendency on surface temperature from previous physical processes (K/s)
     pq(ngrid,nlayer,nq),             &! tracers (../kg_air)
     pdq(ngrid,nlayer,nq),            &! tendency on tracers from previous physical processes
     zdqssed_co2(ngrid),              &! CO2 flux at the surface  (kg.m-2.s-1)
     fluxsurf_sw(ngrid,2),            &! added to calculate flux dependent albedo
     zls,                             &! solar longitude (rad)
     pcondicea_co2microp(ngrid,nlayer) ! tendency due to CO2 condensation (kg/kg.s-1)
!----------------------------------------------------------------------------------------------------------------------!
! In/output arguments:
!---------------------
  real, intent(inout) :: &
     piceco2(ngrid),   &! CO2 ice on the surface (kg.m-2)
     pemisurf(ngrid),  &! emissivity of the surface
     psolaralb(ngrid,2) ! albedo of the surface
!----------------------------------------------------------------------------------------------------------------------!
! Output arguments:
!------------------
  real, intent(out) :: &
     pdtc(ngrid,nlayer),  &! tendency on temperature dT/dt due to cond/sub (K/s)
     pdtsrfc(ngrid),      &! tendency on surface temperature (K/s)
     pdpsrf(ngrid),       &! tendency on surface pressure (Pa/s)
     pduc(ngrid,nlayer),  &! tendency on zonal wind (m.s-2)
     pdvc(ngrid,nlayer),  &! tendency on meridional wind (m.s-2)
     pdqc(ngrid,nlayer,nq) ! tendency on tracers
!----------------------------------------------------------------------------------------------------------------------!
! Local:
!-------
!----1) Parameters:
!------------------
  real, parameter :: &
     latcond = 595594, &! latent heat of solid CO2 ice (J/kg)
     cpice = 1000.,    &! specific heat of CO2 ice (J.kg-1.K-1)
     tcond1mb = 136.27,&! condensation temperature at 1 mbar (K)
     m_co2 = 44.01E-3, &! CO2 molecular mass (kg/mol)
     m_noco2 = 33.37E-3 ! non condensible molecular mass (kg/mol)

  logical, parameter :: &
     improved_ztcond = .true. ! improved_ztcond flag: If set to .true. (AND running with a 'co2' tracer) then
                              !  condensation temperature is computed using partial pressure of CO2
!----------------------------------------------------------------------------------------------------------------------!
!----2) Saved:
!-------------
  real, save :: &
     A,    &! coefficient used to compute mean molecular mass
     B,    &! coefficient used to compute mean molecular mass
     acond,&! coefficient used to compute ztcondsol
     bcond,&! coefficient used to compute ztcondsol
     ccond  ! coefficient used to compute ztcondsol

  logical, save :: &
     firstcall = .true. ! Used to compute saved variables
!----------------------------------------------------------------------------------------------------------------------!
!----3) Variables:
!-----------------
  integer :: &
     l,  &! loop on layers
     ig, &! loop on ngrid points
     iq   ! loop on tracer

  real :: &
     qco2,                   &! effective quantity of CO2, used to compute mean molecular mass
     mmean,                  &! mean molecular mass
     zfallheat,              &! aerodynamical friction and energy used to heat the amount of ice
     zmflux(nlayer+1),       &! mass fluxes through the sigma levels (kg.m-2.s-1)
     zu(nlayer),             &! effective zonal wind
     zv(nlayer),             &! effective meridional wind
     zq(nlayer,nq),          &! effective tracers quantities
     zq1(nlayer),            &! buffer of zq
     ztsrf(ngrid),           &! effective temperature at the surface
     ztm(nlayer+1),          &! temperature fluxes through the sigma levels
     zum(nlayer+1),          &! zonal wind fluxes through the sigma levels
     zvm(nlayer+1),          &! meridional wind fluxes through the sigma levels
     zqm(nlayer+1,nq),       &! quantity of tracers flux through the sigma levels
     zqm1(nlayer),           &! quantity of tracers after Van-Leer scheme
     masse(nlayer),          &! mass layer (kg)
     w(nlayer+1),            &! total mass fluxes through the sigma levels during ptimestep (kg.m-2)
     availco2,               &! available quantity of co2 (kg)
     emisref(ngrid),         &! emissivity reference
     vmr_co2(ngrid,nlayer),  &! CO2 volume mixing ratio
     zt(nlayer),             &! effective temperature in the atmosphere (K)
     ztcond(ngrid,nlayer+1), &! CO2 condensation temperature (atm)
     ztcondsol(ngrid),       &! CO2 condensation temperature (surface)
     zdiceco2(ngrid),        &! tendency on co2ice surface tracer (kg/m2/s)
     zcondicea(ngrid,nlayer),&! condensation rate in layer l (kg/m2/s)
     zcondices(ngrid),       &! condensation rate on the ground (kg/m2/s)
     zfallice(ngrid)          ! amount of ice falling from layer l (kg/m2/s)

  logical :: &
     condsub(ngrid) ! True if there is condensation/sublimation (used for co2snow)


  ! check with co2 cloud parameterisation
  real :: &
     zt_2(ngrid,nlayer), &
     ztcond_2(ngrid,nlayer+1), &
     zfallice_2(ngrid,nlayer+1), &
     pdtc_2(ngrid,nlayer), &
     zfallheat_2, &
     zcondicea_2(ngrid,nlayer), &
     diff_zcondicea(ngrid,nlayer), &
     diff_zfallice(ngrid)
!======================================================================================================================!
! BEGIN ===============================================================================================================!
!======================================================================================================================!
! 1. Initialization
!----------------------------------------------------------------------------------------------------------------------!
  availco2 = 0.
  zfallheat = 0.
  zt(1:nlayer) = 0.
  ztcond(1:ngrid, 1:nlayer+1) = 0.
  ztcondsol(1:ngrid) = 0.
  zmflux(1:nlayer+1) = 0.
  zu(1:nlayer) = 0.
  zv(1:nlayer) = 0.
  zq(1:nlayer, 1:nq) = 0.
  zq1(1:nlayer) = 0.
  ztsrf(1:ngrid) = 0.
  ztm(1:nlayer+1) = 0.
  zum(1:nlayer+1) = 0.
  zvm(1:nlayer+1) = 0.
  zqm(1:nlayer+1, 1:nq) = 0.
  masse(1:nlayer) = 0.
  w(1:nlayer+1) = 0.
  emisref(1:ngrid) = 0.
  vmr_co2(1:ngrid, 1:nlayer) = 0.
  zcondices(1:ngrid) = 0.
  pdtsrfc(1:ngrid) = 0.
  pdpsrf(1:ngrid) = 0.
  zdiceco2(1:ngrid) = 0.
  condsub(1:ngrid) = .false.
  zcondicea(1:ngrid, 1:nlayer) = 0.
  zfallice(1:ngrid) = 0.
  pduc(1:ngrid, 1:nlayer) = 0.
  pdvc(1:ngrid, 1:nlayer) = 0.
  pdqc(1:ngrid, 1:nlayer, 1:nq) = 0.
  pdtc(1:ngrid,1:nlayer) = 0.
!----------------------------------------------------------------------------------------------------------------------!
! 2. Firstcall
!----------------------------------------------------------------------------------------------------------------------!
! AS: firstcall OK absolute
  if (firstcall) then
    firstcall = .false.

    bcond = 1. / tcond1mb
    ccond = cpp / (g*latcond)
    acond = r / latcond

    write(*,*)'CO2condens: improved_ztcond=', improved_ztcond
    write(*,*)'In co2condens:Tcond(P=1mb)=', tcond1mb, ' Lcond=', latcond
    write(*,*)'acond,bcond,ccond', acond, bcond, ccond

!   Prepare Special treatment if one of the tracer is CO2 gas. Compute A and B coefficient use to compute mean molecular
!   mass Mair defined by:
!     1/Mair = q(igcm_co2)/m_co2 + (1-q(igcm_co2))/m_noco2
!     1/Mair = A*q(igcm_co2) + B
    A = (1./m_co2 - 1./m_noco2)
    B = 1./m_noco2
  end if ! of IF (firstcall)
!----------------------------------------------------------------------------------------------------------------------!
! 3. Compute CO2 Volume mixing ratio
!----------------------------------------------------------------------------------------------------------------------!
  if (improved_ztcond.and.(igcm_co2/=0)) then
    do l = 1, nlayer
      do ig = 1, ngrid
        qco2 = pq(ig,l,igcm_co2) + pdq(ig,l,igcm_co2)*ptimestep
!       Mean air molecular mass = 1/(q(igcm_co2)/m_co2 + (1-q(igcm_co2))/m_noco2)
        mmean = 1. / (A*qco2 +B)
        vmr_co2(ig,l) = (qco2*mmean) / m_co2
      end do
    end do
  else
    do l = 1, nlayer
      do ig = 1, ngrid
        vmr_co2(ig,l) = 0.95
      end do
    end do
  end if
!----------------------------------------------------------------------------------------------------------------------!
! 4. Set zcondicea, zfallice from co2clouds condensation rate
!----------------------------------------------------------------------------------------------------------------------!
  do l = nlayer, 1, -1
    do ig = 1, ngrid
      zcondicea(ig,l) = pcondicea_co2microp(ig,l) *  (pplev(ig,l) - pplev(ig,l+1))/g
    end do
  end do

! Only sedimentation falls on the ground !
  do ig = 1, ngrid
      zfallice(ig) = zdqssed_co2(ig)
      piceco2(ig) = piceco2(ig) + zfallice(ig)*ptimestep
  end do

! Compute without microphysics
 diff_zcondicea(1:ngrid, 1:nlayer) = 0.
 diff_zfallice(1:ngrid) = 0.
 do l =1, nlayer
   do ig = 1, ngrid
     zt_2(ig,l) = pt(ig,l) + pdt(ig,l)*ptimestep
   end do
 end do

 do l = 1, nlayer
   do ig = 1, ngrid
     if (pplay(ig,l).ge.1e-4) then
       ztcond_2(ig,l) = 1. / (bcond-acond*log(.01*vmr_co2(ig,l)*pplay(ig,l)))
     else
       ztcond_2(ig,l) = 0.0 !mars Monica
     endif
   end do
 end do

 ztcond_2(:,nlayer+1)=ztcond_2(:,nlayer)
 zfallice_2(:,:) = 0.
 zcondicea_2(:,:) = 0.
 do l = nlayer , 1, -1
   do ig = 1, ngrid
     pdtc_2(ig,l) = 0.
     if ((zt_2(ig,l).lt.ztcond_2(ig,l)).or.(zfallice_2(ig,l+1).gt.0)) then
       if (zfallice_2(ig,l+1).gt.0) then  
         zfallheat_2 = zfallice_2(ig,l+1)*(pphi(ig,l+1)-pphi(ig,l) + cpice*(ztcond_2(ig,l+1)-ztcond_2(ig,l)))/latcond
       else
         zfallheat_2 = 0.
       end if
       pdtc_2(ig,l) = (ztcond_2(ig,l) - zt_2(ig,l))/ptimestep
       zcondicea_2(ig,l) = (pplev(ig,l)-pplev(ig,l+1))*ccond*pdtc_2(ig,l)- zfallheat_2
       ! Case when the ice from above sublimes entirely
       ! """""""""""""""""""""""""""""""""""""""""""""""
       if (zfallice_2(ig,l+1).lt.- zcondicea_2(ig,l)) then
         zcondicea_2(ig,l)= -zfallice_2(ig,l+1)
       end if
       zfallice_2(ig,l) = zcondicea_2(ig,l)+zfallice_2(ig,l+1)
     end if
   end do
 end do
 diff_zcondicea(:,:) = zcondicea_2(:,:) - zcondicea(:,:)
 diff_zfallice(:) = zfallice_2(:,1) - zfallice(:)
 call writediagfi(ngrid, "diff_zfallice", "sans - avec microphysique", "", 2, diff_zfallice)
 call writediagfi(ngrid, "diff_zcondicea", "sans - avec microphysique", "", 3, diff_zcondicea)
!----------------------------------------------------------------------------------------------------------------------!
! 5. Main co2condens
!----------------------------------------------------------------------------------------------------------------------!
  do ig = 1, ngrid
!----------------------------------------------------------------------------------------------------------------------!
!   5.1. Forecast of ground temperature ztsrf and frost temperature ztcondsol
!----------------------------------------------------------------------------------------------------------------------!
    ztcondsol(ig) = 1. / (bcond-acond*log(.01*vmr_co2(ig,1)*pplev(ig,1)))
    ztsrf(ig) = ptsrf(ig) + pdtsrf(ig)*ptimestep
!----------------------------------------------------------------------------------------------------------------------!
!   5.2. Check if we have condensation/sublimation on the ground
!----------------------------------------------------------------------------------------------------------------------!
!        ground condensation       ||   falling snow       ||    ground sublimation
!----------------------------------------------------------------------------------------------------------------------!
    if ((ztsrf(ig)<ztcondsol(ig)) .OR. (zfallice(ig)/=0.) .OR. ((ztsrf(ig)>ztcondsol(ig)) .AND. (piceco2(ig)/=0.))) then
      condsub(ig) = .true.
!----------------------------------------------------------------------------------------------------------------------!
!   5.3. Compute zfallheat
!----------------------------------------------------------------------------------------------------------------------!
      zfallheat = 0.
!----------------------------------------------------------------------------------------------------------------------!
!   5.4. Compute direct condensation/sublimation of CO2 ice
!----------------------------------------------------------------------------------------------------------------------!
      zcondices(ig) = pcapcal(ig) * (ztcondsol(ig)-ztsrf(ig)) / (latcond*ptimestep) - zfallheat
      pdtsrfc(ig) = (ztcondsol(ig) - ztsrf(ig)) / ptimestep
!----------------------------------------------------------------------------------------------------------------------!
!   5.4.a. If there is not enough CO2 tracer in 1st layer to condense
!----------------------------------------------------------------------------------------------------------------------!
!     Available CO2 tracer in layer 1 at end of timestep (kg/m2)
#ifndef MESOSCALE
      availco2 = pq(ig,1,igcm_co2) * ( (ap(1)-ap(2)) + (bp(1)-bp(2)) * (pplev(ig,1)/g - zcondices(ig)*ptimestep) )
#else
      availco2 = pq(ig,1,igcm_co2)
      PRINT*, "MESOSCALE: CO2 tracer AT FIRST LEVEL IS NOT CORRECTED FROM SIGMA LEVELS"
#endif
      if ( zcondices(ig) * ptimestep>availco2 ) then
        zcondices(ig) = availco2/ptimestep
        zdiceco2(ig) = zcondices(ig) + zfallice(ig)
        pdtsrfc(ig) = (latcond/pcapcal(ig)) *  (zcondices(ig)+zfallheat)
      end if
!----------------------------------------------------------------------------------------------------------------------!
!   5.4.b. If the entire CO2 ice layer sublimes (including what has just condensed in the atmosphere)
!----------------------------------------------------------------------------------------------------------------------!
      if ( (piceco2(ig)/ptimestep) <= -zcondices(ig) ) then
        zcondices(ig) = -piceco2(ig)/ptimestep
        pdtsrfc(ig) = (latcond/pcapcal(ig)) *  (zcondices(ig)+zfallheat)
      end if
!----------------------------------------------------------------------------------------------------------------------!
!   5.5. Changing CO2 ice amount and pressure
!----------------------------------------------------------------------------------------------------------------------!
      zdiceco2(ig) = zcondices(ig) + zfallice(ig) + sum(zcondicea(ig,:))
      piceco2(ig) = piceco2(ig) + zcondices(ig)*ptimestep
      pdpsrf(ig) = -zdiceco2(ig) * g

      if (abs(pdpsrf(ig)*ptimestep)>pplev(ig,1)) then
        print *, 'STOP in condens'
        print *, 'condensing more than total mass'
        print *, 'Grid point ', ig
        print *, 'Longitude(degrees): ', longitude_deg(ig)
        print *, 'Latitude(degrees): ', latitude_deg(ig)
        print *, 'Ps = ', pplev(ig,1)
        print *, 'd Ps = ', pdpsrf(ig)
        call abort_physic('co2condens4micro', 'condensing more than total mass', 1)
      end if
!----------------------------------------------------------------------------------------------------------------------!
!   5.6. Surface albedo and emissivity of the surface below the snow (emisref)
!----------------------------------------------------------------------------------------------------------------------!
!   5.6.a. Check that amont of CO2 ice is not problematic
!----------------------------------------------------------------------------------------------------------------------!
      if(.not.piceco2(ig)>=0.) then
        if(piceco2(ig)<=-5.e-8) then
          write(*,*)'WARNING co2condens piceco2(', ig, ')=', piceco2(ig)
          piceco2(ig) = 0.
        end if
      end if
!----------------------------------------------------------------------------------------------------------------------!
!   5.6.c. Set pemisurf to emissiv when there is bare surface (needed for co2snow)
!----------------------------------------------------------------------------------------------------------------------!
      if (piceco2(ig)==0) then
        pemisurf(ig) = emissiv
      end if
!----------------------------------------------------------------------------------------------------------------------!
!  5.7. Correction to account for redistribution between sigma or hybrid layers when changing surface pressure (and
!         warming/cooling of the CO2 currently changing phase).
!----------------------------------------------------------------------------------------------------------------------!
      do l= 1, nlayer
        zt(l) = pt(ig,l) + pdt(ig,l)*ptimestep
        zu(l) = pu(ig,l) + pdu(ig,l)*ptimestep
        zv(l) = pv(ig,l) + pdv(ig,l)*ptimestep
        do iq=1,nq
          zq(l,iq) = pq(ig,l,iq) + pdq(ig,l,iq)*ptimestep
        end do
      end do
!----------------------------------------------------------------------------------------------------------------------!
!   5.7.a. Mass fluxes through the sigma levels (kg.m-2.s-1)  (>0 when up)
!----------------------------------------------------------------------------------------------------------------------!
      zmflux(1) = - zcondices(ig) - zfallice(ig)
      do l = 1, nlayer
#ifndef MESOSCALE
        zmflux(l+1) = zmflux(l) - zcondicea(ig,l) &
                      + (bp(l)-bp(l+1)) * (-pdpsrf(ig)/g)
! zmflux set to 0 if very low to avoid: top layer is disappearing in v1ld
        if (abs(zmflux(l+1))<1E-13.OR.bp(l+1)==0.) then
          zmflux(l+1) = 0.
        end if
#else
        zmflux(l+1) = zmflux(l) - zcondicea(ig,l)
        if (abs(zmflux(l+1))<1E-13) then
          zmflux(l+1) = 0.
        end if
        PRINT*, "MESOSCALE: FLUX THROUGH SIGMA LEVELS from dPS HAVE TO BE IMPLEMENTED"
#endif
      end do
!----------------------------------------------------------------------------------------------------------------------!
!   5.7.b. Mass of each layer at the end of timestep
!----------------------------------------------------------------------------------------------------------------------!
      do l = 1, nlayer
#ifndef MESOSCALE
        masse(l) = (pplev(ig,l) - pplev(ig,l+1) + (bp(l)-bp(l+1))*pdpsrf(ig)*ptimestep)/g
#else
        masse(l) = (pplev(ig,l) - pplev(ig,l+1))/g
        PRINT*, "MESOSCALE: MASS OF EACH LAYER IS NOT CORRECTED AT END OF TIMESTEP (from SIGMA LEVELS and dPS)"
#endif
      end do
!----------------------------------------------------------------------------------------------------------------------!
!   5.7.c. Corresponding fluxes for T, U, V and Q (averaging operator for transport)
!----------------------------------------------------------------------------------------------------------------------!
!   5.7.c.i. Value transfert at the surface interface when condensation/sublimation
!----------------------------------------------------------------------------------------------------------------------!
      ztm(1) = ztsrf(ig) + pdtsrfc(ig)*ptimestep
      zum(1) = 0
      zvm(1) = 0

!     Most tracer do not condense
      do iq = 1, nq
        zqm(1,iq) = 0.
      end do

!     Special case if one of the tracer is CO2 gas
      if (igcm_co2/=0) then
        zqm(1,igcm_co2) = 1. ! flux is 100% CO2
      end if
!----------------------------------------------------------------------------------------------------------------------!
!   5.7.c.ii. Van Leer scheme
!----------------------------------------------------------------------------------------------------------------------!
      do l=1,nlayer
        w(l)=-zmflux(l)*ptimestep
      end do

      call vl1d(nlayer,zt,2.,masse,w,ztm)
      call vl1d(nlayer,zu ,2.,masse,w,zum)
      call vl1d(nlayer,zv ,2.,masse,w,zvm)

      do iq=1, nq
        do l=1, nlayer
          zq1(l) = zq(l,iq)
        end do

        zqm1(1) = zqm(1,iq)
        zqm1(2:nlayer) = 0.

        call vl1d(nlayer,zq1,2.,masse,w,zqm1)

        do l = 2, nlayer
          zq(l,iq) = zq1(l)
          zqm(l,iq) = zqm1(l)
        end do
      end do

!     Surface condensation affects low winds
      if (zmflux(1)<0) then
        zum(1) = zu(1) * (w(1)/masse(1))
        zvm(1) = zv(1) * (w(1)/masse(1))
        if (w(1)>masse(1)) then ! ensure numerical stability
          zum(1) = ((zu(1)-zum(2))*masse(1)/w(1)) + zum(2)
          zvm(1) = ((zv(1)-zvm(2))*masse(1)/w(1)) + zvm(2)
        end if
      end if

      ztm(nlayer+1) = zt(nlayer) ! should not be used, but...
      zum(nlayer+1) = zu(nlayer)  ! should not be used, but...
      zvm(nlayer+1) = zv(nlayer)  ! should not be used, but...

      do iq = 1, nq
        zqm(nlayer+1,iq) = zq(nlayer,iq)
      end do
!----------------------------------------------------------------------------------------------------------------------!
!   5.7.c.iii Compute tendencies on T, U, V, Q
!----------------------------------------------------------------------------------------------------------------------!
#ifdef MESOSCALE
! AS: This part must be commented in the mesoscale model
! AS: ... to avoid instabilities.
! AS: you have to compile with -DMESOSCALE to do so
#else
      do l = 1, nlayer
!       Tendencies on T
        pdtc(ig,l) = (1./masse(l)) * ( zmflux(l)*(ztm(l) - zt(l))  - zmflux(l+1)*(ztm(l+1) - zt(l)) )

!       Tendencies on U
        pduc(ig,l) = (1./masse(l)) * ( zmflux(l)*(zum(l) - zu(l)) - zmflux(l+1)*(zum(l+1) - zu(l)) )

!       Tendencies on V
        pdvc(ig,l) = (1./masse(l)) * ( zmflux(l)*(zvm(l) - zv(l)) - zmflux(l+1)*(zvm(l+1) - zv(l)) )
      end do
#endif
!       Tendencies on Q
      do iq = 1, nq
        if (iq==igcm_co2) then
          do l = 1, nlayer
            pdqc(ig,l,iq) = (1./masse(l)) * (zmflux(l)*(zqm(l,iq) - zq(l,iq))- zmflux(l+1)*(zqm(l+1,iq) - zq(l,iq))&
                            + zcondicea(ig,l)*(zq(l,iq) - 1.))
          end do
        else
          do l = 1, nlayer
            pdqc(ig,l,iq) = (1./masse(l)) * ( zmflux(l)*(zqm(l,iq)-zq(l,iq)) - zmflux(l+1)*(zqm(l+1,iq)-zq(l,iq))&
                            + zcondicea(ig,l)*zq(l,iq) )
          end do
        end if
      end do
    end if ! if
  end do  ! loop on ig
!----------------------------------------------------------------------------------------------------------------------!
!   5.6.b. Set albedo and emissivity of the surface
!----------------------------------------------------------------------------------------------------------------------!
  call albedocaps(zls,ngrid,piceco2,psolaralb,emisref)
!----------------------------------------------------------------------------------------------------------------------!
! 6. CO2 snow / clouds scheme
!----------------------------------------------------------------------------------------------------------------------!
  call co2snow(ngrid, nlayer, ptimestep, emisref, condsub, pplev, zcondicea, zcondices, zdqssed_co2, &
               pemisurf)
!----------------------------------------------------------------------------------------------------------------------!
! 7. Extra special case for surface temperature tendency pdtsrfc:
!      We want to fix the south pole temperature to CO2 condensation temperature.
!----------------------------------------------------------------------------------------------------------------------!
#ifndef MESOSCALE 
  if (caps.and.(obliquit<27.)) then
    ! check if last grid point is the south pole
    if (abs(latitude(ngrid)-(-pi/2.))<1.e-5) then
!     NB: Updated surface pressure, at grid point 'ngrid', is ps(ngrid)=pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep
!     write(*,*)"co2condens: South pole: latitude(ngrid)=", latitude(ngrid)
      ztcondsol(ngrid) = 1./(bcond-acond*log(.01*vmr_co2(ngrid,1) * (pplev(ngrid,1)+pdpsrf(ngrid)*ptimestep)))
      pdtsrfc(ngrid) = (ztcondsol(ngrid)-ztsrf(ngrid))/ptimestep
    end if
  end if
#endif
!======================================================================================================================!
! END =================================================================================================================!
!======================================================================================================================!
  end subroutine co2condens4micro


!**********************************************************************************************************************!
!**********************************************************************************************************************!


!======================================================================================================================!
! SUBROUTINE: Van-Leer scheme =========================================================================================!
!======================================================================================================================!
! Subject:
!---------
!  Operateur de moyenne inter-couche pour calcul de transport type Van-Leer " pseudo amont " dans la verticale
!----------------------------------------------------------------------------------------------------------------------!
! Comments:
!----------
!  q,w are input arguments for the s-pg ....
!----------------------------------------------------------------------------------------------------------------------!
! Paper:
!-------
!  Van-Leer (1977), "Towards the Ultimate Conservative Difference Scheme. IV. A New Approach to Numerical Convection"
!----------------------------------------------------------------------------------------------------------------------!
  subroutine vl1d(nlayer,q,pente_max,masse,w,qm)

  implicit none
!----------------------------------------------------------------------------------------------------------------------!
! VARIABLE DECLARATION
!----------------------------------------------------------------------------------------------------------------------!
! Input arguments:
!-----------------
  integer, intent(in) :: &
     nlayer ! number of layers

  real, intent(in) :: &
     pente_max,     &! coefficient, pente_max = 2 advised
     masse(nlayer), &! masse layer Dp/g (kg)
     q(nlayer)       ! quantity of tracer
!----------------------------------------------------------------------------------------------------------------------!
! In-Output arguments:
!---------------------
  real, intent(inout) :: &
     w(nlayer+1) ! masse d'atm ``transferee'' a chaque pas de temps (kg.m-2)
!----------------------------------------------------------------------------------------------------------------------!
! Output arguments:
!------------------
  real, intent(out) :: &
     qm(nlayer+1) ! quantity of tracer after Van-Leer scheme
!----------------------------------------------------------------------------------------------------------------------!
! Locals variables:
!------------------
  integer :: &
     l, &! loop on nlayer
     m   ! index

  real :: &
     dzqmax,      &! maximum of dzq between two adjacent layers
     sigw,        &!
     Mtot,        &!
     MQtot,       &!
     dzq(nlayer), &!
     dzqw(nlayer),&!
     adzqw(nlayer) !
!======================================================================================================================!
! BEGIN ===============================================================================================================!
!======================================================================================================================!
! 1. On oriente tout dans le sens de la pression: w > 0 when down
!----------------------------------------------------------------------------------------------------------------------!
  do l = 2, nlayer
    dzqw(l) = q(l-1) - q(l)
    adzqw(l) = abs(dzqw(l))
  end do

  do l = 2, nlayer-1
    if(dzqw(l)*dzqw(l+1)>0.) then
      dzq(l) = 0.5 * (dzqw(l)+dzqw(l+1))
    else
      dzq(l) = 0.
    end if

    dzqmax = pente_max * min(adzqw(l), adzqw(l+1))

    dzq(l) = sign(min(abs(dzq(l)),dzqmax), dzq(l))
  end do

  dzq(1)=0.
  dzq(nlayer)=0.

  do l = 1, nlayer-1
!----------------------------------------------------------------------------------------------------------------------!
! 2.1. Regular scheme (transfered mass < layer mass)
!----------------------------------------------------------------------------------------------------------------------!
    if (w(l+1)>0. .and. w(l+1)<=masse(l+1)) then
      sigw = w(l+1) / masse(l+1)
      qm(l+1) = (q(l+1) + 0.5*(1.-sigw)*dzq(l+1))
    else if (w(l+1)<=0. .and. -w(l+1)<=masse(l)) then
      sigw = w(l+1) / masse(l)
      qm(l+1) = (q(l) - 0.5*(1.+sigw)*dzq(l))
!----------------------------------------------------------------------------------------------------------------------!
! 2.2. Extended scheme (transfered mass > layer mass)
!----------------------------------------------------------------------------------------------------------------------!
    else if (w(l+1)>0.) then
      m = l+1
      Mtot = masse(m)
      MQtot = masse(m)*q(m)

      do while ((m<nlayer).and.(w(l+1)>(Mtot+masse(m+1))))
        m = m+1
        Mtot = Mtot + masse(m)
        MQtot = MQtot + masse(m)*q(m)
      end do

      if (m<nlayer) then
        sigw = (w(l+1)-Mtot) / masse(m+1)
        qm(l+1) = (1/w(l+1))*( MQtot + (w(l+1)-Mtot)* (q(m+1)+0.5*(1.-sigw)*dzq(m+1)) )
      else
        w(l+1) = Mtot
        qm(l+1) = Mqtot / Mtot
        call abort_physic('co2condens4micro', 'top layer is disapearing !', 1)
      end if
!----------------------------------------------------------------------------------------------------------------------!
    else ! if(w(l+1).lt.0)
      m = l-1
      Mtot = masse(m+1)
      MQtot = masse(m+1)*q(m+1)
      if (m>0) then ! because some compilers will have problems evaluating masse(0)
        do while ((m>0).and.(-w(l+1)>(Mtot+masse(m))))
          m = m-1
          Mtot = Mtot + masse(m+1)
          MQtot = MQtot + masse(m+1)*q(m+1)
          if (m==0) then
            exit
          end if
        end do
      end if

      if (m>0) then
        sigw = (w(l+1)+Mtot) / masse(m)
        qm(l+1) = (-1/w(l+1)) * ( MQtot + (-w(l+1)-Mtot) * (q(m)-0.5*(1.+sigw)*dzq(m)) )
      else
        qm(l+1) = (-1/w(l+1)) * (MQtot + (-w(l+1)-Mtot)*qm(1))
      end if
    end if
  end do ! l = 1, nlayer-1
!======================================================================================================================!
! END =================================================================================================================!
!======================================================================================================================!
  end subroutine vl1d

end module co2condens_mod4micro
