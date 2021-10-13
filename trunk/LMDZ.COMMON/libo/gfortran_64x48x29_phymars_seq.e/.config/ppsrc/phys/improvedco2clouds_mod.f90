










!======================================================================================================================!
! Module: Scheme of co2 cloud formation ===============================================================================!
!----------------------------------------------------------------------------------------------------------------------!
! Authors: Joaquim Audouard, Constantino Listowski, Anni Määttänen
! Date: 09/2016
!----------------------------------------------------------------------------------------------------------------------!
! Contains subroutines:
!     - improvedco2clouds_mod: nucleation
!======================================================================================================================!
module improvedco2clouds_mod

implicit none

contains
!======================================================================================================================!
! SUBROUTINE: improvedco2clouds =======================================================================================!
!======================================================================================================================!
! Subject:
!---------
!  This routine is used to form CO2 clouds when a parcel of the GCM is saturated.
!----------------------------------------------------------------------------------------------------------------------!
! Comments:
!----------
!  Adaptation for CO2 clouds based on improvedclouds_mod.F
!
!  It includes the ability to have supersaturation, a computation of the nucleation rates, growthrates and the
!  scavenging of dust particles by clouds. It is worth noting that the amount of dust is computed using the dust
!  optical depth computed in aeropacity.F.
!  That's why the variable called "tauscaling" is used to convert pq(dust_mass) and pq(dust_number), which are relative
!  quantities, to absolute and realistic quantities stored in zq. This has to be done to convert the inputs into
!  absolute values, but also to convert the outputs back into relative values which are then used by the sedimentation
!  and advection schemes.
!  CO2 ice particles can nucleate on both dust and water ice particles. When CO2 ice is deposited onto a water ice
!  particles, the particle is removed from the water tracers. Memory of the origin of the co2 particles is kept and
!  thus the water cycle shouldn't be modified by this.
!  There is an energy limit to how much co2 can sublimate/condensate. It is defined by the difference of the GCM
!  temperature with the co2 condensation temperature.
!
!                                                /!\ WARNING /!!  No sedimentation of the water ice origin is performed in the microphysical timestep in co2cloud.F.
!
!  If meteoritic particles are activated and turn into co2 ice particles, then they will be reversed in the dust
!  tracers if the cloud sublimates.
!----------------------------------------------------------------------------------------------------------------------!
! Paper:
!-------
!    see co2cloud.F90
!----------------------------------------------------------------------------------------------------------------------!
! Algorithm:
!-----------
!   0. Firstcall
!     0.1. Bonus: meteoritic component, extract data
!   1. Initialization
!   2. Compute saturation
!   3. Bonus: additional meteoritic particles for nucleation
!   4. Actual microphysics: Main loop over the GCM's grid
!     4.1 Nucleation
!     4.2. Ice growth: scheme for radius evolution
!     4.3 Dust cores releasing if no more co2 ice
!   5. Get cloud tendencies
!======================================================================================================================!
  subroutine improvedCO2clouds(ngrid, nlay, microtimestep, pplay, pplev, pteff, sum_subpdt, pqeff, sum_subpdq, &
                               subpdqcloudco2, subpdtcloudco2, nq, tauscaling, mem_Mccn_co2, mem_Mh2o_co2, &
                               mem_Nccn_co2, rb_cldco2, sigma_iceco2, dev2)

  use comcstfi_h,   only: pi, g, cpp

  use updaterad,    only: updaterice_micro, updaterice_microco2, updaterccnCO2

  use tracer_mod,   only: igcm_dust_mass, igcm_dust_number, rho_dust, igcm_h2o_ice, igcm_ccn_mass, igcm_ccn_number, &
                          nuice_sed, igcm_co2, igcm_co2_ice, igcm_ccnco2_mass, igcm_ccnco2_number, nuiceco2_sed, &
                          nuiceco2_ref

  use conc_mod,     only: mmean

  use datafile_mod, only: datadir

  use density_co2_ice_mod, only: density_co2_ice

  implicit none

  include "callkeys.h"
  include "microphys.h"
!----------------------------------------------------------------------------------------------------------------------!
! VARIABLES DECLARATION
!----------------------------------------------------------------------------------------------------------------------!
! Input arguments:
!-----------------
  integer, intent(in) :: &
     nq,    &! number of tracers
     ngrid, &! number of point grid
     nlay    ! number of layer

  real, intent(in) :: &
     microtimestep,           &! physics time step (s)
     pplay(ngrid,nlay),       &! mid-layer pressure (Pa)
     pplev(ngrid,nlay+1),     &! inter-layer pressure (Pa)
     pteff(ngrid,nlay),       &! temperature at the middle of the layers (K)
     sum_subpdt(ngrid,nlay),  &! tendency on temperature from previous physical parametrizations
     pqeff(ngrid,nlay,nq),    &! tracers (kg/kg)
     tauscaling(ngrid),       &! convertion factor for qdust and Ndust
     sum_subpdq(ngrid,nlay,nq) ! tendencies on tracers before condensation (kg/kg.s-1)

  real,  intent(in) :: &
     sigma_iceco2 ! Variance of the co2 ice and CCN distributions

  double precision, intent(in) :: &
     rb_cldco2(nbinco2_cld+1), & ! boundary values of each rad_cldco2
     dev2
!----------------------------------------------------------------------------------------------------------------------!
! Output arguments:
!------------------
  real, intent(out) :: &
     subpdtcloudco2(ngrid,nlay),  &! tendency on tracers due to CO2 condensation (K/s)
     subpdqcloudco2(ngrid,nlay,nq) ! tendency on tracers due to CO2 condensation (kg/kg.s-1)
!----------------------------------------------------------------------------------------------------------------------!
! Local:
!-------
!----1) Parameters:
!------------------
  integer, parameter :: &
! ---Meteoritic flux input file
     nbin_meteor = 100, &! Dimension 1
     nlev_meteor = 130, &! Dimension 2
     uMeteor = 666,     &! File unit
! ---Latent heat computation
     l0 = 595594d0,     &! coeff from: Azreg-Aïnou (2005)
     l1 = 903.111d0,    &!   Title: "Low-temperature data for carbon dioxide"
     l2 = -11.5959d0,   &!   Pulication: eprint arXiv:1403.4403
     l3 = 0.0528288d0,  &!
     l4 = -0.000103183d0 !

  real, parameter :: &
     threshold = 1e-30 ! limit value

  character(len=20), parameter:: &
     file_meteoritic_flux = 'Meteo_flux_Plane.dat'
!----------------------------------------------------------------------------------------------------------------------!
!----2) Saved:
!-------------
  real, save :: &
     sigma_ice  ! Variance of the h2o ice and CCN distributions

  double precision, save :: &
     meteor(nlev_meteor,nbin_meteor), &! Meteoritic flux read from file uMeteor
     dev3                              ! 1. / ( sqrt(2.) * sigma_ice )

  logical, save :: &
     firstcall = .true. ! Used to compute saved variables
!----------------------------------------------------------------------------------------------------------------------!
!----3) Variables:
!-----------------
  integer :: &
     ig,     &! loop on ngrid
     l,      &! loop on nlay
     i,      &! loop on nbinco2
! ---Variables for meteoritic flux
     ibin,   &! loop on nbin_meteor
     nelem,  &! nb of points during interpolation of meteoritic flux
     lebon1, &! index where P_meteor is the nearest of pplev(ig,l)
     lebon2, &! index where P_meteor is the nearest of pplev(ig,l+1)
     read_ok  ! file uMeteor iostat

  real :: &
    masse(ngrid,nlay),     &! mass layer (kg.m-2)
    rice(ngrid,nlay),      &! water ice mass mean radius (m): used for nucleation of CO2 on ice-coated ccns
    zq(ngrid,nlay,nq),     &! local value of tracers (kg/kg)
    zq0(ngrid,nlay,nq),    &! local init value of tracers (kg/kg)
    zt(ngrid,nlay),        &! local value of temperature (K)
    zqsat(ngrid,nlay),     &! saturation vapor pressure for CO2 (K)
    tcond(ngrid,nlay),     &! condensation temperature of CO2 (K)
    lw,                    &! Latent heat of sublimation (J.kg-1)
    rdust(ngrid,nlay),     &! Dust geometric mean radius (m)
    rhocloud(ngrid,nlay),  &! Cloud density (kg.m-3)
    rhocloudco2(ngrid,nlay) ! Cloud density (kg.m-3)

  real(kind=8) :: &
    derf ! Error function
 
  double precision :: &
    dMice,                             &! mass of condensed ice
    facteurmax,                        &! for energy limit on mass growth
    pco2,                              &! Co2 vapor partial pressure (Pa)
    satu,                              &! Co2 vapor saturation ratio over ice
    Mo,                                &! mass of aerosol particles
    No,                                &! number of aerosol particles
    Rn,                                &! logarithm of rdust/rice
    Rm,                                &! Rn * variance of ice and CCN distribution
    n_derf,                            &! derf( (rb_cldco2(1)+Rn) *dev3)
    m_derf,                            &! derf( (rb_cldco2(1)+Rm) *dev2)
    mem_Mccn_co2(ngrid,nlay),          &! Memory of CCN mass of H2O and dust used by CO2
    mem_Mh2o_co2(ngrid,nlay),          &! Memory of H2O mass integred into CO2 crystal
    mem_Nccn_co2(ngrid,nlay),          &! Memory of CCN number of H2O and dust used by CO2
    n_aer(nbinco2_cld),                &! Radius used by the microphysical scheme (m)
    m_aer(nbinco2_cld),                &! number concentration V-1 of particle/each size bin
    n_aer_h2oice(nbinco2_cld),         &! mass mixing ratio of particle/each size bin
    m_aer_h2oice(nbinco2_cld),         &! Same - for CO2 nucleation
    rad_h2oice(nbinco2_cld),           &! Same - for CO2 nucleation
    Ic_rice,                           &! Mass transfer rate CO2 ice crystal
    ratioh2o_ccn,                      &! 1./(zq(ig,l,igcm_h2o_ice)  + zq(ig,l,igcm_ccn_mass)*tauscaling(ig))
    vo2co2,                            &! volume of co2 ice particle
    dN,                                &! number of particle of dust used as ccn
    dM,                                &! mass of dN
    dNh2o,                             &! number of particle of h2o ice used as ccn
    dMh2o,                             &! mass of dNh2o
    dNN,                               &! min(dN,zq(ig,l,igcm_dust_number))
    dMM,                               &! min(dM,zq(ig,l,igcm_dust_mass))
    dNNh2o,                            &! min(dNNh2o,zq(ig,l,igcm_ccn_number))
    dMh2o_ice,                         &! min(dMh2o*zq(ig,l,igcm_h2o_ice)*ratioh2o_ccn, zq(ig,l,igcm_h2o_ice))
    dMh2o_ccn,                         &! min(dMh2o_ccn,zq(ig,l,igcm_ccn_mass))
    rate(nbinco2_cld),                 &! nucleation rate
    rateh2o(nbinco2_cld),              &! nucleation rate for h2o
    rho_ice_co2T,                      &! density of co2 ice Temperature-dependent
    riceco2(ngrid,nlay),               &! CO2 ice mean radius (m)
    vrat_cld,                          &! Volume ratio
    Proba,                             &! 1.0 - exp(-1.*microtimestep*rate(i))
    Probah2o,                          &! 1.0 - exp(-1.*microtimestep*rateh2o(i))
    mtemp(nbinco2_cld),                &! sum(meteor(lebon1:lebon2,ibin))
    pression_meteor(nlev_meteor),      &! pressure from meteoritic flux input file
    ltemp1(nlev_meteor),               &! abs(pression_meteor(:)-pplev(ig,l))
    ltemp2(nlev_meteor),               &! abs(pression_meteor(:)-pplev(ig,l+1))
    meteor_ccn(ngrid,nlay,nbinco2_cld)  ! nbinco2_cld = 100

  logical :: &
    file_ok ! test if meteoritic input file exists
!======================================================================================================================!
! BEGIN ===============================================================================================================!
!======================================================================================================================!
! 0. Firstcall
!----------------------------------------------------------------------------------------------------------------------!
  if (firstcall) then
    firstcall = .false.

!   Variance of the ice and CCN distributions
    sigma_ice = sqrt(log(1.+nuice_sed))
    dev3 = 1. / ( sqrt(2.) * sigma_ice )
!----------------------------------------------------------------------------------------------------------------------!
! 0.1. Bonus: meteoritic component, extract data
!----------------------------------------------------------------------------------------------------------------------!
    if (meteo_flux) then
      ! Check if file exists
      inquire(file=trim(datadir)//'/'//file_meteoritic_flux, exist=file_ok)
      if (.not. file_ok) then
        call abort_physic("CO2clouds", 'file '//file_meteoritic_flux//' should be in'//trim(datadir), 1)
      end if

      ! open file
      open(unit=uMeteor,file=trim(datadir)//'/'//file_meteoritic_flux, FORM='formatted')

      !skip 1 line
      read(uMeteor,*)

      ! extract pressure_meteor
      do i = 1, nlev_meteor
        read(uMeteor,*)pression_meteor(i)
      end do

      !skip 1 line
      read(uMeteor,*)

      ! extract meteor flux
      do i = 1, nlev_meteor
        ! les mêmes 100 bins size que la distri nuclea : on touche pas
        do ibin = 1, nbin_meteor
          read(uMeteor,'(F12.6)') meteor(i,ibin)
        end do
      end do

      ! close file
      close(uMeteor)
    end if ! of if meteo_flux
  end if ! firstcall
!----------------------------------------------------------------------------------------------------------------------!
! 1. Initialization
!----------------------------------------------------------------------------------------------------------------------!
  rdust(:,:) = 0.
  meteor_ccn(:,:,:) = 0.
  rice(:,:) = 1.e-8
  riceco2(:,:) = 1.e-11

  ! Initialize the tendencies
  subpdqcloudco2(:,:,:) = 0.
  subpdtcloudco2(:,:) = 0.

  ! pteff temperature layer; sum_subpdt dT.s-1
  zt(1:ngrid,1:nlay) = 0.
  zt(:,:) = pteff(:,:) + sum_subpdt(:,:) * microtimestep

  ! pqeff traceur kg/kg; sum_subpdq tendance idem .s-1
  zq(:,:,:) = pqeff(:,:,:) + sum_subpdq(:,:,:) * microtimestep
  where( zq(:,:,:) < threshold ) zq(:,:,:) = threshold

  zq0(:,:,:) = zq(:,:,:)
  zqsat(:,:) = 0.
!----------------------------------------------------------------------------------------------------------------------!
! 2. Compute saturation
!----------------------------------------------------------------------------------------------------------------------!
  call co2sat(ngrid*nlay,zt,zqsat)
  call tcondco2(ngrid,nlay,pplay, zq(:,:,igcm_co2) + zq(:,:,igcm_co2_ice),tcond)
!----------------------------------------------------------------------------------------------------------------------!
! 3. Bonus: additional meteoritic particles for nucleation
!----------------------------------------------------------------------------------------------------------------------!
  ! TODO: instead of intepolation, used only the nearest pplev(ig,l)
  if (meteo_flux) then
    do l = 1, nlay
      do ig = 1, ngrid
        masse(ig,l) = (pplev(ig,l) - pplev(ig,l+1)) / g

        ltemp1 = abs(pression_meteor(:)-pplev(ig,l))
        ltemp2 = abs(pression_meteor(:)-pplev(ig,l+1))

        lebon1 = minloc(ltemp1,DIM=1)
        lebon2 = minloc(ltemp2,DIM=1)

        nelem = lebon2-lebon1+1.

        mtemp(:) = 0d0

        do ibin = 1, nbin_meteor
          mtemp(ibin) = sum(meteor(lebon1:lebon2,ibin))
        end do

        ! Par kg air csi par m carre, x epaisseur/masse pour par kg/air. Check original unit with J. Plane
        meteor_ccn(ig,l,:)=mtemp(:)/nelem/masse(ig,l)
      end do
    end do
  end if
!----------------------------------------------------------------------------------------------------------------------!
! 4. Actual microphysics: Main loop over the GCM's grid
!----------------------------------------------------------------------------------------------------------------------!
  do l = 1, nlay
    do ig = 1, ngrid

      ! Get the partial pressure of co2 vapor and its saturation ratio
      pco2 = zq(ig,l,igcm_co2) * (mmean(ig,l)/(mco2*1e3)) * pplay(ig,l) ! mco2 (kg/mol) => mmean (g/mol)
      satu = pco2 / zqsat(ig,l)

      !T-dependant CO2 ice density
      call density_co2_ice(dble(zt(ig,l)), rho_ice_co2T)

      vo2co2 = m0co2 / rho_ice_co2T
!----------------------------------------------------------------------------------------------------------------------!
! 4.1 Nucleation
!----------------------------------------------------------------------------------------------------------------------!
      ! if there is condensation
      if ( satu >= 1 ) then
        call updaterccnCO2(zq(ig,l,igcm_dust_mass), zq(ig,l,igcm_dust_number), rdust(ig,l), tauscaling(ig))

        ! Expand the dust moments into a binned distribution
        n_aer(:) = 0d0 ! number of aerosol particles
        m_aer(:) = 0d0 ! mass of aerosol particles

        No = zq(ig,l,igcm_dust_number) * tauscaling(ig)

        Mo = (4./3.) * pi * rho_dust * No * rdust(ig,l)**3 *exp(9.*nuiceco2_ref/2.)

        if (No > threshold) then
          Rn = -log(rdust(ig,l))

          Rm = Rn - 3. * sigma_iceco2 * sigma_iceco2

          n_derf = derf( (rb_cldco2(1)+Rn) *dev2)
          m_derf = derf( (rb_cldco2(1)+Rm) *dev2)

          do i = 1, nbinco2_cld
            n_aer(i) = -0.5 * No * n_derf
            m_aer(i) = -0.5 * Mo * m_derf

            n_derf = derf((rb_cldco2(i+1)+Rn) *dev2)
            m_derf = derf((rb_cldco2(i+1)+Rm) *dev2)

            n_aer(i) = n_aer(i) + 0.5 * No * n_derf
            m_aer(i) = m_aer(i) + 0.5 * Mo * m_derf
          end do

          ! Ajout meteor_ccn particles aux particules de poussière background
          if (meteo_flux) then
            do i = 1, nbinco2_cld
              n_aer(i) = n_aer(i) + meteor_ccn(ig,l,i)

              m_aer(i) = m_aer(i) + (4./3.) * pi * rho_dust *meteor_ccn(ig,l,i) * rad_cldco2(i)**3
             end do
          end if

          ! Same but with h2o particles as CCN only if co2useh2o = .true.
          n_aer_h2oice(:)=0.
          m_aer_h2oice(:)=0.

          if (co2useh2o) then
            call updaterice_micro(zq(ig,l,igcm_h2o_ice), zq(ig,l,igcm_ccn_mass), zq(ig,l,igcm_ccn_number), &
                                  tauscaling(ig), rice(ig,l), rhocloud(ig,l))

            Mo = zq(ig,l,igcm_h2o_ice) + zq(ig,l,igcm_ccn_mass) * tauscaling(ig) + threshold

            ! Total mass of H20 crystals,CCN included
            No = zq(ig,l,igcm_ccn_number) * tauscaling(ig) + threshold

            Rn = -log(rice(ig,l))

            Rm = Rn - 3. * sigma_ice * sigma_ice

            n_derf = derf( (rb_cldco2(1)+Rn) *dev3)
            m_derf = derf( (rb_cldco2(1)+Rm) *dev3)

            do i = 1, nbinco2_cld
              n_aer_h2oice(i) = -0.5 * No * n_derf
              m_aer_h2oice(i) = -0.5 * Mo * m_derf

              n_derf = derf( (rb_cldco2(i+1)+Rn) *dev3)
              m_derf = derf( (rb_cldco2(i+1)+Rm) *dev3)

              n_aer_h2oice(i) = n_aer_h2oice(i) + 0.5 * No * n_derf
              m_aer_h2oice(i) = m_aer_h2oice(i) + 0.5 * Mo * m_derf

              rad_h2oice(i) = rad_cldco2(i)
            end do
          end if

          ! Call to nucleation routine
          call nucleaCO2(dble(pco2), zt(ig,l), dble(satu), n_aer, rate, n_aer_h2oice, rad_h2oice, rateh2o, vo2co2)

          dN = 0.
          dM = 0.
          dNh2o = 0.
          dMh2o = 0.

          do i = 1, nbinco2_cld
            Proba = 1.0 - exp(-1.*microtimestep*rate(i))
            dN = dN + n_aer(i) * Proba
            dM = dM + m_aer(i) * Proba
          end do

          if (co2useh2o) then
            do i = 1, nbinco2_cld
              Probah2o = 1.0 - exp(-1.*microtimestep*rateh2o(i))
              dNh2o = dNh2o + n_aer_h2oice(i) * Probah2o
              dMh2o = dMh2o + m_aer_h2oice(i) * Probah2o
            end do
          end if

          ! Now increment CCN tracers and update dust tracers
          dNN = min(dN,zq(ig,l,igcm_dust_number)) ! dNN est devenu DN
          dMM = min(dM,zq(ig,l,igcm_dust_mass))  ! idem dans le min

          zq(ig,l,igcm_ccnco2_mass) = zq(ig,l,igcm_ccnco2_mass) + dMM /tauscaling(ig)

          zq(ig,l,igcm_ccnco2_number) = zq(ig,l,igcm_ccnco2_number) + dNN /tauscaling(ig)

          zq(ig,l,igcm_dust_mass) = zq(ig,l,igcm_dust_mass) - dMM /tauscaling(ig)

          zq(ig,l,igcm_dust_number) = zq(ig,l,igcm_dust_number) - dNN /tauscaling(ig)

          ! Update CCN for CO2 nucleating on H2O CCN : Warning: must keep memory of it
          if (co2useh2o) then
            dNNh2o = dNh2o/tauscaling(ig)
            dNNh2o = min(dNNh2o,zq(ig,l,igcm_ccn_number))

            ratioh2o_ccn = 1./(zq(ig,l,igcm_h2o_ice)  + zq(ig,l,igcm_ccn_mass)*tauscaling(ig))

            dMh2o_ccn = dMh2o * zq(ig,l,igcm_ccn_mass) *  tauscaling(ig) * ratioh2o_ccn
            dMh2o_ccn = dMh2o_ccn/tauscaling(ig)
            dMh2o_ccn = min(dMh2o_ccn,zq(ig,l,igcm_ccn_mass))

            dMh2o_ice = dMh2o * zq(ig,l,igcm_h2o_ice) * ratioh2o_ccn
            dMh2o_ice = min(dMh2o_ice,zq(ig,l,igcm_h2o_ice))

            zq(ig,l,igcm_ccnco2_mass) = zq(ig,l,igcm_ccnco2_mass)  + dMh2o_ice + dMh2o_ccn

            zq(ig,l,igcm_ccnco2_number) =  zq(ig,l,igcm_ccnco2_number) + dNNh2o

            zq(ig,l,igcm_ccn_number) = zq(ig,l,igcm_ccn_number)  - dNNh2o

            zq(ig,l,igcm_h2o_ice) = zq(ig,l,igcm_h2o_ice)  - dMh2o_ice

            zq(ig,l,igcm_ccn_mass) = zq(ig,l,igcm_ccn_mass) - dMh2o_ccn

            mem_Mh2o_co2(ig,l) = mem_Mh2o_co2(ig,l) + dMh2o_ice
            mem_Mccn_co2(ig,l) = mem_Mccn_co2(ig,l) + dMh2o_ccn
            mem_Nccn_co2(ig,l) = mem_Nccn_co2(ig,l) + dNNh2o
          end if ! of if co2useh2o
        end if ! of if No > 1e-30
      end if ! of is satu > 1
!----------------------------------------------------------------------------------------------------------------------!
! 4.2. Ice growth: scheme for radius evolution
!----------------------------------------------------------------------------------------------------------------------!
!    We trigger crystal growth if and only if there is at least one nuclei (N>1). Indeed, if we are supersaturated
!      and still don't have at least one nuclei, we should better wait to avoid unrealistic value for nuclei radius
!      and so on for cases that remain negligible.
!----------------------------------------------------------------------------------------------------------------------!
      ! we trigger crystal growth
      if (zq(ig,l,igcm_ccnco2_number) * tauscaling(ig) + threshold >= 1) then

        call updaterice_microco2(dble(zq(ig,l,igcm_co2_ice)), dble(zq(ig,l,igcm_ccnco2_mass)), &
                                 dble(zq(ig,l,igcm_ccnco2_number)), zt(ig,l), tauscaling(ig), riceco2(ig,l), &
                                 rhocloudco2(ig,l))
        Ic_rice = 0.

        ! J.kg-1
        lw = l0 + l1 * zt(ig,l) + l2 * zt(ig,l)**2 + l3 * zt(ig,l)**3 + l4 * zt(ig,l)**4

        facteurmax = abs(tcond(ig,l)-zt(ig,l)) * (cpp/lw)

        ! call scheme of microphys. mass growth for CO2 (evaporation/condensation)
        call massflowrateCO2(pplay(ig,l), zt(ig,l), satu, riceco2(ig,l), mmean(ig,l), Ic_rice)

        ! Ic_rice Mass transfer rate (kg/s) for a rice particle > 0 si croissance !
        if (isnan(Ic_rice) .or. Ic_rice == 0.) then
          Ic_rice = 0.
          subpdtcloudco2(ig,l) = -sum_subpdt(ig,l)
          dMice = 0
        else
          ! Kg par kg d'air, >0 si croissance !
          ! kg.s-1 par particule * nb particule par kg air*s = kg par kg air
          dMice = zq(ig,l,igcm_ccnco2_number) * Ic_rice * microtimestep * tauscaling(ig)

          ! facteurmax maximum quantity of CO2 that can sublime/condense according to available thermal energy
          ! latent heat release > 0 if growth i.e. if dMice > 0
          dMice = max(dMice,max(-facteurmax,-zq(ig,l,igcm_co2_ice)))
          dMice = min(dMice,min(facteurmax,zq(ig,l,igcm_co2)))

          ! kgco2/kgair* J/kgco2 * 1/(J.kgair-1.K-1)/s = K /s
          subpdtcloudco2(ig,l) = dMice * lw / cpp / microtimestep

          !Now update tracers
          zq(ig,l,igcm_co2_ice) = zq(ig,l,igcm_co2_ice) + dMice
          zq(ig,l,igcm_co2) = zq(ig,l,igcm_co2) - dMice
        end if
      end if ! if zq(ccnco2_number) >= 1
!----------------------------------------------------------------------------------------------------------------------!
! 4.3 Dust cores releasing if no more co2 ice
!----------------------------------------------------------------------------------------------------------------------!
      ! On sublime tout
      if ((zq(ig,l,igcm_co2_ice) <= threshold).or.(zq(ig,l,igcm_ccnco2_number) * tauscaling(ig) < 1.)) then

        if (co2useh2o) then

          if (mem_Mccn_co2(ig,l) > 0) then
            zq(ig,l,igcm_ccn_mass) = zq(ig,l,igcm_ccn_mass) + mem_Mccn_co2(ig,l)
          end if

          if (mem_Mh2o_co2(ig,l) > 0) then
            zq(ig,l,igcm_h2o_ice) = zq(ig,l,igcm_h2o_ice) + mem_Mh2o_co2(ig,l)
          end if

          if (mem_Nccn_co2(ig,l) > 0) then
            zq(ig,l,igcm_ccn_number) = zq(ig,l,igcm_ccn_number) + mem_Nccn_co2(ig,l)
          end if

        end if

        zq(ig,l,igcm_dust_mass) = zq(ig,l,igcm_dust_mass) + zq(ig,l,igcm_ccnco2_mass) - ( mem_Mh2o_co2(ig,l) + &
                                     mem_Mccn_co2(ig,l) )

        zq(ig,l,igcm_dust_number) = zq(ig,l,igcm_dust_number) + zq(ig,l,igcm_ccnco2_number)  - mem_Nccn_co2(ig,l)

        zq(ig,l,igcm_co2) = zq(ig,l,igcm_co2)  + zq(ig,l,igcm_co2_ice)

        zq(ig,l,igcm_co2_ice) = 0.
        zq(ig,l,igcm_ccnco2_mass) = 0.
        zq(ig,l,igcm_ccnco2_number) = 0.
        mem_Nccn_co2(ig,l) = 0.
        mem_Mh2o_co2(ig,l) = 0.
        mem_Mccn_co2(ig,l) = 0.
        riceco2(ig,l) = 0.
      end if !of if co2_ice < threshold or zq(ccnco2_number) < 1
    end do ! of ig loop
  end do ! of nlayer loop
!----------------------------------------------------------------------------------------------------------------------!
! 5. Get cloud tendencies
!----------------------------------------------------------------------------------------------------------------------!
  subpdqcloudco2(:,:,igcm_co2) =   ( zq(:,:,igcm_co2) - zq0(:,:,igcm_co2) ) / microtimestep

  subpdqcloudco2(:,:,igcm_co2_ice) =  ( zq(:,:,igcm_co2_ice) - zq0(:,:,igcm_co2_ice) ) / microtimestep

  subpdqcloudco2(:,:,igcm_ccnco2_mass) =  ( zq(:,:,igcm_ccnco2_mass) - zq0(:,:,igcm_ccnco2_mass) ) / microtimestep

  subpdqcloudco2(:,:,igcm_ccnco2_number) = ( zq(:,:,igcm_ccnco2_number) - zq0(:,:,igcm_ccnco2_number))/microtimestep

  subpdqcloudco2(:,:,igcm_dust_mass) = ( zq(:,:,igcm_dust_mass) - zq0(:,:,igcm_dust_mass) )  / microtimestep

  subpdqcloudco2(:,:,igcm_dust_number) = ( zq(:,:,igcm_dust_number) - zq0(:,:,igcm_dust_number) )   / microtimestep

  if (co2useh2o) then
    subpdqcloudco2(:,:,igcm_h2o_ice) = ( zq(:,:,igcm_h2o_ice) - zq0(:,:,igcm_h2o_ice) )  / microtimestep

    subpdqcloudco2(:,:,igcm_ccn_mass) = ( zq(:,:,igcm_ccn_mass) - zq0(:,:,igcm_ccn_mass) )    / microtimestep

    subpdqcloudco2(:,:,igcm_ccn_number) = ( zq(:,:,igcm_ccn_number) - zq0(:,:,igcm_ccn_number) ) / microtimestep
  end if
!======================================================================================================================!
! END =================================================================================================================!
!======================================================================================================================!
  end subroutine improvedCO2clouds

end module improvedCO2clouds_mod

