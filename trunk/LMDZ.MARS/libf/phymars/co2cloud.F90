!======================================================================================================================!
! Module: CO2 clouds formation ========================================================================================!
!----------------------------------------------------------------------------------------------------------------------!
! Authors: Joachim Audouard, Constantino Listowski, Anni Määttänen
! Date: 09/2016
!----------------------------------------------------------------------------------------------------------------------!
! Contains subroutines:
!      - co2cloud: of co2 cloud microphysics
!
!      - ini_co2cloud: (only used in phys_state_var_init_mod.F90)
!
!      - end_co2cloud: (only used in phys_state_var_init_mod.F90)
!======================================================================================================================!
module co2cloud_mod

implicit none

double precision, allocatable, save :: &
   mem_Mccn_co2(:,:), &! Memory of CCN mass of H2O and dust used by CO2
   mem_Mh2o_co2(:,:), &! Memory of H2O mass integred into CO2 crystal
   mem_Nccn_co2(:,:)   ! Memory of CCN number of H2O and dust used by CO2

contains
!======================================================================================================================!
! SUBROUTINE: co2cloud ================================================================================================!
!======================================================================================================================!
! Subject:
!---------
!   Main subroutine of co2 cloud microphysics
!----------------------------------------------------------------------------------------------------------------------!
! Comments:
!----------
!     - Adaptation of the water ice clouds scheme (with specific microphysics) of Montmessin, Navarro et al.
!
!     - Microphysics subroutine is improvedCO2clouds.F
!
!     - There is a time loop specific to cloud formation due to timescales smaller than the GCM integration timestep
!
!     - The microphysics time step is a fraction of the physical one
!
!     - The CO2 clouds tracers (co2_ice, ccn mass and concentration) are sedimented at each microtimestep. pdqs_sedco2
!       keeps track of the CO2 flux at the surface
!
!     - The subgrid Temperature distribution is modulated (0 or 1) by Spiga et al. (GRL 2012)
!
!     - Saturation Index to account for GW propagation or dissipation upwards
!
!     - 4D and column opacities are computed using Qext values at 1µm
!----------------------------------------------------------------------------------------------------------------------!
! Papers:
!--------
!   "Near-pure vapor condensation in the Martian atmosphere: CO2 ice crystal growth", Listowski et al. (2013), JGRE
!   "Modeling the microphysics of CO2 ice clouds within wave-induced cold pockets in the martian mesosphere", Listowski
!     et al. (2014), Icarus
!   "Global climate modeling of the Martian water cycle with improved microphysics and radiatively active water ice
!     clouds", Navarro et al. (2014), JGRE
!   "Martian GCM with complete CO2 clouds microphysics", Audouard et al. (2017), EPSC abstract
!----------------------------------------------------------------------------------------------------------------------!
! Algorithm:
!-----------
!   0. Firstcall
!     0.1. Initialization of microtimestep from imicroco2
!     0.2. Compute the radius grid of CO2 ice particles (rb_cldco2)
!     0.3. Read file 'optprop_co2ice_1mic.dat' to extract optical properties of CO2 ice at 1 micron (Qext)
!     0.4. Interpole the radius grid (rb_cldco2) to get the corresponding exctinction coefficients (Qext1bins)
!     0.5. Save the radius grid of CO2 particles (rb_cldco2)
!   1. Initialization
!   2. Compute mass and thickness layers
!   3. Define the sub-grid cloud (CLFvaryingCO2)
!     3.1. Representation of sub-grid CO2 ice clouds (CLFvaryingCO2 = True)
!       3.1.a. Saturation index CO2
!       3.1.b. Compute tcond
!       3.1.c. Compute cloud fraction in cells
!     3.2. No sub-grid cloud representation (CLFvaryingCO2 = False)
!   4. Microphysics of CO2 cloud formation
!     4.1. Stepped entry for tendancies: At each micro timestep we add pdt in order to have a stepped entry
!     4.2. Effective tracers quantities in the cloud
!     4.3. Gravitational sedimentation
!       4.3.a. Compute cloud density
!       4.3.b. Save actualized tracer values to compute sedimentation tendancies
!       4.3.c. Sedimentation of co2 ice
!       4.3.d. Sedimentation for other tracers
!       4.3.e. Compute tendencies due to the sedimation process
!     4.4. Main call to the cloud scheme
!     4.5. Updating tendencies after cloud scheme
!   5. Compute final tendencies after time loop
!   6. Update clouds physical values in the cloud (for output)
!     6.1. Update density of co2 ice, riceco2 and opacity
!     6.2. Update rice and rdust
!   7. Correction if a lot of subliming CO2 fills the 1st layer
!   8. Compute water cloud sedimentation radius
!   9. CO2 saturated quantities
!     9.1 Compute CO2 saturation in layers
!     9.2 Compute CO2 saturated quantities in layers
!   10. Everything modified by CO2 microphysics must be wrt co2cloudfrac
!   11. Compute opacity at 1 micron
!   12. Write outputs in diagfi.nc
!======================================================================================================================!
  subroutine co2cloud(ngrid, nlay, ptimestep, pplev, pplay, pdpsrf, pzlay, pt, pdt, pq, pdq, pdqcloudco2, pdtcloudco2, &
                      nq, tau, tauscaling, rdust, rice, riceco2, nuice, rhocloud, rsedcloudco2, rhocloudco2, pzlev, pdqs_sedco2, pdu, &
                      pu, pcondicea, co2ice)

  use ioipsl_getincom, only: getin
  use dimradmars_mod,  only: naerkind
  use comcstfi_h,      only: pi, g, cpp
  use updaterad,       only: updaterice_microco2, updaterice_micro, updaterdust
  use conc_mod,        only: mmean, rnew
  use tracer_mod,      only: igcm_co2, igcm_co2_ice, igcm_dust_mass, igcm_dust_number, igcm_h2o_ice, &
                             igcm_ccn_mass, igcm_ccn_number, igcm_ccnco2_mass, igcm_ccnco2_number, rho_dust, &
                             nuiceco2_sed, nuiceco2_ref, r3n_q, rho_ice, nuice_sed

  use newsedim_mod,    only: newsedim

  use datafile_mod,    only: datadir

  use improvedCO2clouds_mod, only: improvedCO2clouds

#ifndef MESOSCALE
  use vertical_layers_mod, only: ap, bp
#endif

  implicit none

  include "callkeys.h"
  include "microphys.h"
!----------------------------------------------------------------------------------------------------------------------!
! VARIABLES DECLARATION
!----------------------------------------------------------------------------------------------------------------------!
! Input arguments:
!-----------------
  integer, intent(in) ::&
     ngrid, &! Number of grid points
     nlay,  &! Number of layers
     nq      ! Number of tracers

  real, intent(in) :: &
     ptimestep,           &! Physical time step (s)
     pplev(ngrid,nlay+1), &! Inter-layer pressures (Pa)
     pplay(ngrid,nlay),   &! Mid-layer pressures (Pa)
     pdpsrf(ngrid),       &! Tendency on surface pressure
     pzlay(ngrid,nlay),   &! Altitude at the middle of the layers
     pt(ngrid,nlay),      &! Temperature at the middle of the layers (K)
     pdt(ngrid,nlay),     &! Tendency on temperature from other parametrizations
     pq(ngrid,nlay,nq),   &! Tracers (kg/kg)
     pdq(ngrid,nlay,nq),  &! Tendencies before condensation (kg/kg.s-1)
     tau(ngrid,naerkind), &! Column dust optical depth at each point
     tauscaling(ngrid),   &! Convertion factor for dust amount
     pu(ngrid,nlay),      &! Zonal Wind: zu = pu + (pdu * ptimestep)
     pdu(ngrid,nlay),     &! Tendency of zonal wind before condensation
     pzlev(ngrid,nlay+1), &! Altitude at the boundaries of the layers
     nuice(ngrid,nlay),   &! Estimated effective variance of the size distribution
     co2ice(ngrid)         ! Amount of co2 ice at the surface
!----------------------------------------------------------------------------------------------------------------------!
! Output arguments:
!------------------
  real, intent(out) :: &
     rice(ngrid,nlay),          & ! Water Ice mass mean radius (m)
!     rsedcloud(ngrid,nlay),     & ! Water Cloud sedimentation radius
     rhocloud(ngrid,nlay),      & ! Water Cloud density (kg.m-3)
     pdqs_sedco2(ngrid),        & ! CO2 flux at the surface
     pdqcloudco2(ngrid,nlay,nq),& ! Tendency due to CO2 condensation (kg/kg.s-1)
     pcondicea(ngrid,nlay),     & ! Rate of condensation/sublimation of co2 ice in layers
     pdtcloudco2(ngrid,nlay),   & ! Tendency on temperature due to latent heat
     rsedcloudco2(ngrid,nlay)     ! Cloud sedimentation radius

  real, intent(inout) :: &
     rdust(ngrid,nlay) ! Dust geometric mean radius (m)

  double precision, intent(out) :: &
     riceco2(ngrid,nlay)          ! Ice mass mean radius (m) r_c in Montmessin et al. (2004)
!----------------------------------------------------------------------------------------------------------------------!
! Local:
!-------
!-----1) Parameters:
!-------------------
  integer, parameter :: &
     uQext = 555,         &! file_qext unit ID
     var_dim_qext = 10000  ! Exact dimension of radv and qextv1mic from file_qext

  real, parameter :: &
     mincloud = 0.1,  &! Minimum cloud fraction
     beta = 0.85,     &! correction for the shape of the particles (see Murphy et al. JGR 1990 vol.95):
                       !   beta = 1    for spheres
                       !   beta = 0.85 for irregular particles
                       !   beta = 0.5  for disk shaped particles
     threshold = 1e-30 ! limit value

  double precision, parameter :: &
     rmin_cld = 1.e-9,  &! Minimum radius (m)
     rmax_cld = 5.e-6,  &! Maximum radius (m)
     rbmin_cld = 1.e-10,&! Minimum boundary radius (m)
     rbmax_cld = 2.e-4, &! Maximum boundary radius (m)
     Fo = 7.5e-7,       &! for sat index  (J.m-3)
     lambdaH = 150.e3    ! for sat index (km)

  character(len=23), parameter :: &
     file_qext = 'optprop_co2ice_1mic.dat' ! File extinction coefficients of CO2 particles
!----------------------------------------------------------------------------------------------------------------------!
!-----2) Saved:
!--------------
  integer, save :: &
     imicroco2 ! Time subsampling for coupled water microphysics sedimentation microtimestep timeloop for microphysics:
               !   if imicroco2 = 1, subpdt is the same as pdt
  real, save :: &
     sigma_iceco2, &! Variance of the ice and CCN distributions
     microtimestep  ! Integration timestep for coupled water microphysics & sedimentation

  double precision, save :: &
     dev2,                   &! 1. / ( sqrt(2.) * sigma_iceco2 )
     Qext1bins(nbinco2_cld), &! Extinction coefficients for rb_cldco2 radius of CO2 ice particles
     Qextv1mic(var_dim_qext), &
     radv(var_dim_qext),      &    ! radius of CO2 ice at 1 µm (read from file_qext)
     rb_cldco2(nbinco2_cld+1) ! boundary values of each rad_cldco2 bin (m)
  logical, save :: &
     firstcall = .true. ! Used to compute saved variables
!----------------------------------------------------------------------------------------------------------------------!
!-----3) Variables:
!------------------
  integer :: &
     iq,       &! loop on tracers
     ig,       &! loop on grid points
     l,        &! loop on layers
     i,        &! loop on nbinco2_cld
     nelem,    &! number of point between lebon1 and lebon2 => interpolation
     lebon1,   &! bound limit for the interpolation
     lebon2,   &! bound limit for the interpolation
     microstep  ! Time subsampling step variable

  real :: &
! ---Tendency given by clouds inside the micro loop
     subpdqcloudco2(ngrid,nlay,nq), &! On tracers, cf. pdqcloud
     subpdtcloudco2(ngrid,nlay),    &! On temperature, cf. pdtcloud
!  ---Global tendency (clouds+physics)
     sum_subpdq(ngrid,nlay,nq),     &! On tracers, cf. pdqcloud
     sum_subpdt(ngrid,nlay),        &! On temperature, cf. pdtcloud
!  ---Sedimentation
     ztsed(ngrid,nlay),             &! Temperature with real-time value in microtimeloop
     zqsed(ngrid,nlay,nq),          &! Tracers with real-time value in µloop
     zqsed0(ngrid,nlay,nq),         &! For sedimentation tendancy
     subpdqsed(ngrid,nlay,nq),      &! Tendancy due to sedimentation
     sum_subpdqs_sedco2(ngrid),     &! CO2 flux at the surface
!  ---For sub grid T distribution
     zt(ngrid,nlay),                &! Local value of temperature
     zq_co2vap(ngrid, nlay),        &! Local value of CO2 vap
     rhocloudco2t(ngrid, nlay),     &! Cloud density (kg.m-3)
!  ---For Saturation Index computation
     zdelt,                         &! Delta T for the temperature distribution
     co2cloudfrac(ngrid,nlay),      &! Cloud fraction used only with CLFvarying is true
!  ---Misc
     rhocloudco2(ngrid, nlay),      &! Cloud density (kg.m-3)
     Nccnco2,                       &! buffer: number of ccn used for co2 condensation
     Qccnco2,                       &! buffer: mass of ccn used for co2 condensation
     Niceco2,                       &! buffer: mmr co2 ice
     epaisseur(ngrid,nlay),         &! Layer thickness (m)
     masse(ngrid,nlay),             &! Layer  mass (kg.m-2)
     pteff(ngrid, nlay),            &! Effective temperature in the cloud
     pqeff(ngrid, nlay, nq),        &! Effective tracers quantities in the cloud
     wq(ngrid,nlay+1),              &! Displaced tracer mass (kg.m-2) during microtimestep
     satuco2(ngrid,nlay),           &! CO2 satu ratio for output diagfi
     zqsatco2(ngrid,nlay),          &! Saturation co2
     availco2,&
     masslayer, &
     tmp, a,b, &
     new_pdq(ngrid,nlay)
     
  double precision :: &
! ---Extinction coefficients at 1 micron of CO2 particles
     vrat_cld,                      &! Volume ratio
     n_derf,                        &! derf( (rb_cldco2(1)-log(riceco2(ig,l))) *dev2)
     Qtemp,                         &! mean value in the interval during the interpolation
     ltemp1(var_dim_qext),          &! abs(radv(:)-rb_cldco2(i))
     ltemp2(var_dim_qext),          &! abs(radv(:)-rb_cldco2(i+1))
     n_aer(nbinco2_cld),            &! -0.5 * Nccnco2*tauscaling(ig) * n_derf
     tau1mic(ngrid),                &! CO2 ice column opacity at 1µm
     Qext1bins2(ngrid,nlay),        &! CO2 ice opacities
! ---For Saturation Index computation
     rho,                           &! background density
     zu,                            &! absolute value of zonal wind field
     NN,                            &! N^2 static stability
     gradT,                         &! thermal gradient
     SatIndex(ngrid,nlay),          &! Saturation index S in Spiga 2012 paper, assuming like in the paper GW phase speed
                                     !  (stationary waves): c = 0 m.s-1, lambdaH = 150 km, Fo = 7.5e-7 J.m-3
     SatIndexmap(ngrid),            &! maxval(SatIndex(ig,12:26))
! ---Misc
     myT,                           &! Temperature scalar for co2 density computation
     tcond(ngrid,nlay)               ! CO2 condensation temperature

  logical :: &
     file_qext_ok ! Check if file_qext exists
!======================================================================================================================!
! BEGIN ===============================================================================================================!
!======================================================================================================================!
! 0. Firstcall
!----------------------------------------------------------------------------------------------------------------------!
  if (firstcall) then
    firstcall=.false.
!----------------------------------------------------------------------------------------------------------------------!
! 0.1. Initialization of microtimestep from imicroco2
!----------------------------------------------------------------------------------------------------------------------!
#ifdef MESOSCALE
    imicroco2 = 2
#else
    imicroco2 = 30
#endif
    call getin("imicroco2", imicroco2)
    microtimestep = ptimestep/real(imicroco2)
    sigma_iceco2 = sqrt(log(1.+nuiceco2_sed))
    dev2 = 1. / ( sqrt(2.) * sigma_iceco2 )
!----------------------------------------------------------------------------------------------------------------------!
! 0.2. Compute the radius grid of CO2 ice particles (rb_cldco2)
!        > the grid spacing is computed assuming a constant volume ratio between two consecutive bins; i.e. vrat_cld.
!   - rad_cldco2 is the primary radius grid used for microphysics computation.
!   - The grid spacing is computed assuming a constant volume ratio between two consecutive bins; i.e. vrat_cld.
!   - vrat_cld is determined from the boundary values of the size grid: rmin_cld and rmax_cld.
!   - The rb_cldco2 array contains the boundary values of each rad_cldco2 bin.
!----------------------------------------------------------------------------------------------------------------------!
    ! vrat_cld is determined from the boundary values of the size grid: rmin_cld and rmax_cld.
    vrat_cld = exp(log(rmax_cld/rmin_cld) / float(nbinco2_cld-1) * 3.)

    ! rad_cldco2 is the primary radius grid used for microphysics computation.
    rad_cldco2(1) = rmin_cld
    do i = 1, nbinco2_cld-1
      rad_cldco2(i+1) = rad_cldco2(i) * vrat_cld**(1./3.)
    end do

    ! rb_cldco2 array contains the boundary values of each rad_cldco2 bin.
    rb_cldco2(1) = rbmin_cld
    do i = 1, nbinco2_cld
      rb_cldco2(i+1) = ( (2.*vrat_cld) / (vrat_cld+1.) )**(1./3.) * rad_cldco2(i)
    end do
    rb_cldco2(nbinco2_cld+1) = rbmax_cld
!----------------------------------------------------------------------------------------------------------------------!
! 0.3. Read file 'optprop_co2ice_1mic.dat' to extract optical properties of CO2 ice at 1 micron (Qext)
!----------------------------------------------------------------------------------------------------------------------!
    ! get information about file_qext
    inquire(file=trim(datadir)//'/'//file_qext, exist=file_qext_ok)

    ! if file_qext is missing then stop
    if (.not. file_qext_ok) then
      write(*,*)'file'//file_qext//'should be in ', trim(datadir)
      call abort_physic('co2cloud', 'file missing', 1)
    end if

    ! read file_qext
    open(unit=uQext,file=trim(datadir)//'/'//file_qext, form='formatted')

    ! skip 1 line
    read(uQext,*)

    ! extract radv
    do i = 1, var_dim_qext
      read(uQext,'(E12.5)')radv(i)
    end do

    ! skip 1 line
    read(uQext,*)

    ! Qextv1mic
    do i = 1 , var_dim_qext
      read(uQext,'(E12.5)')Qextv1mic(i)
    end do

    ! close file_qext
    close(uQext)
!----------------------------------------------------------------------------------------------------------------------!
! 0.4. Interpole the radius grid (rb_cldco2) to get the corresponding exctinction coefficients (Qext1bins), using
!       file_qext values (radv, Qextv1mic)
!----------------------------------------------------------------------------------------------------------------------!
    do i = 1, nbinco2_cld
      ltemp1 = abs(radv(:)-rb_cldco2(i))
      ltemp2 = abs(radv(:)-rb_cldco2(i+1))
      lebon1 = minloc(ltemp1,DIM=1)
      lebon2 = min(minloc(ltemp2,DIM=1), var_dim_qext)
      nelem = lebon2 - lebon1 + 1.

      ! mean value in the interval
      Qtemp = 0d0
      do l = 0, nelem
        Qtemp = Qtemp + Qextv1mic(min(lebon1+l, var_dim_qext))
      end do

      Qext1bins(i) = Qtemp / nelem
    end do

    Qext1bins(:) = Qext1bins(:) * pi * (rad_cldco2(:)**2)

    ! print result of the interpolation
    write(*,*)'--------------------------------------------'
    write(*,*)'Microphysics co2: size bin-Qext information:'
    write(*,*)'   i, rad_cldco2(i), Qext1bins(i)'
    do i = 1, nbinco2_cld
      write(*,'(i3,3x,3(e13.6,4x))')i, rad_cldco2(i), Qext1bins(i)
    end do
    write(*,*)'--------------------------------------------'
!----------------------------------------------------------------------------------------------------------------------!
! 0.5. Save the radius grid of CO2 particles (rb_cldco2)
!----------------------------------------------------------------------------------------------------------------------!
    do i = 1, nbinco2_cld+1
      rb_cldco2(i) = log(rb_cldco2(i))
    end do
  end if ! of IF (firstcall)
!----------------------------------------------------------------------------------------------------------------------!
! 1. Initialization
!----------------------------------------------------------------------------------------------------------------------!
  sum_subpdq(1:ngrid,1:nlay,1:nq) = 0.
  sum_subpdt(1:ngrid,1:nlay) = 0.

  subpdqcloudco2(1:ngrid,1:nlay,1:nq) = 0.
  subpdtcloudco2(1:ngrid,1:nlay) = 0.

  pdqcloudco2(1:ngrid,1:nlay,1:nq) = 0.
  pdtcloudco2(1:ngrid,1:nlay) = 0.

  ! default value if no ice
  rhocloudco2(1:ngrid,1:nlay) = rho_dust
  rhocloudco2t(1:ngrid,1:nlay) = rho_dust
  epaisseur(1:ngrid,1:nlay) = 0.
  masse(1:ngrid,1:nlay) = 0.
  riceco2(1:ngrid, 1:nlay) = 0.
  zqsed0(1:ngrid,1:nlay,1:nq) = 0.
  sum_subpdqs_sedco2(1:ngrid) = 0.
  subpdqsed(1:ngrid,1:nlay,1:nq) = 0.
!----------------------------------------------------------------------------------------------------------------------!
! 2. Compute mass and thickness layers
!----------------------------------------------------------------------------------------------------------------------!
  do l = 1, nlay
    do ig = 1, ngrid
#ifndef MESOSCALE
      masse(ig,l) = (pplev(ig,l) - pplev(ig,l+1) + (bp(l)-bp(l+1)) ) / g
#else
      masse(ig,l) = (pplev(ig,l) - pplev(ig,l+1)) / g
#endif
      epaisseur(ig,l) = pzlev(ig,l+1) - pzlev(ig,l)
    end do
  end do
!----------------------------------------------------------------------------------------------------------------------!
! 3. Define the sub-grid cloud (CLFvaryingCO2)
!----------------------------------------------------------------------------------------------------------------------!
! 3.1. Representation of sub-grid CO2 ice clouds (CLFvaryingCO2 = True)
!----------------------------------------------------------------------------------------------------------------------!
  if (CLFvaryingCO2) then
    ! effective temperature
    pteff(:,:) = pt(:,:)

    ! min co2cloudfrac when there is ice
    co2cloudfrac(:,:) = mincloud

    ! temperature
    do l=1,nlay
      do ig=1,ngrid
        zt(ig,l) = pt(ig,l) + pdt(ig,l)*ptimestep
      end do
    end do

    ! Quantities of traceurs
    if (igcm_co2 /= 0) then
      do l = 1, nlay
        do ig = 1, ngrid
          zq_co2vap(ig,l) = pq(ig,l,igcm_co2) + pdq(ig,l,igcm_co2)*ptimestep
        end do
      end do
    end if
!----------------------------------------------------------------------------------------------------------------------!
! 3.1.a. Saturation index CO2
!----------------------------------------------------------------------------------------------------------------------!
    ! if saturation index co2 is true
    if (satindexco2) then
      ! layers 12 --> 26 ~ 12->85 km
      do l = 12, 26
        do ig = 1, ngrid
          ! compute N^2 static stability
          gradT = (zt(ig,l+1)-zt(ig,l))/(pzlev(ig,l+1)-pzlev(ig,l))
          NN = sqrt(g/zt(iq,l) * (max(gradT,-g/cpp) + g/cpp))

          ! compute absolute value of zonal wind field
          zu = abs(pu(ig,l) + pdu(ig,l)*ptimestep)

          ! compute background density
          rho = pplay(ig,l) / (rnew(ig,l)*zt(ig,l))

          ! saturation index: Modulate the DeltaT by GW propagation index:
          ! --------------------------------------------------------------
          SatIndex(ig,l) = sqrt(Fo*lambdaH/(2.*pi)*NN / (rho*zu**3) )
        end do
      end do

      ! Then compute Satindex map in layers 12 --> 26 ~ 12->85 km
      do ig = 1, ngrid
        SatIndexmap(ig) = maxval(SatIndex(ig,12:26))
      end do

      ! Write outputs in diagfi.nc
      call WRITEDIAGFI(ngrid, "SatIndex", "SatIndex", " ", 3, SatIndex)

      call WRITEDIAGFI(ngrid, "SatIndexmap", "SatIndexmap", "km", 2, SatIndexmap)
    !------------------------------------------------------------------------------------------------------------------!
    ! if saturation index co2 is false, set saturation index to 0.05
    !------------------------------------------------------------------------------------------------------------------!
    else
      do ig = 1, ngrid
        SatIndexmap(ig)=0.05
      end do
    end if ! of if (satindexco2)
!----------------------------------------------------------------------------------------------------------------------!
! 3.1.b. Compute tcond
!----------------------------------------------------------------------------------------------------------------------!
    call tcondco2(ngrid,nlay,pplay,zq_co2vap,tcond)
!----------------------------------------------------------------------------------------------------------------------!
! 3.1.c. Compute cloud fraction in cells
!----------------------------------------------------------------------------------------------------------------------!
    do ig = 1, ngrid
      if (SatIndexmap(ig) <= 0.1) then
        do l = 1, nlay-1

          ! The entire fraction is saturated
          if (tcond(ig,l) >= (zt(ig,l)+zdelt) .or. tcond(ig,l) <= 0.) then
            pteff(ig,l) = zt(ig,l)
            co2cloudfrac(ig,l) = 1.

          ! No saturation at all
          else if (tcond(ig,l) <= (zt(ig,l)-zdelt)) then
            pteff(ig,l) = zt(ig,l) - zdelt
            co2cloudfrac(ig,l) = mincloud

          ! Mean temperature of the cloud fraction
          else
            pteff(ig,l) = (tcond(ig,l)+zt(ig,l)-zdelt) / 2.
            co2cloudfrac(ig,l) = (tcond(ig,l)-zt(ig,l)+zdelt) / (2.0*zdelt)
          end if

          pteff(ig,l) = pteff(ig,l) - pdt(ig,l)*ptimestep

          ! check boundary values of co2cloudfrac
          if (co2cloudfrac(ig,l) <= mincloud) then
            co2cloudfrac(ig,l) = mincloud
          else if (co2cloudfrac(ig,l)> 1) then
            co2cloudfrac(ig,l) = 1.
          end if
        end do

      ! SatIndex not favorable for GW: leave pt untouched
      else
        pteff(ig,l) = pt(ig,l)
        co2cloudfrac(ig,l) = mincloud
      end if ! of if (SatIndexmap <= 0.1)
    end do ! of ngrid
    ! TODO: Totalcloud frac of the column missing here
!----------------------------------------------------------------------------------------------------------------------!
! 3.2. No sub-grid cloud representation (CLFvarying = False)
!----------------------------------------------------------------------------------------------------------------------!
  else
    do l = 1, nlay
      do ig = 1, ngrid
        pteff(ig,l) = pt(ig,l)
      end do
    end do
  end if ! end if (CLFvaryingco2)
!----------------------------------------------------------------------------------------------------------------------!
! 4. Microphysics of CO2 cloud formation
!----------------------------------------------------------------------------------------------------------------------!
  pqeff(:,:,:) = pq(:,:,:)
  pteff(:,:) = pt(:,:)
!----------------------------------------------------------------------------------------------------------------------!
! 4.2. Effective tracers quantities in the cloud
!----------------------------------------------------------------------------------------------------------------------!
  if (CLFvaryingCO2) then
    pqeff(:,:,igcm_ccnco2_mass) = pq(:,:,igcm_ccnco2_mass) / co2cloudfrac(:,:)

    pqeff(:,:,igcm_ccnco2_number) = pq(:,:,igcm_ccnco2_number) / co2cloudfrac(:,:)

    pqeff(:,:,igcm_co2_ice) = pq(:,:,igcm_co2_ice) / co2cloudfrac(:,:)
  end if
!----------------------------------------------------------------------------------------------------------------------!
  do microstep = 1, imicroco2
!----------------------------------------------------------------------------------------------------------------------!
! 4.1. Stepped entry for tendancies: At each micro timestep we add pdt in order to have a stepped entry
!----------------------------------------------------------------------------------------------------------------------!
   do l = 1, nlay
     do ig = 1, ngrid
       ! on temperature
       sum_subpdt(ig,l) = sum_subpdt(ig,l) + pdt(ig,l)
        
       ! on tracers
       sum_subpdq(ig,l,igcm_dust_mass) = sum_subpdq(ig,l,igcm_dust_mass) + pdq(ig,l,igcm_dust_mass)

       sum_subpdq(ig,l,igcm_dust_number) = sum_subpdq(ig,l,igcm_dust_number) + pdq(ig,l,igcm_dust_number)

       sum_subpdq(ig,l,igcm_ccnco2_mass) = sum_subpdq(ig,l,igcm_ccnco2_mass) + pdq(ig,l,igcm_ccnco2_mass)

       sum_subpdq(ig,l,igcm_ccnco2_number) = sum_subpdq(ig,l,igcm_ccnco2_number) + pdq(ig,l,igcm_ccnco2_number)

       sum_subpdq(ig,l,igcm_co2_ice) = sum_subpdq(ig,l,igcm_co2_ice) + pdq(ig,l,igcm_co2_ice)

       sum_subpdq(ig,l,igcm_co2) = sum_subpdq(ig,l,igcm_co2) + pdq(ig,l,igcm_co2)

       if (co2useh2o) then
         sum_subpdq(ig,l,igcm_h2o_ice) = sum_subpdq(ig,l,igcm_h2o_ice) + pdq(ig,l,igcm_h2o_ice)

         sum_subpdq(ig,l,igcm_ccn_mass) = sum_subpdq(ig,l,igcm_ccn_mass) + pdq(ig,l,igcm_ccn_mass)

         sum_subpdq(ig,l,igcm_ccn_number) = sum_subpdq(ig,l,igcm_ccn_number) + pdq(ig,l,igcm_ccn_number)
       end if
     end do ! ngrid
   end do ! nlay
!----------------------------------------------------------------------------------------------------------------------!
! 4.4. Main call to the cloud scheme
!----------------------------------------------------------------------------------------------------------------------!
   call improvedco2clouds(ngrid, nlay, microtimestep, pplay, pplev, pteff, sum_subpdt, pqeff, sum_subpdq, &
                          subpdqcloudco2, subpdtcloudco2, nq, tauscaling, mem_Mccn_co2, mem_Mh2o_co2, mem_Nccn_co2, &
                          rb_cldco2, sigma_iceco2, dev2)

   do l = 1, nlay
     do ig = 1, ngrid
       if(pq(ig,l,igcm_co2_ice) + microtimestep*(sum_subpdq(ig,l,igcm_co2_ice)+subpdqcloudco2(ig,l,igcm_co2_ice)) &
         .le. 1.e-12) then
         subpdqcloudco2(ig,l,igcm_co2_ice) = -pq(ig,l,igcm_co2_ice)/microtimestep - sum_subpdq(ig,l,igcm_co2_ice)
         subpdqcloudco2(ig,l,igcm_co2) = -subpdqcloudco2(ig,l,igcm_co2_ice)
       end if

       if(pq(ig,l,igcm_co2) + microtimestep*(sum_subpdq(ig,l,igcm_co2)+subpdqcloudco2(ig,l,igcm_co2)) .le. 1.e-12) then
         subpdqcloudco2(ig,l,igcm_co2) =  - pq(ig,l,igcm_co2)/microtimestep - sum_subpdq(ig,l,igcm_co2)
         subpdqcloudco2(ig,l,igcm_co2_ice) = -subpdqcloudco2(ig,l,igcm_co2)
       end if
  ! ccnco2_number and ccnco2_mass
       if (((pq(ig,l,igcm_ccnco2_number)+(sum_subpdq(ig,l,igcm_ccnco2_number)+subpdqcloudco2(ig,l,igcm_ccnco2_number)) &
          *microtimestep).le.1.) .or. &
         (pq(ig,l,igcm_ccnco2_mass)+(sum_subpdq(ig,l,igcm_ccnco2_mass)+subpdqcloudco2(ig,l,igcm_ccnco2_mass)) &
          *microtimestep.le.1e-20)) then
         subpdqcloudco2(ig,l,igcm_ccnco2_number) = - pq(ig,l,igcm_ccnco2_number)/microtimestep + 1. &
                                                  - sum_subpdq(ig,l,igcm_ccnco2_number)
         subpdqcloudco2(ig,l,igcm_dust_number) = - subpdqcloudco2(ig,l,igcm_ccnco2_number)

         subpdqcloudco2(ig,l,igcm_ccnco2_mass) = - pq(ig,l,igcm_ccnco2_mass)/microtimestep + 1e-20 &
                                                - sum_subpdq(ig,l,igcm_ccnco2_mass)
         subpdqcloudco2(ig,l,igcm_dust_mass) = - subpdqcloudco2(ig,l,igcm_ccnco2_mass)
       end if
     end do
   end do
!----------------------------------------------------------------------------------------------------------------------!
! 4.5. Updating tendencies after cloud scheme
!----------------------------------------------------------------------------------------------------------------------!
    do l = 1, nlay
      do ig = 1, ngrid
        sum_subpdt(ig,l) = sum_subpdt(ig,l) + subpdtcloudco2(ig,l)

        sum_subpdq(ig,l,igcm_dust_mass) = sum_subpdq(ig,l,igcm_dust_mass) + subpdqcloudco2(ig,l,igcm_dust_mass)

        sum_subpdq(ig,l,igcm_dust_number) = sum_subpdq(ig,l,igcm_dust_number) + subpdqcloudco2(ig,l,igcm_dust_number)

        sum_subpdq(ig,l,igcm_ccnco2_mass) = sum_subpdq(ig,l,igcm_ccnco2_mass) + subpdqcloudco2(ig,l,igcm_ccnco2_mass)

        sum_subpdq(ig,l,igcm_ccnco2_number) = sum_subpdq(ig,l,igcm_ccnco2_number) + &
                                              subpdqcloudco2(ig,l,igcm_ccnco2_number)

        sum_subpdq(ig,l,igcm_co2_ice) = sum_subpdq(ig,l,igcm_co2_ice) + subpdqcloudco2(ig,l,igcm_co2_ice)

        sum_subpdq(ig,l,igcm_co2) = sum_subpdq(ig,l,igcm_co2) + subpdqcloudco2(ig,l,igcm_co2)

        if (co2useh2o) then
          sum_subpdq(ig,l,igcm_h2o_ice) = sum_subpdq(ig,l,igcm_h2o_ice) + subpdqcloudco2(ig,l,igcm_h2o_ice)

          sum_subpdq(ig,l,igcm_ccn_mass) = sum_subpdq(ig,l,igcm_ccn_mass) + subpdqcloudco2(ig,l,igcm_ccn_mass)

          sum_subpdq(ig,l,igcm_ccn_number) = sum_subpdq(ig,l,igcm_ccn_number) + subpdqcloudco2(ig,l,igcm_ccn_number)
        end if
      end do ! ngrid
    end do ! nlay
!----------------------------------------------------------------------------------------------------------------------!
! 4.3. Gravitational sedimentation
!----------------------------------------------------------------------------------------------------------------------!
    if (sedimentation) then
!----------------------------------------------------------------------------------------------------------------------!
! 4.3.a. Compute cloud density
!----------------------------------------------------------------------------------------------------------------------!
      do l = 1, nlay
        do ig = 1, ngrid
          ! temperature during the sedimentation process
          ztsed(ig,l) = pteff(ig,l) + sum_subpdt(ig,l) * microtimestep

          ! quantities tracers during the sedimentation process
          zqsed(ig,l,:) = pqeff(ig,l,:) + sum_subpdq(ig,l,:) * microtimestep

          ! assure positive value of co2_ice mmr, ccnco2 number, ccnco2 mass
          Niceco2 = max(zqsed(ig,l,igcm_co2_ice), threshold)
          Nccnco2 = max(zqsed(ig,l,igcm_ccnco2_number), threshold)
          Qccnco2 = max(zqsed(ig,l,igcm_ccnco2_mass), threshold)

          ! Get density cloud and co2 ice particle radius
          if (Niceco2.ne.0d0) then
            call updaterice_microco2(dble(Niceco2), dble(Qccnco2), dble(Nccnco2), ztsed(ig,l), tauscaling(ig), &
                                     riceco2(ig,l), rhocloudco2t(ig,l))
          else
            riceco2(ig,l) = 0.
            rhocloudco2t(ig,l) = 0.
          end if
        end do ! ngrid
      end do ! nlay
!----------------------------------------------------------------------------------------------------------------------!
! 4.3.b. Save actualized tracer values to compute sedimentation tendancies
!----------------------------------------------------------------------------------------------------------------------!
      zqsed0(:,:,igcm_co2_ice) = zqsed(:,:,igcm_co2_ice)
      zqsed0(:,:,igcm_ccnco2_mass) = zqsed(:,:,igcm_ccnco2_mass)
      zqsed0(:,:,igcm_ccnco2_number) = zqsed(:,:,igcm_ccnco2_number)
!----------------------------------------------------------------------------------------------------------------------!
! 4.3.c. Sedimentation of co2 ice
!----------------------------------------------------------------------------------------------------------------------!
      do ig = 1, ngrid
        do l = 1, nlay
          rsedcloudco2(ig,l) = max( riceco2(ig,l)*(1.+sigma_iceco2)*(1.+sigma_iceco2)*(1.+sigma_iceco2), &
                                    rdust(ig,l) )
        end do
      end do

      wq(:,:) = 0.
      call newsedim(ngrid, nlay, ngrid*nlay, ngrid*nlay, microtimestep, pplev, masse, epaisseur, ztsed, &
                   rsedcloudco2, rhocloudco2t, zqsed(:,:,igcm_co2_ice), wq, beta)

      do ig = 1, ngrid
        sum_subpdqs_sedco2(ig) = sum_subpdqs_sedco2(ig) + wq(ig,1) / microtimestep !wq est en kg.m-2
      end do
!----------------------------------------------------------------------------------------------------------------------!
! 4.3.d. Sedimentation for other tracers
!----------------------------------------------------------------------------------------------------------------------!
      wq(:,:) = 0.
      ! for ccnco2_mass
      call newsedim(ngrid, nlay, ngrid*nlay, ngrid*nlay, microtimestep, pplev, masse, epaisseur, ztsed, &
                   rsedcloudco2, rhocloudco2t, zqsed(:,:,igcm_ccnco2_mass), wq, beta)
      !TODO: ajouter le calcule de la tendance a la surface comme co2ice

      wq(:,:) = 0.
      ! for ccnco2_number
      call newsedim(ngrid, nlay, ngrid*nlay, ngrid*nlay,microtimestep, pplev, masse, epaisseur, ztsed, &
                   rsedcloudco2, rhocloudco2t, zqsed(:,:,igcm_ccnco2_number), wq, beta)
      !TODO: ajouter le calcule de la tendance a la surface comme co2ice
!----------------------------------------------------------------------------------------------------------------------!
! 4.3.e. Compute tendencies due to the sedimation process
!----------------------------------------------------------------------------------------------------------------------!
      do l = 1, nlay
        do ig = 1, ngrid
          subpdqsed(ig,l,igcm_ccnco2_mass) = ( zqsed(ig,l,igcm_ccnco2_mass) - zqsed0(ig,l,igcm_ccnco2_mass)  ) &
                                             / microtimestep

          subpdqsed(ig,l,igcm_ccnco2_number) = ( zqsed(ig,l,igcm_ccnco2_number) - zqsed0(ig,l,igcm_ccnco2_number) )&
                                               / microtimestep

          subpdqsed(ig,l,igcm_co2_ice) = ( zqsed(ig,l,igcm_co2_ice) - zqsed0(ig,l,igcm_co2_ice) ) / microtimestep
        end do
      end do
      ! update subtimestep tendencies with sedimentation input
      do l = 1, nlay
        do ig = 1, ngrid
          sum_subpdq(ig,l,igcm_ccnco2_mass) = sum_subpdq(ig,l,igcm_ccnco2_mass) + subpdqsed(ig,l,igcm_ccnco2_mass)

          sum_subpdq(ig,l,igcm_ccnco2_number) = sum_subpdq(ig,l,igcm_ccnco2_number) + subpdqsed(ig,l,igcm_ccnco2_number)

          sum_subpdq(ig,l,igcm_co2_ice) = sum_subpdq(ig,l,igcm_co2_ice) + subpdqsed(ig,l,igcm_co2_ice)
        end do
      end do
    end if !(end if sedimentation)
  end do ! of do microstep = 1, imicroco2
!----------------------------------------------------------------------------------------------------------------------!
! 5. Compute final tendencies after time loop
!----------------------------------------------------------------------------------------------------------------------!
! condensation/sublimation rate in the atmosphere (kg.kg-1.s-1)
  do l = nlay, 1, -1
    do ig = 1, ngrid
      pcondicea(ig,l) = sum_subpdq(ig,l,igcm_co2_ice) / real(imicroco2)
    end do
  end do

  ! CO2 flux at surface (kg.m-2.s-1)
  do ig = 1, ngrid
    pdqs_sedco2(ig) = sum_subpdqs_sedco2(ig) / real(imicroco2)
  end do
  ! temperature tendency (T.s-1)
  do l = 1, nlay
    do ig = 1, ngrid
      pdtcloudco2(ig,l) = ( sum_subpdt(ig,l)/real(imicroco2) ) - pdt(ig,l)
    end do
  end do

  ! tracers tendencies
  do l = 1, nlay
    do ig = 1, ngrid
      pdqcloudco2(ig,l,igcm_co2) = 0. ! here is the trick, this tendency is computed in co2condens_mod4micro

      pdqcloudco2(ig,l,igcm_co2_ice) = ( sum_subpdq(ig,l,igcm_co2_ice) / real(imicroco2) ) - pdq(ig,l,igcm_co2_ice)

      pdqcloudco2(ig,l,igcm_ccnco2_mass) = ( sum_subpdq(ig,l,igcm_ccnco2_mass)/real(imicroco2) ) - &
                                           pdq(ig,l,igcm_ccnco2_mass)

      pdqcloudco2(ig,l,igcm_ccnco2_number) = ( sum_subpdq(ig,l,igcm_ccnco2_number) / real(imicroco2) ) - &
                                             pdq(ig,l,igcm_ccnco2_number)

      pdqcloudco2(ig,l,igcm_dust_mass) = ( sum_subpdq(ig,l,igcm_dust_mass) / real(imicroco2) ) - &
                                          pdq(ig,l,igcm_dust_mass)

      pdqcloudco2(ig,l,igcm_dust_number) = ( sum_subpdq(ig,l,igcm_dust_number)/real(imicroco2) ) - &
                                           pdq(ig,l,igcm_dust_number)

      if (co2useh2o) then
        pdqcloudco2(ig,l,igcm_h2o_ice) = ( sum_subpdq(ig,l,igcm_h2o_ice) / real(imicroco2) ) - &
                                          pdq(ig,l,igcm_h2o_ice)

        pdqcloudco2(ig,l,igcm_ccn_mass) = ( sum_subpdq(ig,l,igcm_ccn_mass) / real(imicroco2) ) - &
                                           pdq(ig,l,igcm_ccn_mass)

        pdqcloudco2(ig,l,igcm_ccn_number) = ( sum_subpdq(ig,l,igcm_ccn_number) / real(imicroco2) ) - &
                                            pdq(ig,l,igcm_ccn_number)
      end if
    end do ! ngrid
  end do ! nlay
!----------------------------------------------------------------------------------------------------------------------!
! 6. Update clouds physical values in the cloud (for output)
!----------------------------------------------------------------------------------------------------------------------!
! 6.1. Update density of co2 ice, riceco2 and opacity
!----------------------------------------------------------------------------------------------------------------------!
  do l = 1, nlay
    do ig = 1, ngrid
      Niceco2 = pqeff(ig,l,igcm_co2_ice) + ( pdq(ig,l,igcm_co2_ice) + pdqcloudco2(ig,l,igcm_co2_ice) ) * ptimestep
      Niceco2 = max (Niceco2, threshold)

      Nccnco2 = max( (pqeff(ig,l,igcm_ccnco2_number) + (pdq(ig,l,igcm_ccnco2_number) + &
                       pdqcloudco2(ig,l, igcm_ccnco2_number)) * ptimestep) &
                    , threshold)

      Qccnco2 = max( (pqeff(ig,l,igcm_ccnco2_mass) + (pdq(ig,l,igcm_ccnco2_mass) + &
                       pdqcloudco2(ig,l, igcm_ccnco2_mass)) * ptimestep)&
                    , threshold)

      myT = pteff(ig,l) + (pdt(ig,l)+pdtcloudco2(ig,l))*ptimestep

      ! Compute particle size
      call updaterice_microco2(dble(Niceco2), dble(Qccnco2), dble(Nccnco2), myT, tauscaling(ig), riceco2(ig,l), &
                               rhocloudco2(ig,l))

     ! Compute opacities
      if ( (Niceco2 <= threshold .or. Nccnco2*tauscaling(ig) <= 1.) ) then
        riceco2(ig,l) = 0.
        Qext1bins2(ig,l) = 0.
      else
        n_derf = derf( (rb_cldco2(1)-log(riceco2(ig,l))) *dev2)
        Qext1bins2(ig,l) = 0.

        do i = 1, nbinco2_cld
          n_aer(i) = -0.5 * Nccnco2*tauscaling(ig) * n_derf

          n_derf = derf((rb_cldco2(i+1)-log(riceco2(ig,l))) *dev2)
          n_aer(i) = n_aer(i) + (0.5 * Nccnco2*tauscaling(ig) * n_derf)

          Qext1bins2(ig,l) = Qext1bins2(ig,l) + Qext1bins(i) * n_aer(i)
        end do
      end if
!----------------------------------------------------------------------------------------------------------------------!
! 6.2. Update rice and rdust
!----------------------------------------------------------------------------------------------------------------------!
      ! update rice water only if co2 use h2o ice as CCN
      if (co2useh2o) then
        call updaterice_micro( &
              pqeff(ig,l,igcm_h2o_ice) + (pdq(ig,l,igcm_h2o_ice) + pdqcloudco2(ig,l,igcm_h2o_ice))*ptimestep, &
              pqeff(ig,l,igcm_ccn_mass) + (pdq(ig,l,igcm_ccn_mass) + pdqcloudco2(ig,l,igcm_ccn_mass))*ptimestep, &
              pqeff(ig,l,igcm_ccn_number) + (pdq(ig,l,igcm_ccn_number) + pdqcloudco2(ig,l,igcm_ccn_number))*ptimestep, &
              tauscaling(ig),rice(ig,l),rhocloud(ig,l))
      end if

      ! update rdust
      call updaterdust( &
           pqeff(ig,l,igcm_dust_mass) + (pdq(ig,l,igcm_dust_mass) + pdqcloudco2(ig,l,igcm_dust_mass))*ptimestep, &
           pqeff(ig,l,igcm_dust_number) + (pdq(ig,l,igcm_dust_number) + pdqcloudco2(ig,l,igcm_dust_number))*ptimestep, &
           rdust(ig,l))
    end do ! ngrid
  end do ! nlay
!----------------------------------------------------------------------------------------------------------------------!
! 7. A correction if a lot of subliming CO2 fills the 1st layer FF (04/2005). Then that should not affect the ice
!      particle radius
!----------------------------------------------------------------------------------------------------------------------!
  do ig = 1, ngrid
    if ( pdpsrf(ig)*ptimestep > 0.9*(pplev(ig,1)-pplev(ig,2))) then

      if ( pdpsrf(ig)*ptimestep > 0.9*(pplev(ig,1)-pplev(ig,3)) ) then
        riceco2(ig,2) = riceco2(ig,3)
      end if

      riceco2(ig,1) = riceco2(ig,2)
    end if
  end do
!----------------------------------------------------------------------------------------------------------------------!
! 9. CO2 saturated quantities
!----------------------------------------------------------------------------------------------------------------------!
! 9.1 Compute CO2 saturation in layers
!----------------------------------------------------------------------------------------------------------------------!
  call co2sat(ngrid*nlay, pteff+(pdt+pdtcloudco2)*ptimestep, zqsatco2)
!----------------------------------------------------------------------------------------------------------------------!
! 9.2 Compute CO2 saturated quantities in layers
!----------------------------------------------------------------------------------------------------------------------!
  do l = 1, nlay
    do ig = 1, ngrid
      satuco2(ig,l) = ( pqeff(ig,l,igcm_co2)  + (pdq(ig,l,igcm_co2) + pdqcloudco2(ig,l,igcm_co2))*ptimestep ) *  &
                      (mmean(ig,l)/(mco2*1e3)) * pplay(ig,l) / zqsatco2(ig,l)
    end do
  end do
!----------------------------------------------------------------------------------------------------------------------!
! 10. Everything modified by CO2 microphysics must be wrt co2cloudfrac
!----------------------------------------------------------------------------------------------------------------------!
  if (CLFvaryingCO2) then
    do l = 1, nlay
      do ig = 1, ngrid
        pdqcloudco2(ig,l,igcm_ccnco2_mass) = pdqcloudco2(ig,l,igcm_ccnco2_mass) * co2cloudfrac(ig,l)

        pdqcloudco2(ig,l,igcm_ccnco2_number) = pdqcloudco2(ig,l,igcm_ccnco2_number) * co2cloudfrac(ig,l)

        pdqcloudco2(ig,l,igcm_dust_mass) = pdqcloudco2(ig,l,igcm_dust_mass) * co2cloudfrac(ig,l)

        pdqcloudco2(ig,l,igcm_dust_number) = pdqcloudco2(ig,l,igcm_dust_number) * co2cloudfrac(ig,l)

        pdqcloudco2(ig,l,igcm_co2_ice) = pdqcloudco2(ig,l,igcm_co2_ice) * co2cloudfrac(ig,l)

        pdqcloudco2(ig,l,igcm_co2) = pdqcloudco2(ig,l,igcm_co2) * co2cloudfrac(ig,l)

        pdtcloudco2(ig,l) = pdtcloudco2(ig,l) * co2cloudfrac(ig,l)

        Qext1bins2(ig,l) = Qext1bins2(ig,l) * co2cloudfrac(ig,l)

        if (co2useh2o) then
          pdqcloudco2(ig,l,igcm_h2o_ice) = pdqcloudco2(ig,l,igcm_h2o_ice) * co2cloudfrac(ig,l)

          pdqcloudco2(ig,l,igcm_ccn_mass) = pdqcloudco2(ig,l,igcm_ccn_mass) * co2cloudfrac(ig,l)

          pdqcloudco2(ig,l,igcm_ccn_number) = pdqcloudco2(ig,l,igcm_ccn_number) * co2cloudfrac(ig,l)
        end if
      end do ! ngrid
    end do ! nlay
  end if ! if CLFvaryingCO2 is true
!----------------------------------------------------------------------------------------------------------------------!
! 11. Compute opacity at 1 micron: Opacity in mesh ig is the sum over l of Qext1bins2. Is this true ?
!----------------------------------------------------------------------------------------------------------------------!
  tau1mic(:)=0.
  do l = 1, nlay
    do ig = 1, ngrid
      tau1mic(ig) = tau1mic(ig) + Qext1bins2(ig,l)
    end do
  end do
!----------------------------------------------------------------------------------------------------------------------!
! 12. Write outputs in diagfi.nc
!----------------------------------------------------------------------------------------------------------------------!
  call WRITEDIAGFI(ngrid, "satuco2", "vap in satu", " ", 3, satuco2)

  call WRITEDIAGFI(ngrid, "precip_co2_ice", "surface deposition of co2 ice", "kg.m-2", 2,  pdqs_sedco2(:)*ptimestep)

  call WRITEDIAGFI(ngrid, "precip_co2_ice_rate", "surface deposition rate of co2 ice", "kg.m-2.s-1", 2, pdqs_sedco2(:))

  call WRITEDIAGFI(ngrid, "co2ice_cond_rate", "CO2 condensation rate in atm layers", "kg.kg-1.s-1", 3, pcondicea)

  call WRITEDIAGFI(ngrid, "pdtcloudco2", "temperature variation of CO2 latent heat", "K.s-1", 3, pdtcloudco2)

  call writediagfi(ngrid, "riceco2", "ice radius", "m", 3, riceco2)

  call WRITEDIAGFI(ngrid, "Tau3D1mic", " co2 ice opacities", " ", 3, Qext1bins2)

  call WRITEDIAGFI(ngrid, "tau1mic", "co2 ice opacity 1 micron", " ", 2, tau1mic)

  call WRITEDIAGFI(ngrid, "mem_Nccn_co2", "CCN number used by CO2", "kg.kg-1", 3, mem_Nccn_co2)

  call WRITEDIAGFI(ngrid, "mem_Mccn_co2", "CCN mass used by CO2",  "kg.kg-1", 3, mem_Mccn_co2)

  call WRITEDIAGFI(ngrid, "mem_Mh2o_co2", "H2O mass in CO2 crystal", "kg.kg-1", 3, mem_Mh2o_co2)
  if (CLFvaryingCO2) then
    call WRITEDIAGFI(ngrid, "co2cloudfrac", "co2 cloud fraction", " ", 3, co2cloudfrac)
  end if
!======================================================================================================================!
! END =================================================================================================================!
!======================================================================================================================!
  end subroutine co2cloud


!**********************************************************************************************************************!
!**********************************************************************************************************************!


!======================================================================================================================!
! SUBROUTINE: ini_co2cloud ============================================================================================!
!======================================================================================================================!
! Subject:
!---------
!   Allocate arrays used for co2 microphysics
!======================================================================================================================!
  subroutine ini_co2cloud(ngrid,nlayer)

  implicit none

  integer, intent(in) :: &
     ngrid, &! number of atmospheric columns
     nlayer  ! number of atmospheric layers

  allocate(mem_Nccn_co2(ngrid,nlayer))
  allocate(mem_Mccn_co2(ngrid,nlayer))
  allocate(mem_Mh2o_co2(ngrid,nlayer))

  end subroutine ini_co2cloud

!**********************************************************************************************************************!
!**********************************************************************************************************************!

!======================================================================================================================!
! SUBROUTINE: end_co2cloud ============================================================================================!
!======================================================================================================================!
! Subject:
!---------
!   Deallocate arrays used for co2 microphysics
!======================================================================================================================!
  subroutine end_co2cloud

  implicit none

  if (allocated(mem_Nccn_co2)) deallocate(mem_Nccn_co2)
  if (allocated(mem_Mccn_co2)) deallocate(mem_Mccn_co2)
  if (allocated(mem_Mh2o_co2)) deallocate(mem_Mh2o_co2)

  end subroutine end_co2cloud

end module co2cloud_mod
