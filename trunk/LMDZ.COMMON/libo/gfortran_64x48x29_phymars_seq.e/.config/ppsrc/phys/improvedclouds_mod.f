










      MODULE improvedclouds_mod

      IMPLICIT NONE

      CONTAINS

      subroutine improvedclouds(ngrid,nlay,microtimestep,
     &             pplay,pteff,sum_subpdt,
     &             pqeff,sum_subpdq,subpdqcloud,subpdtcloud,
     &             nq,tauscaling)
      USE updaterad, ONLY: updaterice_micro, updaterccn
      USE watersat_mod, ONLY: watersat
      use tracer_mod, only: rho_ice, nuice_sed, igcm_h2o_vap,
     &                      igcm_h2o_ice, igcm_dust_mass,
     &                      igcm_dust_number, igcm_ccn_mass,
     &                      igcm_ccn_number,
     &                      igcm_hdo_vap,igcm_hdo_ice,
     &                      qperemin
      use conc_mod, only: mmean
      use comcstfi_h, only: pi, cpp
      implicit none
      
      
c------------------------------------------------------------------
c  This routine is used to form clouds when a parcel of the GCM is
c    saturated. It includes the ability to have supersaturation, a
c    computation of the nucleation rates, growthrates and the
c    scavenging of dust particles by clouds.
c  It is worth noting that the amount of dust is computed using the
c    dust optical depth computed in aeropacity.F. That's why
c    the variable called "tauscaling" is used to convert
c    pq(dust_mass) and pq(dust_number), which are relative
c    quantities, to absolute and realistic quantities stored in zq.
c    This has to be done to convert the inputs into absolute
c    values, but also to convert the outputs back into relative
c    values which are then used by the sedimentation and advection
c    schemes.

c  Authors: J.-B. Madeleine, based on the work by Franck Montmessin
c           (October 2011)
c           T. Navarro, debug,correction, new scheme (October-April 2011)
c           A. Spiga, optimization (February 2012)
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
!-----------------------------------------------------------------------
! INCLUDE 'microphys.h'
! Parameters and physical constants used by the microphysal scheme;
! Parameters for CO2 microphysics are also in this file
!-----------------------------------------------------------------------

!     Number of bins
      INTEGER, PARAMETER :: nbin_cld = 5

!     Reference temperature, T=273.15 K
      REAL, PARAMETER :: To = 273.15
!     Avogadro number
      DOUBLE PRECISION, PARAMETER :: nav = 6.023d23
!     Perfect gas constant
      DOUBLE PRECISION, PARAMETER :: rgp = 8.3143
!     Boltzman constant
      DOUBLE PRECISION, PARAMETER :: kbz = 1.381d-23
!     Molecular weight of H2O (kg.mol-1)
      DOUBLE PRECISION, PARAMETER :: mh2o = 18.d-3
!     Molecular weight of HDO (kg.mol-1)
      DOUBLE PRECISION, PARAMETER :: mhdo = 19.d-3
!     Molecular weight of CO2 (kg.mol-1)
      DOUBLE PRECISION, PARAMETER :: mco2 = 44.d-3
!     Molecular weight of N2 (kg.mol-1)
      DOUBLE PRECISION, PARAMETER :: mn2 = 28.01d-3
!     Effective CO2 gas molecular radius (m)
  !    bachnar 2016 value :1.97d-10   ! old value = 2.2d-10
      DOUBLE PRECISION, PARAMETER :: molco2 = 1.97d-10
!     Effective H2O gas molecular radius (m)
      DOUBLE PRECISION, PARAMETER :: molh2o = 1.2d-10
!     Effective HDO gas molecular radius (m)
      DOUBLE PRECISION, PARAMETER :: molhdo = 1.2d-10
!     Surface tension of ice/vapor (N.m)
      DOUBLE PRECISION, PARAMETER :: sigh2o = 0.12
!     Activation energy for desorption of
!       water on a dust-like substrate
!       (J/molecule)
      DOUBLE PRECISION, PARAMETER :: desorp = 0.288e-19
!     Jump frequency of a water molecule (s-1)
      DOUBLE PRECISION, PARAMETER :: nus = 1.e+13
!     Estimated activation energy for
!       surface diffusion of water molecules
!       (J/molecule)
      DOUBLE PRECISION, PARAMETER :: surfdif = desorp / 10.
!     Weight of a water molecule (kg)
      DOUBLE PRECISION, PARAMETER :: m0 = mh2o / nav

!     Contact parameter ( m=cos(theta) )
!       (initialized in improvedclouds.F)
      REAL mteta

!     Volume of a water molecule (m3)
      DOUBLE PRECISION vo1
!     Radius used by the microphysical scheme (m)
      DOUBLE PRECISION rad_cld(nbin_cld)




!CO2 part
!      number of bins for nucleation
      INTEGER, PARAMETER :: nbinco2_cld=100
!     Surface tension of ice/vapor (J.m-2)
      DOUBLE PRECISION, PARAMETER :: sigco2 = 0.08
!     Activation energy for desorption of
!       water on a dust-like substrate
!       (J/molecule)
      DOUBLE PRECISION, PARAMETER ::desorpco2=3.07e-20
!     bachnar 2016 value 3.07d-20 
!old value 3.20e-20
!     Jump frequency of a co2 molecule (s-1)
      DOUBLE PRECISION, PARAMETER :: nusco2 =  2.9e+12
!     Estimated activation energy for
!       surface diffusion of co2 molecules
!       (J/molecule)
      DOUBLE PRECISION, PARAMETER :: surfdifco2 = desorpco2 / 10.
!     Weight of a co2 molecule (kg)
      DOUBLE PRECISION, PARAMETER :: m0co2 = mco2 / nav
!     Contact parameter ( m=cos(theta) )
!       (initialized in improvedCO2clouds.F)
!    bachnar 2016 value :0.78 
!old value 0.95
      REAL, parameter :: mtetaco2 = 0.95
!     Volume of a co2 molecule (m3)
       DOUBLE PRECISION :: vo1co2
!     Radius used by the microphysical scheme (m)
      DOUBLE PRECISION :: rad_cldco2(nbinco2_cld)
       REAL, parameter :: threshJA = 1
!     COMMON/microphys/vo1co2,rad_cldco2

! NB: to keep commons aligned: 
!     split them in groups (reals, integers and characters)

      COMMON/microphys/rad_cld,vo1,rad_cldco2,vo1co2
		  COMMON/microphys_2/mteta
      
!     EXAMPLE:
!     COMMON/tracer/radius,rho_q,alpha_lift,alpha_devil,mmol,           &
!    & varian,r3n_q,rho_dust,rho_ice,nuice_ref,nuice_sed,               &
!    & ref_r0,dryness
!-----------------------------------------------------------------------
c------------------------------------------------------------------
c     Inputs/outputs:

      INTEGER, INTENT(IN) :: ngrid,nlay
      INTEGER, INTENT(IN) :: nq                 ! nombre de traceurs
      REAL, INTENT(IN) :: microtimestep         ! pas de temps physique (s)
      REAL, INTENT(IN) :: pplay(ngrid,nlay)     ! pression au milieu des couches (Pa)            
      REAL, INTENT(IN) :: pteff(ngrid,nlay)     ! temperature at the middle of the
                                                ! layers (K)
      REAL, INTENT(IN) :: sum_subpdt(ngrid,nlay)! tendance temperature des autres
                                                !   param.
      REAL, INTENT(IN) :: pqeff(ngrid,nlay,nq)  ! traceur (kg/kg)
      REAL, INTENT(IN) :: sum_subpdq(ngrid,nlay,nq)    ! tendance avant condensation
                                                       !   (kg/kg.s-1)
      REAL, INTENT(IN) :: tauscaling(ngrid)     ! Convertion factor for qdust and Ndust
      
      REAL, INTENT(OUT) :: subpdqcloud(ngrid,nlay,nq)  ! tendance de la condensation
                                                       !   H2O(kg/kg.s-1)
      REAL, INTENT(OUT) :: subpdtcloud(ngrid,nlay)     ! tendance temperature due
                                                       !   a la chaleur latente

c------------------------------------------------------------------
c     Local variables:

      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall

      REAL*8   derf ! Error function
      !external derf
     
      INTEGER ig,l,i
      
      REAL zq(ngrid,nlay,nq)  ! local value of tracers
      REAL zq0(ngrid,nlay,nq) ! local initial value of tracers
      REAL zt(ngrid,nlay)       ! local value of temperature
      REAL zqsat(ngrid,nlay)    ! saturation
      REAL lw                         !Latent heat of sublimation (J.kg-1) 
      REAL cste
      REAL dMice           ! mass of condensed ice
      REAL dMice_hdo       ! mass of condensed HDO ice
      REAL alpha(ngrid,nlay) ! HDO equilibrium fractionation coefficient (Saturation=1)
      REAL alpha_c(ngrid,nlay) ! HDO real fractionation coefficient
!      REAL sumcheck
      REAL*8 ph2o          ! Water vapor partial pressure (Pa)
      REAL*8 satu          ! Water vapor saturation ratio over ice
      REAL*8 Mo,No
      REAL*8 Rn, Rm, dev2, n_derf, m_derf
      REAL*8 n_aer(nbin_cld) ! number conc. of particle/each size bin
      REAL*8 m_aer(nbin_cld) ! mass mixing ratio of particle/each size bin

      REAL*8 sig      ! Water-ice/air surface tension  (N.m)
      EXTERNAL sig

      REAL dN,dM
      REAL rate(nbin_cld)  ! nucleation rate
      REAL seq

      REAL rice(ngrid,nlay)      ! Ice mass mean radius (m)
                                 ! (r_c in montmessin_2004)
      REAL rhocloud(ngrid,nlay)  ! Cloud density (kg.m-3)
      REAL rdust(ngrid,nlay) ! Dust geometric mean radius (m)

      REAL res      ! Resistance growth
      REAL Dv,Dv_hdo ! Water/HDO vapor diffusion coefficient
      

c     Parameters of the size discretization
c       used by the microphysical scheme
      DOUBLE PRECISION, PARAMETER :: rmin_cld = 0.1e-6 ! Minimum radius (m)
      DOUBLE PRECISION, PARAMETER :: rmax_cld = 10.e-6 ! Maximum radius (m)
      DOUBLE PRECISION, PARAMETER :: rbmin_cld = 0.0001e-6
                                           ! Minimum boundary radius (m)
      DOUBLE PRECISION, PARAMETER :: rbmax_cld = 1.e-2 ! Maximum boundary radius (m)
      DOUBLE PRECISION vrat_cld ! Volume ratio
      DOUBLE PRECISION rb_cld(nbin_cld+1)! boundary values of each rad_cld bin (m)
      SAVE rb_cld
      DOUBLE PRECISION dr_cld(nbin_cld)   ! width of each rad_cld bin (m)
      DOUBLE PRECISION vol_cld(nbin_cld)  ! particle volume for each bin (m3)


      REAL sigma_ice ! Variance of the ice and CCN distributions
      SAVE sigma_ice


      
c----------------------------------      
c TESTS

      INTEGER countcells
      
      LOGICAL test_flag    ! flag for test/debuging outputs
      SAVE    test_flag    


      REAL satubf(ngrid,nlay),satuaf(ngrid,nlay) 
      REAL res_out(ngrid,nlay)
 

c------------------------------------------------------------------

      ! AS: firstcall OK absolute
      IF (firstcall) THEN
!=============================================================
! 0. Definition of the size grid
!=============================================================
c       rad_cld is the primary radius grid used for microphysics computation.
c       The grid spacing is computed assuming a constant volume ratio
c       between two consecutive bins; i.e. vrat_cld.
c       vrat_cld is determined from the boundary values of the size grid: 
c       rmin_cld and rmax_cld.
c       The rb_cld array contains the boundary values of each rad_cld bin.
c       dr_cld is the width of each rad_cld bin.

c       Volume ratio between two adjacent bins
!        vrat_cld = log(rmax_cld/rmin_cld) / float(nbin_cld-1) *3.
!        vrat_cld = exp(vrat_cld)
        vrat_cld = log(rmax_cld/rmin_cld) / float(nbin_cld-1) *3.
        vrat_cld = exp(vrat_cld)
        write(*,*) "vrat_cld", vrat_cld

        rb_cld(1)  = rbmin_cld
        rad_cld(1) = rmin_cld
        vol_cld(1) = 4./3. * dble(pi) * rmin_cld*rmin_cld*rmin_cld
!        vol_cld(1) = 4./3. * pi * rmin_cld*rmin_cld*rmin_cld

        do i=1,nbin_cld-1
          rad_cld(i+1)  = rad_cld(i) * vrat_cld**(1./3.)
          vol_cld(i+1)  = vol_cld(i) * vrat_cld
        enddo
        
        do i=1,nbin_cld
          rb_cld(i+1)= ( (2.*vrat_cld) / (vrat_cld+1.) )**(1./3.) *
     &      rad_cld(i)
          dr_cld(i)  = rb_cld(i+1) - rb_cld(i)
        enddo
        rb_cld(nbin_cld+1) = rbmax_cld
        dr_cld(nbin_cld)   = rb_cld(nbin_cld+1) - rb_cld(nbin_cld)

        print*, ' '
        print*,'Microphysics: size bin information:'
        print*,'i,rb_cld(i), rad_cld(i),dr_cld(i)'
        print*,'-----------------------------------'
        do i=1,nbin_cld
          write(*,'(i2,3x,3(e12.6,4x))') i,rb_cld(i), rad_cld(i),
     &      dr_cld(i)
        enddo
        write(*,'(i2,3x,e12.6)') nbin_cld+1,rb_cld(nbin_cld+1)
        print*,'-----------------------------------'

        do i=1,nbin_cld+1
!           rb_cld(i) = log(rb_cld(i))  
            rb_cld(i) = log(rb_cld(i))  !! we save that so that it is not computed
                                         !! at each timestep and gridpoint
        enddo

c       Contact parameter of water ice on dust ( m=cos(theta) )
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       mteta  = 0.95
        write(*,*) 'water_param contact parameter:', mteta

c       Volume of a water molecule (m3)
        vo1 = mh2o / dble(rho_ice)
c       Variance of the ice and CCN distributions
        sigma_ice = sqrt(log(1.+nuice_sed))
       
 
        write(*,*) 'Variance of ice & CCN distribs :', sigma_ice
        write(*,*) 'nuice for sedimentation:', nuice_sed
        write(*,*) 'Volume of a water molecule:', vo1


        test_flag = .false.
 
        firstcall=.false.
      END IF


!=============================================================
! 1. Initialisation
!=============================================================
      cste = 4*pi*rho_ice*microtimestep

      res_out(:,:) = 0
      rice(:,:) = 1.e-8

c     Initialize the tendencies
      subpdqcloud(1:ngrid,1:nlay,1:nq)=0
      subpdtcloud(1:ngrid,1:nlay)=0
    
      
      zt(1:ngrid,1:nlay) = 
     &      pteff(1:ngrid,1:nlay) + 
     &      sum_subpdt(1:ngrid,1:nlay) * microtimestep

      zq(1:ngrid,1:nlay,1:nq) = 
     &      pqeff(1:ngrid,1:nlay,1:nq) + 
     &      sum_subpdq(1:ngrid,1:nlay,1:nq) * microtimestep
      
      
      WHERE( zq(1:ngrid,1:nlay,1:nq) < 1.e-30 )
     &       zq(1:ngrid,1:nlay,1:nq) = 1.e-30

      zq0(1:ngrid,1:nlay,1:nq) = zq(1:ngrid,1:nlay,1:nq)
      
!=============================================================
! 2. Compute saturation
!=============================================================


      dev2 = 1. / ( sqrt(2.) * sigma_ice )

      call watersat(ngrid*nlay,zt,pplay,zqsat)
            
      countcells = 0

c     Main loop over the GCM's grid
      DO l=1,nlay
        DO ig=1,ngrid

c       Get the partial pressure of water vapor and its saturation ratio
        ph2o = zq(ig,l,igcm_h2o_vap) * (mmean(ig,l)/18.) * pplay(ig,l)
        satu = zq(ig,l,igcm_h2o_vap) / zqsat(ig,l)

!=============================================================
! 3. Nucleation
!=============================================================

       IF ( satu .ge. 1. ) THEN         ! if there is condensation

        call updaterccn(zq(ig,l,igcm_dust_mass),
     &          zq(ig,l,igcm_dust_number),rdust(ig,l),tauscaling(ig))


c       Expand the dust moments into a binned distribution
        Mo = zq(ig,l,igcm_dust_mass)* tauscaling(ig)   + 1.e-30
        No = zq(ig,l,igcm_dust_number)* tauscaling(ig) + 1.e-30
        Rn = rdust(ig,l)
        Rn = -log(Rn) 
        Rm = Rn - 3. * sigma_ice*sigma_ice  
        n_derf = derf( (rb_cld(1)+Rn) *dev2)
        m_derf = derf( (rb_cld(1)+Rm) *dev2)
        do i = 1, nbin_cld
          n_aer(i) = -0.5 * No * n_derf !! this ith previously computed
          m_aer(i) = -0.5 * Mo * m_derf !! this ith previously computed
          n_derf = derf( (rb_cld(i+1)+Rn) *dev2)
          m_derf = derf( (rb_cld(i+1)+Rm) *dev2)
          n_aer(i) = n_aer(i) + 0.5 * No * n_derf
          m_aer(i) = m_aer(i) + 0.5 * Mo * m_derf
        enddo
        
!        sumcheck = 0
!        do i = 1, nbin_cld
!          sumcheck = sumcheck + n_aer(i)
!        enddo
!        sumcheck = abs(sumcheck/No - 1)
!        if ((sumcheck .gt. 1e-5).and. (1./Rn .gt. rmin_cld)) then
!          print*, "WARNING, No sumcheck PROBLEM"
!          print*, "sumcheck, No",sumcheck, No
!          print*, "min radius, Rn, ig, l", rmin_cld, 1./Rn, ig, l
!          print*, "Dust binned distribution", n_aer
!        endif
!        
!        sumcheck = 0
!        do i = 1, nbin_cld
!          sumcheck = sumcheck + m_aer(i)
!        enddo
!        sumcheck = abs(sumcheck/Mo - 1)
!        if ((sumcheck .gt. 1e-5) .and.  (1./Rn .gt. rmin_cld)) then
!          print*, "WARNING, Mo sumcheck PROBLEM"
!          print*, "sumcheck, Mo",sumcheck, Mo
!          print*, "min radius, Rm, ig, l", rmin_cld, 1./Rm, ig, l
!          print*, "Dust binned distribution", m_aer
!        endif

  
c       Get the rates of nucleation
        call nuclea(ph2o,zt(ig,l),satu,n_aer,rate)
        
        dN = 0.
        dM = 0.
        do i = 1, nbin_cld
          dN       = dN + n_aer(i)*(exp(-rate(i)*microtimestep)-1.)
          dM       = dM + m_aer(i)*(exp(-rate(i)*microtimestep)-1.)
        enddo


c       Update Dust particles
        zq(ig,l,igcm_dust_mass)   = 
     &  zq(ig,l,igcm_dust_mass)   + dM/ tauscaling(ig) !max(tauscaling(ig),1.e-10) 
        zq(ig,l,igcm_dust_number) = 
     &  zq(ig,l,igcm_dust_number) + dN/ tauscaling(ig) !max(tauscaling(ig),1.e-10)
c       Update CCNs
        zq(ig,l,igcm_ccn_mass)   = 
     &  zq(ig,l,igcm_ccn_mass)   - dM/ tauscaling(ig) !max(tauscaling(ig),1.e-10)
        zq(ig,l,igcm_ccn_number) = 
     &  zq(ig,l,igcm_ccn_number) - dN/ tauscaling(ig) !max(tauscaling(ig),1.e-10)
        
        ENDIF ! of is satu >1

!=============================================================
! 4. Ice growth: scheme for radius evolution
!=============================================================

c We trigger crystal growth if and only if there is at least one nuclei (N>1).
c Indeed, if we are supersaturated and still don't have at least one nuclei, we should better wait
c to avoid unrealistic value for nuclei radius and so on for cases that remain negligible.

       IF ( zq(ig,l,igcm_ccn_number)*tauscaling(ig).ge. 1.) THEN ! we trigger crystal growth 

     
        call updaterice_micro(zq(ig,l,igcm_h2o_ice),
     &    zq(ig,l,igcm_ccn_mass),zq(ig,l,igcm_ccn_number),
     &    tauscaling(ig),rice(ig,l),rhocloud(ig,l))

        No   = zq(ig,l,igcm_ccn_number)* tauscaling(ig) + 1.e-30

c       saturation at equilibrium
c       rice should not be too small, otherwise seq value is not valid
        seq  = exp(2.*sig(zt(ig,l))*mh2o / (rho_ice*rgp*zt(ig,l)*
     &             max(rice(ig,l),1.e-7)))
        
c       get resistance growth
        call growthrate(zt(ig,l),pplay(ig,l),
     &             real(ph2o/satu),rice(ig,l),res,Dv)

        res_out(ig,l) = res

ccccccc  implicit scheme of mass growth

        dMice =
     & (zq(ig,l,igcm_h2o_vap)-seq*zqsat(ig,l))
     &   /(res*zqsat(ig,l)/(cste*No*rice(ig,l)) + 1.)


! With the above scheme, dMice cannot be bigger than vapor, 
! but can be bigger than all available ice.
       dMice = max(dMice,-zq(ig,l,igcm_h2o_ice))
       dMice = min(dMice,zq(ig,l,igcm_h2o_vap)) ! this should be useless...

       zq(ig,l,igcm_h2o_ice) = zq(ig,l,igcm_h2o_ice)+dMice
       zq(ig,l,igcm_h2o_vap) = zq(ig,l,igcm_h2o_vap)-dMice


       countcells = countcells + 1 
       
       ! latent heat release
       lw=(2834.3-0.28*(zt(ig,l)-To)-
     &     0.004*(zt(ig,l)-To)*(zt(ig,l)-To))*1.e+3
       subpdtcloud(ig,l)= dMice*lw/cpp/microtimestep
          
c Special case of the isotope of water HDO    
       if (hdo) then
         !! condensation
         if (dMice.gt.0.0) then
         !! do we use fractionation?
           if (hdofrac) then
             !! Calculation of the HDO vapor coefficient
             Dv_hdo = 1./3. * sqrt( 8*kbz*zt(ig,l)/(pi*mhdo/nav) )
     &          * kbz * zt(ig,l) / 
     &          ( pi * pplay(ig,l) * (molco2+molhdo)*(molco2+molhdo) 
     &          * sqrt(1.+mhdo/mco2) )
             !! Calculation of the fractionnation coefficient at equilibrium
             alpha(ig,l) = exp(16288./zt(ig,l)**2.-9.34d-2)
c             alpha = exp(13525./zt(ig,l)**2.-5.59d-2)  !Lamb
             !! Calculation of the 'real' fractionnation coefficient
             alpha_c(ig,l) = (alpha(ig,l)*satu)/
     &          ( (alpha(ig,l)*(Dv/Dv_hdo)*(satu-1.)) + 1.)
c             alpha_c(ig,l) = alpha(ig,l) ! to test without the effect of cinetics
           else
             alpha_c(ig,l) = 1.d0
           endif
           if (zq0(ig,l,igcm_h2o_vap).gt.qperemin) then
              dMice_hdo= 
     &          dMice*alpha_c(ig,l)*
     &     ( zq0(ig,l,igcm_hdo_vap)
     &             /zq0(ig,l,igcm_h2o_vap) )
           else
             dMice_hdo=0.
           endif
         !! sublimation
         else
           if (zq0(ig,l,igcm_h2o_ice).gt.qperemin) then
             dMice_hdo= 
     &            dMice*
     &     ( zq0(ig,l,igcm_hdo_ice)
     &             /zq0(ig,l,igcm_h2o_ice) )
           else
             dMice_hdo=0.
           endif
         endif !if (dMice.gt.0.0)

       dMice_hdo = max(dMice_hdo,-zq(ig,l,igcm_hdo_ice))
       dMice_hdo = min(dMice_hdo,zq(ig,l,igcm_hdo_vap))

       zq(ig,l,igcm_hdo_ice) = zq(ig,l,igcm_hdo_ice)+dMice_hdo
       zq(ig,l,igcm_hdo_vap) = zq(ig,l,igcm_hdo_vap)-dMice_hdo

       endif ! if (hdo)        
     
!=============================================================
! 5. Dust cores released, tendancies, latent heat, etc ...
!=============================================================

c         If all the ice particles sublimate, all the condensation
c         nuclei are released:
          if (zq(ig,l,igcm_h2o_ice).le.1.e-28) then
          
c           Water
            zq(ig,l,igcm_h2o_vap) = zq(ig,l,igcm_h2o_vap) 
     &                            + zq(ig,l,igcm_h2o_ice)
            zq(ig,l,igcm_h2o_ice) = 0.
            if (hdo) then
              zq(ig,l,igcm_hdo_vap) = zq(ig,l,igcm_hdo_vap) 
     &                            + zq(ig,l,igcm_hdo_ice)
              zq(ig,l,igcm_hdo_ice) = 0.
            endif
c           Dust particles
            zq(ig,l,igcm_dust_mass) = zq(ig,l,igcm_dust_mass)
     &                              + zq(ig,l,igcm_ccn_mass)
            zq(ig,l,igcm_dust_number) = zq(ig,l,igcm_dust_number)
     &                                + zq(ig,l,igcm_ccn_number)
c           CCNs
            zq(ig,l,igcm_ccn_mass) = 0.
            zq(ig,l,igcm_ccn_number) = 0.

          endif
         
          ENDIF !of if Nccn>1
         
        ENDDO ! of ig loop
      ENDDO ! of nlayer loop
      
      
      ! Get cloud tendencies
        subpdqcloud(1:ngrid,1:nlay,igcm_h2o_vap) =
     &   (zq(1:ngrid,1:nlay,igcm_h2o_vap) - 
     &    zq0(1:ngrid,1:nlay,igcm_h2o_vap))/microtimestep
        subpdqcloud(1:ngrid,1:nlay,igcm_h2o_ice) =
     &   (zq(1:ngrid,1:nlay,igcm_h2o_ice) -
     &    zq0(1:ngrid,1:nlay,igcm_h2o_ice))/microtimestep
        if (hdo) then
          subpdqcloud(1:ngrid,1:nlay,igcm_hdo_vap) =
     &     (zq(1:ngrid,1:nlay,igcm_hdo_vap) - 
     &      zq0(1:ngrid,1:nlay,igcm_hdo_vap))/microtimestep
          subpdqcloud(1:ngrid,1:nlay,igcm_hdo_ice) =
     &     (zq(1:ngrid,1:nlay,igcm_hdo_ice) -
     &      zq0(1:ngrid,1:nlay,igcm_hdo_ice))/microtimestep
        endif
        subpdqcloud(1:ngrid,1:nlay,igcm_ccn_mass) =
     &   (zq(1:ngrid,1:nlay,igcm_ccn_mass) -
     &    zq0(1:ngrid,1:nlay,igcm_ccn_mass))/microtimestep
        subpdqcloud(1:ngrid,1:nlay,igcm_ccn_number) =
     &   (zq(1:ngrid,1:nlay,igcm_ccn_number) -
     &    zq0(1:ngrid,1:nlay,igcm_ccn_number))/microtimestep
     
      if (scavenging) then
      
        subpdqcloud(1:ngrid,1:nlay,igcm_dust_mass) =
     &   (zq(1:ngrid,1:nlay,igcm_dust_mass) -
     &    zq0(1:ngrid,1:nlay,igcm_dust_mass))/microtimestep
        subpdqcloud(1:ngrid,1:nlay,igcm_dust_number) =
     &   (zq(1:ngrid,1:nlay,igcm_dust_number) -
     &    zq0(1:ngrid,1:nlay,igcm_dust_number))/microtimestep
          
      endif
     
     
     
!!!!!!!!!!!!!! TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS 
!!!!!!!!!!!!!! TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS 
!!!!!!!!!!!!!! TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS
      IF (test_flag) then
      
!       error2d(:) = 0.
       DO l=1,nlay
       DO ig=1,ngrid
!         error2d(ig) = max(abs(error_out(ig,l)),error2d(ig))
          satubf(ig,l) = zq0(ig,l,igcm_h2o_vap)/zqsat(ig,l) 
          satuaf(ig,l) = zq(ig,l,igcm_h2o_vap)/zqsat(ig,l) 
       ENDDO
       ENDDO

      print*, 'count is ',countcells, ' i.e. ',
     &     countcells*100/(nlay*ngrid), '% for microphys computation'

!      IF (ngrid.ne.1) THEN ! 3D
!         call WRITEDIAGFI(ngrid,"satu","ratio saturation","",3,
!     &                    satu_out)
!         call WRITEDIAGFI(ngrid,"dM","ccn variation","kg/kg",3,
!     &                    dM_out)
!         call WRITEDIAGFI(ngrid,"dN","ccn variation","#",3,
!     &                    dN_out)
!         call WRITEDIAGFI(ngrid,"error","dichotomy max error","%",2,
!     &                    error2d)
!         call WRITEDIAGFI(ngrid,"zqsat","zqsat","kg",3,
!     &                    zqsat)
!      ENDIF

!      IF (ngrid.eq.1) THEN ! 1D
!         call WRITEDIAGFI(ngrid,"error","incertitude sur glace","%",1,
!     &                    error_out)
         call WRITEdiagfi(ngrid,"resist","resistance","s/m2",1,
     &                    res_out)
         call WRITEdiagfi(ngrid,"satu_bf","satu before","kg/kg",1,
     &                    satubf)
         call WRITEdiagfi(ngrid,"satu_af","satu after","kg/kg",1,
     &                    satuaf)
         call WRITEdiagfi(ngrid,"vapbf","h2ovap before","kg/kg",1,
     &                    zq0(1,1,igcm_h2o_vap))
         call WRITEdiagfi(ngrid,"vapaf","h2ovap after","kg/kg",1,
     &                    zq(1,1,igcm_h2o_vap))
         call WRITEdiagfi(ngrid,"icebf","h2oice before","kg/kg",1,
     &                    zq0(1,1,igcm_h2o_ice))
         call WRITEdiagfi(ngrid,"iceaf","h2oice after","kg/kg",1,
     &                    zq(1,1,igcm_h2o_ice))
         call WRITEdiagfi(ngrid,"ccnbf","ccn before","/kg",1,
     &                    zq0(1,1,igcm_ccn_number))
         call WRITEdiagfi(ngrid,"ccnaf","ccn after","/kg",1,
     &                    zq(1,1,igcm_ccn_number))
c         call WRITEDIAGFI(ngrid,"growthrate","growth rate","m^2/s",1,
c     &                    gr_out)
c         call WRITEDIAGFI(ngrid,"nuclearate","nucleation rate","",1,
c     &                    rate_out)
c         call WRITEDIAGFI(ngrid,"dM","ccn variation","kg",1,
c     &                    dM_out)
c         call WRITEDIAGFI(ngrid,"dN","ccn variation","#",1,
c     &                    dN_out)
         call WRITEdiagfi(ngrid,"zqsat","p vap sat","kg/kg",1,
     &                    zqsat)
!         call WRITEDIAGFI(ngrid,"satu","ratio saturation","",1,
!     &                    satu_out)
         call WRITEdiagfi(ngrid,"rice","ice radius","m",1,
     &                    rice)
!         call WRITEDIAGFI(ngrid,"rdust_sca","rdust","m",1,
!     &                    rdust)
!         call WRITEDIAGFI(ngrid,"rsedcloud","rsedcloud","m",1,
!     &                    rsedcloud)
!         call WRITEDIAGFI(ngrid,"rhocloud","rhocloud","kg.m-3",1,
!     &                    rhocloud)
!      ENDIF
      
      ENDIF ! endif test_flag
!!!!!!!!!!!!!! TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS 
!!!!!!!!!!!!!! TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS 
!!!!!!!!!!!!!! TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS TESTS OUTPUTS 

      return

      
      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c The so -called "phi" function is such as phi(r) - phi(r0) = t - t0
c It is an analytical solution to the ice radius growth equation, 
c with the approximation of a constant 'reduced' cunningham correction factor 
c (lambda in growthrate.F) taken at radius req instead of rice    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      subroutine phi(rice,req,coeff1,coeff2,time)
c      
c      implicit none
c      
c      ! inputs
c      real rice ! ice radius
c      real req  ! ice radius at equilibirum
c      real coeff1  ! coeff for the log
c      real coeff2  ! coeff for the arctan
c
c      ! output      
c      real time
c      
c      !local
c      real var
c      
c      ! 1.73205 is sqrt(3)
c      
c      var = max(
c     &  abs(rice-req) / sqrt(rice*rice + rice*req  + req*req),1e-30)
c            
c       time = 
c     &   coeff1 * 
c     &   log( var )
c     & + coeff2 * 1.73205 *
c     &   atan( (2*rice+req) / (1.73205*req) )
c      
c      return
c      end
      
      
      
      END SUBROUTINE improvedclouds
      
      END MODULE improvedclouds_mod
