










      subroutine growthrate(temp,pmid,psat,rcrystal,res,Dv)

      use tracer_mod, only: rho_ice
      USE comcstfi_h
      IMPLICIT NONE

c=======================================================================
c
c     Determination of the water ice crystal growth rate
c
c     Authors: F. Montmessin
c       Adapted for the LMD/GCM by J.-B. Madeleine (October 2011)
c       Use of resistances in the analytical function 
c            instead of growth rate - T. Navarro (2012)
c     
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

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

c
c   arguments:
c   ----------

c     Input
      REAL temp     ! temperature in the middle of the layer (K)
      REAL pmid     ! pressure in the middle of the layer (K)
      REAL psat   ! water vapor saturation pressure (Pa) 
      REAL rcrystal ! crystal radius before condensation (m)

c     Output
      REAL res      ! growth resistance (res=Rk+Rd)
      REAL Dv       ! water vapor diffusion coefficient

c   local:
c   ------

      REAL k,Lv                 
      REAL knudsen           ! Knudsen number (gas mean free path/particle radius)
      REAL afactor,lambda       ! Intermediate computations for growth rate
      REAL Rk,Rd
      
      

c-----------------------------------------------------------------------
c      Ice particle growth rate by diffusion/impegement of water molecules
c                r.dr/dt = (S-Seq) / (Seq*Rk+Rd)
c        with r the crystal radius, Rk and Rd the resistances due to 
c        latent heat release and to vapor diffusion respectively 
c----------------------------------------------------------------------- 

c     - Equilibrium saturation accounting for KeLvin Effect
c      seq=exp(2*sigh2o*mh2o/(rho_ice*rgp*t*r))
c      (already computed in improvedcloud.F)

c     - Thermal conductibility of CO2
      k  = (0.17913 * temp - 13.9789) * 4.184e-4
c     - Latent heat of h2o (J.kg-1)
      Lv = (2834.3 
     &        - 0.28  * (temp-To) 
     &        - 0.004 * (temp-To) * (temp-To) ) * 1.e+3

c     - Constant to compute gas mean free path
c     l= (T/P)*a, with a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))
      afactor = 0.707*rgp/(4 * pi * molco2 * molco2 * nav)

c     - Compute Dv, water vapor diffusion coefficient
c       accounting for both kinetic and continuum regime of diffusion,
c       the nature of which depending on the Knudsen number.

      Dv = 1./3. * sqrt( 8*kbz*temp/(pi*mh2o/nav) )* kbz * temp / 
     &   ( pi * pmid * (molco2+molh2o)*(molco2+molh2o) 
     &        * sqrt(1.+mh2o/mco2) )
      
      knudsen = temp / pmid * afactor / rcrystal
      lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen)
      
c      Dv is not corrected. Instead, we use below coefficients coeff1, coeff2
c      Dv      = Dv / (1. + lambda * knudsen)

c     - Compute Rk
      Rk = Lv*Lv* rho_ice * mh2o / (k*rgp*temp*temp)
c     - Compute Rd
      Rd = rgp * temp *rho_ice / (Dv*psat*mh2o)
      
      
      res = Rk + Rd*(1. + lambda * knudsen)
      
      !coeff1 = real(Rk + Rd*(1. + lambda * knudsen))
      !coeff2 = real(Rk + Rd*(1. - lambda * knudsen))
      
c Below are growth rate used for other schemes
c     - Compute growth=rdr/dt, then r(t+1)= sqrt(r(t)**2.+2.*growth*dt)
c      growth = 1. / (seq*Rk+Rd)
c      growth = (ph2o/psat-seq) / (seq*Rk+Rd)
c      rf   = sqrt( max( r**2.+2.*growth*timestep , 0. ) )
c      dr   = rf-r

      RETURN
      END

