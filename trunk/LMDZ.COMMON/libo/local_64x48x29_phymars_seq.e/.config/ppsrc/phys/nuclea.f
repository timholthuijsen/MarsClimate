










*******************************************************
*                                                     *
      subroutine nuclea(ph2o,temp,sat,n_ccn,nucrate)
      USE comcstfi_h
      implicit none
*                                                     *
*   This subroutine computes the nucleation rate      *
*   as given in Pruppacher & Klett (1978) in the      *
*   case of water ice forming on a solid substrate.   *
*     Definition refined by Keese (jgr,1989)          *
*   Authors: F. Montmessin                            *
*     Adapted for the LMD/GCM by J.-B. Madeleine      *
*     (October 2011)                                  *
*     Optimisation by A. Spiga (February 2012)        *  
*******************************************************

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
      include "callkeys.h"

c     Inputs
      DOUBLE PRECISION ph2o,sat
      DOUBLE PRECISION n_ccn(nbin_cld)
      REAL temp

c     Output
   !   DOUBLE PRECISION nucrate(nbin_cld)
      REAL nucrate(nbin_cld)

c     Local variables
      DOUBLE PRECISION nh2o
      DOUBLE PRECISION sig      ! Water-ice/air surface tension  (N.m)
      external sig
      DOUBLE PRECISION rstar    ! Radius of the critical germ (m)
      DOUBLE PRECISION gstar    ! # of molecules forming a critical embryo
      DOUBLE PRECISION fistar   ! Activation energy required to form a critical embryo (J)
!      DOUBLE PRECISION zeldov   ! Zeldovitch factor (no dim)
      DOUBLE PRECISION fshape   ! function defined at the end of the file
      DOUBLE PRECISION deltaf

c     Ratio rstar/radius of the nucleating dust particle
c     double precision xratio
      
      double precision mtetalocal ! local mteta in double precision

      double precision fshapesimple,zefshape


      integer i
      
      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall

c     *************************************************

      mtetalocal = mteta  !! use mtetalocal for better performance

      IF (temp_dependant_m) THEN
         mtetalocal = min(0.0044*temp + 0.1831,0.97)
      ENDIF ! (temp_dependant_m) THEN
cccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc ESSAIS TN MTETA = F (T) cccccccccccccc
c      if (temp .gt. 200) then
c         mtetalocal = mtetalocal
c      else if (temp .lt. 190) then
c         mtetalocal = mtetalocal-0.05
c      else
c         mtetalocal = mtetalocal - (190-temp)*0.005
c      endif
c----------------exp law, see Trainer 2008, J. Phys. Chem. C 2009, 113, 2036\u20132040
       !mtetalocal = max(mtetalocal - 6005*exp(-0.065*temp),0.1)
       !mtetalocal = max(mtetalocal - 6005*exp(-0.068*temp),0.1)
               !print*, mtetalocal, temp
cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc 
      IF (firstcall.and.temp_dependant_m) THEN
          print*, ' '  
          print*, 'dear user, please keep in mind that'
          print*, 'contact parameter IS NOT constant ;'
          print*, 'Using the following linear fit from'
          print*, 'Maattanen et al. 2014 (SM linear fit) :'
          print*, 'min(0.0044*temp + 0.1831,0.97)'
          print*, ' '  
         firstcall=.false.
      ELSE IF (firstcall.and.(.not.(temp_dependant_m))) THEN
          print*, ' '  
          print*, 'dear user, please keep in mind that'
          print*, 'contact parameter IS constant'
          print*, ' '  
         firstcall=.false.
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccc
   

      if (sat .gt. 1.) then    ! minimum condition to activate nucleation

        nh2o   = ph2o / kbz / temp
        rstar  = 2. * sig(temp) * vo1 / (rgp*temp*log(sat))
        gstar  = 4. * nav * pi * (rstar * rstar * rstar) / (3.*vo1)
        
        fshapesimple = (2.+mtetalocal)*(1.-mtetalocal)*(1.-mtetalocal)
     &                   / 4.

c       Loop over size bins
        do 200 i=1,nbin_cld

          if ( n_ccn(i) .lt. 1e-10 ) then
c           no dust, no need to compute nucleation!
            nucrate(i)=0.
            goto 200
          endif

          if (rad_cld(i).gt.3000.*rstar) then
            zefshape = fshapesimple
          else
            zefshape = fshape(mtetalocal,rad_cld(i)/rstar)
          endif

          fistar = (4./3.*pi) * sig(temp) * (rstar * rstar) * 
     &             zefshape
          deltaf = (2.*desorp-surfdif-fistar)/
     &             (kbz*temp)
          deltaf = min( max(deltaf, -100.d0), 100.d0)

          if (deltaf.eq.-100.) then
            nucrate(i) = 0.
          else
            nucrate(i)= real(sqrt ( fistar /
     &               (3.*pi*kbz*temp*(gstar*gstar)) )
     &                  * kbz * temp * rstar
     &                  * rstar * 4. * pi
     &                  * ( nh2o*rad_cld(i) )
     &                  * ( nh2o*rad_cld(i) )
     &                  / ( zefshape * nus * m0 )
     &                  * exp (deltaf))
          endif

200     continue

      else

        do i=1,nbin_cld
          nucrate(i) = 0.
        enddo

      endif

      return
      end

*********************************************************
      double precision function fshape(cost,rap)
      implicit none
*        function computing the f(m,x) factor           *
* related to energy required to form a critical embryo  *
*********************************************************

      double precision cost,rap
      double precision yeah

          !! PHI
          yeah = sqrt( 1. - 2.*cost*rap + rap*rap )
          !! FSHAPE = TERM A
          fshape = (1.-cost*rap) / yeah
          fshape = fshape * fshape * fshape
          fshape = 1. + fshape
          !! ... + TERM B
          yeah = (rap-cost)/yeah
          fshape = fshape + 
     & rap*rap*rap*(2.-3.*yeah+yeah*yeah*yeah)
          !! ... + TERM C 
          fshape = fshape + 3. * cost * rap * rap * (yeah-1.)
          !! FACTOR 1/2
          fshape = 0.5*fshape

      return 
      end
