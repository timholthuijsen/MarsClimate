










c**************************************************************************
c
      subroutine nltecool(ngrid,nlayer,nq,pplay,pt,pq,dtnlte)
c
c This code was designed as a delivery for the "Martian Environment Models" 
c project ( ESA contract 11369/95/nl/jg CCN2 )
c Computes non-LTE heating rates from CO2 emission at 15 um
c in the Martian upper atmosphere.
c Uses a simplified model consisting of two excited levels with two
c emission bands, one of them stronger than the other, which correspond
c to the behaviours of the 626 fundamental band and the isotopic fund.bands.
c It uses a cool-to-space approximation with tabulated escape functions.
c These escape functions have been precomputed for the strong and weak bands,
c and are given as a function of pressure in separate files.
c The output values are the heating rates (actually, cooling, since they
c are always negative) for the two bands, i.e., the total cooling is the
c sum of them.
c Miguel A. Lopez Valverde
c Instituto de Astrofisica de Andalucia (CSIC), Granada, Spain
c
c Version 1b.  See description above.  22-March-2000.
c Adapted as a subroutine for use in GCM -- PLR/SRL 6/2000
c Version 1c.  Inclusion of VMR in the tabulation of escape functions. 
c              Table contains now only 1 input file -- Miguel 11/Jul/2000
c Version 1d  data contained in original input file "nlte_escape.dat"
c now stored in include file "nltedata.h" Y.Wanherdrick 09/2000

c       jul 2011 fgg   Modified to allow variable O
c     
c***************************************************************************

      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_n2, mmol
      use conc_mod, only: mmean
      implicit none

c INCLUDE nlte_escape.h
       !
       ! Escape functions and VMRs from tabulated values.
       ! Origin: nlte_escape.dat file form Miguel L.
       ! (Y. Wanherdrick, 09/2000)

       integer     np
       parameter ( np = 68 )                ! # data points in tabulation
       real        pnb(np)                  ! Pressure in tabulation
       real        ef1(np)                  ! Esc.funct.#1, tabulated
       real        ef2(np)                  ! Esc.funct.#2, tabulated
       real        co2vmr(np)               ! CO2 VMR tabulated
       real        o3pvmr(np)               ! CO2 VMR tabulated
       real        n2covmr(np)              ! N2+CO VMR tabulated
 
 
      DATA pnb/
     &    12.0000,    11.0000,    10.8000,
     &    10.6000,   10.40000,   10.20000,
     &   10.00000,    9.80000,    9.60000,
     &    9.40000,    9.20000,    9.00000,
     &    8.80000,    8.60000,    8.40000,
     &    8.20000,    8.00000,    7.80000,
     &    7.60000,    7.40000,    7.20000,
     &    7.00000,    6.80000,    6.60000,
     &    6.40000,    6.20000,    6.00000,
     &    5.80000,    5.60000,    5.40000,
     &    5.20000,    5.00000,    4.80000,
     &    4.60000,    4.40000,    4.20000,
     &    4.00000,    3.80000,    3.60000,
     &    3.40000,    3.20000,    3.00000,
     &    2.80000,    2.60000,    2.40000,
     &    2.20000,    2.00000,    1.80000,
     &    1.60000,    1.40000,    1.20000,
     &    1.00000,   0.800000,   0.599999,
     &   0.400000,   0.200000,  0.,
     &  -0.200000,  -0.400001,  -0.600000,
     &  -0.800000,   -1.00000,   -1.20000,
     &   -1.40000,   -1.60000,   -1.80000,
     &   -2.00000,   -3.00000/
 
      DATA ef1/
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.58112E-04,    4.58112E-04,
     &    4.58112E-04,    4.61707E-04,    4.76886E-04,
     &    4.95638E-04,    5.20935E-04,    5.55511E-04,
     &    6.01219E-04,    6.63734E-04,    7.50691E-04,
     &    8.63474E-04,    1.00900E-03,    1.19642E-03,
     &    1.42690E-03,    1.71398E-03,    2.06663E-03,
     &    2.48974E-03,    3.01578E-03,    3.64350E-03,
     &    4.40323E-03,    5.32066E-03,    6.40456E-03,
     &    7.72069E-03,    9.25684E-03,    1.10905E-02,
     &    1.32374E-02,    1.57643E-02,    1.87388E-02,
     &    2.22072E-02,    2.63099E-02,    3.10614E-02,
     &    3.66948E-02,    4.32373E-02,    5.15022E-02,
     &    6.21455E-02,    7.77212E-02,    9.92027E-02,
     &   0.131155,   0.179470,   0.258913,
     &   0.380549,   0.530450,   0.643180,
     &   0.741061,   0.826336,   0.922787,
     &   0.997203,    1.00000/
 
      DATA ef2/
     &    1.98992E-03,    1.98992E-03,    1.98992E-03,
     &    1.98992E-03,    1.98992E-03,    1.98992E-03,
     &    1.98992E-03,    1.98992E-03,    1.98992E-03,
     &    1.98992E-03,    1.98992E-03,    2.01376E-03,
     &    2.09450E-03,    2.22993E-03,    2.42056E-03,
     &    2.68018E-03,    3.04398E-03,    3.43896E-03,
     &    3.80282E-03,    4.20622E-03,    4.76121E-03,
     &    8.01698E-03,    1.19947E-02,    1.69149E-02,
     &    2.24497E-02,    2.85244E-02,    3.54813E-02,
     &    4.39264E-02,    5.46248E-02,    6.75367E-02,
     &    8.29931E-02,    1.01717E-01,   0.123422,
     &   0.148468,   0.177096,   0.208816,
     &   0.244003,   0.282013,   0.322559,
     &   0.365542,   0.410518,   0.457384,
     &   0.505358,   0.553627,   0.600472,
     &   0.644807,   0.687185,   0.727429,
     &   0.764734,   0.798562,   0.828699,
     &   0.854797,   0.877717,   0.897874,
     &   0.915258,   0.929904,   0.942381,
     &   0.952906,   0.962173,   0.970191,
     &   0.976437,   0.981501,   0.985406,
     &   0.988560,   0.991111,   0.993653,
     &   0.995561,    1.00000/
 
      DATA co2vmr/
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.950000,   0.950000,   0.950000,
     &   0.949619,   0.947694,   0.945830,
     &   0.944016,   0.940557,   0.937068,
     &   0.932366,   0.893661/
 
      DATA o3pvmr/
     &    5.06756E-08,    9.16539E-07,    1.68217E-06,
     &    3.00843E-06,    5.03151E-06,    8.07489E-06,
     &    1.23137E-05,    1.79029E-05,    2.45308E-05,
     &    3.27431E-05,    4.26692E-05,    5.44396E-05,
     &    6.78865E-05,    8.33147E-05,    1.00148E-04,
     &    1.18846E-04,    1.39681E-04,    1.64909E-04,
     &    1.93617E-04,    2.25161E-04,    2.60834E-04,
     &    3.01501E-04,    3.44953E-04,    3.91011E-04,
     &    4.40377E-04,    4.90820E-04,    5.43200E-04,
     &    5.95335E-04,    6.45420E-04,    6.93166E-04,
     &    7.43729E-04,    7.93710E-04,    8.44394E-04,
     &    8.94318E-04,    9.44732E-04,    9.94964E-04,
     &    1.04901E-03,    1.10008E-03,    1.16302E-03,
     &    1.22989E-03,    1.30026E-03,    1.37131E-03,
     &    1.45556E-03,    1.55186E-03,    1.66328E-03,
     &    1.77802E-03,    1.91546E-03,    2.07503E-03,
     &    2.24903E-03,    2.47117E-03,    2.71728E-03,
     &    2.99739E-03,    3.33582E-03,    3.73507E-03,
     &    4.20819E-03,    4.76887E-03,    5.42558E-03,
     &    6.20815E-03,    7.14473E-03,    8.28545E-03,
     &    9.51779E-03,    1.08140E-02,    1.22359E-02,
     &    1.36870E-02,    1.51495E-02,    1.67196E-02,
     &    1.85485E-02,    3.36252E-02/
 
      DATA n2covmr/
     &    2.71412E-02,    2.71464E-02,    2.71490E-02,
     &    2.71523E-02,    2.71558E-02,    2.71617E-02,
     &    2.71672E-02,    2.71749E-02,    2.71837E-02,
     &    2.71943E-02,    2.72058E-02,    2.72189E-02,
     &    2.72326E-02,    2.72483E-02,    2.72661E-02,
     &    2.72848E-02,    2.73054E-02,    2.73279E-02,
     &    2.73514E-02,    2.73775E-02,    2.74048E-02,
     &    2.74345E-02,    2.74672E-02,    2.75021E-02,
     &    2.75404E-02,    2.75826E-02,    2.76340E-02,
     &    2.77013E-02,    2.78220E-02,    2.79707E-02,
     &    2.81759E-02,    2.84339E-02,    2.87587E-02,
     &    2.91600E-02,    2.96561E-02,    3.02558E-02,
     &    3.09922E-02,    3.18062E-02,    3.27010E-02,
     &    3.35635E-02,    3.44388E-02,    3.53327E-02,
     &    3.62143E-02,    3.70941E-02,    3.79315E-02,
     &    3.87626E-02,    3.95524E-02,    4.03071E-02,
     &    4.10071E-02,    4.16229E-02,    4.21231E-02,
     &    4.25167E-02,    4.27964E-02,    4.29773E-02,
     &    4.30488E-02,    4.29638E-02,    4.28049E-02,
     &    4.26788E-02,    4.26822E-02,    4.29426E-02,
     &    4.34634E-02,    4.42559E-02,    4.53038E-02,
     &    4.65879E-02,    4.80262E-02,    4.96303E-02,
     &    5.14885E-02,    6.91651E-02/
 
!--------------------------------------------
!     data for photochemistry
!--------------------------------------------

!--------------------------------------------
!     dimensions of photolysis lookup table
!--------------------------------------------

      integer, parameter :: nd    = 13  ! species
      integer, parameter :: nz    = 143 ! altitude
      integer, parameter :: nozo  = 7   ! ozone
      integer, parameter :: nsza  = 27  ! solar zenith angle
      integer, parameter :: ntemp = 4   ! temperature
      integer, parameter :: ntau  = 8   ! dust

!--------------------------------------------

      common/chimiedata/jphot,colairtab,table_ozo

      real jphot(ntemp,nsza,nz,nozo,ntau,nd)
      real colairtab(nz)
      real szatab(nsza)
      real table_ozo(nozo)
      real tautab(ntau)

      data szatab/0.,  5., 10., 15., 20., 25.,                          &
     &            30., 35., 40., 45., 50., 55.,                         &
     &            60., 65., 70., 75., 80., 82.,                         &
     &            84., 86., 88., 90., 91., 92.,                         &
     &            93., 94., 95./

      data tautab/0., 0.2, 0.4, 0.6, 0.8, 1., 2., 4./
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

c Input and output variables
c
      integer     ngrid                           ! no. of horiz. gridpoints
      integer     nlayer                          ! no. of atmospheric layers
      integer     nq                              ! no. of tracers
      real        pplay(ngrid,nlayer)             ! input pressure grid
      real        pt(ngrid,nlayer)                ! input temperatures
      real        pq(ngrid,nlayer,nq)                ! input mmrs
      real        dtnlte(ngrid,nlayer)            ! output temp. tendencies

c
c Standard atmosphere variables
c
      real        nt                              ! number density [cm-3]
      real        co2(nlayer)                     !  "   of CO2
      real        o3p(nlayer)                     !  "   of atomic oxygen
      real        n2co(nlayer)                    !  "  of N2 + CO
      real        pyy(nlayer)                     ! auxiliary pressure grid

c
c Vectors and indexes for the tabulation of escape functions and VMR
c
c                 np                              ! # data points in tabulation
c                 pnb(np)                         ! Pressure in tabulation
c                 ef1(np)                         ! Esc.funct.#1, tabulated
c                 ef2(np)                         ! Esc.funct.#2, tabulated
c                 co2vmr(np)                      ! CO2 VMR tabulated
c                 o3pvmr(np)                      ! CO2 VMR tabulated
c                 n2covmr(np)                     ! N2+CO VMR tabulated
      real        escf1(nlayer)                   ! Esc.funct.#1, interpolated
      real        escf2(nlayer)                   ! Esc.funct.#2, interpolated


c
c Local Constants
c
      real       nu1, nu2                         ! freq. of energy levels
      real       imr1, imr2                       ! isotopic abundances
      real       hplanck, gamma, vlight           ! physical constants
      real       ee
      real       rfvt                             ! collisional rate
      real       rfvto3p                          !     "
      real       rfvv                             !     "

c
c Local variables for the main loop
c
      real       n1, n2, co2t                     ! ground populations
      real       l1, p1, p12                      ! prod & losses
      real       l2, p2, p21
      real       tt                               ! dummy variable
      real       c1, c2                           ! molecular constants
      real       ae1, ae2                         ! einstein spontaneous emission
      real       a1, a2, a12, a21
      real       pl1, pl2
      real       el1, el2
      real       hr1, hr2                         ! heating rate due to each band
      real       hr(nlayer)                       ! total heating rate

c
c Indexes
c
      integer    i
      integer    j,ii

c
c Rate coefficients
c
      real       k19xca, k19xcb
      real       k19cap1, k19cap2
      real       k19cbp1, k19cbp2
      real       d19c, d19cp1, d19cp2
      real       k20xc, k20cp1, k20cp2
      real       k21xc, k21cp2

      logical    firstcall
      data       firstcall/.true./
      save       firstcall,ef1,ef2,co2vmr,n2covmr,o3pvmr,pnb

c
c Data
c
      data       nu1, nu2, hplanck, gamma, vlight, ee/
     1     667.38, 662.3734, 6.6261e-27, 1.191e-5, 3.e10, 1.438769/

c*************************************************************************
c       PROGRAM  STARTS
c*************************************************************************

      imr1 = 0.987
      imr2 = 0.00408 + 0.0112
      rfvt = 0.1
      rfvto3p = 1.0
      rfvv = 0.1

      if(firstcall) then

         do i=1,np
            pnb(i)=1.0e-4*exp(pnb(i)) ! p into Pa
         end do

         firstcall = .false.

      endif
c
c MAIN LOOP, for each gridpoint and altitude:
c
      do j=1,ngrid  ! loop over grid points
c
c set up local pressure grid
c
         do ii=1,nlayer
            pyy(ii)=pplay(j,ii)
         enddo
!
! Interpolate escape functions and VMR to the desired grid
!
         call interp1(escf2,pyy,nlayer,ef2,pnb,np)
         call interp1(escf1,pyy,nlayer,ef1,pnb,np)
         if(nltemodel.eq.0) then
            call interp3(co2,o3p,n2co,pyy,nlayer,
     &           co2vmr,o3pvmr,n2covmr,pnb,np)
         endif
        
         do i=1,nlayer  ! loop over layers
C
C test if p lies outside range (p > 3.5 Pa)
C changed to 1 Pa since transition will always be higher than this
C
            if(pyy(i) .gt. 1.0 .or. pyy(i) .lt. 4.0e-6) then 
               hr(i)=0.0
               dtnlte(j,i)=0.0
            else
c
c           if(pt(j,i).lt.1.0)print*,pt(j,i)
               nt = pyy(i)/(1.381e-17*pt(j,i)) ! nt in cm-3
               !Dynamic composition
               if(nltemodel.eq.1) then
                  co2(i)=pq(j,i,igcm_co2)*mmean(j,i)/mmol(igcm_co2)
                  o3p(i)=pq(j,i,igcm_o)*mmean(j,i)/mmol(igcm_o)
                  n2co(i)=pq(j,i,igcm_co)*mmean(j,i)/mmol(igcm_co) +
     $                 pq(j,i,igcm_n2)*mmean(j,i)/mmol(igcm_n2)
               endif

               !Mixing ratio to density
               co2(i)=co2(i)*nt                ! CO2 density in cm-3
               o3p(i)=o3p(i)*nt                ! O3p density in cm-3
               n2co(i)=n2co(i)*nt              ! N2+CO in cm-3
c molecular populations
               n1 = co2(i) * imr1
               n2 = co2(i) * imr2
               co2t = n1 + n2

c intermediate collisional rates
               tt = pt(j,i)*pt(j,i)

               if (pt(j,i).le.175.0) then
                  k19xca = 3.3e-15
                  k19xcb = 7.6e-16
               else
                  k19xca = 4.2e-12 * exp( -2988.0/pt(j,i) + 303930.0/tt)
                  k19xcb = 2.1e-12 * exp( -2659.0/pt(j,i) + 223052.0/tt)
               endif
               k19xca = k19xca * rfvt
               k19xcb = k19xcb * rfvt
               k19cap1 = k19xca * 2.0 * exp( -ee*nu1/pt(j,i) )
               k19cap2 = k19xca * 2.0 * exp( -ee*nu2/pt(j,i) )
               k19cbp1 = k19xcb * 2.0 * exp( -ee*nu1/pt(j,i) )
               k19cbp2 = k19xcb * 2.0 * exp( -ee*nu2/pt(j,i) )
               d19c = k19xca*co2t + k19xcb*n2co(i)
               d19cp1 = k19cap1*co2t + k19cbp1*n2co(i)
               d19cp2 = k19cap2*co2t + k19cbp2*n2co(i)
                                !
               k20xc = 3.e-12 * rfvto3p
               k20cp1 = k20xc * 2.0 * exp( -ee/pt(j,i) * nu1 )
               k20cp2 = k20xc * 2.0 * exp( -ee/pt(j,i) * nu2 )
                                !
               k21xc = 2.49e-11 * 0.5 * rfvv
               k21cp2 = k21xc * exp( - ee/pt(j,i) * (nu2-nu1) )
                                !
               l1 = d19c + k20xc*o3p(i) + k21cp2*n2
               p1 = ( d19cp1 + k20cp1*o3p(i) ) * n1
               p12 = k21xc*n1
                                !
               l2 = d19c + k20xc*o3p(i) + k21xc*n1
               p2 = ( d19cp2 + k20cp2*o3p(i) ) * n2
               p21 = k21cp2*n2

c radiative rates
               ae1 = 1.3546 * 1.66 / 4.0 * escf1(i)
               ae2 = ( 1.3452 + 1.1878 ) * 1.66 / 4.0 * escf2(i)
               l1 = l1 + ae1
               l2 = l2 + ae2

c solving the system
               c1 = gamma*nu1**3. * 0.5
               c2 = gamma*nu2**3. * 0.5
               a1 = c1 * p1 / (n1*l1)
               a2 = c2 * p2 / (n2*l2)
               a12 = (nu1/nu2)**3. * n2/n1 * p12/l1
               a21 = (nu2/nu1)**3. * n1/n2 * p21/l2
               el2 = (a2 + a21 * a1 ) / ( 1.0 - a21 * a12 )
               el1 = a1 + a12 * el2
               pl1 = el1 * n1 / c1
               pl2 = el2 * n2 / c2

c  heating rate
               hr1 = - hplanck*vlight * nu1 * ae1 * pl1
               hr2 = - hplanck*vlight * nu2 * ae2 * pl2
               hr(i) = hr1 + hr2
               dtnlte(j,i)=0.1*hr(i)*pt(j,i)/(4.4*pyy(i)) ! dtnlte in K s-1
c              write(7,25)pxx(i),hr1,hr2,hr(i),qt
c  25         format(' ',1p5e12.4)

            endif

         enddo  ! end loop over layers
      enddo     ! end loop over grid points
c     close(7)
c
        return
        end

c***********************************************************************

      subroutine interp1(escout,p,nlayer,escin,pin,nl)
C
C subroutine to perform linear interpolation in pressure from 1D profile 
C escin(nl) sampled on pressure grid pin(nl) to profile
C escout(nlayer) on pressure grid p(nlayer).
C
      real escout(nlayer),p(nlayer)
      real escin(nl),pin(nl),wm,wp
      integer nl,nlayer,n1,n,nm,np
      do n1=1,nlayer
         if(p(n1) .gt. 3.5 .or. p(n1) .lt. 4.0e-6) then
            escout(n1) = 0.0
         else
            do n = 1,nl-1
               if (p(n1).le.pin(n).and.p(n1).ge.pin(n+1)) then
                  nm=n
                  np=n+1
                  wm=abs(pin(np)-p(n1))/(pin(nm)-pin(np))
                  wp=1.0 - wm
               endif
            enddo
            escout(n1) = escin(nm)*wm + escin(np)*wp
         endif
      enddo
      return
      end

c***********************************************************************

      subroutine interp3(esco1,esco2,esco3,p,nlayer,
     1     esci1,esci2,esci3,pin,nl)
C
C subroutine to perform 3 simultaneous linear interpolations in pressure from 
C 1D profiles esci1-3(nl) sampled on pressure grid pin(nl) to 1D profiles 
C esco1-3(nlayer) on pressure grid p(ngrid,nlayer).
C
      real esco1(nlayer),esco2(nlayer),esco3(nlayer),p(nlayer)
      real esci1(nl),    esci2(nl),    esci3(nl),    pin(nl),wm,wp
      integer nl,nlayer,n1,n,nm,np
      do n1=1,nlayer
         if (p(n1).gt. 3.5 .or. p(n1) .lt. 4.0e-6) then
            esco1(n1)=0.0
            esco2(n1)=0.0
            esco3(n1)=0.0
         else 
            do n = 1,nl-1
               if (p(n1).le.pin(n).and.p(n1).ge.pin(n+1)) then
                  nm=n
                  np=n+1
                  wm=abs(pin(np)-p(n1))/(pin(nm)-pin(np))
                  wp=1.0 - wm
               endif
            enddo
            esco1(n1) = esci1(nm)*wm + esci1(np)*wp
            esco2(n1) = esci2(nm)*wm + esci2(np)*wp
            esco3(n1) = esci3(nm)*wm + esci3(np)*wp
         endif
      enddo
      return
      end
