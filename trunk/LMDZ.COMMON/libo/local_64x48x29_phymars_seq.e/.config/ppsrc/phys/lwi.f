










      subroutine lwi (ig0,kdlon,kflev
     .                ,psi,zdblay,pdp
     .                ,newpcolc )


      use dimradmars_mod, only: ndlo2, ndlon, nflev, nir
      use yomlw_h, only: gcp, nlaylte, xi
      USE comcstfi_h, ONLY: g, cpp
      USE time_phylmdz_mod, ONLY: dtphys
      implicit none

c.......................................................................
c  le COMMON pour GRADS-1D
c  (Utilise pour les sorties format Grads dans la version 1D du modele)
c
c  on peut se dire : "on ne sauvera pas plus de 1000 variables ... hein ?"
c
      INTEGER g1d_nvarmx
      PARAMETER(g1d_nvarmx=1000)
c
c         * g1d_nlayer     ---> nombre de couches verticales
c         * g1d_nomfich    ---> nom du fichier grads
c         * g1d_unitfich   ---> code du fichier grads
c         * g1d_nomctl     ---> nom du fichier ctl
c         * g1d_unitctl    ---> code du fichier ctl
c         * g1d_premier    ---> variable logique pour dire si le fichier
c                               est deja ouvert
c         * g1d_irec       ---> indice de derniere ecriture
c         * g1d_nvar       ---> nombre de variables deja definies a la
c                               derniere ecriture
c         * g1d_nomvar     ---> noms des vecteurs existants
c         * g1d_dimvar     ---> taille des vecteurs
c         * g1d_titrevar   ---> titres des vecteurs
c         * g1d_tmp1       ---> caractere 
c         * g1d_tmp2       ---> caractere 
c
      INTEGER g1d_nlayer
      CHARACTER*100 g1d_nomfich
      INTEGER g1d_unitfich
      CHARACTER*100 g1d_nomctl
      INTEGER g1d_unitctl
      LOGICAL g1d_premier
      LOGICAL g2d_premier
      INTEGER g1d_irec
      INTEGER g2d_irec
      INTEGER g2d_appel
      INTEGER g1d_nvar
      CHARACTER*100 g1d_nomvar
      INTEGER g1d_dimvar
      CHARACTER*100 g1d_titrevar
      CHARACTER*100 g1d_tmp1,g1d_tmp2
c
      COMMON/COMG1DI/g1d_nlayer
     &             ,g1d_unitfich
     &             ,g1d_unitctl
     &             ,g1d_irec
     &             ,g2d_irec
     &             ,g2d_appel
     &             ,g1d_nvar
      COMMON/COMG1DC/g1d_dimvar(0:g1d_nvarmx)
     &             ,g1d_nomfich
     &             ,g1d_nomctl
     &             ,g1d_nomvar(0:g1d_nvarmx)
     &             ,g1d_titrevar(0:g1d_nvarmx)
     &             ,g1d_tmp1
     &             ,g1d_tmp2
      COMMON/COMG1DL/g1d_premier
     &             ,g2d_premier
c
c.......................................................................
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
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                             -   lwi    -    
C
C     PURPOSE:       Shema semi - implicite 
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c************************************************************************
c
c        0.    Declarations
c              ------------
c
c-------------------------------------------------------------------------
c        0.1   Arguments
c              ---------
c
 
      integer ig0,kdlon,kflev

      real    psi(ndlo2,kflev)
     .     ,  zdblay(ndlo2,nir,kflev)
     .     ,  pdp(ndlo2,kflev)


      real    newpcolc(ndlo2,kflev)

c-------------------------------------------------------------------------
c        0.2   local arrays
c              ------------
c
      real    di(ndlon,nflev) 
     .      , hi(ndlon,nflev) 
     .      , bi(ndlon,nflev) 

      real    ci(ndlon,nflev) 
     .      , ai(ndlon,nflev) 
      real    deltat

      real   semit, denom

      integer i, jl

c************************************************************************
c
c        1.    Initialisations
c              ---------------
c
c-----------------------------------------------------------------------
 
        deltat = dtphys * iradia
c       print*,'SEMI = ',semi, '(expl:0  semi-implicite:0.5  impl:1)'
        semit = semi * deltat
c       semi = 0.

c       print*,'dtphys,iradia,deltat,semit:',dtphys,iradia,deltat,semit
c       print*,'g,cpp',g,cpp


c************************************************************************
c
c        2.    
c              ---------------
c
c-------------------------------------------------------------------------
c        2.1   Calcul des di
c              -------------
c


      do i = 1 , nlaylte-1
        do jl = 1 , kdlon
c     -------------------
      di(jl,i) =  1 + semit * (g / pdp(jl,i) / cpp) * (
     .    ( xi(ig0+jl,1,i,nlaylte+1)
     .    + xi(ig0+jl,1,i,i+1)
     .    + xi(ig0+jl,1,i,i-1) )
     .    *    zdblay(jl,1,i)
     .  + ( xi(ig0+jl,2,i,nlaylte+1)
     .    + xi(ig0+jl,2,i,i+1)
     .    + xi(ig0+jl,2,i,i-1) )
     .    *    zdblay(jl,2,i)
     .     )
c     -------------------
        enddo
      enddo

c couche nlaylte
c ------------
c      , on enleve i,i+1 sinon on a 2 fois le cooling2space

      do jl = 1 , kdlon
c     -------------------
      di(jl,nlaylte) =  1 + semit * (g / pdp(jl,nlaylte) / cpp) * (
     .    ( xi(ig0+jl,1,nlaylte,nlaylte+1)
     .    + xi(ig0+jl,1,nlaylte,nlaylte-1) )
     .    *    zdblay(jl,1,nlaylte)
     .  + ( xi(ig0+jl,2,nlaylte,nlaylte+1)
     .    + xi(ig0+jl,2,nlaylte,nlaylte-1) )
     .    *    zdblay(jl,2,nlaylte)
     .     )
c     -------------------
      enddo

c-------------------------------------------------------------------------
c        2.2   Calcul des hi
c              -------------
c

      do i = 1 , nlaylte-1
        do jl = 1 , kdlon
c     -------------------
      hi(jl,i) =    - semit * (g / pdp(jl,i) / cpp) *
     .            (    xi(ig0+jl,1,i,i+1) * zdblay(jl,1,i+1)   
     .               + xi(ig0+jl,2,i,i+1) * zdblay(jl,2,i+1)   )
c     -------------------
        enddo
      enddo

c-------------------------------------------------------------------------
c        2.3   Calcul des bi
c              -------------
c


      do i = 2 , nlaylte
        do jl = 1 , kdlon
c     -------------------
      bi(jl,i) =   - semit * (g / pdp(jl,i) / cpp) * 
     .           (     xi(ig0+jl,1,i,i-1) * zdblay(jl,1,i-1)   
     .               + xi(ig0+jl,2,i,i-1) * zdblay(jl,2,i-1)   )
c     -------------------
        enddo
      enddo


c couche 1
c --------
c  tant qu'on a pas un calcul propre de zdblay(0) qui tienne compte de 
c    la discontinuite de temperature au sol , on met  b1 = 0


      do jl = 1 , kdlon
        bi(jl,1) = 0 
      enddo

c-------------------------------------------------------------------------
c        2.4   
c              -------------
c

c couche nlaylte
c ------------

      do jl = 1 , kdlon
c     -------------------
      ci(jl,nlaylte) = (gcp * psi(jl,nlaylte) / pdp(jl,nlaylte))
     .                   / di(jl,nlaylte)

      ai(jl,nlaylte) = - bi(jl,nlaylte) / di(jl,nlaylte)
c     -------------------
      enddo



      do i = nlaylte-1 , 1 , -1
        do jl = 1 , kdlon
c     -------------------
      denom = di(jl,i) + hi(jl,i) * ai(jl,i+1)

      ci(jl,i) = (  gcp * psi(jl,i) / pdp(jl,i)
     .             - hi(jl,i) * ci(jl,i+1)  )  / denom
 
      ai(jl,i) = -bi(jl,i) / denom
c     -------------------
        enddo
      enddo


c-------------------------------------------------------------------------
c        2.5   
c              -------------
c

c couche 1
c -------
      do jl = 1 , kdlon
        newpcolc(jl,1) = ci(jl,1)
      enddo


      do i = 2 , nlaylte
        do jl = 1 , kdlon
           newpcolc(jl,i) = ci(jl,i) + ai(jl,i) * newpcolc(jl,i-1)
        enddo
      enddo



c-------------------------------------------------------------------------
      RETURN
      END
