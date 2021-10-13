










      subroutine lwdiff (kdlon,kflev
     .         ,pbsur,pbtop,pdbsl
     .         ,tautotal,omegtotal,gtotal
     .         ,pemis,pfluc)

      use dimradmars_mod, only: nir, npademx, nabsmx, nflev, ndlon,
     &                          ndlo2
      use yomlw_h, only: nlaylte
      USE comcstfi_h
      IMPLICIT NONE
 
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

C-----------------------------------------------------------------------
C
c  ABSORPTION ET DIFUSION HORS DE LA BANDE A 15 MICRONS :
c! 1) Dans la bande a 15 micron (CO2), les poussieres 
c! n'interviennent que comme un milieu gris non diffusif avec 
c!                      Q=Qext*(1-Omega)
c! cette bande est decoupee en deux sous bandes (indices 1 et 2)
c! pour lesquelles les parametres optiques des poussieres sont
c! identiques
c! 2)  le reste est decoupe en "nir-2" bandes : une bande qui recouvre toutes
c! les longueurs d'onde inferieures a celles de la bande a 15 microns
c! (indice 3) et nir-3 bandes pour les grandes longueurs d'onde
c! (indices 4...nir) sur chacune de ces  bandes, les poussieres
c! sont supposees diffusantes grises.
c!
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      integer kdlon,kflev
      REAL PBSUR(NDLO2,nir), PBTOP(NDLO2,nir)
     S  ,  PDBSL(NDLO2,nir,KFLEV*2), PEMIS(NDLO2)

      real PFLUC(NDLO2,2,KFLEV+1)
      real tautotal(ndlon,nflev,nir)
      real omegtotal(ndlon,nflev,nir), gtotal(ndlon,nflev,nir)

C
C
C-------------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
C
      integer jl, jk, ndd, indd, iir, j1
      integer  j2, j2dd2, j2dd1,j2bot,j2top, j2dd
      REAL ZADJD(NDLON,NFLEV+1), ZADJU(NDLON,NFLEV+1)
     S  ,  ZDBDT(NDLON,nir,NFLEV)
     S  ,  ZDISD(NDLON,NFLEV+1), ZDISU(NDLON,NFLEV+1)
     S  ,  ZFD(NDLON), ZFDN(NDLON,NFLEV+1), ZFU(NDLON)
     S  ,  ZFUP(NDLON,NFLEV+1),ZGLAYD(NDLON),ZGLAYU(NDLON)
     S  ,  ZOMEGADD(NDLON,NFLEV*2),ZGDD(NDLON,NFLEV*2)
     S  ,  ZTAUDD(NDLON,NFLEV*2)
     S  ,  ZBHDD(NDLON,NFLEV*2+1),ZBSDD(NDLON)
     S  ,  ZZBHDD(NDLON,NFLEV*2+1),ZZBSDD(NDLON)
     S  ,  ZFAHDD(NDLON,NFLEV*2+1),ZFDHDD(NDLON,NFLEV*2+1)
     S  ,  ZZFAHDD(NDLON,NFLEV*2+1),ZZFDHDD(NDLON,NFLEV*2+1)
C
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
 110  CONTINUE
C

      do jl = 1 , kdlon
        do jk = 1 , nlaylte
          PFLUC(jl,1,jk) = 0.
          PFLUC(jl,2,jk) = 0.
        enddo
      enddo

      DO 112 JK = 1 , nlaylte+1
        DO 111 JL = 1 , KDLON
          ZADJD(JL,JK) = 0.
          ZADJU(JL,JK) = 0.
          ZDISD(JL,JK) = 0.
          ZDISU(JL,JK) = 0.
 111    CONTINUE
 112  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
C     ------------------------------------------------------------------
C
C
C  ==================================================================
C*         2.0     contribution des bandes "hors co2"
C  ==================================================================
C
 200  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.0.1   preparation des Planck a chaque hauteur
C                  ----------------------------------
C
c!
c! le nombre de couche pour la diffusion sera le nombre de layer * 2
c! soit NDD=nlaylte*2, donc la taille du vecteur des Planck sera
c! nlaylte*2 + 1. la taille des vecteurs omega / g / tau sera
c! par contre nlaylte*2 (voir dans FLUSV.F).
c!
      NDD=nlaylte*2
      DO indd=1,ndd+1
                                            do jl=1,kdlon
         ZFAHDD(jl,indd)=0.
         ZFDHDD(jl,indd)=0.
         ZBHDD(jl,indd)=0.
                                            enddo
      ENDDO
                                            do jl=1,kdlon
      ZBSDD(jl)=0.
                                            enddo
c!
c! boucle sur les  bandes hors CO2
c!
      DO 10001 iir=3,nir
c!
                                            do jl=1,kdlon
        ZZBHDD(JL,1)=PBTOP(JL,iir)/pi
                                            enddo
        DO J1=2,NDD+1
                                            do jl=1,kdlon
           ZZBHDD(JL,J1)=
     &     ZZBHDD(JL,J1-1)-PDBSL(JL,iir,NDD-J1+2)/pi
                                            enddo
        ENDDO
                                            do jl=1,kdlon
        ZZBSDD(JL)=PBSUR(JL,iir)/pi
                                            enddo

C
C     ------------------------------------------------------------------
C
C*         2.0.2   preparation des coefficients de diffusion
C                  -----------------------------------------
C
c! les omega, g, tau ... boucle de bas en haut
        DO J2=1,nlaylte-1
          J2DD2=(nlaylte-J2+1)*2
          J2DD1=J2DD2-1
          J2BOT=3*J2-2
          J2TOP=3*J2+1
          do jl=1,kdlon
            ZTAUDD(JL,J2DD1)=tautotal(jl,J2,iir)*0.5
            ZTAUDD(JL,J2DD2)=ZTAUDD(JL,J2DD1)
            ZOMEGADD(JL,J2DD1)=omegtotal(jl,J2,iir)
            ZOMEGADD(JL,J2DD2)=omegtotal(jl,J2,iir)
            ZGDD(JL,J2DD1)=gtotal(jl,J2,iir)
            ZGDD(JL,J2DD2)=gtotal(jl,J2,iir)
          enddo
        ENDDO
        J2=nlaylte
        J2DD2=2
        J2DD1=1
        J2BOT=3*J2-2
                                            do jl=1,kdlon
        ZTAUDD(JL,J2DD1)= tautotal(jl,J2,iir)*0.5
        ZTAUDD(JL,J2DD2)= tautotal(jl,J2,iir)*0.5
        ZOMEGADD(JL,J2DD1)= omegtotal(jl,J2,iir)
        ZOMEGADD(JL,J2DD2)= omegtotal(jl,J2,iir)
        ZGDD(JL,J2DD1)= gtotal(jl,J2,iir)
        ZGDD(JL,J2DD2)= gtotal(jl,J2,iir)
                                            enddo
C
C     ------------------------------------------------------------------
C
C*         2.0.3   calcul de la diffusion
C                  ----------------------
C

c-----------------------------------------------------------------------
        CALL flusv(KDLON,0
     &  ,NDD,ZOMEGADD,ZGDD,ZTAUDD,PEMIS
     &  ,ZZBHDD,ZZBSDD
     &  ,ZZFAHDD,ZZFDHDD)
c!
c!  Cumul des flux sur le spectre hors bande du CO2
c!
        DO indd=1,ndd+1
           do jl=1,kdlon
             ZFAHDD(jl,indd)=ZFAHDD(jl,indd)+ZZFAHDD(jl,indd)
             ZFDHDD(jl,indd)=ZFDHDD(jl,indd)+ZZFDHDD(jl,indd)
           enddo
        ENDDO
10001 CONTINUE

      DO J2=1,nlaylte+1
        J2DD=(nlaylte-J2+1)*2+1
        do jl=1,kdlon
          PFLUC(JL,1,J2)=PFLUC(JL,1,J2)+ZFAHDD(JL,J2DD)
          PFLUC(JL,2,J2)=PFLUC(JL,2,J2)-ZFDHDD(JL,J2DD)
        enddo
      ENDDO


      END
