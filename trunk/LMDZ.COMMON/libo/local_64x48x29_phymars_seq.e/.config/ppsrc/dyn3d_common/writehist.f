










!
! $Id: writehist.F 1403 2010-07-01 09:02:53Z fairhead $
!
      subroutine writehist(time,vcov,ucov,teta,phi,q,masse,ps,phis)

      USE ioipsl
      USE infotrac, ONLY : nqtot, ttext
      use com_io_dyn_mod, only : histid,histvid,histuid
      USE temps_mod, ONLY: itau_dyn
      implicit none

C
C   Ecriture du fichier histoire au format IOIPSL
C
C   Appels succesifs des routines: histwrite
C
C   Entree:
C      time: temps de l'ecriture
C      vcov: vents v covariants
C      ucov: vents u covariants
C      teta: temperature potentielle
C      phi : geopotentiel instantane
C      q   : traceurs
C      masse: masse
C      ps   :pression au sol
C      phis : geopotentiel au sol
C      
C
C   L. Fairhead, LMD, 03/99
C
C =====================================================================
C
C   Declarations
!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 64,jjm=48,llm=29,ndm=1)

!-----------------------------------------------------------------------
!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!-----------------------------------------------------------------------
!   INCLUDE 'paramet.h'

      INTEGER  iip1,iip2,iip3,jjp1,llmp1,llmp2,llmm1
      INTEGER  kftd,ip1jm,ip1jmp1,ip1jmi1,ijp1llm
      INTEGER  ijmllm,mvar
      INTEGER jcfil,jcfllm

      PARAMETER( iip1= iim+1,iip2=iim+2,iip3=iim+3                       &
     &    ,jjp1=jjm+1-1/jjm)
      PARAMETER( llmp1 = llm+1,  llmp2 = llm+2, llmm1 = llm-1 )
      PARAMETER( kftd  = iim/2 -ndm )
      PARAMETER( ip1jm  = iip1*jjm,  ip1jmp1= iip1*jjp1 )
      PARAMETER( ip1jmi1= ip1jm - iip1 )
      PARAMETER( ijp1llm= ip1jmp1 * llm, ijmllm= ip1jm * llm )
      PARAMETER( mvar= ip1jmp1*( 2*llm+1) + ijmllm )
      PARAMETER( jcfil=jjm/2+5, jcfllm=jcfil*llm )

!-----------------------------------------------------------------------
!
! $Header$
!
!CDK comgeom
      COMMON/comgeom/                                                   &
     & cu(ip1jmp1),cv(ip1jm),unscu2(ip1jmp1),unscv2(ip1jm),             &
     & aire(ip1jmp1),airesurg(ip1jmp1),aireu(ip1jmp1),                  &
     & airev(ip1jm),unsaire(ip1jmp1),apoln,apols,                       &
     & unsairez(ip1jm),airuscv2(ip1jm),airvscu2(ip1jm),                 &
     & aireij1(ip1jmp1),aireij2(ip1jmp1),aireij3(ip1jmp1),              &
     & aireij4(ip1jmp1),alpha1(ip1jmp1),alpha2(ip1jmp1),                &
     & alpha3(ip1jmp1),alpha4(ip1jmp1),alpha1p2(ip1jmp1),               &
     & alpha1p4(ip1jmp1),alpha2p3(ip1jmp1),alpha3p4(ip1jmp1),           &
     & fext(ip1jm),constang(ip1jmp1),rlatu(jjp1),rlatv(jjm),            &
     & rlonu(iip1),rlonv(iip1),cuvsurcv(ip1jm),cvsurcuv(ip1jm),         &
     & cvusurcu(ip1jmp1),cusurcvu(ip1jmp1),cuvscvgam1(ip1jm),           &
     & cuvscvgam2(ip1jm),cvuscugam1(ip1jmp1),                           &
     & cvuscugam2(ip1jmp1),cvscuvgam(ip1jm),cuscvugam(ip1jmp1),         &
     & unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2,                 &
     & unsair_gam1(ip1jmp1),unsair_gam2(ip1jmp1),unsairz_gam(ip1jm),    &
     & aivscu2gam(ip1jm),aiuscv2gam(ip1jm),xprimu(iip1),xprimv(iip1)

!
        REAL                                                            &
     & cu,cv,unscu2,unscv2,aire,airesurg,aireu,airev,unsaire,apoln     ,&
     & apols,unsairez,airuscv2,airvscu2,aireij1,aireij2,aireij3,aireij4,&
     & alpha1,alpha2,alpha3,alpha4,alpha1p2,alpha1p4,alpha2p3,alpha3p4 ,&
     & fext,constang,rlatu,rlatv,rlonu,rlonv,cuvscvgam1,cuvscvgam2     ,&
     & cvuscugam1,cvuscugam2,cvscuvgam,cuscvugam,unsapolnga1,unsapolnga2&
     & ,unsapolsga1,unsapolsga2,unsair_gam1,unsair_gam2,unsairz_gam    ,&
     & aivscu2gam ,aiuscv2gam,cuvsurcv,cvsurcuv,cvusurcu,cusurcvu,xprimu&
     & , xprimv
!
!
! $Header$
!
!
! gestion des impressions de sorties et de débogage
! lunout:    unité du fichier dans lequel se font les sorties 
!                           (par defaut 6, la sortie standard)
! prt_level: niveau d'impression souhaité (0 = minimum)
!
      INTEGER lunout, prt_level
      COMMON /comprint/ lunout, prt_level

C
C   Arguments
C

      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm),phi(ip1jmp1,llm)                   
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL phis(ip1jmp1)                  
      REAL q(ip1jmp1,llm,nqtot)
      integer time


! This routine needs IOIPSL to work
C   Variables locales
C
      integer iq, ii, ll
      integer ndexu(ip1jmp1*llm),ndexv(ip1jm*llm),ndex2d(ip1jmp1)
      logical ok_sync
      integer itau_w
      REAL vnat(ip1jm,llm),unat(ip1jmp1,llm)

C
C  Initialisations
C
      ndexu = 0
      ndexv = 0
      ndex2d = 0
      ok_sync =.TRUE.
      itau_w = itau_dyn + time
!  Passage aux composantes naturelles du vent
      call covnat(llm, ucov, vcov, unat, vnat)
C
C  Appels a histwrite pour l'ecriture des variables a sauvegarder
C
C  Vents U
C
      call histwrite(histuid, 'u', itau_w, unat, 
     .               iip1*jjp1*llm, ndexu)
C
C  Vents V
C
      call histwrite(histvid, 'v', itau_w, vnat, 
     .               iip1*jjm*llm, ndexv)

C
C  Temperature potentielle
C
      call histwrite(histid, 'teta', itau_w, teta, 
     .                iip1*jjp1*llm, ndexu)
C
C  Geopotentiel
C
      call histwrite(histid, 'phi', itau_w, phi, 
     .                iip1*jjp1*llm, ndexu)
C
C  Traceurs
C
!        DO iq=1,nqtot
!          call histwrite(histid, ttext(iq), itau_w, q(:,:,iq), 
!     .                   iip1*jjp1*llm, ndexu)
!        enddo
!C
C  Masse
C
      call histwrite(histid,'masse',itau_w, masse,iip1*jjp1*llm,ndexu)
C
C  Pression au sol
C
      call histwrite(histid, 'ps', itau_w, ps, iip1*jjp1, ndex2d)
C
C  Geopotentiel au sol
C
!      call histwrite(histid, 'phis', itau_w, phis, iip1*jjp1, ndex2d)
C
C  Fin
C
      if (ok_sync) then
        call histsync(histid)
        call histsync(histvid)
        call histsync(histuid)
      endif
! of #ifdef 1
      return
      end
