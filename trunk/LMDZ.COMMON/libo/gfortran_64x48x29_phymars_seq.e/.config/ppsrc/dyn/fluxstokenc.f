










!
! $Id: fluxstokenc.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE fluxstokenc(pbaru,pbarv,masse,teta,phi,phis,
     . time_step,itau )
! This routine is designed to work with ioipsl

       USE IOIPSL
c
c     Auteur :  F. Hourdin
c
c
ccc   ..   Modif. P. Le Van  ( 20/12/97 )  ...
c
      IMPLICIT NONE
c
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
      common /tracstoke/istdyn,istphy,unittrac
      integer istdyn,istphy,unittrac
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

      REAL time_step,t_wrt, t_ops
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL masse(ip1jmp1,llm),teta(ip1jmp1,llm),phi(ip1jmp1,llm)
      REAL phis(ip1jmp1)

      REAL pbaruc(ip1jmp1,llm),pbarvc(ip1jm,llm)
      REAL massem(ip1jmp1,llm),tetac(ip1jmp1,llm),phic(ip1jmp1,llm)

      REAL pbarug(ip1jmp1,llm),pbarvg(iip1,jjm,llm),wg(ip1jmp1,llm)

      REAL pbarvst(iip1,jjp1,llm),zistdyn
	real dtcum

      INTEGER iadvtr,ndex(1) 
      integer nscal
      real tst(1),ist(1),istp(1)
      INTEGER ij,l,irec,i,j,itau
      INTEGER, SAVE :: fluxid, fluxvid,fluxdid
 
      SAVE iadvtr, massem,pbaruc,pbarvc,irec
      SAVE phic,tetac
      logical first
      save first
      data first/.true./
      DATA iadvtr/0/


c AC initialisations
      pbarug(:,:)   = 0.
      pbarvg(:,:,:) = 0.
      wg(:,:)       = 0.
      

      if(first) then

	CALL initfluxsto( 'fluxstoke',
     .  time_step,istdyn* time_step,istdyn* time_step,
     .  fluxid,fluxvid,fluxdid) 
	
	ndex(1) = 0
        call histwrite(fluxid, 'phis', 1, phis, iip1*jjp1, ndex)
        call histwrite(fluxid, 'aire', 1, aire, iip1*jjp1, ndex)
	
	ndex(1) = 0
        nscal = 1
        tst(1) = time_step
        call histwrite(fluxdid, 'dtvr', 1, tst, nscal, ndex)
        ist(1)=istdyn
        call histwrite(fluxdid, 'istdyn', 1, ist, nscal, ndex)
        istp(1)= istphy
        call histwrite(fluxdid, 'istphy', 1, istp, nscal, ndex)
	
	first = .false.

      endif


      IF(iadvtr.EQ.0) THEN
         phic(:,:)=0
         tetac(:,:)=0
         pbaruc(:,:)=0
         pbarvc(:,:)=0
      ENDIF

c   accumulation des flux de masse horizontaux
      DO l=1,llm
         DO ij = 1,ip1jmp1
            pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
            tetac(ij,l) = tetac(ij,l) + teta(ij,l)
            phic(ij,l) = phic(ij,l) + phi(ij,l)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
         ENDDO
      ENDDO

c   selection de la masse instantannee des mailles avant le transport.
      IF(iadvtr.EQ.0) THEN
         CALL SCOPY(ip1jmp1*llm,masse,1,massem,1)
      ENDIF

      iadvtr   = iadvtr+1


c   Test pour savoir si on advecte a ce pas de temps
      IF ( iadvtr.EQ.istdyn ) THEN
c    normalisation
      DO l=1,llm
         DO ij = 1,ip1jmp1
            pbaruc(ij,l) = pbaruc(ij,l)/REAL(istdyn)
            tetac(ij,l) = tetac(ij,l)/REAL(istdyn)
            phic(ij,l) = phic(ij,l)/REAL(istdyn)
         ENDDO
         DO ij = 1,ip1jm
            pbarvc(ij,l) = pbarvc(ij,l)/REAL(istdyn)
         ENDDO
      ENDDO

c   traitement des flux de masse avant advection.
c     1. calcul de w
c     2. groupement des mailles pres du pole.

        CALL groupe( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )

        do l=1,llm
           do j=1,jjm
              do i=1,iip1
                 pbarvst(i,j,l)=pbarvg(i,j,l)
              enddo
           enddo
           do i=1,iip1
              pbarvst(i,jjp1,l)=0.
           enddo
        enddo

         iadvtr=0
	write(lunout,*)'ITAU auquel on stoke les fluxmasses',itau
	
	call histwrite(fluxid, 'masse', itau, massem,
     .               iip1*jjp1*llm, ndex)
	
	call histwrite(fluxid, 'pbaru', itau, pbarug,
     .               iip1*jjp1*llm, ndex)
	
	call histwrite(fluxvid, 'pbarv', itau, pbarvg,
     .               iip1*jjm*llm, ndex)
	
        call histwrite(fluxid, 'w' ,itau, wg, 
     .             iip1*jjp1*llm, ndex) 
	
	call histwrite(fluxid, 'teta' ,itau, tetac, 
     .             iip1*jjp1*llm, ndex) 
	
	call histwrite(fluxid, 'phi' ,itau, phic, 
     .             iip1*jjp1*llm, ndex) 
	
C

      ENDIF ! if iadvtr.EQ.istdyn

! of #ifdef 1
      RETURN
      END
