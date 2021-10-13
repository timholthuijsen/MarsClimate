










!
! $Id: initfluxsto.F 1279 2009-12-10 09:02:56Z fairhead $
!
      subroutine initfluxsto
     .  (infile,tstep,t_ops,t_wrt,
     .                    fileid,filevid,filedid)

       USE IOIPSL
       USE comvert_mod, ONLY: nivsigs
       USE comconst_mod, ONLY: pi
       USE temps_mod, ONLY: annee_ref,day_ref,itau_dyn
      implicit none

C
C   Routine d'initialisation des ecritures des fichiers histoires LMDZ
C   au format IOIPSL
C
C   Appels succesifs des routines: histbeg
C                                  histhori
C                                  histver
C                                  histdef
C                                  histend
C
C   Entree:
C
C      infile: nom du fichier histoire a creer
C      day0,anne0: date de reference
C      tstep: duree du pas de temps en seconde
C      t_ops: frequence de l'operation pour IOIPSL
C      t_wrt: frequence d'ecriture sur le fichier
C
C   Sortie:
C      fileid: ID du fichier netcdf cree
C      filevid:ID du fichier netcdf pour la grille v
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

C   Arguments
C
      character*(*) infile
      real tstep, t_ops, t_wrt
      integer fileid, filevid,filedid

! This routine needs IOIPSL to work
C   Variables locales
C
      real nivd(1)
      integer tau0
      real zjulian
      character*3 str
      character*10 ctrac
      integer iq
      real rlong(iip1,jjp1), rlat(iip1,jjp1),rl(1,1)
      integer uhoriid, vhoriid, thoriid, zvertiid,dhoriid,dvertiid
      integer ii,jj
      integer zan, idayref
      logical ok_sync
C
C  Initialisations
C
      pi = 4. * atan (1.)
      str='q  '
      ctrac = 'traceur   '
      ok_sync = .true.
C
C  Appel a histbeg: creation du fichier netcdf et initialisations diverses
C         

      zan = annee_ref
      idayref = day_ref
      CALL ymds2ju(zan, 1, idayref, 0.0, zjulian)
      tau0 = itau_dyn
	
	do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonu(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo
 
      call histbeg(infile, iip1, rlong(:,1), jjp1, rlat(1,:),
     .             1, iip1, 1, jjp1,
     .             tau0, zjulian, tstep, uhoriid, fileid)
C
C  Creation du fichier histoire pour la grille en V (oblige pour l'instant,
C  IOIPSL ne permet pas de grilles avec des nombres de point differents dans 
C  un meme fichier)


      do jj = 1, jjm
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatv(jj) * 180. / pi
        enddo
      enddo

      call histbeg('fluxstokev.nc', iip1, rlong(:,1), jjm, rlat(1,:),
     .             1, iip1, 1, jjm,
     .             tau0, zjulian, tstep, vhoriid, filevid)
	
	rl(1,1) = 1.	
      call histbeg('defstoke.nc', 1, rl, 1, rl,
     .             1, 1, 1, 1,
     .             tau0, zjulian, tstep, dhoriid, filedid)

C
C  Appel a histhori pour rajouter les autres grilles horizontales
C
      do jj = 1, jjp1
        do ii = 1, iip1
          rlong(ii,jj) = rlonv(ii) * 180. / pi
          rlat(ii,jj) = rlatu(jj) * 180. / pi
        enddo
      enddo

      call histhori(fileid, iip1, rlong, jjp1, rlat, 'scalar',
     .              'Grille points scalaires', thoriid)
	
C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'sig_s', 'Niveaux sigma',
     . 'sigma_level',
     .              llm, nivsigs, zvertiid)
C Pour le fichier V
      call histvert(filevid, 'sig_s', 'Niveaux sigma',
     .  'sigma_level',
     .              llm, nivsigs, zvertiid)
c pour le fichier def
      nivd(1) = 1
      call histvert(filedid, 'sig_s', 'Niveaux sigma',
     .  'sigma_level',
     .              1, nivd, dvertiid)

C
C  Appels a histdef pour la definition des variables a sauvegarder
	
	CALL histdef(fileid, "phis", "Surface geop. height", "-",
     .                iip1,jjp1,thoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)

         CALL histdef(fileid, "aire", "Grid area", "-",
     .                iip1,jjp1,thoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
	
	CALL histdef(filedid, "dtvr", "tps dyn", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
        
         CALL histdef(filedid, "istdyn", "tps stock", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)
         
         CALL histdef(filedid, "istphy", "tps stock phy", "s",
     .                1,1,dhoriid, 1,1,1, -99, 32,
     .                "once", t_ops, t_wrt)


C
C Masse 
C
      call histdef(fileid, 'masse', 'Masse', 'kg',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Pbaru 
C
      call histdef(fileid, 'pbaru', 'flx de masse zonal', 'kg m/s',
     .             iip1, jjp1, uhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Pbarv 
C
      call histdef(filevid, 'pbarv', 'flx de masse mer', 'kg m/s',
     .             iip1, jjm, vhoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  w 
C
      call histdef(fileid, 'w', 'flx de masse vert', 'kg m/s',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)

C
C  Temperature potentielle
C
      call histdef(fileid, 'teta', 'temperature potentielle', '-',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C

C
C Geopotentiel 
C
      call histdef(fileid, 'phi', 'geopotentiel instantane', '-',
     .             iip1, jjp1, thoriid, llm, 1, llm, zvertiid,
     .             32, 'inst(X)', t_ops, t_wrt)
C
C  Fin
C
      call histend(fileid)
      call histend(filevid)
      call histend(filedid)
      if (ok_sync) then
        call histsync(fileid)
        call histsync(filevid)
        call histsync(filedid)
      endif
	
! of #ifdef 1
      return
      end
