










      SUBROUTINE writeg1d(ngrid,nx,x,nom,titre)
      IMPLICIT NONE
c.......................................................................
c
c  ecriture de x pour GRADS-1D
c
c  in :
c         * ngrid      ---> pour controler que l'on est bien en 1D
c         * nx         ---> taille du vecteur a stocker
c                             "1" pour une variable de surface
c                             "nlayer" pour une variable de centre de couche
c                             "nlayer+1" pour une variable d'interface
c         * x          ---> variable a stocker
c         * nom        ---> nom "pour grads"
c         * titre      ---> titre "pour grads"
c
c.......................................................................
c
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

c
c.......................................................................
c  declaration des arguments 
c
      INTEGER ngrid,nx,i
      REAL*4 xr4(1000)
      REAL x(nx)
      CHARACTER*(*) nom,titre
c
c  declaration des arguments 
c....................................................................... 
c  declaration des variables locales
c
      INTEGER ilayer,ivar
      LOGICAL test 
c
c  declaration des variables locales
c.......................................................................
c  controle 1D
c
c     print*,'ngrid=',ngrid
      IF (ngrid.NE.1) return
c
c  controle 1D
c.......................................................................
c  copy => force en reel 4 pour l'ecriture dans grads1d.dat

      do i=1,nx
        xr4(i) = x(i)
      enddo

c  copy => force en reel 4 pour l'ecriture dans grads1d.dat
c.......................................................................
c  ouverture du fichier au premier appel
c
      IF (g1d_premier) THEN
        OPEN (g1d_unitfich,FILE=g1d_nomfich
     &       ,FORM='unformatted',ACCESS='direct',RECL=4)
        g1d_irec=0
        g1d_nvar=0
        g1d_premier=.false.
      ENDIF
c
c  ouverture du fichier au premier appel
c.......................................................................
c  pour l'ecriture du fichier ctl
c
      test=.true.
      DO ivar=1,g1d_nvar
        IF (nom.EQ.g1d_nomvar(ivar)) test=.false.
        IF (nx .GT. 1000) then
          print*,'ERROR:  nx > 1000 dans writeg1d.F' 
          print*,'Changer la dimension de xr4'
          call exit(1)
        ENDIF
      ENDDO
      IF (test) THEN
        g1d_nvar=g1d_nvar+1
        g1d_nomvar(g1d_nvar)=nom
        g1d_titrevar(g1d_nvar)=titre
        IF (nx.EQ.1) THEN
           g1d_dimvar(g1d_nvar)=0
        ELSEIF (nx.EQ.g1d_nlayer) THEN
           g1d_dimvar(g1d_nvar)=g1d_nlayer
        ELSEIF (nx.EQ.g1d_nlayer+1) THEN
           g1d_dimvar(g1d_nvar)=g1d_nlayer+1
        ELSE
           PRINT *,'._. probleme de dimension dans GRADS-1D ._.'
           print*,'NX = ',nx
           print*,'g1d_nlayer = ',g1d_nlayer
        ENDIF
      ENDIF
c
c  pour l'ecriture du fichier ctl
c.......................................................................
c  ecriture
c
      IF (nx.EQ.1) THEN
        g1d_irec=g1d_irec+1
        WRITE(g1d_unitfich,REC=g1d_irec) xr4(1)
      ELSE
        DO ilayer=1,nx
          g1d_irec=g1d_irec+1
          WRITE(g1d_unitfich,REC=g1d_irec) xr4(ilayer)
        ENDDO
      ENDIF
c
c  ecriture
c.......................................................................
c
10001 CONTINUE
c
c.......................................................................
c
      RETURN
      END





c *********************************************************************
c *********************************************************************

      SUBROUTINE endg1d(ngrid,nlayer,zlayer,ndt)
      USE time_phylmdz_mod, ONLY: dtphys, daysec
      IMPLICIT NONE
c.......................................................................
c
c  ecriture du fichier de controle pour GRADS-1D
c
c  in :
c         * ngrid      ---> pour controler que l'on est bien en 1D
c         * nlayer     ---> nombre de couches
c         * zlayer     ---> altitude au centre de chaque couche (km)
c         * ndt        ---> nombre de pas de temps
c
c.......................................................................
c
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
c
c.......................................................................
c  declaration des arguments 
c
      INTEGER ngrid,nlayer
      REAL zlayer(nlayer)
      INTEGER ndt
c
c  declaration des arguments 
c....................................................................... 
c  declaration des variables locales
c
      INTEGER ivar,ilayer
c
c  declaration des variables locales
c.......................................................................
c  contole 1D
c
      IF (ngrid.NE.1) GOTO 10001
c
c  contole 1D
c.......................................................................
c
      IF (nlayer.ne.g1d_nlayer) 
     &PRINT *,'._. probleme de dimension dans GRADS-1D (endg1d.F) '
c
c.......................................................................
c
      CLOSE (g1d_unitfich)
c
c.......................................................................
c
      OPEN (g1d_unitctl,FILE=g1d_nomctl,FORM='formatted',RECL=4*100)
      WRITE (g1d_unitctl,'(a4,2x,a1,a20)') 'DSET','^',g1d_nomfich
      WRITE (g1d_unitctl,'(a5,2x,a20)') 'UNDEF ','1.E+30'
      WRITE (g1d_unitctl,'(a11)') 'FORMAT YREV'
      WRITE (g1d_unitctl,'(a5,2x,a30)') 'TITLE ','champs 1D'
      WRITE (g1d_unitctl,'(a5,i4,a20)') 'XDEF ',1,' LINEAR 0 1'
      WRITE (g1d_unitctl,'(a5,i4,a20)') 'YDEF ',1,' LINEAR 0 1'
      WRITE (g1d_unitctl,'(a5,i4,a20)') 'ZDEF ',g1d_nlayer,' LEVELS'
      WRITE (g1d_unitctl,'(5(1x,f13.5))')
     &      (zlayer(ilayer),ilayer=1,g1d_nlayer)

c     Writing true timestep in g1d.ctl (in planet "minutes"= sol/(60*24)
      ivar =min( max(int(1440.*dtphys/daysec +0.5),1) , 99)   
      WRITE (g1d_unitctl,'(a4,2x,i10,a19,i2,a2)')
     &      'TDEF ',ndt,' LINEAR 01JAN2000 ', ivar,'MN '

      WRITE (g1d_unitctl,'(a5,i5)') 'VARS ',g1d_nvar
      DO ivar=1,g1d_nvar
      WRITE (g1d_unitctl,'(a9,3x,i4,i3,1x,a39)') 
     &       g1d_nomvar(ivar),g1d_dimvar(ivar),99,g1d_titrevar(ivar)
      ENDDO
      WRITE (g1d_unitctl,'(a7)') 'ENDVARS'
      CLOSE (g1d_unitctl)
c
c.......................................................................
c
10001 CONTINUE
c
c.......................................................................
c
      RETURN
      END
