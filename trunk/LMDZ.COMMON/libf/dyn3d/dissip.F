!
! $Id: $
!
      SUBROUTINE dissip( vcov,ucov,teta,p, dv,du,dh )
c
      USE comconst_mod, ONLY: dtvr,dtdiss

      IMPLICIT NONE


c ..  Avec nouveaux operateurs star :  gradiv2 , divgrad2, nxgraro2  ...
c                                 (  10/01/98  )

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   Dissipation horizontale
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comdissnew.h"
#include "comdissipn.h"

c   Arguments:
c   ----------

      REAL,INTENT(IN) :: vcov(ip1jm,llm) ! covariant meridional wind
      REAL,INTENT(IN) :: ucov(ip1jmp1,llm) ! covariant zonal wind
      REAL,INTENT(IN) :: teta(ip1jmp1,llm) ! potential temperature
      REAL,INTENT(IN) :: p(ip1jmp1,llmp1) ! pressure
      ! tendencies (.../s) on covariant winds and potential temperature
      REAL,INTENT(OUT) :: dv(ip1jm,llm)
      REAL,INTENT(OUT) :: du(ip1jmp1,llm)
      REAL,INTENT(OUT) :: dh(ip1jmp1,llm)

c   Local:
c   ------

      REAL gdx(ip1jmp1,llm),gdy(ip1jm,llm)
      REAL grx(ip1jmp1,llm),gry(ip1jm,llm)
      REAL te1dt(llm),te2dt(llm),te3dt(llm)
      REAL deltapres(ip1jmp1,llm)

      INTEGER l,ij

      REAL  SSUM

c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

      DO l=1,llm
         te1dt(l) = tetaudiv(l) * dtdiss
         te2dt(l) = tetaurot(l) * dtdiss
         te3dt(l) = tetah(l)    * dtdiss
      ENDDO
      du=0.
      dv=0.
      dh=0.

c-----------------------------------------------------------------------
c   Calcul de la dissipation:
c   -------------------------

c   Calcul de la partie   grad  ( div ) :
c   -------------------------------------


      IF(lstardis) THEN
         CALL gradiv2( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ELSE
         CALL gradiv ( llm,ucov,vcov,nitergdiv,gdx,gdy )
      ENDIF

      DO l=1,llm

         DO ij = 1, iip1
            gdx(     ij ,l) = 0.
            gdx(ij+ip1jm,l) = 0.
         ENDDO

         DO ij = iip2,ip1jm
            du(ij,l) = du(ij,l) - te1dt(l) *gdx(ij,l)
         ENDDO
         DO ij = 1,ip1jm
            dv(ij,l) = dv(ij,l) - te1dt(l) *gdy(ij,l)
         ENDDO

       ENDDO

c   calcul de la partie   n X grad ( rot ):
c   ---------------------------------------

      IF(lstardis) THEN
         CALL nxgraro2( llm,ucov, vcov, nitergrot,grx,gry )
      ELSE
         CALL nxgrarot( llm,ucov, vcov, nitergrot,grx,gry )
      ENDIF


      DO l=1,llm
         DO ij = 1, iip1
            grx(ij,l) = 0.
         ENDDO

         DO ij = iip2,ip1jm
            du(ij,l) = du(ij,l) - te2dt(l) * grx(ij,l)
         ENDDO
         DO ij =  1, ip1jm
            dv(ij,l) = dv(ij,l) - te2dt(l) * gry(ij,l)
         ENDDO
      ENDDO

c   calcul de la partie   div ( grad ):
c   -----------------------------------

        
      IF(lstardis) THEN

       DO l = 1, llm
          DO ij = 1, ip1jmp1
            deltapres(ij,l) = AMAX1( 0.,  p(ij,l) - p(ij,l+1) )
          ENDDO
       ENDDO

         CALL divgrad2( llm,teta, deltapres  ,niterh, gdx )
      ELSE
         CALL divgrad ( llm,teta, niterh, gdx        )
      ENDIF

      DO l = 1,llm
         DO ij = 1,ip1jmp1
            dh( ij,l ) = dh( ij,l ) - te3dt(l) * gdx( ij,l )
         ENDDO
      ENDDO

      RETURN
      END
