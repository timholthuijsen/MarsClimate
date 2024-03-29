!
! $Header$
!
      SUBROUTINE divgrad (klevel,h, lh, divgra )
      IMPLICIT NONE
c
c=======================================================================
c
c  Auteur :   P. Le Van
c  ----------
c
c                              lh
c      calcul de  (div( grad ))   de h  .....
c      h  et lh  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
c=======================================================================
c
c   declarations:
c   -------------
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comdissipn.h"
c
      INTEGER klevel
      REAL h( ip1jmp1,klevel ), divgra( ip1jmp1,klevel )
c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)

      INTEGER  l,ij,iter,lh
c
c
c
      CALL SCOPY ( ip1jmp1*klevel,h,1,divgra,1 )
c
      DO 10 iter = 1,lh

      CALL filtreg ( divgra,jjp1,klevel,2,1,.true.,1  )

      CALL    grad (klevel,divgra, ghx  , ghy          )
      CALL  diverg (klevel,  ghx , ghy  , divgra       )

      CALL filtreg ( divgra,jjp1,klevel,2,1,.true.,1)

      DO 5 l = 1,klevel
      DO 4  ij = 1, ip1jmp1
      divgra( ij,l ) = - cdivh * divgra( ij,l )
   4  CONTINUE
   5  CONTINUE
c
  10  CONTINUE
      RETURN
      END
