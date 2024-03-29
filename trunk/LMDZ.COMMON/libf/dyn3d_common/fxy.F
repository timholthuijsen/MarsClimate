!
! $Id: fxy.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE fxy (rlatu,yprimu,rlatv,yprimv,rlatu1,yprimu1,
     ,                    rlatu2,yprimu2,
     , rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)

      USE comconst_mod, ONLY: pi
      USE serre_mod, ONLY: pxo,pyo,alphax,alphay,transx,transy

      IMPLICIT NONE

c     Auteur  :  P. Le Van
c
c     Calcul  des longitudes et des latitudes  pour une fonction f(x,y)
c           a tangente sinusoidale et eventuellement avec zoom  .
c
c
#include "dimensions.h"
#include "paramet.h"

       INTEGER i,j

       REAL rlatu(jjp1), yprimu(jjp1),rlatv(jjm), yprimv(jjm),
     , rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
       REAL rlonu(iip1),xprimu(iip1),rlonv(iip1),xprimv(iip1),
     , rlonm025(iip1),xprimm025(iip1), rlonp025(iip1),xprimp025(iip1)

#include "fxy_new.h"


c    ......  calcul  des  latitudes  et de y'   .....
c
       DO j = 1, jjm + 1 
          rlatu(j) = fy    ( REAL( j )        )
         yprimu(j) = fyprim( REAL( j )        )
       ENDDO


       DO j = 1, jjm

         rlatv(j)  = fy    ( REAL( j ) + 0.5  )
         rlatu1(j) = fy    ( REAL( j ) + 0.25 ) 
         rlatu2(j) = fy    ( REAL( j ) + 0.75 ) 

        yprimv(j)  = fyprim( REAL( j ) + 0.5  ) 
        yprimu1(j) = fyprim( REAL( j ) + 0.25 )
        yprimu2(j) = fyprim( REAL( j ) + 0.75 )

       ENDDO

c
c     .....  calcul   des  longitudes et de  x'   .....
c
       DO i = 1, iim + 1
           rlonv(i)     = fx    (   REAL( i )          )
           rlonu(i)     = fx    (   REAL( i ) + 0.5    )
        rlonm025(i)     = fx    (   REAL( i ) - 0.25  )
        rlonp025(i)     = fx    (   REAL( i ) + 0.25  )

         xprimv  (i)    = fxprim (  REAL( i )          )
         xprimu  (i)    = fxprim (  REAL( i ) + 0.5    )
        xprimm025(i)    = fxprim (  REAL( i ) - 0.25   )
        xprimp025(i)    = fxprim (  REAL( i ) + 0.25   )
       ENDDO

c
       RETURN
       END

