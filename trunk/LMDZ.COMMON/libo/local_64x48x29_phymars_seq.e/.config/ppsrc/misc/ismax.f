










!
! $Header$
!
      function ismax(n,sx,incx)
c
      IMPLICIT NONE
c
      INTEGER n,i,incx,ismax,ix
      real sx((n-1)*incx+1),sxmax
c
      ix=1
      ismax=1
      sxmax=sx(1)
      do 10 i=1,n-1
       ix=ix+incx
       if(sx(ix).gt.sxmax) then
         sxmax=sx(ix)
         ismax=i+1
       endif
10    continue
c
      return
      end

