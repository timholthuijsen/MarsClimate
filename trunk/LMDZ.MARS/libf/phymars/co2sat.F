c=======================================================================
c     SUBROUTINE co2sat
c-----------------------------------------------------------------------
c     Aim:
c     ----
c     Compute saturated steam pressure (from James et al, 1992)
c=======================================================================
      subroutine co2sat(naersize, t, psat)
     
      implicit none
c-----------------------------------------------------------------------
c     VARIABLES
c-----------------------------------------------------------------------
c     Inputs:
c     ------- 
      integer, intent(in) ::
     &   naersize ! dimension of tables t and psat
     
      real, intent(in) ::
     &   t(naersize) ! temperature table

c     Output:
c     -------
      real, intent(out) ::
     &   psat(naersize) ! Saturated steam pressure (Pa)

c     Local:
c     ------
      integer ::
     &   i ! loop on naersize
c=======================================================================
c===== BEGIN
c=======================================================================
      do i = 1, naersize
        psat(i) = 1.382 * 1e12 * exp(-3182.48/t(i))
      end do
c=======================================================================
c===== END
c=======================================================================
      RETURN
      END

