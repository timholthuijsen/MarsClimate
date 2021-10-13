










! $Id: $
      SUBROUTINE abort_physic(modname, message, ierr)
     
      USE IOIPSL
      USE mod_phys_lmdz_para
      USE print_control_mod, ONLY: lunout
      IMPLICIT NONE
!
! Stops the simulation cleanly, closing files and printing various
! comments
!
!  Input: modname = name of calling program
!         message = stuff to print
!         ierr    = severity of situation ( = 0 normal )

      character(len=*), intent(in):: modname
      integer ierr, ierror_mpi
      character(len=*), intent(in):: message

      write(lunout,*) 'in abort_physic'
!$OMP MASTER
      call histclo
      call restclo
      if (mpi_rank .eq. 0) then
         call getin_dump
      endif
!$OMP END MASTER

      write(lunout,*) 'Stopping in ', modname
      write(lunout,*) 'Reason = ',message
      if (ierr .eq. 0) then
        write(lunout,*) 'Everything is cool'
      else
        write(lunout,*) 'Houston, we have a problem, ierr = ', ierr
        stop 1
      endif
      END
