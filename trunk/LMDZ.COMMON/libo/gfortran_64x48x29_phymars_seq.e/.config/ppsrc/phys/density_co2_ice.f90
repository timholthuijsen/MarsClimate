










!======================================================================================================================!
! SUBROUTINE: density_co2_ice =========================================================================================!
!======================================================================================================================!
! Subject:
!---------
!   Compute co2 ice particles density
!----------------------------------------------------------------------------------------------------------------------!
! Reference:
!-----------
!   Mangan et al. (2017), 'CO2 ice structure and density under Martian atmospheric conditions', Icarus
!   Valid for: P = 1e-2 mbar, 80 K <= T <= 195 K
!======================================================================================================================!
module density_co2_ice_mod

  implicit none
  contains

  subroutine density_co2_ice(temperature, density)

  implicit none

  double precision, intent(in) :: temperature
  double precision, intent(out) :: density

  density = 1000. * (1.72391 - 2.53e-4*temperature - 2.87e-6*temperature*temperature)

  end subroutine density_co2_ice

end module density_co2_ice_mod
