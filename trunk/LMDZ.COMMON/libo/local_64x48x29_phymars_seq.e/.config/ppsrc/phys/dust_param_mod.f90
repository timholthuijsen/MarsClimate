










module dust_param_mod
! This module contains flags and saved variables for the dust cycle 
implicit none
  LOGICAL,SAVE :: active ! is dust radiatively active?
  LOGICAL,SAVE :: doubleq ! use 2-moment schmeme for dust? (dustbin must then be >=2)
  LOGICAL,SAVE :: submicron ! include a secondary dust distribution (submicron particules)
  LOGICAL,SAVE :: lifting ! flag to activate injection of dust from the surface
  LOGICAL,SAVE :: freedust ! if true: no rescaling (via tauscaling) of the dust mass and number
  LOGICAL,SAVE :: callddevil ! flag to activate dust devil (dust lifing/injection) parametrization
  
  INTEGER,SAVE :: dustbin ! number of bins of dust tracers

  REAL,PARAMETER :: odpref = 610. ! Reference pressure (Pa) of
                     ! DOD (Dust optical Depth) tau_pref_*
  
  REAL,SAVE,ALLOCATABLE :: tauscaling(:)   ! Convertion factor for qdust and Ndust
  INTEGER,SAVE :: dustscaling_mode ! dust scaling modes
                  ! =0, no rescaling (freedust)
                  ! =1, prescribed scaling GCM5.3 style (using tauscaling)
                  ! =2, only radiative scaling (using dust_rad_adjust)
  REAL,SAVE,ALLOCATABLE :: dust_rad_adjust(:) ! radiative scaling for dust
  REAL,PARAMETER :: t_scenario_sol=14/24. ! time of day (sol) at which
                    ! tau_pref_scenario is deemed exact

contains

  subroutine ini_dust_param_mod(ngrid)
    implicit none
    integer,intent(in) :: ngrid ! number of atmospheric columns
  
    allocate(tauscaling(ngrid))
    allocate(dust_rad_adjust(ngrid))
  
  end subroutine ini_dust_param_mod

  subroutine end_dust_param_mod
    implicit none
    
    if (allocated(tauscaling))  deallocate(tauscaling)
    if (allocated(dust_rad_adjust)) deallocate(dust_rad_adjust)

  end subroutine end_dust_param_mod
end module dust_param_mod
