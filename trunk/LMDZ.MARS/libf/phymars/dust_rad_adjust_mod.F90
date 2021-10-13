module dust_rad_adjust_mod

implicit none

real,save,allocatable :: dust_rad_adjust_prev(:) ! adjustment coefficient
                         ! computed when at current t_scenario
real,save,allocatable :: dust_rad_adjust_next(:) ! adjustment coefficient
                         ! computed for t_scenario of the next sol

contains

  subroutine compute_dust_rad_adjust(ngrid,nlayer,zday,pplev, &
                                     taudust,dust_rad_adjust)
 
  use geometry_mod, only: longitude_deg
  use time_phylmdz_mod, only: dtphys, daysec
  use dust_param_mod, only: odpref, t_scenario_sol
  
  implicit none
 
  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nlayer ! number of atmospheric levels
  real,intent(in) :: zday ! tim (in sols and fraction thereof)
  real,intent(in) :: pplev(ngrid,nlayer+1) ! pressure (Pa) at layer boundaries
  real,intent(in) :: taudust(ngrid) ! visible dust columns opacity in the GCM
  real,intent(out) :: dust_rad_adjust(ngrid) ! radiative adjustment coefficient 
                      ! for dust
  
  real,allocatable,save :: local_time(:) ! LT at current physics time step
  real,allocatable,save :: local_time_prevdt(:) ! LT at previous physics time step
  real :: zday_prevdt !value of zday at previous physics time step
  real,save :: zday_scenario ! to fetch dod values from the scenario
  real,save :: zday_scenario_next ! to fetch dod values from the scenario the next day
  logical,save :: firstcall=.true.
  integer :: ig
!  real,allocatable,save :: tau_pref_scenario(:)
  real,allocatable,save :: tau_pref_scenario_next(:)
  real :: weight ! interpolation weight
  real,save :: zday_prev_call=-666. ! stored value of zday from previous call
  
  ! 0. preliminary stuff
  ! NB: this routine may be called multiple times per physics
  ! so we have to save some arrays to store the information and not 
  ! recompute it for each call
  
  if (firstcall) then
    write(*,*) "compute_dust_rad_adjust: dust scenario assumed exact at", &
               " time(sol)=",t_scenario_sol
    allocate(local_time(ngrid))
    allocate(local_time_prevdt(ngrid))
!    allocate(tau_pref_scenario(ngrid))
    allocate(tau_pref_scenario_next(ngrid))
    firstcall=.false.
  endif ! of if firstcall
  
  ! 1. Compute local times (in sol fraction), if not already done
  if (zday/=zday_prev_call) then
    local_time(1:ngrid)=modulo(1.+(zday-INT(zday)) + &
                         (longitude_deg(1:ngrid)/15)/24,1.)
    zday_prevdt=zday-dtphys/daysec
    local_time_prevdt(1:ngrid)=modulo(1.+(zday_prevdt-INT(zday_prevdt)) + &
                         (longitude_deg(1:ngrid)/15)/24,1.)
    
    zday_scenario=zday-modulo(zday,1.) ! integer value of the day: the scenario
    ! opacity is assumed to be measured at 2pm but stored at nidnight
    zday_scenario_next=zday_scenario+1
  endif ! of if (zday/=zday_prev_call)
  
  ! 2. Load dust opacities for zday_scenario and zday_scenario_next
  !    if not already done
  if (zday/=zday_prev_call) then
!    call read_dust_scenario(ngrid,nlayer,zday_scenario,pplev,     &
!                                         tau_pref_scenario)
    call read_dust_scenario(ngrid,nlayer,zday_scenario_next,pplev, &
                                         tau_pref_scenario_next)
  endif ! of if (zday/=zday_prev_call)

  ! 3. Update dust_rad_adjust_* for grid points which just reached 2pm
  ! but only when this routine is called for the first time
  ! during this time step
  if (zday/=zday_prev_call) then
   do ig=1,ngrid
    if ((local_time(ig).ge.t_scenario_sol).and. &
                 (local_time_prevdt(ig).lt.(t_scenario_sol))) then
      ! store previous "next" as "prev" (NB we could also decide to recompute
      ! it using the current taudust...)
      dust_rad_adjust_prev(ig)=dust_rad_adjust_next(ig)
      ! compute new target based on current dust opacity
      dust_rad_adjust_next(ig)=tau_pref_scenario_next(ig)* &
                               pplev(ig,1)/odpref/taudust(ig)
    endif
   enddo
  endif ! of if (zday/=zday_prev_call)
  
  ! 4. Compute dust_rad_adjust using linear interpolation
  ! between dust_rad_adjust_prev and dust_rad_adjust_next
  do ig=1,ngrid
    ! prev and next are separated by a sol exactly
    ! we just need the distance (in sol) between current local time
    ! and 2pm the day before
    if (local_time(ig).ge.t_scenario_sol) then
      ! we are between t_scenario_sol and midnight
      weight=local_time(ig)-t_scenario_sol
    else
      ! we are between midnight and t_scenario_sol of the next day
      weight=(1.-t_scenario_sol)+local_time(ig)
    endif
    dust_rad_adjust(ig)=dust_rad_adjust_prev(ig)+ &
                        weight* &
                        (dust_rad_adjust_next(ig)-dust_rad_adjust_prev(ig))
  enddo! of do=ig=1,ngrid
  
  ! update zday_prev_call
  zday_prev_call=zday
  
  end subroutine compute_dust_rad_adjust

!=======================================================================
! Initialization of the module variables

  subroutine ini_dust_rad_adjust_mod(ngrid)
       
  implicit none
       
  integer, intent(in) :: ngrid
       
  allocate(dust_rad_adjust_prev(ngrid))
  allocate(dust_rad_adjust_next(ngrid))
       
  end subroutine ini_dust_rad_adjust_mod
       
  subroutine end_dust_rad_adjust_mod
       
  implicit none
       
  if (allocated(dust_rad_adjust_prev)) deallocate(dust_rad_adjust_prev)
  if (allocated(dust_rad_adjust_next)) deallocate(dust_rad_adjust_next)

  end subroutine end_dust_rad_adjust_mod

end module dust_rad_adjust_mod
