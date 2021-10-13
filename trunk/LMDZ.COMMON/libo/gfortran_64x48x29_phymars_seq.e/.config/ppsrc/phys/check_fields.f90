










module check_fields_mod

  real,parameter :: default_temp_min=50.  ! minimum reasonable temperature (K)
  real,parameter :: default_temp_max=350. ! maximum reasonable temperature (K)

  real,parameter :: default_wind_max=500. ! maximum reasonable wind magnitude (m/s)

  real,parameter :: default_ps_min=80.  ! minimum reasonable surface pressure (Pa)
  real,parameter :: default_ps_max=2000. ! maximum reasonable surface pressure (Pa)

contains

subroutine check_physics_fields(message,temp,u,v,pplev)
use dimphy, only: klon, klev
implicit none
character(len=*),intent(in):: message 
real,intent(in) :: temp(klon,klev)
real,intent(in) :: u(klon,klev) ! zonal wind (m/s)
real,intent(in) :: v(klon,klev) ! meridional wind (m/s)
real,intent(in) :: pplev(klon,klev+1) ! pressure at level interfaces (Pa)

character(len=50) :: name="check_physics_fields"
logical ok_t,ok_w,ok_ps

! 1. Initialisations
ok_t=.true.
ok_w=.true.
ok_ps=.true.

! 2. Check temperature, winds and surface pressure
call check_temperature(temp,ok_t)
call check_winds(u,v,ok_w)
call check_ps(pplev(:,1),ok_ps)

if ((.not.ok_t).or.(.not.ok_w).or.(.not.ok_ps)) then
  ! Something is wrong, might as well stop here
  call abort_physic(trim(name),trim(message)//" Invalid field values",1)
endif

end subroutine check_physics_fields


subroutine check_temperature(temp,ok,temp_min,temp_max)
use dimphy, only: klon, klev
implicit none
real,intent(in) :: temp(klon,klev)
logical,intent(out) :: ok ! returns .true. if everything OK, .false. otherwise
real,intent(in),optional :: temp_min ! user provided acceptable minimum
real,intent(in),optional :: temp_max ! user provided acceptable maximum

character(len=50) :: name="check_temperature"
real :: tmin,tmax
integer :: i,k

! 0. Check optional inputs
if (present(temp_min)) then
  tmin=temp_min
else
  tmin=default_temp_min
endif

if (present(temp_max)) then
  tmax=temp_max
else
  tmax=default_temp_max
endif

! 1. initializations
ok=.true.

! 2. do the checking
do i=1,klon
  do k=1,klev
    ! look for NaN
    if (temp(i,k).ne.temp(i,k)) then
      ok=.false.
      write(*,*)trim(name)//" temp(i,k)=",temp(i,k)," for i=",i," k=",k
    endif
    ! check for temperatures too low
    if (temp(i,k).lt.tmin) then
      ok=.false.
      write(*,*)trim(name)//" temp(i,k)=",temp(i,k)," for i=",i," k=",k,&
      "<",tmin
    endif
    ! check for temperatures too high
    if (temp(i,k).gt.tmax) then
      ok=.false.
      write(*,*)trim(name)//" temp(i,k)=",temp(i,k)," for i=",i," k=",k,&
      ">",tmax
    endif
  enddo
enddo

end subroutine check_temperature

subroutine check_winds(u,v,ok,wind_max)
use dimphy, only: klon, klev
implicit none
real,intent(in) :: u(klon,klev) ! zonal wind (m/s)
real,intent(in) :: v(klon,klev) ! meridional wind (m/s)
logical,intent(out) :: ok ! returns .true. if everything OK, .false. otherwise
real,intent(in),optional :: wind_max ! user provided acceptable maximum magnitude

character(len=50) :: name="check_winds"
real :: wmax
integer :: i,k

! 0. Check optional inputs

if (present(wind_max)) then
  wmax=wind_max
else
  wmax=default_wind_max
endif

! 1. initializations
ok=.true.

! 2. do the checking
do i=1,klon
  do k=1,klev
    ! look for NaN
    if (u(i,k).ne.u(i,k)) then
      ok=.false.
      write(*,*)trim(name)//" u(i,k)=",u(i,k)," for i=",i," k=",k
    endif
    if (v(i,k).ne.v(i,k)) then
      ok=.false.
      write(*,*)trim(name)//" v(i,k)=",v(i,k)," for i=",i," k=",k
    endif
    ! check for magnitudes too high
    if (abs(u(i,k)).gt.wmax) then
      ok=.false.
      write(*,*)trim(name)//" u(i,k)=",u(i,k)," for i=",i," k=",k,&
      ">",wmax
    endif
    if (abs(v(i,k)).gt.wmax) then
      ok=.false.
      write(*,*)trim(name)//" v(i,k)=",v(i,k)," for i=",i," k=",k,&
      ">",wmax
    endif
  enddo
enddo

end subroutine check_winds

subroutine check_ps(ps,ok,ps_min,ps_max)
use dimphy, only: klon
implicit none
real,intent(in) :: ps(klon)
logical,intent(out) :: ok ! returns .true. if everything OK, .false. otherwise
real,intent(in),optional :: ps_min ! user provided acceptable minimum
real,intent(in),optional :: ps_max ! user provided acceptable maximum

character(len=50) :: name="check_ps"
real :: pmin,pmax
integer :: i

! 0. Check optional inputs
if (present(ps_min)) then
  pmin=ps_min
else
  pmin=default_ps_min
endif

if (present(ps_max)) then
  pmax=ps_max
else
  pmax=default_ps_max
endif

! 1. initializations
ok=.true.

! 2. do the checking
do i=1,klon
  ! look for NaN
  if (ps(i).ne.ps(i)) then
    ok=.false.
    write(*,*)trim(name)//" ps(i)=",ps(i)," for i=",i
  endif
  ! check for pressures too low
  if (ps(i).lt.pmin) then
    ok=.false.
    write(*,*)trim(name)//" ps(i)=",ps(i)," for i=",i,"<",pmin
  endif
  ! check for pressures too high
  if (ps(i).gt.pmax) then
    ok=.false.
    write(*,*)trim(name)//" ps(i)=",ps(i)," for i=",i,">",pmax
  endif
enddo

end subroutine check_ps

end module check_fields_mod
