










subroutine read_dust_scenario(ngrid,nlayer,zday,pplev,tau_pref_scenario)

! Reading of the dust scenario file

use netcdf
use geometry_mod, only: latitude, longitude ! in radians
use datafile_mod, only: datadir
use dust_param_mod, only: odpref, dustscaling_mode
use planete_h, only: year_day
implicit none

include "callkeys.h"

integer, intent(in) :: ngrid,nlayer
real, intent(in) :: zday ! date in martian days
real, dimension(ngrid,nlayer+1), intent(in) :: pplev
real, dimension(ngrid), intent(out) :: tau_pref_scenario ! visible dust column
                                       ! opacity at odpref from scenario

! Local variables

real :: realday
integer nid,nvarid,ierr
integer ltloop,lsloop,iloop,jloop,varloop,ig
real, dimension(2) :: taubuf
real tau1(4),tau
real alt(4)
integer latp(4),lonp(4)
integer yinf,ysup,xinf,xsup,tinf,tsup
real latinf,latsup,loninf,lonsup
real latintmp,lonintmp,latdeg,londeg
real colat,dlat,dlon,colattmp
logical, save :: firstcall=.true.
logical :: timeflag
real,save :: radeg,pi
integer :: timedim,londim,latdim
real, dimension(:), allocatable, save :: lat,lon,time
real, dimension(:,:,:), allocatable, save :: tautes
integer, save :: timelen,lonlen,latlen
character(len=33),save :: filename

realday=mod(zday,669.)

! AS: firstcall OK absolute
if (firstcall) then
   firstcall=.false.

   pi=acos(-1.)
   radeg=180/pi
   
   ! assimilated dust file: (NB: iaervar is a common in "callkeys.h")
   ! iaervar=4 means read dust_tes.nc file
   ! iaervar=6 means read dust_cold.nc file
   ! iaervar=7 means read dust_warm.nc file
   ! iaervar=8 means read dust_clim.nc file
   ! iaervar=24 means read dust_MY24.nc file
   ! iaervar=25 means read dust_MY25.nc file
   ! iaervar=26 means read dust_MY26.nc file, etc.
   if (iaervar.eq.4) then
   ! NB: 4: old TES assimilated MY24 dust scenarios (at 700Pa ref pressure!)
     filename="dust_tes.nc"
   else if (iaervar.eq.6) then
     filename="dust_cold.nc"
   else if (iaervar.eq.7) then
     filename="dust_warm.nc"
   else if (iaervar.eq.8) then
     filename="dust_clim.nc"
   else if ((iaervar.ge.24).and.(iaervar.le.35)) then
     write(filename,'(a7,i2,a3)')"dust_MY",iaervar,".nc"
   ! 124,125,126: old TES assimilated dust scenarios (at 700Pa ref pressure!)
   else if  (iaervar.eq.124) then
     filename="dust_tes_MY24.nc" 
   else if (iaervar.eq.125) then
     filename="dust_tes_MY25.nc" 
   else if (iaervar.eq.126) then
     filename="dust_tes_MY26.nc"
   endif
   
   ierr=nf90_open(trim(datadir)//"/"//trim(filename),NF90_NOWRITE,nid)
   IF (ierr.NE.nf90_noerr) THEN
     write(*,*)'Problem opening ',trim(filename),' (in phymars/read_dust_scenario.F90)'
     write(*,*)'It should be in :',trim(datadir),'/'
     write(*,*)'1) You can change this directory address in callfis.def with'
     write(*,*)'   datadir=/path/to/datafiles'
     write(*,*)'2) If necessary, dust*.nc (and other datafiles)'
     write(*,*)'   can be obtained online on:'
     write(*,*)'   http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
     CALL abort_physic("read_dust_scenario","Cannot find file",1)
   ENDIF

   ierr=nf90_inq_dimid(nid,"Time",timedim)
   if (ierr.ne.nf90_noerr) then
    ! 'Time' dimension not found, look for 'time'
    ierr=nf90_inq_dimid(nid,"time",timedim)
    if (ierr.ne.nf90_noerr) then
      write(*,*)"Error: read_dust_scenario <time> not found"
    endif
   endif
   ierr=nf90_inquire_dimension(nid,timedim,len=timelen)
   ierr=nf90_inq_dimid(nid,"latitude",latdim)
   ierr=nf90_inquire_dimension(nid,latdim,len=latlen)
   ierr=nf90_inq_dimid(nid,"longitude",londim)
   ierr=nf90_inquire_dimension(nid,londim,len=lonlen)


   if (.not.allocated(tautes)) allocate(tautes(lonlen,latlen,timelen))
   if (.not.allocated(lat)) allocate(lat(latlen), lon(lonlen), time(timelen))

   ! "dustop" if loading visible extinction opacity
   ! "cdod" if loading IR absorption opacity
   ierr=nf90_inq_varid(nid,"dustop",nvarid)
   if (ierr.eq.nf90_noerr) then
     ierr=nf90_get_var(nid,nvarid,tautes)
     IF (ierr .NE. nf90_noerr) THEN
        PRINT*, "Error: read_dust_scenario <dustop> not found"
        write(*,*)trim(nf90_strerror(ierr))
        call abort_physic("read_dust_scenario","dustop not found",1)
     ENDIF
   else
     ! did not find "dustop" , look for "cdod"
     ierr=nf90_inq_varid(nid,"cdod",nvarid)
     ierr=nf90_get_var(nid,nvarid,tautes)
     IF (ierr .NE. nf90_noerr) THEN
        PRINT*, "Error: read_dust_scenario <cdod> not found"
        write(*,*)trim(nf90_strerror(ierr))
        call abort_physic("read_dust_scenario","cdod not found",1)
     ENDIF
     ! and multiply by 2*1.3=2.6 to convert from IR absorption
     ! to visible extinction opacity
     tautes(:,:,:)=2.6*tautes(:,:,:)
   endif

   ierr=nf90_inq_varid(nid,"Time",nvarid)
   ierr=nf90_get_var(nid,nvarid,time)
   IF (ierr .NE. nf90_noerr) THEN
      PRINT*, "Error: read_dust_scenario <Time> not found"
      write(*,*)trim(nf90_strerror(ierr))
      call abort_physic("read_dust_scenario","Time not found",1)
   ENDIF

   ! Adapt time axis values if using dustinjection scheme and/or
   ! dustscaling_mode == 2
   IF ((dustinjection>0).or.(dustscaling_mode==2)) THEN
     ! Sanity check: we expect to have one map per sol
     if (timelen/=year_day) then
       write(*,*)"read_dust_scenario error: timelen=",timelen
       write(*,*)" whereas it was exepceted to be year_day=",year_day
       call abort_physic("read_dust_scenario","timelen/=year_day",1)
     endif
     ! fill time() with intergers from 0 to timelen-1
     write(*,*)"read_dust_scenario: adapting time axis to integer values"
     do iloop=1,timelen
       time(iloop)=iloop-1
     enddo
   ENDIF ! of IF ((dustinjection>0).or.(dustscaling_mode==2))

   ierr=nf90_inq_varid(nid,"latitude",nvarid)
   ierr=nf90_get_var(nid,nvarid,lat)
   IF (ierr .NE. nf90_noerr) THEN
      PRINT*, "Error: read_dust_scenario <latitude> not found"
      write(*,*)trim(nf90_strerror(ierr))
      call abort_physic("read_dust_scenario","latitude not found",1)
   ENDIF

   ierr=nf90_inq_varid(nid,"longitude",nvarid)
   ierr=nf90_get_var(nid,nvarid,lon)
   IF (ierr .NE. nf90_noerr) THEN
      PRINT*, "Error: read_dust_scenario <longitude> not found"
      write(*,*)trim(nf90_strerror(ierr))
      call abort_physic("read_dust_scenario","longitude not found",1)
   ENDIF

   ierr=nf90_close(nid)

endif ! of if (firstcall)

do ig=1,ngrid

! Find the four nearest points, arranged as follows:
!                               1 2
!                               3 4

   latdeg=latitude(ig)*radeg  ! latitude, in degrees
   londeg=longitude(ig)*radeg  ! longitude, in degrees east
   colat=90-latdeg        ! colatitude, in degrees

! Ehouarn: rounding effects and/or specific compiler issues
!          sometime cause londeg to be sligthly below -180 ...
  if (londeg.lt.-180) then
    ! check if it is by a large amount
    if ((londeg+180).lt.-1.e-3) then
      write(*,*) ' ig=',ig,' londeg=',londeg
      call abort_physic("read_dust_scenario","longitude way out of range",1)
    else
      londeg=-180
    endif
  endif

 ! Find encompassing latitudes
   if (colat<(90-lat(1))) then ! between north pole and lat(1)
      ysup=1
      yinf=1
   else if (colat>=90-(lat(latlen))) then ! between south pole and lat(laten)
      ysup=latlen
      yinf=latlen
   else ! general case
      do iloop=2,latlen
         if(colat<(90-lat(iloop))) then
           ysup=iloop-1
           yinf=iloop
           exit
         endif
      enddo
   endif

   latinf=lat(yinf)
   latsup=lat(ysup)


 ! Find encompassing longitudes
   ! Note: in input file, lon(1)=-180.
   if (londeg>lon(lonlen)) then
      xsup=1
      xinf=lonlen
      loninf=lon(xsup)
      lonsup=180.0 ! since lon(1)=-180.0
   else
      do iloop=2,lonlen
         if(londeg<lon(iloop)) then
              xsup=iloop
              xinf=iloop-1
              exit
         endif
      enddo
      loninf=lon(xinf)
      lonsup=lon(xsup)
   endif

   if ((xsup.gt.lonlen).OR.(yinf.gt.latlen).OR.(xinf.lt.1)&
     .OR.(ysup.lt.1)) then
      write (*,*) "read_dust_scenario: SYSTEM ERROR on x or y in ts_gcm"
      write (*,*) "xinf: ",xinf
      write (*,*) "xsup: ",xsup
      write (*,*) "yinf: ",yinf
      write (*,*) "ysup: ",ysup
      call abort_physic("read_dust_scenario","invalid coordinates",1)
   endif

!   loninf=lon(xinf)
!   lonsup=lon(xsup)
!   latinf=lat(yinf)
!   latsup=lat(ysup)

! The four neighbouring points are arranged as follows:
!                               1 2
!                               3 4

   latp(1)=ysup
   latp(2)=ysup
   latp(3)=yinf
   latp(4)=yinf

   lonp(1)=xinf
   lonp(2)=xsup
   lonp(3)=xinf
   lonp(4)=xsup

! Linear interpolation on time, for all four neighbouring points

   if ((realday<time(1)).or.(realday>time(timelen))) then
      tinf=timelen
      tsup=1
      timeflag=.true.
   else 
      timeflag=.false.
      do iloop=2,timelen
         if (realday<=time(iloop)) then
            tinf=iloop-1
            tsup=iloop
            exit
         endif
      enddo
   endif

! Bilinear interpolation on the four nearest points

   if ((colat<(90-lat(1))).OR.(colat>(90-lat(latlen))).OR.(latsup==latinf)) then
      dlat=0
   else
      dlat=((90-latinf)-colat)/((90-latinf)-(90-latsup))
   endif

   if (lonsup==loninf) then
      dlon=0
   else
      dlon=(londeg-loninf)/(lonsup-loninf)
   endif

   do iloop=1,4
      taubuf(1)=tautes(lonp(iloop),latp(iloop),tinf)
      taubuf(2)=tautes(lonp(iloop),latp(iloop),tsup)
      if (timeflag) then
         if (realday>time(timelen)) then
            tau1(iloop)=taubuf(1)+(taubuf(2)-taubuf(1))*(realday-time(tinf))/(time(tsup)+(669-time(tinf))) 
         else
            tau1(iloop)=taubuf(1)+(taubuf(2)-taubuf(1))*realday/(time(tsup)+(669-time(tinf)))
         endif
      else
         tau1(iloop)=taubuf(1)+(taubuf(2)-taubuf(1))*(realday-time(tinf))/(time(tsup)-time(tinf))
      endif
      if (tau1(iloop)<0) then
          write (*,*) "read_dust_scenario: SYSTEM ERROR on tau"
          write (*,*) "utime ",realday
          write (*,*) "time(tinf) ",time(tinf)
          write (*,*) "time(tsup) ",time(tsup)
          write (*,*) "tau1 ",taubuf(1)
          write (*,*) "tau2 ",taubuf(2)
          write (*,*) "tau ",tau1(iloop)
          call abort_physic("read_dust_scenario","negative tau!",1)
      endif
   enddo

   if ((dlat>1).OR.(dlon>1) .OR. (dlat<0) .OR. (dlon<0)) then
      write (*,*) "read_dust_scenario: SYSTEM ERROR on dlat or dlon in ts_gcm"
      write (*,*) "dlat: ",dlat
      write (*,*) "lat: ",latdeg
      write (*,*) "dlon: ",dlon
      write (*,*) "lon: ",londeg
      call abort_physic("read_dust_scenario","negative dlat or dlon!",1)
   endif

   tau= dlat*(dlon*(tau1(2)+tau1(3)-tau1(1)-tau1(4))+tau1(1)-tau1(3)) +dlon*(tau1(4)-tau1(3))+tau1(3)

   tau_pref_scenario(ig)=tau
! 
enddo ! of ig=1,ngrid

if (filename(1:8)=="dust_tes") then
  ! when using old TES files, correction for:
  ! - the reference pressure was 700Pa (unlike 610Pa now)
  ! - the 1.3 conversion factor from IR absorbtion opacity to
  !   IR extinction opacity
  tau_pref_scenario(:)=tau_pref_scenario(:)*1.3*(odpref/700.)
endif

if (swrtype.eq.1) then ! Fouquart (NB: swrtype is set in callkeys.h)
 ! when using old radiative transfer (like in MCD 4.x)
 ! needed to decrease opacity (*0.825) to compensate overestimation of
 ! heating rates
  tau_pref_scenario(:)=tau_pref_scenario(:)*0.825/1.3
endif

end
