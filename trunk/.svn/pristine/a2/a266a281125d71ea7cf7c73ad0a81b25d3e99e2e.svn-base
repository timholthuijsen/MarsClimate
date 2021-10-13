

program lslin

! Program to interpolate data in Solar Longitude linear time coordinate
! (usable with grads) from Netcdf diagfi or concatnc  files
! author Y. Wanherdrick, K. Dassas 2004
! Modified by F.Forget 04/2005
! More modifications by Ehouarn Millour 12/2005
! Modified by Ehouarn Millour 10/2007 (changed evaluation of 'start_var'
! from hard-coded values to a computed value) 
! Read controle field, if available TN, October 2013
! Allow to chose a starting Ls and to average within Ls timestep a
! (great to compare with TES and MCS data. F. Forget december 2013)


implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=50) :: infile
! infile(): input file name
character (len=50), dimension(:), allocatable :: var
! var(): names of the variables
character (len=50) :: long_name,units
! long_name(),units(): [netcdf] "long_name" and "units" attributes
character (len=50) :: altlong_name,altunits,altpositive
! altlong_name(): [netcdf] altitude "long_name" attribute
! altunits(): [netcdf] altitude "units" attribute
! altpositive(): [netcdf] altitude "positive" attribute

character (len=100) :: outfile
! outfile(): output file name
character (len=100) :: vartmp
! vartmp(): used for temporary storage of various strings
character (len=1) :: answer
! answer: user's answers at multiple questions (especially, answer='y'/'Y')
character (len=1) :: average
! average: the character 'y' to average within the Ls timestep
integer :: nid,varid,ierr,miss,validr
! nid: [netcdf] ID # of input file
! varid: [netcdf] ID # of a variable
! ierr: [netcdf] error code (returned by subroutines)
integer :: nout,varidout
! nout: [netcdf] ID # of output file
! varidout: [netcdf] ID # of a variable (to write in the output file)
integer :: i,j,k,l,x,y,n
! counters for various loops
integer :: start_var
! starting index/ID # from which physical variables are to be found
integer :: reptime ! Ehouarn: integer or real ? 
! rep_time: starting date/time of the run (given by user)
integer :: day_ini ! Ehouarn: integer or real ?
! day_ini: starting date/time of the run stored in input file
real, dimension(:), allocatable:: lat,lon,alt,time,lsgcm,timels,ctl
! lat(): latitude coordinates (read from input file)
! lon(): longitude coordinates (read from input file)
! alt(): altitude coordinates (read from input file)
! time(): time coordinates (in "sol", read from input file)
! lsgcm(): time coordinate (in unevenly spaced "Ls")
! timels(): new time coordinates (evenly spaced "Ls"; written to output file)
! ctl(): array, stores controle array
integer :: latlen,lonlen,altlen,timelen,Nls,Sls,ctllen
! latlen: # of elements of lat() array
! lonlen: # of elements of lon() array
! altvar: # of elements of alt() array
! timelen: # of elements of time() and lsgcm() arrays
! Nls: # of elements of timels() array
integer :: nbvarfile,nbvar,ndim !,nbfile
! nbvar: # of time-dependent variables
! nbvarfile: total # of variables (in netcdf file)
! ndim: [netcdf] # of dimensions (3 or 4) of a variable
integer :: latdim,londim,altdim,timedim,ctldim
! latdim: [netcdf] "latitude" dim ID
! londim: [netcdf] "longitude" dim ID
! altdim: [netcdf] "altdim" dim ID
! timedim: [netcdf] "timedim" dim ID
! ctldim: [netcdf] "controle" dim ID
integer :: latvar,lonvar,altvar,timevar,ctlvar
! latvar: [netcdf] ID of "latitude" variable
! lonvar: [netcdf] ID of "longitude" variable
! altvar: [netcdf] ID of "altitude" variable
! timevar: [netcdf] ID of "Time" variable
integer :: latdimout,londimout,altdimout,timedimout,timevarout
! latdimout: [netcdf] output latitude (dimension) ID
! londimout: [netcdf] output longitude (dimension) ID
! altdimout: [netcdf] output altitude (dimension) ID
! timedimout: [netcdf] output time (dimension) ID
! timevarout: [netcdf] ID of output "Time" variable
integer, dimension(4) :: corner,edges,dim
! corner(4): [netcdf] start indexes (where block of data will be written)
! edges(4): [netcdf] lenght (along dimensions) of block of data to write
! dim(4): [netcdf] lat, long, alt and time dimensions
real, dimension(:,:,:,:), allocatable :: var3d,var3dls
! var3d(,,,): 4D array to store a variable (on initial lat/lon/alt/sol grid)
! var3dls(,,,): 4D array to store a variable (on new lat/lon/alt/Ls grid)
real, dimension(:), allocatable :: var3dxy
! var3dxy(): to store the temporal evolution of a variable (at fixed lat/lon/alt)
real :: deltatsol,deltals,resultat,ls0
! deltatsol: time step (in sols) of input file data
! deltals: time step (in Ls) for the data sent to output
! ls0: first Ls value for the data sent to output
! resultat: to temporarily store the result of the interpolation
character (len=3) :: mon
! mon(3): to store (and write to file) the 3 first letters of a month
real :: date
! date: used to compute/build a fake date (for grads output)
real :: missing
! to handle missing value when reading /writing files
real, dimension(2) :: valid_range
! valid range

!==============================================================================
! 1. Initialisation step
!==============================================================================

!==============================================================================
! 1.1. Get input file name and 'type' (to initialize start_var and reptime)
!==============================================================================

write(*,*) "which file do you want to use?"
read(*,'(a50)') infile 

!==============================================================================
! 1.2. Open input file and read/list the variables it contains
!==============================================================================

write(*,*) "opening "//trim(infile)//"..."
ierr = NF_OPEN(infile,NF_NOWRITE,nid)
if (ierr.NE.NF_NOERR) then
   write(*,*) 'Failed to open file '//infile
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! Compute 'start_var', the index from which variables are lon-lat-time
! and/or lon-lat-alt-time
! WARNING: We assume here that the ordering of variables in the NetCDF
! file is such that 0D, 1D and 2D variables are placed BEFORE 3D and 4D
! variables

i=1 ! dummy initialization to enter loop below
start_var=0 ! initialization
do while (i.lt.3)
  start_var=start_var+1
  ! request number of dims of variable number 'start_var'
  ierr=NF_INQ_VARNDIMS(nid,start_var,i)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Failed to get number of dims for variable number:",start_var
    write(*,*) NF_STRERROR(ierr)
    stop ""
  endif
enddo
write(*,*) "start_var=",start_var


ierr=NF_INQ_NVARS(nid,nbvarfile)
! nbvarfile is now equal to the (total) number of variables in input file
                                       
! Get the variables' names from the input file (and store them in var())
write(*,*) ""
nbvar=nbvarfile-(start_var-1)
allocate(var(nbvar))
do i=1,nbvar
   ierr=NF_INQ_VARNAME(nid,i+(start_var-1),var(i))
   write(*,'(a9,1x,i2,1x,a1,1x,a50)') "variable ",i,":",var(i)
enddo

!==============================================================================
! 1.3. Output file name
!==============================================================================
outfile=infile(1:len_trim(infile)-3)//"_Ls.nc"
write(*,*) ""
write(*,*) "Output file : "//trim(outfile)


!==============================================================================
! 2. Work: read input, build new time coordinate and write it to output
!==============================================================================

!==============================================================================
! 2.1. Read (and check) latitude, longitude and altitude from input file
!==============================================================================

   ierr=NF_INQ_DIMID(nid,"latitude",latdim)
   ierr=NF_INQ_VARID(nid,"latitude",latvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <latitude> is missing'
      write(*,*) NF_STRERROR(ierr)
      stop ""  
   endif
   ierr=NF_INQ_DIMLEN(nid,latdim,latlen)
!  write(*,*) "latlen: ",latlen

   ierr=NF_INQ_DIMID(nid,"longitude",londim)
   ierr=NF_INQ_VARID(nid,"longitude",lonvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <longitude> is missing'
      write(*,*) NF_STRERROR(ierr)
      stop "" 
   endif
   ierr=NF_INQ_DIMLEN(nid,londim,lonlen)
!  write(*,*) "lonlen: ",lonlen

   ierr=NF_INQ_DIMID(nid,"altitude",altdim)
   ierr=NF_INQ_VARID(nid,"altitude",altvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <altitude> is missing'
      write(*,*) NF_STRERROR(ierr)
!     stop ""
       altlen = 1 
   else
      ierr=NF_INQ_DIMLEN(nid,altdim,altlen)
      
    ! Get altitude attributes to handle files with any altitude type
      ierr=nf_get_att_text(nid,altvar,'long_name',altlong_name)
      ierr=nf_get_att_text(nid,altvar,'units',altunits)
      ierr=nf_get_att_text(nid,altvar,'positive',altpositive)
   endif
!   write(*,*) "altlen: ",altlen


! Allocate lat(), lon() and alt()
   allocate(lat(latlen))
   allocate(lon(lonlen))
   allocate(alt(altlen))

! Read lat(),lon() and alt() from input file
   ierr = NF_GET_VAR_REAL(nid,latvar,lat)
   ierr = NF_GET_VAR_REAL(nid,lonvar,lon)
   ierr = NF_GET_VAR_REAL(nid,altvar,alt)
   if (ierr.NE.NF_NOERR) then
       if (altlen.eq.1) alt(1)=0.
   end if


   ierr=NF_INQ_DIMID(nid,"index",ctldim)
   if (ierr.NE.NF_NOERR) then
      write(*,*) ' Dimension <index> is missing in file '//trim(infile)
      ctllen=0
      !stop ""  
   endif
   ierr=NF_INQ_VARID(nid,"controle",ctlvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'Field <controle> is missing in file '//trim(infile)
      ctllen=0
      !stop ""
   else
      ierr=NF_INQ_DIMLEN(nid,ctldim,ctllen)
   endif

   allocate(ctl(ctllen),stat=ierr)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to allocate ctl(ctllen)"
     stop
   endif

   if (ctllen .ne. 0) then
     ierr = NF_GET_VAR_REAL(nid,ctlvar,ctl)
     if (ierr.ne.0) then
       write(*,*) "Error: failed to load controle"
       write(*,*) NF_STRERROR(ierr)
       stop
     endif
   endif ! of if (ctllen .ne. 0)



!==============================================================================
! 2.2. Create (and initialize) latitude, longitude and altitude in output file
!==============================================================================

  ! Initialize output file's lat,lon,alt and time dimensions
   call initiate(outfile,lat,lon,alt,ctl,latlen,lonlen,altlen,ctllen,&
         altlong_name,altunits,altpositive,&
         nout,latdimout,londimout,altdimout,timedimout,timevarout)

  ! Initialize output file's aps,bps and phisinit variables
   call init2(nid,lonlen,latlen,altlen,altdim,&
              nout,londimout,latdimout,altdimout)

   write(*,*)""
!==============================================================================
! 2.3. Read time from input file
!==============================================================================

ierr=NF_INQ_DIMID(nid,"Time",timedim)
if (ierr.NE.NF_NOERR) then
   write(*,*) 'ERROR: Dimension <Time> is missing'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif
ierr=NF_INQ_VARID(nid,"Time",timevar)
if (ierr.NE.NF_NOERR) then
   write(*,*) 'ERROR: Field <Time> is missing'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

ierr=NF_INQ_DIMLEN(nid,timedim,timelen)
write(*,*) "timelen: ",timelen

allocate(time(timelen))
allocate(lsgcm(timelen))

ierr = NF_GET_VAR_REAL(nid,timevar,time)

!==============================================================================
! 2.4. Initialize day_ini (starting day of the run)
!==============================================================================

if (ctllen .ne. 0) then ! controle is present in input file
  write(*,*) "Do you want to specify the beginning day of the file (y/n) ?"
  write(*,*) "If 'n', the value stored in the controle field will be used."
  read(*,*) answer
else ! controle is not present in the input file
  write(*,*) "There is no field controle in file "//trim(infile)
  write(*,*) "You have to specify the beginning day of the file."
  answer="y"
endif

if ((answer=="y").or.(answer=="Y")) then
   write(*,*) "Beginning day of the file?"
   read(*,*) reptime
endif

write(*,*) 'answer : ',answer
if ((answer/="y").and.(answer/="Y")) then
   ! the starting date of the run
   ! is stored in the "controle" array
   
   day_ini = nint(ctl(4))                                               
   day_ini = modulo(day_ini,669)                                                    
   write(*,*) 'day_ini', day_ini
else
   day_ini= reptime
   write(*,*) 'day_ini', day_ini
endif

!==============================================================================
! 2.5. Build timels() (i.e.: time, but in evenly spaced "Ls" time steps)
!==============================================================================

! compute time step, in sols, of input dataset
deltatsol=time(2)-time(1)
write(*,*) 'deltatsol',deltatsol

! compute time/dates in "Ls"; result stored in lsgcm()
do i=1,timelen
   call sol2ls(day_ini+time(i),lsgcm(i)) 
  if (lsgcm(i).lt.lsgcm(1)) lsgcm(i)=lsgcm(i) + 360.
enddo

write(*,*) 'Input data LS range :', lsgcm(1),lsgcm(timelen)

! generate the time step of the new Ls axis for the output
write(*,*) "Automatic Ls timestep (y/n)?"
read(*,*) answer
if ((answer=="y").or.(answer=="Y")) then
!   *********************************************
!   Trick of the trade:
!   Use the value of deltatsol to determine
!   a suitable value for deltals
!   *********************************************
    deltals=1./12.
    if (0.6*deltatsol.ge.1/6.) deltals=1./6.
    if (0.6*deltatsol.ge.0.25) deltals=0.25
    if (0.6*deltatsol.ge.0.5) deltals=0.5
    if (0.6*deltatsol.ge.0.75) deltals=0.75
    if (0.6*deltatsol.ge.1) deltals=1.
    if (0.6*deltatsol.ge.2) deltals=2.
    if (0.6*deltatsol.ge.3) deltals=3.
    if (0.6*deltatsol.ge.5) deltals=5.
    ls0=lsgcm(1) ! First Ls date
else
    write(*,*) "New timestep in Ls (deg) :"
    read(*,*) deltals
    write(*,*) "First Ls (deg) :"
    read(*,*) ls0
    if (ls0.lt.lsgcm(1)) then 
       write(*,*) 'with this file, the earliest Ls is ',lsgcm(1),'let s use it'
       ls0=lsgcm(1)
    end if
    write(*,*) "First Ls date (deg) = ", ls0  ! FF2013
endif

NLs=int(int((lsgcm(timelen)-ls0)/deltals) +1)   ! FF2013
allocate(timels(Nls))

! Build timels()
timels(1) = 0.01*nint(100.*ls0) ! 2 decimals
do k=2,Nls
   timels(k) = timels(k-1) + deltals
end do

write(*,*) 'New timestep in Ls (deg) ', deltals
write(*,*) 'Output data LS range : ', timels(1),timels(Nls)

! create averaged bins or interpolate at exact Ls ?
write(*,*) "Average data within the Ls timestep (y/n) or interpolate ? "
read(*,*) average
write(*,*)""

!==============================================================================
! 2.6. Write timels() to output file
!==============================================================================

do k=1,Nls
  ierr= NF_PUT_VARA_REAL(nout,timevarout,k,1,timels(k))
enddo

!==============================================================================
! 3. Read variables, interpolate them along the new time coordinate
!    and send the result to output
!==============================================================================

do j=1,nbvar ! loop on the variables to read/interpolate/write

!==============================================================================
! 3.1 Check that variable exists, and get some of its attributes
!==============================================================================
   write(*,*) "variable ",var(j)
   ! Get the variable's ID
   ierr=NF_INQ_VARID(nid,var(j),varid)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <',var(j),'> not found'
      write(*,*) NF_STRERROR(ierr)
      stop "Better stop now..."
   endif
   
   ! Get the value of 'ndim' for this varriable
   ierr=NF_INQ_VARNDIMS(nid,varid,ndim)
   write(*,*) 'ndim', ndim
   
   ! Check that it is a 3D or 4D variable
   if (ndim.lt.3) then
     write(*,*) "Error:",trim(var(j))," is neither a 3D nor a 4D variable"
     write(*,*) "We'll skip that variable..."
     CYCLE ! go directly to the next variable
   endif

!==============================================================================
! 3.2 Prepare a few things in order to interpolate/write
!==============================================================================

   if (ndim==3) then
      allocate(var3d(lonlen,latlen,timelen,1))
      allocate(var3dls(lonlen,latlen,Nls,1))
      allocate(var3dxy(timelen))

      dim(1)=londimout
      dim(2)=latdimout
      dim(3)=timedimout

      corner(1)=1
      corner(2)=1
      corner(3)=1
      corner(4)=1

      edges(1)=lonlen 
      edges(2)=latlen 
      edges(3)=Nls
      edges(4)=1

   else if (ndim==4) then
      allocate(var3d(lonlen,latlen,altlen,timelen))
      allocate(var3dls(lonlen,latlen,altlen,Nls))
      allocate(var3dxy(timelen))

      dim(1)=londimout
      dim(2)=latdimout
      dim(3)=altdimout
      dim(4)=timedimout

      corner(1)=1
      corner(2)=1
      corner(3)=1
      corner(4)=1

      edges(1)=lonlen
      edges(2)=latlen
      edges(3)=altlen
      edges(4)=Nls
   endif

!==============================================================================
! 3.3 Write this variable's definition and attributes to the output file
!==============================================================================
   units=" "
   long_name=" "
   ierr=nf_get_att_text(nid,varid,"long_name",long_name)
   if (ierr.ne.nf_noerr) then
   ! if no attribute "long_name", try "title"
     ierr=nf_get_att_text(nid,varid,"title",long_name)
   endif
   ierr=nf_get_att_text(nid,varid,"units",units)

   call def_var(nout,var(j),long_name,units,ndim,dim,varidout,ierr)

!==============================================================================
! 3.4 Read variable
!==============================================================================

   ierr = NF_GET_VAR_REAL(nid,varid,var3d)
   miss=NF_GET_ATT_REAL(nid,varid,"missing_value",missing)
   validr=NF_GET_ATT_REAL(nid,varid,"valid_range",valid_range)

!==============================================================================
! 3.6 interpolate or average
!==============================================================================

! 2D variable  (lon, lat, time)    
! interpolation of var at timels
   if (ndim==3) then
      do x=1,lonlen
       do y=1,latlen
!        write(*,*) 'd: x, y', x, y
        do l=1,timelen
          var3dxy(l)=var3d(x,y,l,1)
        enddo
        do n=1,Nls
          if(average.eq.'y') then
            resultat=0.
            Sls=0 ! (gcm data counter within each Ls timestep) 
            do l=1,timelen
               if((lsgcm(l).ge.(timels(n)-deltals/2.)).and.   &
               (lsgcm(l).lt.(timels(n)+deltals/2.))) then
                  if(var3dxy(l) .ne. missing) then
                    Sls= Sls +1 
                    resultat = resultat + var3dxy(l)
                  end if 
               end if 
               if (Sls.ne.0) then
                  var3dls(x,y,n,1)=resultat/float(Sls)
               else
                  var3dls(x,y,n,1)=missing
               endif 
            enddo
          else   ! average = 'n'
            call interpolf(timels(n),resultat,missing,lsgcm,var3dxy,timelen)
            var3dls(x,y,n,1)=resultat
          endif
        enddo
       enddo
      enddo
! 3D variable (lon, lat, alt, time)      
! interpolation of var at timels
   else if (ndim==4) then
      do x=1,lonlen
       do y=1,latlen
        do k=1,altlen
        do l=1,timelen
        var3dxy(l)=var3d(x,y,k,l)
        enddo
        do n=1,Nls
          if(average.eq.'y') then
            resultat=0.
            Sls=0 ! (gcm data counter within each Ls timestep) 
            do l=1,timelen
               if((lsgcm(l).ge.(timels(n)-deltals/2.)).and.   &
               (lsgcm(l).lt.(timels(n)+deltals/2.))) then
                  if(var3dxy(l) .ne. missing) then
                    Sls= Sls +1 
                    resultat = resultat + var3dxy(l)
                  end if 
               end if 
               if (Sls.ne.0) then
                  var3dls(x,y,k,n)=resultat/float(Sls)
               else
                  var3dls(x,y,k,n)=missing
               endif 
            enddo
          else   ! average = 'n'
             call interpolf(timels(n),resultat,missing,lsgcm,var3dxy,timelen)
             var3dls(x,y,k,n)=resultat
          endif
         enddo
        enddo
       enddo
      enddo
   endif

!==============================================================================
! 3.7 Write variable to output file
!==============================================================================

   ierr= NF_PUT_VARA_REAL(nout,varidout,corner,edges,var3dls)

   if (ierr.ne.NF_NOERR) then
     write(*,*) 'PUT_VAR ERROR: ',NF_STRERROR(ierr)
     stop ""
   endif

! In case there is a "missing value" attribute and "valid range"
   if (miss.eq.NF_NOERR) then
     call missing_value(nout,varidout,missing)  
   endif

   deallocate(var3d)
   deallocate(var3dls)
   deallocate(var3dxy)

enddo ! of do j=1,nbvar loop

deallocate(time)
deallocate(lsgcm)

! close input and output files
ierr=nf_close(nid)
ierr=NF_CLOSE(nout)


!==============================================================================
! 4. Build a companion file 'lslin.ctl', so that output file can be 
!    processed by Grads
!==============================================================================

!  ----------------------------------------------------
!  Writing ctl file to directly read Ls coordinate in grads
!  (because of bug in grads that refuse to read date like 0089 in .nc files)
write(*,*)""
write(*,*) "Writing a controle file to directly read Ls coordinate with Grads..."

 open(33,file=infile(1:len_trim(infile)-3)//"_Ls.ctl")
 date= (timels(1)-int(timels(1)))*365.
 mon='jan'
 if(date.ge.32) mon='feb'
 if(date.ge.60) mon='mar'
 if(date.ge.91) mon='apr'
 if(date.ge.121) mon='may'
 if(date.ge.152) mon='jun'
 if(date.ge.182) mon='jul'
 if(date.ge.213) mon='aug'
 if(date.ge.244) mon='sep'
 if(date.ge.274) mon='oct'
 if(date.ge.305) mon='nov'
 if(date.ge.335) mon='dec'
write(33,98) "^"//outfile
98 format("DSET ",a)
  write(33,99) Nls, 15,mon, int(timels(1)),nint(deltals*12),'mo'
99 format("TDEF Time ",i5," LINEAR ", i2,a3,i4.4,1x,i2,a2)
    
deallocate(timels)
   
write(*,*)""   
write(*,*) "Program completed !"

contains

!******************************************************************************

subroutine initiate(filename,lat,lon,alt,ctl,latlen,lonlen,altlen,ctllen,&
                     altlong_name,altunits,altpositive,&
                     nout,latdimout,londimout,altdimout,timedimout,timevarout)
!==============================================================================
! Purpose:
! Create and initialize a data file (NetCDF format)
!==============================================================================
! Remarks:
! The NetCDF file (created in this subroutine) remains open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
character (len=*), intent(in):: filename
! filename(): the file's name
integer,intent(in) ::latlen,lonlen,altlen,ctllen
real, intent(in):: lat(latlen)
! lat(): latitude
real, intent(in):: lon(lonlen)
! lon(): longitude
real, intent(in):: alt(altlen)
! alt(): altitude
real, intent(in):: ctl(ctllen)
! ctl(): controle
character (len=*), intent(in) :: altlong_name
! altlong_name(): [netcdf] altitude "long_name" attribute
character (len=*), intent(in) :: altunits
! altunits(): [netcdf] altitude "units" attribute
character (len=*), intent(in) :: altpositive
! altpositive(): [netcdf] altitude "positive" attribute
integer, intent(out):: nout
! nout: [netcdf] file ID
integer, intent(out):: latdimout
! latdimout: [netcdf] lat() (i.e.: latitude)  ID
integer, intent(out):: londimout
! londimout: [netcdf] lon()  ID
integer, intent(out):: altdimout
! altdimout: [netcdf] alt()  ID
integer, intent(out):: timedimout
! timedimout: [netcdf] "Time"  ID
integer, intent(out):: timevarout
! timevarout: [netcdf] "Time" (considered as a variable) ID

!==============================================================================
! Local variables:
!==============================================================================
!integer :: latdim,londim,altdim,timedim
integer :: nvarid,ierr
integer :: ctldimout
! nvarid: [netcdf] ID of a variable
! ierr: [netcdf]  return error code (from called subroutines)

!==============================================================================
! 1. Create (and open) output file
!==============================================================================
write(*,*) "creating "//trim(adjustl(filename))//'...'
ierr = NF_CREATE(filename,IOR(NF_CLOBBER,NF_64BIT_OFFSET),nout)
! NB: setting NF_CLOBBER mode means that it's OK to overwrite an existing file
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'ERROR: Impossible to create the file ',trim(filename)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

!==============================================================================
! 2. Define/write "dimensions" and get their IDs
!==============================================================================

ierr = NF_DEF_DIM(nout, "latitude", latlen, latdimout)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed to define dimension <latitude>'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif
ierr = NF_DEF_DIM(nout, "longitude", lonlen, londimout)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed to define dimension <longitude>'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif
ierr = NF_DEF_DIM(nout, "altitude", altlen, altdimout)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed to define dimension <altitude>'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif
if (size(ctl).ne.0) then
  ierr = NF_DEF_DIM(nout, "index", ctllen, ctldimout)
  if (ierr.NE.NF_NOERR) then
    WRITE(*,*)'initiate: error failed to define dimension <index>'
    write(*,*) NF_STRERROR(ierr)
    stop ""
  endif
endif
ierr = NF_DEF_DIM(nout, "Time", NF_UNLIMITED, timedimout)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed to define dimension <Time>'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! End netcdf define mode
ierr = NF_ENDDEF(nout)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed to switch out of define mode'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

!==============================================================================
! 3. Write "Time" and its attributes
!==============================================================================

call def_var(nout,"Time","Ls (Solar Longitude)","degrees",1,&
             (/timedimout/),timevarout,ierr)
! Switch to netcdf define mode
ierr = NF_REDEF (nout)
ierr = NF_ENDDEF(nout)

!==============================================================================
! 4. Write "latitude" (data and attributes)
!==============================================================================

call def_var(nout,"latitude","latitude","degrees_north",1,&
             (/latdimout/),nvarid,ierr)

ierr = NF_PUT_VAR_REAL (nout,nvarid,lat)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed writing variable <latitude>'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

!==============================================================================
! 4. Write "longitude" (data and attributes)
!==============================================================================

call def_var(nout,"longitude","East longitude","degrees_east",1,&
             (/londimout/),nvarid,ierr)

ierr = NF_PUT_VAR_REAL (nout,nvarid,lon)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed writing variable <longitude>'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

!==============================================================================
! 5. Write "altitude" (data and attributes)
!==============================================================================

! Switch to netcdf define mode
ierr = NF_REDEF (nout)

ierr = NF_DEF_VAR (nout,"altitude",NF_FLOAT,1,altdimout,nvarid)

ierr = NF_PUT_ATT_TEXT (nout,nvarid,'long_name',len_trim(adjustl(altlong_name)),adjustl(altlong_name))
ierr = NF_PUT_ATT_TEXT (nout,nvarid,'units',len_trim(adjustl(altunits)),adjustl(altunits))
ierr = NF_PUT_ATT_TEXT (nout,nvarid,'positive',len_trim(adjustl(altpositive)),adjustl(altpositive))

! End netcdf define mode
ierr = NF_ENDDEF(nout)

ierr = NF_PUT_VAR_REAL (nout,nvarid,alt)
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'initiate: error failed writing variable <altitude>'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

!==============================================================================
! 6. Write "controle" (data and attributes)
!==============================================================================

if (size(ctl).ne.0) then
   ! Switch to netcdf define mode
   ierr = NF_REDEF (nout)

   ierr = NF_DEF_VAR (nout,"controle",NF_FLOAT,1,ctldimout,nvarid)

   ierr = NF_PUT_ATT_TEXT (nout,nvarid,"long_name",18,"Control parameters")

   ! End netcdf define mode
   ierr = NF_ENDDEF(nout)

   ierr = NF_PUT_VAR_REAL (nout,nvarid,ctl)
   if (ierr.NE.NF_NOERR) then
      WRITE(*,*)'initiate: error failed writing variable <controle>'
      write(*,*) NF_STRERROR(ierr)
      stop ""
   endif
endif

end Subroutine initiate

!******************************************************************************
subroutine init2(infid,lonlen,latlen,altlen,altdimid, &
                 outfid,londimout,latdimout,altdimout)
!==============================================================================
! Purpose:
! Copy aps(), bps() and phisinit() from input file to output file
!==============================================================================
! Remarks:
! The NetCDF files must be open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: infid  ! NetCDF output file ID
integer, intent(in) :: lonlen ! # of grid points along longitude
integer, intent(in) :: latlen ! # of grid points along latitude
integer, intent(in) :: altlen ! # of grid points along altitude
integer, intent(in) :: altdimid ! ID of altitude dimension
integer, intent(in) :: outfid ! NetCDF output file ID
integer, intent(in) :: londimout ! longitude dimension ID
integer, intent(in) :: latdimout ! latitude dimension ID
integer, intent(in) :: altdimout ! altitude dimension ID
!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:,:),allocatable :: phisinit ! Ground geopotential
integer :: apsid,bpsid,phisinitid
integer :: ierr
integer :: tmpdimid ! temporary dimension ID
integer :: tmpvarid ! temporary variable ID
integer :: tmplen ! temporary variable length
logical :: phis, aps_ok, bps_ok ! are "phisinit" "aps" "bps" available ?


!==============================================================================
! 1. Read data from input file
!==============================================================================

! hybrid coordinate aps
  allocate(aps(altlen),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "init2: failed to allocate aps(altlen)"
    stop
  endif

ierr=NF_INQ_VARID(infid,"aps",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "oops Failed to get aps ID. OK"
  aps_ok=.false.
else
  ! Check that aps() is the right size (it most likely won't be if
  ! zrecast has be used to generate the input file)
  ierr=NF_INQ_VARDIMID(infid,tmpvarid,tmpdimid)
  ierr=NF_INQ_DIMLEN(infid,tmpdimid,tmplen)
  if (tmplen.ne.altlen) then
    write(*,*) "tmplen,altlen", tmpdimid, altdimid
    write(*,*) "init2: wrong dimension size for aps(), skipping it ..."
    aps_ok=.false.
  else
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,aps)
    if (ierr.ne.NF_NOERR) then
     stop "init2 error: Failed reading aps"
     aps_ok=.false.
    endif
    aps_ok=.true.
  endif
endif

! hybrid coordinate bps
  allocate(bps(altlen),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "init2: failed to allocate bps(altlen)"
    stop
  endif

ierr=NF_INQ_VARID(infid,"bps",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "oops: Failed to get bps ID. OK"
  bps_ok=.false.
else
  ! Check that bps() is the right size
  ierr=NF_INQ_VARDIMID(infid,tmpvarid,tmpdimid)
  ierr=NF_INQ_DIMLEN(infid,tmpdimid,tmplen)
  if (tmplen.ne.altlen) then
    write(*,*) "init2: wrong dimension size for bps(), skipping it ..."
    bps_ok=.false.
  else
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,bps)
    if (ierr.ne.NF_NOERR) then
      stop "init2 Error: Failed reading bps"
      bps_ok=.false.
    endif
    bps_ok=.true.
  endif
endif

! ground geopotential phisinit
allocate(phisinit(lonlen,latlen),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate phisinit(lonlen,latlen)"
  stop
endif
ierr=NF_INQ_VARID(infid,"phisinit",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Failed to get phisinit ID. OK"
  phisinit = 0.
  phis = .false.
else
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,phisinit)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed reading phisinit"
  endif
  phis = .true.
endif


!==============================================================================
! 2. Write
!==============================================================================

!==============================================================================
! 2.2. Hybrid coordinates aps() and bps()
!==============================================================================

IF (aps_ok) then 
  ! define aps
  call def_var(outfid,"aps","hybrid pressure at midlayers"," ",1,&
               (/altdimout/),apsid,ierr)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to def_var aps"
  endif

  ! write aps
  ierr=NF_PUT_VAR_REAL(outfid,apsid,aps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to write aps"
  endif
ENDIF ! of IF (aps_ok)

IF (bps_ok) then 
  ! define bps
  call def_var(outfid,"bps","hybrid sigma at midlayers"," ",1,&
               (/altdimout/),bpsid,ierr)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to def_var bps"
  endif

  ! write bps
  ierr=NF_PUT_VAR_REAL(outfid,bpsid,bps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to write bps"
  endif
ENDIF ! of IF (bps_ok)

!==============================================================================
! 2.2. phisinit()
!==============================================================================

IF (phis) THEN

  !define phisinit
  call def_var(outfid,"phisinit","Ground level geopotential"," ",2,&
              (/londimout,latdimout/),phisinitid,ierr)
  if (ierr.ne.NF_NOERR) then
     stop "Error: Failed to def_var phisinit"
  endif

  ! write phisinit
  ierr=NF_PUT_VAR_REAL(outfid,phisinitid,phisinit)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to write phisinit"
  endif

ENDIF ! of IF (phis)


! Cleanup
deallocate(aps)
deallocate(bps)
deallocate(phisinit)

end subroutine init2

!******************************************************************************
subroutine def_var(nid,name,long_name,units,nbdim,dim,nvarid,ierr)
!==============================================================================
! Purpose: Write a variable (i.e: add a variable to a dataset)
! called "name"; along with its attributes "long_name", "units"...
! to a file (following the NetCDF format)
!==============================================================================
! Remarks:
! The NetCDF file must be open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: nid
! nid: [netcdf] file ID #
character (len=*), intent(in) :: name
! name(): [netcdf] variable's name
character (len=*), intent(in) :: long_name
! long_name(): [netcdf] variable's "long_name" attribute
character (len=*), intent(in) :: units
! unit(): [netcdf] variable's "units" attribute
integer, intent(in) :: nbdim
! nbdim: number of dimensions of the variable
integer, dimension(nbdim), intent(in) :: dim
! dim(nbdim): [netcdf] dimension(s) ID(s)
integer, intent(out) :: nvarid
! nvarid: [netcdf] ID # of the variable
integer, intent(out) :: ierr
! ierr: [netcdf] subroutines returned error code

! Switch to netcdf define mode
ierr=NF_REDEF(nid)

! Insert the definition of the variable
ierr=NF_DEF_VAR(nid,adjustl(name),NF_FLOAT,nbdim,dim,nvarid)

! Write the attributes
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"long_name",len_trim(adjustl(long_name)),adjustl(long_name))
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"units",len_trim(adjustl(units)),adjustl(units))

! End netcdf define mode
ierr=NF_ENDDEF(nid)

end subroutine def_var

!******************************************************************************
subroutine  missing_value(nout,nvarid,missing)
!==============================================================================
! Purpose:
! Write "valid_range" and "missing_value" attributes (of a given
! variable) to a netcdf file
!==============================================================================
! Remarks:
! NetCDF file must be open
! Variable (nvarid) ID must be set
!==============================================================================

implicit none

include "netcdf.inc"  ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: nout
! nout: [netcdf] file ID #
integer, intent(in) :: nvarid
! varid: [netcdf] variable ID #
!real, dimension(2), intent(in) :: valid_range
! valid_range(2): [netcdf] "valid_range" attribute (min and max)
real, intent(in) :: missing
! missing: [netcdf] "missing_value" attribute

!==============================================================================
! Local variables:
!==============================================================================
integer :: ierr
! ierr: [netcdf] subroutine returned error code
!      INTEGER nout,nvarid,ierr
!      REAL missing
!      REAL valid_range(2)

! Switch to netcdf dataset define mode
ierr = NF_REDEF (nout)
if (ierr.ne.NF_NOERR) then
   print*,'missing_value: '
   print*, NF_STRERROR(ierr)
endif


!********* valid range not used in Lslin ****************
! Write valid_range() attribute
!ierr=NF_PUT_ATT_REAL(nout,nvarid,'valid_range',NF_FLOAT,2,valid_range)
!
!if (ierr.ne.NF_NOERR) then
!   print*,'missing_value: valid_range attribution failed'
!   print*, NF_STRERROR(ierr)
!   stop ""
!endif
!*********************************************************

! Write "missing_value" attribute
ierr= NF_PUT_ATT_REAL(nout,nvarid,'missing_value',NF_FLOAT,1,missing)
if (ierr.NE.NF_NOERR) then
   print*, 'missing_value: missing value attribution failed'
   print*, NF_STRERROR(ierr)
   stop ""
endif

! End netcdf dataset define mode
ierr = NF_ENDDEF(nout)

end subroutine  missing_value

!******************************************************************************
subroutine sol2ls(sol,Ls)
!==============================================================================
! Purpose: 
! Convert a date/time, given in sol (martian day),
! into solar longitude date/time, in Ls (in degrees),
! where sol=0 is (by definition) the northern hemisphere
!  spring equinox (where Ls=0).
!==============================================================================
! Notes:
! Even though "Ls" is cyclic, if "sol" is greater than N (martian) year,
! "Ls" will be increased by N*360
! Won't work as expected if sol is negative (then again,
! why would that ever happen?)
!==============================================================================

implicit none

!==============================================================================
! Arguments:
!==============================================================================
real,intent(in) :: sol
real,intent(out) :: Ls

!==============================================================================
! Local variables:
!==============================================================================
real year_day,peri_day,timeperi,e_elips,twopi,degrad
data year_day /669./            ! # of sols in a martian year
data peri_day /485.0/           
data timeperi /1.9082314/! True anomaly at vernal equinox = 2*PI-Lsperi
data e_elips  /0.093358/
data twopi       /6.2831853/    ! 2.*pi
data degrad   /57.2957795/      ! pi/180

real zanom,xref,zx0,zdx,zteta,zz

integer count_years
integer iter

!==============================================================================
! 1. Compute Ls
!==============================================================================

zz=(sol-peri_day)/year_day
zanom=twopi*(zz-nint(zz))
xref=abs(zanom)

!  The equation zx0 - e * sin (zx0) = xref, solved by Newton
zx0=xref+e_elips*sin(xref)
do iter=1,20 ! typically, 2 or 3 iterations are enough
   zdx=-(zx0-e_elips*sin(zx0)-xref)/(1.-e_elips*cos(zx0))
   zx0=zx0+zdx
   if(abs(zdx).le.(1.e-7)) then
!      write(*,*)'iter:',iter,'     |zdx|:',abs(zdx)
      exit
   endif
enddo

if(zanom.lt.0.) zx0=-zx0

zteta=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
Ls=zteta-timeperi

if(Ls.lt.0.) then
   Ls=Ls+twopi
else
   if(Ls.gt.twopi) then
      Ls=Ls-twopi
   endif
endif

Ls=degrad*Ls
! Ls is now in degrees

!==============================================================================
! 1. Account for (eventual) years included in input date/time sol
!==============================================================================

count_years=0 ! initialize
zz=sol  ! use "zz" to store (and work on) the value of sol
do while (zz.ge.year_day)
   count_years=count_years+1
   zz=zz-year_day
enddo

! Add 360 degrees to Ls for every year
if (count_years.ne.0) then
   Ls=Ls+360.*count_years
endif


end subroutine sol2ls

!******************************************************************************
subroutine interpolf(x,y,missing,xd,yd,nd)
!==============================================================================
! Purpose:
! Yield y=f(x), where xd() end yd() are arrays of known values,
! using linear interpolation
! If x is not included in the interval spaned by xd(), then y is set
! to a default value 'missing'
! Note:
! Array xd() should contain ordered (either increasing or decreasing) abscissas
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
real,intent(in) :: x ! abscissa at which interpolation is to be done
real,intent(in) :: missing ! missing value (if no interpolation is performed)
integer :: nd ! size of arrays
real,dimension(nd),intent(in) :: xd ! array of known absissas
real,dimension(nd),intent(in) :: yd ! array of correponding values

real,intent(out) :: y ! interpolated value
!==============================================================================
! Local variables:
!==============================================================================
integer :: i

! default: set y to 'missing'
y=missing

   do i=1,nd-1
     if (((x.ge.xd(i)).and.(x.le.xd(i+1))).or.&
          ((x.le.xd(i)).and.(x.ge.xd(i+1)))) then
        if ((yd(i).eq.missing).or.(yd(i+1).eq.missing)) then
          ! cannot perform the interpolation if an encompasing value
          ! is already set to 'missing'
        else
          !linear interpolation based on encompassing values
          y=yd(i)+(x-xd(i))*(yd(i+1)-yd(i))/(xd(i+1)-xd(i))
        endif
        exit
     endif
   enddo


end subroutine interpolf

!******************************************************************************

end program lslin
