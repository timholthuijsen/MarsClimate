program localtime

! ----------------------------------------------------------------------------
! Program to redistribute and interpolate the variable a the same
! local times everywhere
! input : diagfi.nc  / concat.nc / stats.nc kind of files
! author: F. Forget
! ----------------------------------------------------------------------------

implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=50)  file
! file(): input file(s) names(s)
character (len=30), dimension(16) :: notprocessed
! notprocessed(): names of the (16) variables that won't be processed
character (len=50), dimension(:), allocatable :: var
! var(): name(s) of variable(s) that will be processed
character (len=100) :: tmpvar,long_name,units
! tmpvar(): used to temporarily store a variable name
! long_name(): [netcdf] long_name attribute
! units(): [netcdf] units attribute
character (len=100) :: filename,vartmp
! filename(): output file name
! vartmp(): temporary variable name (read from netcdf input file)
!character (len=1) :: ccopy
! ccpy: 'y' or 'n' answer

character (len=50) :: altlong_name,altunits,altpositive
! altlong_name(): [netcdf] altitude "long_name" attribute
! altunits(): [netcdf] altitude "units" attribute
! altpositive(): [netcdf] altitude "positive" attribute

integer :: nid,ierr,miss,validr
! nid: [netcdf] file ID #
! ierr: [netcdf] subroutine returned error code
! miss: [netcdf] subroutine returned error code
integer :: i,j,k,inter
! for various loops
integer :: varid
! varid: [netcdf] variable ID #
real, dimension(:), allocatable:: lat,lon,alt,ctl,time
! lat(): array, stores latitude coordinates
! lon(): array, stores longitude coordinates
! alt(): array, stores altitude coordinates
! ctl(): array, stores controle array
! time(): array, stores time coordinates
integer :: nbvar,nbvarfile,ndim
! nbvar: # of variables to concatenate
! nbfile: # number of input file(s)
! nbvarfile: total # of variables in an input file
! ndim: [netcdf] # (3 or 4) of dimensions (for variables)
integer :: latdim,londim,altdim,ctldim,timedim
! latdim: [netcdf] "latitude" dim ID
! londim: [netcdf] "longitude" dim ID
! altdim: [netcdf] "altdim" dim ID
! ctldim: [netcdf] "controle" dim ID
! timedim: [netcdf] "timedim" dim ID
integer :: gcmlayerdim ! NetCDF dimension ID for # of layers in GCM
integer :: latvar,lonvar,altvar,ctlvar,timevar
! latvar: [netcdf] ID of "latitude" variable
! lonvar: [netcdf] ID of "longitude" variable
! altvar: [netcdf] ID of "altitude" variable
! ctlvar: [netcdf] ID of "controle" variable
! timevar: [netcdf] ID of "Time" variable
integer :: latlen,lonlen,altlen,ctllen,timelen,timelen_lt,timelen_tot
integer :: ilat,ilon,ialt,it
! latlen: # of elements of lat() array
! lonlen: # of elements of lon() array
! altlen: # of elements of alt() array
! ctllen: # of elements of ctl() array
! timelen: # of elemnets of time() array
! timelen_tot:# =timelen or timelen+1 (if 1 more time to interpolate needed)
! timelen_lt: # of elemnets of time() array in output
integer :: nhour,nsol
integer :: GCM_layers ! number of GCM atmospheric layers (may not be
! same as altlen if file is an output of zrecast)
integer :: nout,latdimout,londimout,altdimout,timedimout,timevarout
! nout: [netcdf] output file ID
! latdimout: [netcdf] output latitude (dimension) ID
! londimout: [netcdf] output longitude (dimension) ID
! altdimout: [netcdf] output altitude (dimension) ID
! timedimout: [netcdf] output time (dimension) ID
! timevarout: [netcdf] ID of output "Time" variable
integer :: layerdimout ! NetCDF dimension ID for # of layers in GCM
integer :: varidout 
! varidout: [netcdf] variable ID # (of a variable to write to the output file)
integer :: Nnotprocessed,var_ok
! Nnotprocessed: # of (leading)variables that won't be concatenated
! var_ok: flag (0 or 1)
integer, dimension(4) :: corner,edges,dim
! corner: [netcdf]
! edges: [netcdf]
! dim: [netcdf]
real, dimension(:,:,:,:), allocatable :: var3d
! var3D(,,,): 4D array to store a field

real, dimension(:,:,:,:), allocatable :: var3d_lt
! var3D_lt(,,,): 4D array to store a field in local time coordinate
real, dimension(:), allocatable :: lt_gcm,lt_hour
real, dimension(:), allocatable :: lt_out,lt_outc
real, dimension(:), allocatable :: var_gcm

real :: missing
!PARAMETER(missing=1E+20)
! missing: [netcdf] to handle "missing" values when reading/writing files
real, dimension(2) :: valid_range
! valid_range(2): [netcdf] interval in which a value is considered valid
logical :: stats   ! stats=T when reading a "stats" kind of ile


!==============================================================================
! 1.1. Get input file name(s)
!==============================================================================
write(*,*) 
write(*,*) "which file do you want to use?  (diagfi... stats...  concat...)"

read(*,'(a50)') file

if (len_trim(file).eq.0) then
   write(*,*) "no file... game over"
   stop ""
endif

!==============================================================================
! 1.3. Open the input file
!==============================================================================

ierr = NF_OPEN(file,NF_NOWRITE,nid)
if (ierr.NE.NF_NOERR) then
   write(*,*) 'ERROR: Pb opening file '//trim(file)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

ierr=NF_INQ_NVARS(nid,nbvarfile)
! nbvarfile now set to be the (total) number of variables in file
if (ierr.NE.NF_NOERR) then
   write(*,*) 'ERROR: Pb with NF_INQ_NVARS'
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

!==============================================================================
! 1.4. List of variables that should not be processed
!==============================================================================

notprocessed(1)='Time'
notprocessed(2)='controle'
notprocessed(3)='rlonu'
notprocessed(4)='latitude'
notprocessed(5)='longitude'
notprocessed(6)='altitude'
notprocessed(7)='rlatv'
notprocessed(8)='aps'
notprocessed(9)='bps'
notprocessed(10)='ap'
notprocessed(11)='bp'
notprocessed(12)='soildepth'
notprocessed(13)='cu'
notprocessed(14)='cv'
notprocessed(15)='aire'
notprocessed(16)='phisinit'


!==============================================================================
! 1.5. Get (and check) list of variables to process
!==============================================================================
write(*,*)
   Nnotprocessed=0
do i=1,nbvarfile
   ierr=NF_INQ_VARNAME(nid,i,vartmp)
   ! vartmp now contains the "name" of variable of ID # i
   var_ok=0
   do inter=1,16
      if (vartmp.eq.notprocessed(inter)) then
         var_ok=1 
         Nnotprocessed=Nnotprocessed+1
      endif 
   enddo        
   if (var_ok.eq.0)  write(*,*) trim(vartmp)
enddo

! Nnotprocessed: # of variables that won't be processed
! nbvarfile: total # of variables in file
allocate(var(nbvarfile-Nnotprocessed),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error: failed allocation of var(nbvarfile-Nnotprocessed)"
  write(*,*) "  nbvarfile=",nbvarfile
  write(*,*) "  Nnotprocessed=",Nnotprocessed
  stop
endif


write(*,*)
write(*,*) "which variables do you want to redistribute ?"
write(*,*) "all / list of <variables> (separated by <Enter>s)"
write(*,*) "(an empty line , i.e: just <Enter>, implies end of list)"
nbvar=0
read(*,'(a50)') tmpvar
do while ((tmpvar/=' ').AND.(trim(tmpvar)/='all'))
   nbvar=nbvar+1
   var(nbvar)=tmpvar
   read(*,'(a50)') tmpvar
enddo

if (tmpvar=="all") then
         nbvar=nbvarfile-Nnotprocessed
         do j=Nnotprocessed+1,nbvarfile
            ierr=nf_inq_varname(nid,j,var(j-Nnotprocessed))
         enddo
! Variables names from the file are stored in var()
   nbvar=nbvarfile-Nnotprocessed
   do i=1,nbvar
      ierr=nf_inq_varname(nid,i+Nnotprocessed,var(i))
      write(*,'(a9,1x,i2,1x,a1,1x,a50)') "variable ",i,":",var(i)
   enddo
else if(nbvar==0) then
   write(*,*) "no variable... game over"
   stop ""
endif ! of if (tmpvar=="all")

!==============================================================================
! 1.6. Get output file name
!==============================================================================
 filename=file(1:len_trim(file)-3)//"_LT.nc"
 write(*,*) "Output file :"//trim(filename)

!==============================================================================
! 2.1. Open input file
!==============================================================================

   write(*,*) 
   write(*,*) "opening "//trim(file)//"..."
   ierr = NF_OPEN(file,NF_NOWRITE,nid)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Pb opening file '//trim(file)
      write(*,*) NF_STRERROR(ierr)
      stop ""
   endif

!==============================================================================
! 2.2. Read (and check) dimensions of variables from input file
!==============================================================================

   ierr=NF_INQ_DIMID(nid,"latitude",latdim)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Dimension <latitude> is missing in file '//trim(file)
      stop ""  
   endif
   ierr=NF_INQ_VARID(nid,"latitude",latvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <latitude> is missing in file '//trim(file)
      stop ""  
   endif
   ierr=NF_INQ_DIMLEN(nid,latdim,latlen)
!  write(*,*) "latlen: ",latlen

   ierr=NF_INQ_DIMID(nid,"longitude",londim)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Dimension <longitude> is missing in file '//trim(file)
      stop ""  
   endif
   ierr=NF_INQ_VARID(nid,"longitude",lonvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <longitude> is missing in file '//trim(file)
      stop "" 
   endif
   ierr=NF_INQ_DIMLEN(nid,londim,lonlen)
!  write(*,*) "lonlen: ",lonlen

   ierr=NF_INQ_DIMID(nid,"altitude",altdim)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Dimension <altitude> is missing in file '//trim(file)
      stop ""  
   endif
   ierr=NF_INQ_VARID(nid,"altitude",altvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <altitude> is missing in file '//trim(file)
      stop ""
   endif
   ierr=NF_INQ_DIMLEN(nid,altdim,altlen)
!  write(*,*) "altlen: ",altlen

! load size of aps() or sigma() (in case it is not altlen)
   ! default is that GCM_layers=altlen
   ! but for outputs of zrecast, it may be a different value
   ierr=NF_INQ_DIMID(nid,"GCM_layers",gcmlayerdim)
   if (ierr.ne.NF_NOERR) then
     ! didn't find a GCM_layers dimension; therefore we have:
     GCM_layers=altlen
   else
     ! load value of GCM_layers
     ierr=NF_INQ_DIMLEN(nid,gcmlayerdim,GCM_layers)
   endif
!   write(*,*) "GCM_layers=",GCM_layers

   ierr=NF_INQ_DIMID(nid,"index",ctldim)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'Dimension <index> is missing in file '//trim(file)
      ctllen=0
      !stop ""  
   endif
   ierr=NF_INQ_VARID(nid,"controle",ctlvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'Field <controle> is missing in file '//trim(file)
      ctllen=0
      !stop ""
   else
      ierr=NF_INQ_DIMLEN(nid,ctldim,ctllen)
   endif
!  write(*,*) "controle: ",controle

!==============================================================================
! 2.3. Read (and check compatibility of) dimensions of
!       variables from input file
!==============================================================================

   allocate(lat(latlen),stat=ierr)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to allocate lat(latlen)"
     stop
   endif
   allocate(lon(lonlen),stat=ierr)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to allocate lon(lonlen)"
     stop
   endif
   allocate(alt(altlen),stat=ierr)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to allocate alt(altlen)"
     stop
   endif
   allocate(ctl(ctllen),stat=ierr)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to allocate ctl(ctllen)"
     stop
   endif

   ierr = NF_GET_VAR_REAL(nid,latvar,lat)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to load latitude"
     write(*,*) NF_STRERROR(ierr)
     stop
   endif

   ierr = NF_GET_VAR_REAL(nid,lonvar,lon)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to load longitude"
     write(*,*) NF_STRERROR(ierr)
     stop
   endif

   ierr = NF_GET_VAR_REAL(nid,altvar,alt)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to load altitude"
     write(*,*) NF_STRERROR(ierr)
     stop
   endif
! Get altitude attributes to handle files with any altitude type
   ierr=nf_get_att_text(nid,altvar,'long_name',altlong_name)
   if (ierr.ne.nf_noerr) then
   ! if no attribute "long_name", try "title"
     ierr=nf_get_att_text(nid,timevar,"title",long_name)
   endif
   ierr=nf_get_att_text(nid,altvar,'units',altunits)
   ierr=nf_get_att_text(nid,altvar,'positive',altpositive)

   if (ctllen .ne. 0) then
     ierr = NF_GET_VAR_REAL(nid,ctlvar,ctl)
     if (ierr.ne.0) then
       write(*,*) "Error: failed to load controle"
       write(*,*) NF_STRERROR(ierr)
       stop
     endif
   endif ! of if (ctllen .ne. 0)

!==============================================================================
! 2.4. Handle "Time" dimension from input file
!==============================================================================

!==============================================================================
! 2.4.0 Read "Time" dimension from input file
!==============================================================================
   ierr=NF_INQ_DIMID(nid,"Time",timedim)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Dimension <Time> is missing in file'//trim(file)
      stop ""
   endif
   ierr=NF_INQ_VARID(nid,"Time",timevar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <Time> is missing in file'//trim(file)
      stop ""
   endif
   ierr=NF_INQ_DIMLEN(nid,timedim,timelen)
!  write(*,*) "timelen: ",timelen

   ! Sanity Check: Does "Time" has a "long_name" attribute
   ! and is it "Solar longitude" ?
   ierr=nf_get_att_text(nid,timevar,"long_name",long_name)
   if (ierr.ne.nf_noerr) then
   ! if no attribute "long_name", try "title"
     ierr=nf_get_att_text(nid,timevar,"title",long_name)
   endif
   if ((ierr.EQ.NF_NOERR).and.(index(long_name,"Solar").ne.0)) then
     ! if long_name contains "Solar"
     write(*,*) "ERROR: Time axis in input file is in Solar Longitude!"
     write(*,*) "       localtime requires sols as time axis!"
     write(*,*) " Might as well stop here."
     stop
   endif

   ! allocate time() array and fill it with values from input file
   allocate(time(timelen),stat=ierr)
   if (ierr.ne.0) then
     write(*,*) "Error: failed to allocate time(timelen)"
     stop
   endif

   ierr = NF_GET_VAR_REAL(nid,timevar,time)
   if (ierr.NE.NF_NOERR) then
     write(*,*) "Error , failed to load Time"
     write(*,*) NF_STRERROR(ierr)
     stop
   endif


!==============================================================================
! 2.4.1 Select local times to be saved "Time" in output file
!==============================================================================
      write(*,*) 'number of local time to be used ?'
      read(*,*) nhour
      allocate(lt_hour(nhour),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: failed to allocate lt_hour(nhour)"
        stop
      endif
      write(*,*) 'list of selected local time (0<t<24)'
      do it=1,nhour
        read(*,*) lt_hour(it)
      end do

      if ((timelen.eq.12).and.(time(1).eq.2).and.(time(timelen).eq.24))then
        write(*,*) 'We detect a ""stats"" file'
        stats=.true.
        timelen_lt=nhour
        allocate(lt_out(timelen_lt),stat=ierr)
        if (ierr.ne.0) then
          write(*,*) "Error: failed to allocate lt_hour(nhour)"
          stop
        endif
        do it=1,nhour
          lt_out(it)=lt_hour(it)/24.
        end do
!        We rewrite the time from "stats" from 0 to 1 sol...
        do it=1,timelen  ! loop temps in
              time(it) = time(it)/24.
        end do
        nsol =1
      else   ! case of a diagfi or concatnc file 
        stats=.false.
!       Number of sol in input file
        nsol = int(time(timelen)) - int(time(1))
        timelen_lt=nhour*nsol
        allocate(lt_out(timelen_lt),stat=ierr)
        if (ierr.ne.0) then
           write(*,*) "Error: failed to allocate lt_hour(nhour)"
           stop
         endif
        i=0
        do k=1,nsol
          do it=1,nhour
            i=i+1
            lt_out(i)=int(time(1)) + k-1  + lt_hour(it)/24.
          end do
        end do
      end if

      if (nsol.gt.1) then
         timelen_tot=timelen
      else
!        if only 1 sol available, we must add 1 timestep for the interpolation
         timelen_tot=timelen+1
      endif      
      allocate(lt_gcm(timelen_tot),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: failed to allocate lt_gcm(timelen_tot)"
        stop
      endif
      allocate(var_gcm(timelen_tot),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: failed to allocate var_gcm(timelen_tot)"
        stop
      endif
      allocate(lt_outc(timelen_lt),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: failed to allocate lt_outc(timelen_lt)"
        stop
      endif

!==============================================================================
! 2.4.1.5 Initiate dimension in output file
!==============================================================================


   ! Initialize output file's lat,lon,alt and time dimensions
     call initiate(filename,lat,lon,alt,ctl,latlen,lonlen,altlen,ctllen,GCM_layers,&
           altlong_name,altunits,altpositive,&
           nout,latdimout,londimout,altdimout,timedimout,&
           layerdimout,timevarout)

   ! Initialize output file's aps,bps and phisinit variables
     call init2(nid,lonlen,latlen,altlen,GCM_layers,&
                nout,londimout,latdimout,altdimout,&
                layerdimout)


!==============================================================================
! 2.4.2 Write/extend "Time" dimension/values in output file
!==============================================================================

   if (stats) then 

     do it=1,timelen_lt
        ierr= NF_PUT_VARA_REAL(nout,timevarout,it,1,lt_hour(it))
     enddo
   else
     do it=1,timelen_lt
        ierr= NF_PUT_VARA_REAL(nout,timevarout,it,1,lt_out(it))
        if (ierr.NE.NF_NOERR) then
          write(*,*) "Error , failed to write Time"
          stop
        endif
     enddo
  end if

!==============================================================================
! 2.5 Read/write variables
!==============================================================================

   do j=1,nbvar ! loop on variables to read/write

!==============================================================================
! 2.5.1 Check that variable to be read existe in input file
!==============================================================================

      write(*,*) "variable ",var(j)
      ierr=nf_inq_varid(nid,var(j),varid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Field <',var(j),'> not found in file'//file
         stop ""
      endif
      ierr=nf_inq_varndims(nid,varid,ndim)

      ! Check that it is a 3D or 4D variable
      if (ndim.lt.3) then
        write(*,*) "Error:",trim(var(j))," is nor a 3D neither a 4D variable"
        write(*,*) "We'll skip that variable..."
        CYCLE ! go directly to the next variable
      endif
!==============================================================================
! 2.5.2 Prepare things in order to read/write the variable
!==============================================================================

      ! build dim(),corner() and edges() arrays
      ! and allocate var3d() array
      if (ndim==3) then
         allocate(var3d(lonlen,latlen,timelen,1))
         allocate(var3d_lt(lonlen,latlen,timelen_lt,1))
         dim(1)=londimout
         dim(2)=latdimout
         dim(3)=timedimout

         ! start indexes (where data values will be written)
         corner(1)=1
         corner(2)=1
         corner(3)=1
         corner(4)=1

         ! length (along dimensions) of block of data to be written
         edges(1)=lonlen 
         edges(2)=latlen 
         edges(3)=timelen_lt
         edges(4)=1
   
      else if (ndim==4) then
         allocate(var3d(lonlen,latlen,altlen,timelen))
         allocate(var3d_lt(lonlen,latlen,altlen,timelen_lt))
         dim(1)=londimout
         dim(2)=latdimout
         dim(3)=altdimout
         dim(4)=timedimout

         ! start indexes (where data values will be written)
         corner(1)=1
         corner(2)=1
         corner(3)=1
         corner(4)=1

         ! length (along dimensions) of block of data to be written
         edges(1)=lonlen
         edges(2)=latlen
         edges(3)=altlen
         edges(4)=timelen_lt
      endif

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
! 2.5.3 Read from input file 
!==============================================================================

      ierr = NF_GET_VAR_REAL(nid,varid,var3d)
      miss=NF_GET_ATT_REAL(nid,varid,"missing_value",missing)
      validr=NF_GET_ATT_REAL(nid,varid,"valid_range",valid_range)

!==============================================================================
! 2.5.4 Read from input file 
!==============================================================================
! Interpolation in local time :

        do ilon=1,lonlen
!          write(*,*) 'processing longitude =', lon(ilon)
!          Local time at each longitude :
           do it=1,timelen  ! loop temps in
               lt_gcm(it) = time(it) +0.5*lon(ilon)/180.
           end do
           if (nsol.eq.1) lt_gcm(timelen+1) = lt_gcm(1)+1

           do it=1,timelen_lt ! loop time local time
              lt_outc(it)=lt_out(it)
!             Correction to use same local time previous or following
!             day where data are missing :

              if(lt_outc(it).gt.lt_gcm(timelen_tot))lt_outc(it)=lt_out(it)-1
              if(lt_outc(it).lt.lt_gcm(1))lt_outc(it)=lt_out(it)+1
           end do

           if (ndim==3) then
             do ilat=1,latlen
               do it=1,timelen  ! loop temps in
                  var_gcm(it) = var3d(ilon,ilat,it,1)
               end do
               if (nsol.eq.1) var_gcm(timelen+1) = var_gcm(1)
               do it=1,timelen_lt ! loop time local time
                  call interpolf(lt_outc(it),var3d_lt(ilon,ilat,it,1),&
                                 missing,lt_gcm,var_gcm,timelen_tot)
               end do    
             enddo

           else if  (ndim==4) then
             do ilat=1,latlen
               do ialt=1,altlen
                 do it=1,timelen ! loop temps in
                    var_gcm(it) = var3d(ilon,ilat,ialt,it)
                 end do
                 if (nsol.eq.1) var_gcm(timelen+1) = var_gcm(1)
                 do it=1,timelen_lt ! loop time local time
                    call interpolf(lt_outc(it),var3d_lt(ilon,ilat,ialt,it),&
                                   missing,lt_gcm,var_gcm,timelen_tot)
                 end do    
               enddo
             enddo
           end if

        end do


!==============================================================================
! 2.5.5 write (append) to the output file
!==============================================================================

      ierr= NF_PUT_VARA_REAL(nout,varidout,corner,edges,var3d_lt)

      if (ierr.ne.NF_NOERR) then
         write(*,*) 'PUT_VAR ERROR: ',NF_STRERROR(ierr)
         stop ""
      endif

! In case there is a "valid_range" attribute
      ! Write "valid_range" and "missing_value" attributes in output file
      if (miss.eq.NF_NOERR) then
         call missing_value(nout,varidout,validr,miss,valid_range,missing)
      endif

      ! free var3d() array
      deallocate(var3d)
      deallocate(var3d_lt)

   enddo ! of do j=1,nbvar

   deallocate(time)
   deallocate(lt_gcm)
   deallocate(lt_out)
   deallocate(lt_outc)
   deallocate(var_gcm)

   ! Close input file
   ierr=nf_close(nid)

! Close output file
ierr=NF_CLOSE(nout)

write(*,*) "Program completed !"

contains

!******************************************************************************
Subroutine initiate(filename,lat,lon,alt,ctl,latlen,lonlen,altlen,ctllen,GCM_layers,&
                     altlong_name,altunits,altpositive,&
                     nout,latdimout,londimout,altdimout,timedimout,&
                     layerdimout,timevarout)
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
integer,intent(in) :: GCM_layers ! number of GCM layers
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
integer,intent(out) :: layerdimout
! layerdimout: [netcdf] "GCM_layers" ID
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

ierr = NF_DEF_DIM(nout, "GCM_layers", GCM_layers, layerdimout)

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

call def_var(nout,"Time","Time","years since 0000-00-0 00:00:00",1,&
             (/timedimout/),timevarout,ierr)

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
subroutine init2(infid,lonlen,latlen,altlen,GCM_layers,&
                 outfid,londimout,latdimout,altdimout,&
                 layerdimout)
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
integer, intent(in) :: GCM_layers ! # of GCM atmospheric layers
integer, intent(in) :: outfid ! NetCDF output file ID
integer, intent(in) :: londimout ! longitude dimension ID
integer, intent(in) :: latdimout ! latitude dimension ID
integer, intent(in) :: altdimout ! altitude dimension ID
integer, intent(in) :: layerdimout ! GCM_layers dimension ID
!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:,:),allocatable :: phisinit ! Ground geopotential
integer :: apsid,bpsid,phisinitid
integer :: ierr
integer :: tmpdimid ! temporary dimension ID
integer :: tmpvarid ! temporary variable ID
logical :: phis, hybrid ! are "phisinit" and "aps"/"bps" available ?


!==============================================================================
! 1. Read data from input file
!==============================================================================

! hybrid coordinate aps
!write(*,*) "aps: altlen=",altlen," GCM_layers=",GCM_layers
allocate(aps(GCM_layers),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate aps!"
  stop
endif
ierr=NF_INQ_VARID(infid,"aps",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Ooops. Failed to get aps ID. OK."
  hybrid=.false.
else
  hybrid=.true.
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,aps)
  hybrid=.true.
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed reading aps"
  endif

  ! hybrid coordinate bps
!  write(*,*) "bps: altlen=",altlen," GCM_layers=",GCM_layers
  allocate(bps(GCM_layers),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "init2: failed to allocate bps!"
    stop
  endif
  ierr=NF_INQ_VARID(infid,"bps",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Ooops. Failed to get bps ID. OK."
    hybrid=.false.
  else
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,bps)
    if (ierr.ne.NF_NOERR) then
      stop "init2 Error: Failed reading bps"
    endif
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

IF (hybrid) then 
  ! define aps
  call def_var(outfid,"aps","hybrid pressure at midlayers"," ",1,&
               (/layerdimout/),apsid,ierr)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to def_var aps"
  endif

  ! write aps
  ierr=NF_PUT_VAR_REAL(outfid,apsid,aps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to write aps"
  endif
  
  ! define bps
  call def_var(outfid,"bps","hybrid sigma at midlayers"," ",1,&
               (/layerdimout/),bpsid,ierr)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to def_var bps"
  endif

  ! write bps
  ierr=NF_PUT_VAR_REAL(outfid,bpsid,bps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to write bps"
  endif
  
  deallocate(aps)
  deallocate(bps)
ENDIF ! of IF (hybrid)

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

  deallocate(phisinit)
ENDIF ! of IF (phis)


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
subroutine  missing_value(nout,nvarid,validr,miss,valid_range,missing)
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
integer, intent(in) :: validr
! validr : [netcdf] routines return code for "valid_range" attribute
integer, intent(in) :: miss
! miss : [netcdf] routines return code for "missing_value" attribute
real, dimension(2), intent(in) :: valid_range
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

! Write valid_range() attribute
if (validr.eq.NF_NOERR) then
  ierr=NF_PUT_ATT_REAL(nout,nvarid,'valid_range',NF_FLOAT,2,valid_range)

  if (ierr.ne.NF_NOERR) then
     print*,'missing_value: valid_range attribution failed'
     print*, NF_STRERROR(ierr)
     !write(*,*) 'NF_NOERR', NF_NOERR
     !CALL abort
     stop ""
  endif
endif
if (miss.eq.NF_NOERR) then
! Write "missing_value" attribute
  ierr= NF_PUT_ATT_REAL(nout,nvarid,'missing_value',NF_FLOAT,1,missing)

  if (ierr.NE.NF_NOERR) then
     print*, 'missing_value: missing value attribution failed'
     print*, NF_STRERROR(ierr)
  !    WRITE(*,*) 'NF_NOERR', NF_NOERR
  !    CALL abort
     stop ""
  endif
endif

! End netcdf dataset define mode
ierr = NF_ENDDEF(nout)

end subroutine  missing_value

!*****************************************************************************
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

end program localtime

