program simu_MCS
! program for a MCS observer simulator of GCM data
! author : Antoine Bierjon, 2019-2020
! contact : antoine.bierjon@lmd.jussieu.fr
!
!===================================================================================================
!     PREFACE
!===================================================================================================
! This program loads on one hand, a GCM file, which can be :
!       - a GCM diagfi.nc or concat.nc that covers the observation period,
!       - a GCM stats file (stats.nc) which gets extended to cover the observation period
!
! On the other hand, it loads a MRO/MCS data file (binned by Luca Montabone).  
! This MCS file serves then as a reference for the interpolation and the binning of the
! GCM variables contained in the GCM file and given as inputs in the simu_MCS.def
!
! Since the MCS data is binned by intervals of Ls, the GCM simulation's data is interpolated at 
! the MCS spatial coordinates, at every GCM sol value contained in each interval of Ls, before  
! being averaged to create output bins for each interval. The binning in sols also takes into
! account the variability of Local Time in the bin and tries to represent it (with a number of LT
! values in each sol bin equal to the number of values in each associated MCS bin : numbintemp/dust/wice,
! and centered around the MCS bin LT average : timeave).
!
! There are also specific GCM variables that the program looks for in order to make them
! comparable with their equivalent in MCS files. These variables are :
!                                     GCM                         MCS
!                           # temp                    --->   dtemp,ntemp
!                           # dso(/dsodust/qdust)+rho --->   ddust,ndust (dust opacity/km)
!                           # h2o_ice+rho             --->   dwice,nwice (water ice opacity/km)
!
! Eventually, the program outputs a netcdf file, filled with the GCM binned data for dayside
! and nightside and edited with the same format than the MCS file (and with the same 
! dimensions' declaration order).
!
! Minimal requirements :
!     - the MCS file must contain the variables :
!                             dtimeave,dtimemax,dtimemin,
!                             dtemp,dnumbintemp,
!                             ddust,dnumbindust,
!                             dwice,dnumbinwice,
!                             ntimeave,ntimemax,ntimemin,
!                             ntemp,nnumbintemp,
!                             ndust,nnumbindust,
!                             nwice,nnumbinwice
!     - the GCM file must :
!                           # have the same altitude type as the MCS file ;
!                           # cover completely the observation period
!
! See also NOTA BENE in section 2.2
!
!
!
! Algorithm :
!  0. Variable declarations
!  1. OPENING OF THE FILES AND INITIALIZATION OF THE DIMENSIONS
!    1.1 MCS data file : obsfile
!      1.1.1 Open the Observer data file
!      1.1.2 Get dimensions lon,lat,alt,time from the observation file
!    1.2. GCM simulation file : gcmfile
!      1.2.1 Open the GCM simulation file
!      1.2.2 Get dimensions lon,lat,alt,time from the GCM file
!    1.3 Create the output file & initialize the coordinates
!  2. VARIABLES MANAGEMENT
!    2.1 List of the GCM variables to be interpolated
!      2.1.1 Read the GCM variables
!      2.1.2 Handle dust and wice opacities (first set-ups)
!    2.2 Definition of LT variables from obsfile to outfile
![day/night loop begins...
!      2.2.1 Average of Local Time in the OBS bins
!      2.2.2 Maximum of Local Time in the OBS bins
!      2.2.3 Minimum of Local Time in the OBS bins
!    2.3 Definition of numbin variables from obsfile to outfile
!      2.3.1 Number of values in the OBS "temp" bins
!      2.3.2 Number of values in the OBS "dust" bins
!      2.3.3 Number of values in the OBS "wice" bins
!    2.4 Opening of the GCM variables
! [var loop begins...
!      2.4.1 Generic reading of the variable
!      2.4.2 Handle dust and wice opacities (second part)
!    2.5 Opening of the associated MCS variables
!      2.5.1 MCS reference variable (for the missing values)
!      2.5.2 Number of values in the OBS bin (for the sol binning)
!    2.6 Definition of GCM variables in outfile
!  3. EXTRACTION OF THE VARIABLE
!  [coordinates loop begins...
!    3.1 Do some checks and preparations before the extraction
!    3.2 Compute GCM sol date corresponding to Observer Ls (via m_(min/max)sol)
!        and LT (via OBSLT(min/max))
!    3.3 Do the interpolation and binning for the given location
!  ..coordinates loop ends]
!    3.4 Write the data in the netcdf output file
! ..var loop ends]
!  4. END OF THE DAY/NIGHT LOOP
!..day/night loop ends]
!  5. CLOSE THE FILES
!
!  Subroutines
!    extraction
!    inidim
!    ls2sol
!    gen_sol_list
!    status_check
!    LTmod
!===================================================================================================


use netcdf


!===================================================================================
! 0. Variable declarations
!===================================================================================

implicit none ! for no implicitly typed variables

!------------------------
! Files:
character(len=256) :: gcmfile ! GCM simulation file
logical :: is_stats = .false. ! to check if the GCM file is a stats.nc or a diagfi.nc file
character(len=256) :: obsfile ! observation file = MCS/Observer data file
character(len=256) :: outfile ! output file

!------------------------
! NetCDF stuff
integer :: status                                                    ! NetCDF routines return code
character (len=256) :: error_text                                    ! to store the error text
integer :: gcmfid                                                    ! NetCDF gcm file ID
integer :: obsfid                                                    ! NetCDF observation file ID
integer :: outfid                                                    ! NetCDF output file ID
integer :: GCMvarid, OBSvarid, LT_id, numbin_id, outvarid            ! to store the ID of a variable
integer :: lat_dimid_obs,lon_dimid_obs,alt_dimid_obs,time_dimid_obs  ! dimensions' ID in OBS and output files
integer :: lat_dimid_gcm,lon_dimid_gcm,alt_dimid_gcm,time_dimid_gcm  ! dimensions' ID in GCM file
integer :: GCMvarshape(4), OBSvarshape(4), LTshape(3),numbinshape(4) ! to store a variable's coordinates order

!------------------------
! Dimensions
real,dimension(:),allocatable :: GCMlon, OBSlon  ! longitude in the GCM & the Observer files
integer GCMlonlen, OBSlonlen                     ! # of grid points along GCMlon & OBSlon
real,dimension(:),allocatable :: GCMlat, OBSlat  ! latitude in the GCM & the Observer files
integer GCMlatlen, OBSlatlen                     ! # of grid points along GCMlat & OBSlat
real,dimension(:),allocatable :: GCMalt, OBSalt  ! altitude/pressure in the GCM & the Observer files
integer GCMaltlen, OBSaltlen                     ! # of grid point along GCMalt & OBSalt
character(len=1) :: GCMalttype, OBSalttype       ! altitude coord. type:'z' (altitude, m) 'p' (pressure, Pa)
real,dimension(:),allocatable :: GCMtime, OBSLs  ! time in the GCM diagfi (sols) & the Observer files (Ls)
real,dimension(:),allocatable :: GCMstatstime    ! time in the GCM stats file (LT at lon 0°)
integer :: GCMtimelen, GCMstatstimelen, OBSLslen ! # of points along GCMtime, GCMstatstime, OBSLs
real :: starttimeoffset=0.                       ! offset (in sols) wrt Ls=0 of sol 0 in GCM file

!------------------------
! Variables
character(len=10) :: dayornight                        ! are we in the "dayside" or "nightside" loop?
character(len=64),dimension(:),allocatable :: gcm_vars ! list of GCM variables to interpolate
character(len=10),dimension(15) :: notprocessed        ! names of the (15) variables that won't be processed    
integer :: nbvarfile,Nnotprocessed                     ! nbs of variables to deal with the non-processed ones
integer :: nbvar                                       ! nb of variables that will be processed
logical :: var_ok                                      ! is this variable to be processed?
logical :: dustok1,dustok2,dustok3,wiceok              ! is it possible to compute opacities and how?
character(len=64) :: GCMvarname,OBSvarname,outvarname  ! name of the variables
integer :: nbdim                                       ! nb of dimensions of a variable
real,dimension(:,:,:,:),allocatable :: GCM_var,OBS_var ! the 4D variable extracted from GCM & OBS files
real,dimension(:,:,:,:),allocatable :: rho             ! atmospheric density for opacities' computation
character(len=64) :: long_name,units,comment           ! netcdf attributes of the variable

character(len=64) :: OBSLTave_name,OBSLTmax_name,OBSLTmin_name
                                                       ! names of the average, max and min of LT
real,dimension(:,:,:),allocatable :: OBSLT,OBSLTmax,OBSLTmin
                                                       ! 3D variables extracted from obsfile (ave, max and min of LT in a bin)

character (len=64) :: numbin_name                      ! name of the nb of values in an OBS bin
real,dimension(:,:,:,:),allocatable :: numbin          ! nb of values in an OBS temp/dust/wice bin

real :: GCMmiss_val, OBSmiss_val, LTmiss_val           ! value to denote non-existant data
real :: extr_value                                     ! result of the extraction subroutine
real, dimension(:,:,:,:), allocatable :: outvar        ! outvar(,,,): 4D array to store the output variable's data

!------------------------
! Time binning management
real :: OBSdeltaLs            ! difference of Ls between each observation bin
real :: sol, maxsol, minsol   ! sol date corresponding to Observed Ls and LT
integer :: m_maxsol, m_minsol ! indexes of the maximum and minimum GCM sol taken in a bin for interpolation

external LTmod     ! declaration of the function LTmod
real :: LTmod      ! declaration of the type of the function LTmod
integer :: LTcount ! nb of LT samples on which the interpolation is performed,
                   ! for every LT interval (= numbin[lon,lat,alt,Ls])

integer :: solcount                           ! number of GCM sol integer values in one Ls interval
real,dimension(:),allocatable :: int_sol_list ! list of the integer values of GCM sol
real,dimension(:),allocatable :: sol_list     ! list of the sol values used for the interpolation
integer :: solerrcount                        ! nb of GCM missing values during interpolation, removed for the binning 
integer :: errcount = 0                       ! total number of GCM missing values
real :: solbinned_value                       ! extracted value averaged on sol samples, which is finally put in the output bin

!------------------------
! Extraction & loop indices
real :: lon_val, lat_val, alt_val, Ls_val, LT_val ! where and when the output file is written at

integer :: i,j,k,l ! loop iteration indices
integer :: m       ! sol binning loops iteration index
integer :: v,vnot  ! variable loops indices


!===================================================================================
! 1. OPENING OF THE FILES AND INITIALIZATION OF THE DIMENSIONS
!===================================================================================
write(*,*) "Welcome in the MRO/MCS Observer Simulator program !"
!===============================================================================
! 1.1 MCS data file : obsfile
!===============================================================================
!================================================================
! 1.1.1 Open the Observer data file
!================================================================
! Ask the user to give a netcdf observation file
WRITE(*,*) "-> Enter observation file name :"
READ(*,*) obsfile

status=NF90_OPEN(obsfile,nf90_nowrite,obsfid) ! nowrite mode=the program can only read the opened file
error_text="Error: could not open file "//trim(obsfile)
call status_check(status,error_text)

!================================================================
! 1.1.2 Get dimensions lon,lat,alt,time from the observation file
!================================================================
! OBS Latitude
!--------------
status=nf90_inq_dimid(obsfid,"latitude",lat_dimid_obs)
error_text="Failed to find Observer latitude dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(obsfid,lat_dimid_obs,len=OBSlatlen)
error_text="Failed to find Observer latitude length"
call status_check(status,error_text)
allocate(OBSlat(OBSlatlen))

status=nf90_inq_varid(obsfid,"latitude",OBSvarid)
error_text="Failed to find Observer latitude ID"
call status_check(status,error_text)

! Read OBSlat
status=NF90_GET_VAR(obsfid,OBSvarid,OBSlat)
error_text="Failed to load OBSlat"
call status_check(status,error_text)


! OBS Longitude
!--------------
status=nf90_inq_dimid(obsfid,"longitude",lon_dimid_obs)
error_text="Failed to find Observer longitude dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(obsfid,lon_dimid_obs,len=OBSlonlen)
error_text="Failed to find Observer longitude length"
call status_check(status,error_text)
allocate(OBSlon(OBSlonlen))

status=nf90_inq_varid(obsfid,"longitude",OBSvarid)
error_text="Failed to find Observer longitude ID"
call status_check(status,error_text)

! Read OBSlon
status=NF90_GET_VAR(obsfid,OBSvarid,OBSlon)
error_text="Failed to load OBSlon"
call status_check(status,error_text)


! OBS Time (Ls)
!--------------
status=nf90_inq_dimid(obsfid,"time",time_dimid_obs)
error_text="Failed to find Observer time (Ls) dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(obsfid,time_dimid_obs,len=OBSLslen)
error_text="Failed to find Observer time (Ls) length"
call status_check(status,error_text)
allocate(OBSLs(OBSLslen))

status=nf90_inq_varid(obsfid,"time",OBSvarid)
error_text="Failed to find Observer time (Ls) ID"
call status_check(status,error_text)

! Read OBSLs
status=NF90_GET_VAR(obsfid,OBSvarid,OBSLs)
error_text="Failed to load OBSLs"
call status_check(status,error_text)

! Get the observation timestep between bins
OBSdeltaLs=OBSLs(2)-OBSLs(1)


! OBS Altitude
!--------------
status=nf90_inq_dimid(obsfid,"altitude",alt_dimid_obs)
error_text="Failed to find Observer altitude dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(obsfid,alt_dimid_obs,len=OBSaltlen)
error_text="Failed to find Observer altitude length"
call status_check(status,error_text)
allocate(OBSalt(OBSaltlen))

status=nf90_inq_varid(obsfid,"altitude",OBSvarid)
error_text="Failed to find Observer altitude ID"
call status_check(status,error_text)

! Read OBSalt
status=NF90_GET_VAR(obsfid,OBSvarid,OBSalt)
error_text="Failed to load OBSalt"
call status_check(status,error_text)

! Check altitude attribute "units" to find out altitude type and compare with the GCM file
status=nf90_get_att(obsfid,OBSvarid,"units",units)
error_text="Failed to load Observer altitude units attribute"
call status_check(status,error_text)
! an unknown and invisible character is placed just after the unit's 
! characters in the Observer file so we only take the first characters
! corresponding to the sought unit
if (trim(units(1:2)).eq."Pa") then
  units="Pa"
  OBSalttype='p' ! pressure coordinate
else if (trim(units(1:2)).eq."m") then
  units="m"
  OBSalttype='z' ! altitude coordinate
else
  write(*,*)" I do not understand this unit ",trim(units)," for Observer altitude!"
  stop
endif

!===============================================================================
! 1.2. GCM simulation file : gcmfile
!===============================================================================
!================================================================
! 1.2.1 Open the GCM simulation file
!================================================================
! Ask the user to give a netcdf input file
write(*,*)"";WRITE(*,*) "-> Enter input file name (GCM simulation) :"
READ(*,*) gcmfile

! Open GCM file
status=NF90_OPEN(gcmfile,nf90_nowrite,gcmfid)
! nowrite mode=the program can only read the opened file
error_text="Failed to open datafile "//trim(gcmfile)
call status_check(status,error_text)

!================================================================
! 1.2.2 Get dimensions lon,lat,alt,time from the GCM file
!================================================================
! GCM Latitude
!--------------
status=nf90_inq_dimid(gcmfid,"latitude",lat_dimid_gcm)
error_text="Failed to find GCM latitude dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(gcmfid,lat_dimid_gcm,len=GCMlatlen)
error_text="Failed to find GCM latitude length"
call status_check(status,error_text)
allocate(GCMlat(GCMlatlen))

status=nf90_inq_varid(gcmfid,"latitude",GCMvarid)
error_text="Failed to find GCM latitude ID"
call status_check(status,error_text)

! Read GCMlat
status=NF90_GET_VAR(gcmfid,GCMvarid,GCMlat)
error_text="Failed to load GCMlat"
call status_check(status,error_text)


! GCM Longitude
!--------------
status=nf90_inq_dimid(gcmfid,"longitude",lon_dimid_gcm)
error_text="Failed to find GCM longitude dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(gcmfid,lon_dimid_gcm,len=GCMlonlen)
error_text="Failed to find GCM longitude length"
call status_check(status,error_text)
allocate(GCMlon(GCMlonlen))

status=nf90_inq_varid(gcmfid,"longitude",GCMvarid)
error_text="Failed to find GCM longitude ID"
call status_check(status,error_text)

! Read GCMlon
status=NF90_GET_VAR(gcmfid,GCMvarid,GCMlon)
error_text="Failed to load GCMlon"
call status_check(status,error_text)


! GCM Time
!--------------
status=nf90_inq_dimid(gcmfid,"Time",time_dimid_gcm)
error_text="Failed to find GCM time dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(gcmfid,time_dimid_gcm,len=GCMtimelen)
error_text="Failed to find GCM time length"
call status_check(status,error_text)
allocate(GCMtime(GCMtimelen))

status=nf90_inq_varid(gcmfid,"Time",GCMvarid)
error_text="Failed to find GCM time ID"
call status_check(status,error_text)

status=NF90_GET_VAR(gcmfid,GCMvarid,GCMtime)
error_text="Failed to load GCMtime"
call status_check(status,error_text)

! is_stats ?
IF ((GCMtimelen.eq.12).and.(GCMtime(1).eq.2.).and.(GCMtime(GCMtimelen).eq.24.)) then 
  ! if GCM file is a stats, time is in LT at longitude 0° and not in sols at longitude 0°
  write(*,*)"The GCM file is recognized as a stats file."
  is_stats = .true.
  deallocate(GCMtime)
  GCMstatstimelen = GCMtimelen
  allocate(GCMstatstime(GCMstatstimelen))
  status=NF90_GET_VAR(gcmfid,GCMvarid,GCMstatstime)
  error_text="Failed to load GCMstatstime (LT at lon 0)"
  call status_check(status,error_text)
ELSE
  write(*,*)"The GCM file is recognized as a diagfi/concatnc file."
ENDIF

! Simulation time offset management
WRITE(*,*) "Beginning date of the simulation file?"
WRITE(*,*) "(i.e. number of sols since Ls=0 at the Time=0.0 in the GCM file)"
READ(*,*) starttimeoffset
if (.not.is_stats) then
  ! Add the offset to GCMtime(:) if the file is not a stats file
  GCMtime(:)=GCMtime(:)+starttimeoffset
endif 

! Check of temporal coherence between gcmfile & obsfile
call ls2sol(OBSLs(OBSLslen),maxsol) ! maximum date considered
call ls2sol(OBSLs(1),minsol) ! minimum date considered

IF (.not.is_stats) then ! if it is a diagfi, we check the time coherence between the 2 files
  if ((maxsol.gt.maxval(GCMtime)).or.(minsol.lt.minval(GCMtime))) then 
    write(*,*)"Error : obsfile temporal bounds exceed the GCM simulation bounds."
    write(*,*)"Please use a GCM file whose time interval contains the observation period."
    stop
  else
    write(*,*)"Both files are temporally coherent. The program continues..."
  endif

ELSE ! if it is a stats, we create the array GCMtime array (in sols) covering the observation period
     ! and filled with the mean GCM day stored in stats.nc

  GCMtimelen = ((ceiling(maxsol)-floor(minsol)+1)+2) ! we add 2 days in the beginning and the end
                                                     ! to be sure we cover the observation period
  allocate(GCMtime(GCMstatstimelen * GCMtimelen))
  do l=1,GCMtimelen
    do m=1,GCMstatstimelen
      GCMtime(m+(l-1)*GCMstatstimelen) = (floor(minsol)-1) + (l-1) + GCMstatstime(m)/24.
    enddo
  enddo
  GCMtimelen = GCMstatstimelen * GCMtimelen
  write(*,*)"GCMtime has been created from the stats.nc time and the observation period. The program continues..."
ENDIF


! GCM Altitude
!--------------
status=nf90_inq_dimid(gcmfid,"altitude",alt_dimid_gcm)
error_text="Failed to find GCM altitude dimension"
call status_check(status,error_text)

status=nf90_inquire_dimension(gcmfid,alt_dimid_gcm,len=GCMaltlen)
error_text="Failed to find GCM altitude length"
call status_check(status,error_text)
allocate(GCMalt(GCMaltlen))

status=nf90_inq_varid(gcmfid,"altitude",GCMvarid)
error_text="Failed to find GCM altitude ID"
call status_check(status,error_text)

! Read GCMalt
status=NF90_GET_VAR(gcmfid,GCMvarid,GCMalt)
error_text="Failed to load GCMalt"
call status_check(status,error_text)

! Check altitude attribute "units" to find out altitude type
status=nf90_get_att(gcmfid,GCMvarid,"units",units)
error_text="Failed to load GCM altitude units attribute"
call status_check(status,error_text)
if (trim(units).eq."Pa") then
  GCMalttype='p' ! pressure coordinate
else if (trim(units).eq."m") then
  GCMalttype='z' ! altitude coordinate
else
  write(*,*)"I do not understand this unit ",trim(units)," for GCM altitude!"
  if (OBSalttype.eq.'p') then
    write(*,*)"Please use zrecast to put the altitude in the same type as the MCS file (pressure in Pa)"
  else if (OBSalttype.eq.'z') then
    write(*,*)"Please use zrecast to put the altitude in the same type as the MCS file (altitude in m)"
  endif
  stop
endif
IF(OBSalttype.ne.GCMalttype) then
  write(*,*)"Observer altitude type (", OBSalttype,") and ", &
            "GCM altitude type (",GCMalttype,") don't match!"
  stop
ENDIF

!===============================================================================
! 1.3 Create the output file & initialize the coordinates
!===============================================================================
! Name of the outfile
IF (.not.is_stats) then
  outfile=obsfile(1:index(obsfile, ".nc")-1)//"_GCMdiagfi.nc"
ELSE
  outfile=obsfile(1:index(obsfile, ".nc")-1)//"_GCMstats.nc"
ENDIF

! Creation of the outfile
status=NF90_CREATE(outfile,nf90_clobber,outfid)!NB: clobber mode=overwrite any dataset with the same file name !
error_text="Error: could not create file "//trim(outfile)
call status_check(status,error_text)
write(*,*)"";WRITE(*,*)"-> Output file is: ",trim(outfile)

! Creation of the dimensions
call inidim(outfid,OBSlonlen,OBSlatlen,OBSaltlen,OBSLslen,OBSlon,OBSlat,OBSalt,OBSLs,OBSalttype,&
            lon_dimid_obs,lat_dimid_obs,alt_dimid_obs,time_dimid_obs)

!write(*,*)"Dimensions ID in the outfile are :"           
!write(*,*)"lon_dimid=",lon_dimid_obs
!write(*,*)"lat_dimid=",lat_dimid_obs
!write(*,*)"alt_dimid=",alt_dimid_obs
!write(*,*)"time_dimid=",time_dimid_obs



!===================================================================================
! 2. VARIABLES MANAGEMENT
!===================================================================================
!===============================================================================
! 2.1 List of the GCM variables to be interpolated
!===============================================================================
!================================================================
! 2.1.1 Read the GCM variables
!================================================================
! Initialize logicals
dustok1 = .false.
dustok2 = .false.
dustok3 = .false.
wiceok = .false.

! Get nbvarfile (total number of variables in the GCM file)
status=NF90_INQUIRE(gcmfid,nVariables=nbvarfile)
error_text="Error : Pb with nf90_inquire(gcmfid,nVariables=nbvarfile)"
call status_check(status,error_text)

! List of variables that should not be processed
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
notprocessed(12)='cu'
notprocessed(13)='cv'
notprocessed(14)='aire'
notprocessed(15)='phisinit'


! List of variables in the GCM file
write(*,*)""
write(*,*)"List of variables in the GCM file :"
Nnotprocessed=0
do v=1,nbvarfile
  status=NF90_INQUIRE_VARIABLE(gcmfid,v,name=GCMvarname)
  ! GCMvarname now contains the "name" of variable of ID # v
  var_ok=.true.
  do vnot=1,15
    if (GCMvarname.eq.notprocessed(vnot)) then
      var_ok=.false.
      Nnotprocessed=Nnotprocessed+1
    endif 
  enddo        
  if (var_ok) write(*,*) trim(GCMvarname)
  
  ! Detect if we can compute dust and wice opacities
  if (trim(GCMvarname).eq."dso") then
    dustok1 = .true.
  else if (trim(GCMvarname).eq."dsodust") then
    dustok2 = .true.
  else if (trim(GCMvarname).eq."dustq") then
    dustok3 = .true.
  else if (trim(GCMvarname).eq."h2o_ice") then
    wiceok = .true.
  endif
enddo

! Nnotprocessed: # of variables that won't be processed
! nbvarfile: total # of variables in file
! +2: the dust and wice opacities
allocate(gcm_vars(nbvarfile-Nnotprocessed+2),stat=status)
if (status.ne.0) then
  write(*,*) "Error: failed allocation of gcm_vars(nbvarfile-Nnotprocessed+2)"
  write(*,*) "  nbvarfile=",nbvarfile
  write(*,*) "  Nnotprocessed=",Nnotprocessed
  stop
endif

! List of variables to process
write(*,*)
write(*,*) "Which variables do you want to redistribute ?"
write(*,*) "list of <variables> (separated by <Return>s)"
write(*,*) "(an empty line , i.e: just <Return>, implies end of list)"
write(*,*) "NB: this program handles only 4D netcdf variables for now"
nbvar=0
read(*,'(a50)') GCMvarname
do while ((GCMvarname.ne.' ').AND.(trim(GCMvarname).ne."all"))
  nbvar=nbvar+1
  gcm_vars(nbvar)=GCMvarname
  read(*,'(a50)') GCMvarname
enddo

if (GCMvarname.eq."all") then
  nbvar=nbvarfile-Nnotprocessed
  do v=Nnotprocessed+1,nbvarfile
    status=nf90_inquire_variable(gcmfid,v,name=gcm_vars(v-Nnotprocessed))
  enddo
! Variables names from the file are stored in gcm_vars()
  nbvar=nbvarfile-Nnotprocessed
  do v=1,nbvar
    status=nf90_inquire_variable(gcmfid,v+Nnotprocessed,name=gcm_vars(v))
    write(*,'(a9,1x,i2,1x,a1,1x,a64)') "variable ",v,":",gcm_vars(v)
  enddo
else if(nbvar==0) then
  write(*,*) "No variables to process in the GCM file... program stopped"
  stop
endif ! of if (GCMvarname.eq."all")

!================================================================
! 2.1.2 Handle dust and wice opacities (first set-ups)
!================================================================
! 2nd part is in section 2.4.2
write(*,*)
! Load atmospheric density "rho"
if (dustok1.or.dustok2.or.dustok3.or.wiceok) then
  ! Check that the GCM file contains that variable
  status=nf90_inq_varid(gcmfid,"rho",GCMvarid)
  if (status.ne.nf90_noerr) then
    write(*,*) "Failed to find variable rho in "//trim(gcmfile)
    write(*,*) "No computation of opacities..."
    dustok1 =.false.
    dustok2 =.false.
    dustok3 =.false.
    wiceok  =.false.
  else
    ! Length allocation for each dimension of the 4D variable
    allocate(rho(GCMlonlen,GCMlatlen,GCMaltlen,GCMtimelen))

    ! Load datasets
    if (.not.is_stats) then
      status=NF90_GET_VAR(gcmfid,GCMvarid,rho)
      error_text="Failed to load rho"
      call status_check(status,error_text)
    else
      ! if it is a stats file, we load only the first sol, and then copy it to all the other sols
      status=NF90_GET_VAR(gcmfid,GCMvarid,rho(:,:,:,1:GCMstatstimelen))
      error_text="Failed to load rho"
      call status_check(status,error_text)
    !  write(*,*)"GCMstatstimelen = ", GCMstatstimelen
    !  write(*,*)"GCMtimelen = ", GCMtimelen
      do l=(GCMstatstimelen+1),GCMtimelen
        if (modulo(l,GCMstatstimelen).ne.0) then
          rho(:,:,:,l) = rho(:,:,:,modulo(l,GCMstatstimelen))
        else ! if l is a multiple of GCMstatstimelen, since the index modulo(l,GCMstatstimelen)=0
             ! doesn't exist, we make a special case
          rho(:,:,:,l) = rho(:,:,:,GCMstatstimelen)
        endif
      enddo
    endif
    write(*,*) "Variable rho loaded from the GCM file"
  endif
endif ! dustok1.or.dustok2.or.dustok3.or.wiceok

! Dust and wice opacity booleans
if (dustok1.or.dustok2.or.dustok3) then
  nbvar=nbvar+1
  gcm_vars(nbvar)="dust"
endif

if (wiceok) then
  nbvar=nbvar+1
  gcm_vars(nbvar)="wice"
endif

!write(*,*) "gcm_vars retrieved : ",gcm_vars(1:nbvar)

!===============================================================================
! 2.2 Definition of LT variables from obsfile to outfile
!===============================================================================
! --> the day/night loop begins here

!******************** NOTA BENE (cf sections 2.2 and 4)*************************
! We execute the program a first time with the daytime values, and then a second
! time with the nighttime values. 
!*******************************************************************************

dayornight = "dayside" ! we begin with daytime temperature
write(*,*)"" ; write(*,*) "Beginning the 1st loop, on daytime values"; write(*,*)""
DAY_OR_NIGHT: DO ! (the end of the loop is in section 4.)

  SELECT CASE (dayornight)
  CASE ("dayside")
    OBSLTave_name = "dtimeave"
    OBSLTmax_name = "dtimemax"
    OBSLTmin_name = "dtimemin"
  CASE ("nightside")
    OBSLTave_name = "ntimeave"
    OBSLTmax_name = "ntimemax"
    OBSLTmin_name = "ntimemin"
  END SELECT

!================================================================
! 2.2.1 Average of Local Time in the OBS bins
!================================================================
  ! Read the OBS file
  !------------------
  status=nf90_inq_varid(obsfid,trim(OBSLTave_name),LT_id)
  error_text="Failed to find Observer local time ("//trim(OBSLTave_name)//") ID in "//trim(obsfile)
  call status_check(status,error_text)
  status=nf90_inquire_variable(obsfid,LT_id,dimids=LTshape)
  error_text="Failed to get the dim shape of variable "//trim(OBSLTave_name)
  call status_check(status,error_text)
  
  ! Length allocation for each dimension of the 3D variable
  allocate(OBSLT(OBSlonlen,OBSlatlen,OBSLslen))
  
  ! Load datasets
  status=NF90_GET_VAR(obsfid,LT_id,OBSLT)
  error_text="Failed to load "//trim(OBSLTave_name)//" from the obsfile"
  call status_check(status,error_text)
  write(*,*) trim(OBSLTave_name)," loaded from the obsfile"
  
  ! Get LT missing_value attribute
  status=nf90_get_att(obsfid,LT_id,"_FillValue",LTmiss_val)
  error_text="Failed to load missing_value attribute"
  call status_check(status,error_text)

  ! Create the variable in the outfile
  !-----------------------------------
  ! Switch to netcdf define mode
  status=nf90_redef(outfid)
  error_text="Error: could not switch to define mode in the outfile"
  call status_check(status,error_text)
  
  ! Definition of the variable
  status=NF90_DEF_VAR(outfid,trim(OBSLTave_name),nf90_float,LTshape,LT_id)
  error_text="Error: could not define the variable "//trim(OBSLTave_name)//" in the outfile"
  call status_check(status,error_text)
  
  ! Write the attributes
  select case (dayornight)
  case ("dayside")
    status=nf90_put_att(outfid,LT_id,"long_name","Average local time in bin - day side [6h, 18h]")
  case ("nightside")
    status=nf90_put_att(outfid,LT_id,"long_name","Average local time in bin - night side [18h, 6h]")
  end select
  status=nf90_put_att(outfid,LT_id,"units","hours")
  status=nf90_put_att(outfid,LT_id,"_FillValue",LTmiss_val)
  
  ! End the netcdf define mode (and thus enter the "data writing" mode)
  status=nf90_enddef(outfid)
  error_text="Error: could not close the define mode of the outfile"
  call status_check(status,error_text)
  
  ! Write the data in the output file
  status = NF90_PUT_VAR(outfid, LT_id, OBSLT) ! write the MCS d/ntimeave as the output Local Time
  error_text="Error: could not write "//trim(OBSLTave_name)//" data in the outfile"
  call status_check(status,error_text)

  write(*,*)"Local Time (",trim(OBSLTave_name),") has been created in the outfile"
  write(*,'("  with missing_value attribute : ",1pe12.5)')LTmiss_val

!================================================================
! 2.2.2 Maximum of Local Time in the OBS bins
!================================================================
  ! Read the OBS file
  !------------------
  status=nf90_inq_varid(obsfid,trim(OBSLTmax_name),LT_id)
  error_text="Failed to find Observer max local time ("//trim(OBSLTmax_name)//") ID in "//trim(obsfile)
  call status_check(status,error_text)
  status=nf90_inquire_variable(obsfid,LT_id,dimids=LTshape)
  error_text="Failed to get the dim shape of variable "//trim(OBSLTmax_name)
  call status_check(status,error_text)
  
  ! Length allocation for each dimension of the 3D variable
  allocate(OBSLTmax(OBSlonlen,OBSlatlen,OBSLslen))
  
  ! Load datasets
  status=NF90_GET_VAR(obsfid,LT_id,OBSLTmax)
  error_text="Failed to load "//trim(OBSLTmax_name)//" from the obsfile"
  call status_check(status,error_text)
  write(*,*) trim(OBSLTmax_name)," loaded from the obsfile"

!================================================================
! 2.2.3 Minimum of Local Time in the OBS bins
!================================================================
  ! Read the OBS file
  !------------------
  status=nf90_inq_varid(obsfid,trim(OBSLTmin_name),LT_id)
  error_text="Failed to find Observer min local time ("//trim(OBSLTmin_name)//") ID in "//trim(obsfile)
  call status_check(status,error_text)
  status=nf90_inquire_variable(obsfid,LT_id,dimids=LTshape)
  error_text="Failed to obtain information on variable "//trim(OBSLTmin_name)
  call status_check(status,error_text)

  ! Length allocation for each dimension of the 3D variable
  allocate(OBSLTmin(OBSlonlen,OBSlatlen,OBSLslen))
  
  ! Load datasets
  status=NF90_GET_VAR(obsfid,LT_id,OBSLTmin)
  error_text="Failed to load "//trim(OBSLTmin_name)//" from the obsfile"
  call status_check(status,error_text)
  write(*,*) trim(OBSLTmin_name)," loaded from the obsfile"
  write(*,*)""

!===============================================================================
! 2.3 Definition of numbin variables from obsfile to outfile
!===============================================================================
!================================================================
! 2.3.1 Number of values in the OBS "temp" bins
!================================================================
  SELECT CASE (dayornight)
  CASE ("dayside")
    numbin_name = "dnumbintemp"
  CASE ("nightside")
    numbin_name = "nnumbintemp"
  END SELECT
  
  ! Read the OBS file
  !------------------
  status=nf90_inq_varid(obsfid,trim(numbin_name),numbin_id)
  error_text="Failed to find Observer nb of temp values in bin ("//trim(numbin_name)//")'s ID in "//trim(obsfile)
  call status_check(status,error_text)
  status=nf90_inquire_variable(obsfid,numbin_id,dimids=numbinshape)
  error_text="Failed to obtain information on variable "//trim(numbin_name)
  call status_check(status,error_text)
  
  ! Length allocation for each dimension of the 4D variable
  allocate(numbin(OBSlonlen,OBSlatlen,OBSaltlen,OBSLslen))
  
  ! Load datasets
  status=NF90_GET_VAR(obsfid,numbin_id,numbin)
  error_text="Failed to load "//trim(numbin_name)//" from the obsfile"
  call status_check(status,error_text)
  write(*,*) trim(numbin_name)," loaded from the obsfile"

  ! Create the variable in the outfile
  !-----------------------------------
  ! Switch to netcdf define mode
  status=nf90_redef(outfid)
  error_text="Error: could not switch to define mode in the outfile"
  call status_check(status,error_text)
  
  ! Definition of the variable
  status=NF90_DEF_VAR(outfid,trim(numbin_name),nf90_float,numbinshape,numbin_id)
  error_text="Error: could not define the variable "//trim(numbin_name)//" in the outfile"
  call status_check(status,error_text)
  
  ! Write the attributes
  select case (dayornight)
  case ("dayside")
    status=nf90_put_att(outfid,numbin_id,"long_name","Number of temp values in bin - day side")
  case ("nightside")
    status=nf90_put_att(outfid,numbin_id,"long_name","Number of temp values in bin - night side")
  end select

  ! End the netcdf define mode (and thus enter the "data writing" mode)
  status=nf90_enddef(outfid)
  error_text="Error: could not close the define mode of the outfile"
  call status_check(status,error_text)
  
  ! Write the data in the output file
  status = NF90_PUT_VAR(outfid, numbin_id, numbin)
  error_text="Error: could not write "//trim(numbin_name)//" data in the outfile"
  call status_check(status,error_text)
  
  write(*,*)"Number of temp values in bin (",trim(numbin_name),") has been created in the outfile"
  write(*,*)""
  deallocate(numbin)

!================================================================
! 2.3.2 Number of values in the OBS "dust" bins
!================================================================
  SELECT CASE (dayornight)
  CASE ("dayside")
    numbin_name = "dnumbindust"
  CASE ("nightside")
    numbin_name = "nnumbindust"
  END SELECT
  
  ! Read the OBS file
  !------------------
  status=nf90_inq_varid(obsfid,trim(numbin_name),numbin_id)
  error_text="Failed to find Observer nb of dust values in bin ("//trim(numbin_name)//")'s ID in "//trim(obsfile)
  call status_check(status,error_text)
  status=nf90_inquire_variable(obsfid,numbin_id,dimids=numbinshape)
  error_text="Failed to obtain information on variable "//trim(numbin_name)
  call status_check(status,error_text)
  
  ! Length allocation for each dimension of the 4D variable
  allocate(numbin(OBSlonlen,OBSlatlen,OBSaltlen,OBSLslen))
  
  ! Load datasets
  status=NF90_GET_VAR(obsfid,numbin_id,numbin)
  error_text="Failed to load "//trim(numbin_name)//" from the obsfile"
  call status_check(status,error_text)
  write(*,*) trim(numbin_name)," loaded from the obsfile"

  ! Create the variable in the outfile
  !-----------------------------------
  ! Switch to netcdf define mode
  status=nf90_redef(outfid)
  error_text="Error: could not switch to define mode in the outfile"
  call status_check(status,error_text)
  
  ! Definition of the variable
  status=NF90_DEF_VAR(outfid,trim(numbin_name),nf90_float,numbinshape,numbin_id)
  error_text="Error: could not define the variable "//trim(numbin_name)//" in the outfile"
  call status_check(status,error_text)
  
  ! Write the attributes
  select case (dayornight)
  case ("dayside")
    status=nf90_put_att(outfid,numbin_id,"long_name","Number of dust values in bin - day side")
  case ("nightside")
    status=nf90_put_att(outfid,numbin_id,"long_name","Number of dust values in bin - night side")
  end select

  ! End the netcdf define mode (and thus enter the "data writing" mode)
  status=nf90_enddef(outfid)
  error_text="Error: could not close the define mode of the outfile"
  call status_check(status,error_text)
  
  ! Write the data in the output file
  status = NF90_PUT_VAR(outfid, numbin_id, numbin)
  error_text="Error: could not write "//trim(numbin_name)//" data in the outfile"
  call status_check(status,error_text)
  
  write(*,*)"Number of dust values in bin (",trim(numbin_name),") has been created in the outfile"
  write(*,*)""
  deallocate(numbin)
  
!================================================================
! 2.3.3 Number of values in the OBS "wice" bins
!================================================================
  SELECT CASE (dayornight)
  CASE ("dayside")
    numbin_name = "dnumbinwice"
  CASE ("nightside")
    numbin_name = "nnumbinwice"
  END SELECT
  
  ! Read the OBS file
  !------------------
  status=nf90_inq_varid(obsfid,trim(numbin_name),numbin_id)
  error_text="Failed to find Observer nb of wice values in bin ("//trim(numbin_name)//")'s ID in "//trim(obsfile)
  call status_check(status,error_text)
  status=nf90_inquire_variable(obsfid,numbin_id,dimids=numbinshape)
  error_text="Failed to obtain information on variable "//trim(numbin_name)
  call status_check(status,error_text)
  
  ! Length allocation for each dimension of the 4D variable
  allocate(numbin(OBSlonlen,OBSlatlen,OBSaltlen,OBSLslen))
  
  ! Load datasets
  status=NF90_GET_VAR(obsfid,numbin_id,numbin)
  error_text="Failed to load "//trim(numbin_name)//" from the obsfile"
  call status_check(status,error_text)
  write(*,*) trim(numbin_name)," loaded from the obsfile"

  ! Create the variable in the outfile
  !-----------------------------------
  ! Switch to netcdf define mode
  status=nf90_redef(outfid)
  error_text="Error: could not switch to define mode in the outfile"
  call status_check(status,error_text)
  
  ! Definition of the variable
  status=NF90_DEF_VAR(outfid,trim(numbin_name),nf90_float,numbinshape,numbin_id)
  error_text="Error: could not define the variable "//trim(numbin_name)//" in the outfile"
  call status_check(status,error_text)
  
  ! Write the attributes
  select case (dayornight)
  case ("dayside")
    status=nf90_put_att(outfid,numbin_id,"long_name","Number of wice values in bin - day side")
  case ("nightside")
    status=nf90_put_att(outfid,numbin_id,"long_name","Number of wice values in bin - night side")
  end select

  ! End the netcdf define mode (and thus enter the "data writing" mode)
  status=nf90_enddef(outfid)
  error_text="Error: could not close the define mode of the outfile"
  call status_check(status,error_text)
  
  ! Write the data in the output file
  status = NF90_PUT_VAR(outfid, numbin_id, numbin)
  error_text="Error: could not write "//trim(numbin_name)//" data in the outfile"
  call status_check(status,error_text)
  
  write(*,*)"Number of wice values in bin (",trim(numbin_name),") has been created in the outfile"
  write(*,*)""
  deallocate(numbin)


!===============================================================================
! 2.4 Opening of the GCM variables
!===============================================================================
! --> the var loop begins here

  VAR: DO v=1,nbvar ! LOOP ON ALL THE GCM VARIABLES TO PROCESS
                    ! (the end of the loop is in section 3.4)
    
    GCMvarname = gcm_vars(v)
    
    ! Detect the dust and wice opacities special cases
    if (trim(GCMvarname).eq."dust") then
      if (dustok1) then ! "dso" is detected in gcmfile
        GCMvarname="dso"
        dustok2=.false.
        dustok3=.false.
      else if (dustok2) then ! "dsodust" is detected in gcmfile
        GCMvarname="dsodust"
        dustok3=.false.
      else if (dustok3) then ! "dustq" is detected in gcmfile
        GCMvarname="dustq"
      endif
      write(*,*) "Computing dust opacity..."
    endif
    if (trim(GCMvarname).eq."wice") then ! "h2o_ice" detected in gcmfile
      GCMvarname="h2o_ice"
      write(*,*) "Computing water ice opacity..."
    endif

!================================================================
! 2.4.1 Generic reading of the variable
!================================================================
    ! Check that the GCM file contains that variable
    status=nf90_inq_varid(gcmfid,trim(GCMvarname),GCMvarid)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed to find variable "//trim(GCMvarname)//" in "//trim(gcmfile)
      write(*,*) "We'll skip that variable..."
      CYCLE VAR ! go directly to the next variable
    endif

    ! Sanity checks on the variable
    status=nf90_inquire_variable(gcmfid,GCMvarid,ndims=nbdim,dimids=GCMvarshape)
    error_text="Failed to obtain information on variable "//trim(GCMvarname)
    call status_check(status,error_text)

    ! Check that it is a 4D variable
    if (nbdim.ne.4) then
      write(*,*) "Error:",trim(GCMvarname)," is not a 4D variable"
      write(*,*) "We'll skip that variable...";write(*,*)""
      CYCLE VAR ! go directly to the next variable
    endif
    ! Check that its dimensions are indeed lon,lat,alt,time (in the right order)
    if (GCMvarshape(1).ne.lon_dimid_gcm) then
      write(*,*) "Error, expected first dimension of ",trim(GCMvarname)," to be longitude!"
      write(*,*) "We'll skip that variable..."
      CYCLE VAR ! go directly to the next variable
    endif
    if (GCMvarshape(2).ne.lat_dimid_gcm) then
      write(*,*) "Error, expected second dimension of ",trim(GCMvarname)," to be latitude!"
      write(*,*) "We'll skip that variable..."
      CYCLE VAR ! go directly to the next variable
    endif
    if (GCMvarshape(3).ne.alt_dimid_gcm) then
      write(*,*) "Error, expected third dimension of ",trim(GCMvarname)," to be altitude!"
      write(*,*) "We'll skip that variable..."
      CYCLE VAR ! go directly to the next variable
    endif
    if (GCMvarshape(4).ne.time_dimid_gcm) then
      write(*,*) "Error, expected fourth dimension of ",trim(GCMvarname)," to be time!"
      write(*,*) "We'll skip that variable..."
      CYCLE VAR ! go directly to the next variable
    endif

    ! Length allocation for each dimension of the 4D variable
    allocate(GCM_var(GCMlonlen,GCMlatlen,GCMaltlen,GCMtimelen))

    ! Load datasets
    if (.not.is_stats) then
      status=NF90_GET_VAR(gcmfid,GCMvarid,GCM_var)
      error_text="Failed to load "//trim(GCMvarname)
      call status_check(status,error_text)
    else
      ! if it is a stats file, we load only the first sol, and then copy it to all the other sols
      status=NF90_GET_VAR(gcmfid,GCMvarid,GCM_var(:,:,:,1:GCMstatstimelen))
      error_text="Failed to load "//trim(GCMvarname)
      call status_check(status,error_text)
    !  write(*,*)"GCMstatstimelen = ", GCMstatstimelen
    !  write(*,*)"GCMtimelen = ", GCMtimelen
      do l=(GCMstatstimelen+1),GCMtimelen
        if (modulo(l,GCMstatstimelen).ne.0) then
          GCM_var(:,:,:,l) = GCM_var(:,:,:,modulo(l,GCMstatstimelen))
        else ! if l is a multiple of GCMstatstimelen, since the index modulo(l,GCMstatstimelen)=0
             ! doesn't exist, we make a special case
          GCM_var(:,:,:,l) = GCM_var(:,:,:,GCMstatstimelen)
        endif
      enddo
    endif
    write(*,*) "Variable ",trim(GCMvarname)," loaded from the GCM file"

    ! Get dataset's missing_value attribute
    status=nf90_get_att(gcmfid,GCMvarid,"missing_value",GCMmiss_val)
    error_text="Failed to load missing_value attribute"
    call status_check(status,error_text)

    ! Get other variable's attributes
    status=nf90_get_att(gcmfid,GCMvarid,"long_name",long_name)
    if (status.ne.nf90_noerr) then
    ! if no attribute "long_name", try "title"
      status=nf90_get_att(gcmfid,GCMvarid,"title",long_name)
    endif
    status=nf90_get_att(gcmfid,GCMvarid,"units",units)

!================================================================
! 2.4.2 Handle dust and wice opacities (second part)
!================================================================
    ! DUST
    !-----
    if (trim(gcm_vars(v)).eq."dust") then

     IF (dustok1.or.dustok2) THEN
     ! Dust opacity computed from its density-scaled opacity
!       write(*,*)long_name(index(long_name,"(")+1:index(long_name,")")-1)
       do i=1,GCMlonlen
        do j=1,GCMlatlen
         do k=1,GCMaltlen
          do l=1,GCMtimelen
            if (GCM_var(i,j,k,l).ne.GCMmiss_val) then
              ! Multiply by rho to have opacity [1/km]
              GCM_var(i,j,k,l) = GCM_var(i,j,k,l) * rho(i,j,k,l) *1000.

              if (long_name(index(long_name,"(")+1:index(long_name,")")-1).eq."TES") then
               ! The density-scaled opacity was calibrated on TES wavelength (9.3um) 
               ! so we must recalibrate it to MCS wavelength (21.6um) using recalibration 
               ! coefficients from Montabone et al. 2015, section 2.3
                GCM_var(i,j,k,l) = 1.3/2.7 * GCM_var(i,j,k,l)
              endif
            endif
          enddo
         enddo
        enddo
       enddo

       long_name = "IR Dust opacity (from DSO)"

     ELSE IF (dustok3) THEN
     ! Dust opacity computed from its mass mixing ratio
       do i=1,GCMlonlen
        do j=1,GCMlatlen
         do k=1,GCMaltlen
          do l=1,GCMtimelen
            if (GCM_var(i,j,k,l).ne.GCMmiss_val) then
              ! Opacity is computed from the equation of Heavens et al. 2014, section 2.3,
              ! assuming a rho_dust=3000 kg/m3 and an effective radius reff=1.06 microns
              ! + the opacity is here in 1/km and the MMR in kg/kg
              GCM_var(i,j,k,l) = GCM_var(i,j,k,l) * rho(i,j,k,l) / 0.012 * 1000 
            endif
          enddo
         enddo
        enddo
       enddo

       long_name = "IR Dust opacity (from MMR)"

     ENDIF
    
     GCMvarname = gcm_vars(v) ! reput the right name in GCMvarname
     units = "opacity/km"
    endif ! trim(gcm_vars(v)).eq."dust"
    
    
    ! WICE
    !-----
    if (trim(gcm_vars(v)).eq."wice") then
    ! Water ice opacity computed from its mass mixing ratio
      do i=1,GCMlonlen
       do j=1,GCMlatlen
        do k=1,GCMaltlen
         do l=1,GCMtimelen
           if (GCM_var(i,j,k,l).ne.GCMmiss_val) then
             ! Opacity at MCS wavelength (11.9um) is computed from an equation
             ! similar to the one of Heavens et al. 2014, section 2.3.
             ! We assume a rho_wice=920 kg/m3, an effective radius reff=3um,
             ! an extinction coefficient Qext(wvl,reff)=1.54471
             GCM_var(i,j,k,l) = 750*1.54471* GCM_var(i,j,k,l) * rho(i,j,k,l) / (920*3e-6)
           endif
         enddo
        enddo
       enddo
      enddo

      long_name = "IR Water ice opacity (from MMR)"
    
     GCMvarname = gcm_vars(v) ! reput the right name in GCMvarname
     units = "opacity/km"
    endif ! trim(gcm_vars(v)).eq."wice"

!===============================================================================
! 2.5 Opening of the associated MCS variables
!===============================================================================
    ! Observer variables to extract :
    IF ((index(GCMvarname,"dust").ne.0).or.(index(GCMvarname,"dso").ne.0)) THEN
      ! if the variable name contains "dust" or "dso". Especially for the targeted variables :
      ! dustq,dustN,dsodust,reffdust,opadust, and their equivalents for stormdust & topdust
      OBSvarname  = "dust"
      numbin_name = "numbindust"
    ELSE IF ((trim(GCMvarname).eq."h2o_ice").or.(trim(GCMvarname).eq."wice") &
           .or.(trim(GCMvarname).eq."reffice").or.(trim(GCMvarname).eq."opawice")) THEN
      OBSvarname  = "wice"
      numbin_name = "numbinwice"
    ELSE ! default case is temp binning, since it contains the most values
      OBSvarname  = "temp"
      numbin_name = "numbintemp"
    ENDIF
    
    SELECT CASE (dayornight)
    CASE ("dayside")
      OBSvarname  = "d"//OBSvarname
      numbin_name = "d"//numbin_name
    CASE ("nightside")
      OBSvarname  = "n"//OBSvarname
      numbin_name = "n"//numbin_name
    END SELECT
    
!================================================================
! 2.5.1 MCS reference variable (for the missing values)
!================================================================    
    ! Check that the observation file contains that variable
    status=nf90_inq_varid(obsfid,trim(OBSvarname),OBSvarid)
    error_text="Failed to find variable "//trim(OBSvarname)//" in "//trim(obsfile)
    call status_check(status,error_text)

    ! Sanity checks on the variable
    status=nf90_inquire_variable(obsfid,OBSvarid,ndims=nbdim,dimids=OBSvarshape)
    error_text="Failed to obtain information on variable "//trim(OBSvarname)
    call status_check(status,error_text)

    ! Check that it is a 4D variable
    if (nbdim.ne.4) then
      write(*,*) "Error, expected a 4D (lon-lat-alt-time) variable for ",trim(OBSvarname)
      stop
    endif
    ! Check that its dimensions are indeed lon,lat,alt,time (in the right order)
    if (OBSvarshape(1).ne.lon_dimid_obs) then
      write(*,*) "Error, expected first dimension of ",trim(OBSvarname)," to be longitude!"
      stop
    endif
    if (OBSvarshape(2).ne.lat_dimid_obs) then
      write(*,*) "Error, expected second dimension of ",trim(OBSvarname)," to be latitude!"
      stop
    endif
    if (OBSvarshape(3).ne.alt_dimid_obs) then
      write(*,*) "Error, expected third dimension of ",trim(OBSvarname)," to be altitude!"
      stop
    endif
    if (OBSvarshape(4).ne.time_dimid_obs) then
      write(*,*) "Error, expected fourth dimension of ",trim(OBSvarname)," to be time!"
      stop
    endif

    ! Length allocation for each dimension of the 4D variable
    allocate(OBS_var(OBSlonlen,OBSlatlen,OBSaltlen,OBSLslen))

    ! Load datasets
    status=NF90_GET_VAR(obsfid,OBSvarid,OBS_var)
    error_text="Failed to load "//trim(OBSvarname)//" from the obsfile"
    call status_check(status,error_text)
    write(*,*) trim(OBSvarname)," loaded from the obsfile as reference variable"

    ! Get OBS_var missing_value attribute
    status=nf90_get_att(obsfid,OBSvarid,"_FillValue",OBSmiss_val)
    error_text="Failed to load missing_value attribute"
    call status_check(status,error_text)

!================================================================
! 2.5.2 Number of values in the OBS bin (for the sol binning)
!================================================================
    ! Check that the observation file contains that variable
    status=nf90_inq_varid(obsfid,trim(numbin_name),numbin_id)

    ! Checks have already been done in section 2.3
    
    ! Length allocation for each dimension of the 4D variable
    allocate(numbin(OBSlonlen,OBSlatlen,OBSaltlen,OBSLslen))

    ! Load datasets
    status=NF90_GET_VAR(obsfid,numbin_id,numbin)
    
!===============================================================================
! 2.6 Definition of GCM variables in outfile
!===============================================================================
    ! Switch to netcdf define mode
    status=nf90_redef(outfid)
    error_text="Error: could not switch to define mode in the outfile"
    call status_check(status,error_text)

    ! Definition of the variable
    SELECT CASE (dayornight)
    CASE ("dayside")
      outvarname = "d"//GCMvarname
    CASE ("nightside")
      outvarname = "n"//GCMvarname
    END SELECT
    status=NF90_DEF_VAR(outfid,trim(outvarname),nf90_float,OBSvarshape,outvarid)
    error_text="Error: could not define the variable "//trim(outvarname)//" in the outfile"
    call status_check(status,error_text)

    ! Write the attributes
    SELECT CASE (dayornight)
    CASE ("dayside")
      long_name = trim(long_name)//" - day side"
      status=nf90_put_att(outfid,outvarid,"long_name",long_name)
    CASE ("nightside")
      long_name = trim(long_name)//" - night side"
      status=nf90_put_att(outfid,outvarid,"long_name",long_name)
    END SELECT
    status=nf90_put_att(outfid,outvarid,"units",units)
    status=nf90_put_att(outfid,outvarid,"_FillValue",OBSmiss_val)
    comment = "Reference numbin: "//trim(numbin_name)
    status=nf90_put_att(outfid,outvarid,"comment",comment)
    
    write(*,*)trim(outvarname)," has been created in the outfile"
    write(*,'("  with missing_value attribute : ",1pe12.5)')OBSmiss_val
    write(*,*)""
    
    ! End the netcdf define mode (and thus enter the "data writing" mode)
    status=nf90_enddef(outfid)
    error_text="Error: could not close the define mode of the outfile"
    call status_check(status,error_text)
    
    allocate(outvar(OBSlonlen,OBSlatlen,OBSaltlen,OBSLslen))


!===================================================================================
! 3. EXTRACTION OF THE VARIABLE
!===================================================================================
!===============================================================================
! 3.1 Do some checks and preparations before the extraction
!===============================================================================
    Ls: do l=1,OBSLslen ! loop on all observed Solar longitudes
      Ls_val=OBSLs(l) ! Ls_val=center of the output bin
      if ((Ls_val.lt.0.).or.(Ls_val.gt.360.)) then
        write(*,*) "Unexpected value for OBSLs: ",Ls_val
        stop
      endif

      ! Convert the Ls bin into a sol interval on which the binning is done :
      !----------------------------------------------------------------------
      !-> get the index of the maximum value of GCM sol (m_maxsol) that is lower than Ls bin's superior bound (maxsol)
      call ls2sol(Ls_val+OBSdeltaLs/2.,maxsol)
      m_maxsol=1
      do while (((GCMtime(m_maxsol)+0.5).lt.maxsol).AND.(m_maxsol.le.(GCMtimelen-1)))
    !! The +0.5 is there to take into account the whole planet (lon=[-180°;180°]) and not just the lon=0° from the GCM
        m_maxsol=m_maxsol+1
      enddo
      !-> get the index of the minimum value of GCM sol (m_minsol) that is greater than Ls bin's inferior bound (minsol)
      call ls2sol(Ls_val-OBSdeltaLs/2.,minsol)
      m_minsol=1
      do while (((GCMtime(m_minsol)-0.5).le.minsol).AND.(m_minsol.le.(GCMtimelen-1)))
    !! Same comment for the -0.5
        m_minsol=m_minsol+1
      enddo
      if (m_minsol.gt.m_maxsol) then
        write(*,*) "No value in gcmfile between sol=",minsol," and sol=",maxsol," (Ls=",Ls_val,"Â°)"
        ! Write a missing_value to output
        outvar(:,:,:,l) = OBSmiss_val
        CYCLE Ls ! go directly to the next Ls
      endif  
      ! Get all the integer values of GCM sols that fit in this interval
      solcount=floor(GCMtime(m_maxsol))-ceiling(GCMtime(m_minsol))+1
      ! sols that are not fully in the interval are not counted
      allocate(int_sol_list(solcount))
!      write(*,*) "GCMminsol=", GCMtime(m_minsol)
!      write(*,*)"GCMmaxsol=", GCMtime(m_maxsol)
      do m=1,solcount
        int_sol_list(m)=ceiling(GCMtime(m_minsol)) + m-1
      enddo
!      write(*,*)"int_sol_list=",int_sol_list


      latitude: do j=1,OBSlatlen ! loop on all observed latitudes
        lat_val=OBSlat(j)
        if ((lat_val.lt.-90.).or.(lat_val.gt.90.)) then
          write(*,*) "Unexpected value for OBSlat: ",lat_val
          stop
        endif

        longitude: do i=1,OBSlonlen ! loop on all observed longitudes
          lon_val=OBSlon(i)
          if ((lon_val.lt.-360.).or.(lon_val.gt.360.)) then
            write(*,*) "Unexpected value for lon_val: ",lon_val
            stop
          endif
          ! We want lon_val in [-180:180] for the subroutine extraction
          if (lon_val.lt.-180.) lon_val=lon_val+360.
          if (lon_val.gt.180.) lon_val=lon_val-360.

          LT_val=OBSLT(i,j,l) ! find the Observer average LT value at bin(lon_val, lat_val, Ls_val)

          if ((LT_val.lt.0.).or.(LT_val.gt.24.)) then
            if (LT_val.eq.LTmiss_val) then
!              write(*,*) "Missing value in obsfile for LT_val"
              ! Write a missing_value to output
              outvar(i,j,:,l) = OBSmiss_val
              CYCLE longitude ! go directly to the next longitude
            else
              write(*,*) "Unexpected value for LT_val: ",LT_val
              stop
            endif
          endif

          altitude: do k=1,OBSaltlen ! loop on all observed altitudes
            alt_val=OBSalt(k)
            if (OBS_var(i,j,k,l).eq.OBSmiss_val) then
!              write(*,*) "Missing value in obsfile for ",OBSvarname
              ! Write a missing_value to output
              outvar(i,j,k,l) = OBSmiss_val
              CYCLE altitude ! go directly to the next altitude
            endif

!===============================================================================
! 3.2 Compute GCM sol date corresponding to Observer Ls (via m_(min/max)sol)
!       and LT (via OBSLT(min/max))
!===============================================================================
            LTcount=floor(numbin(i,j,k,l)) ! find the Observer number of temp values
                                           ! at bin(lon_val,lat_val,alt_val,Ls_val)
            if (LTcount.eq.0.) then
              ! Write a missing_value to output
              outvar(i,j,k,l) = OBSmiss_val
              CYCLE altitude ! go directly to the next altitude
            endif
            if (LTcount.lt.0.) then
              write(*,*) "Unexpected value for LTcount: ",LTcount
              stop
            endif

            ! Generate the sol list for the interpolation
            allocate(sol_list(solcount*LTcount))
            call gen_sol_list(solcount,int_sol_list,LTcount,LT_val,OBSLTmax(i,j,l),&
                              OBSLTmin(i,j,l),lon_val,LTmod,dayornight,&
                              sol_list) 

            solerrcount=0
            solbinned_value=0
            sol_bin: do m=1,solcount*LTcount ! loop on all GCM sols of the bin
              sol=sol_list(m)
!              write(*,*)"sol=",sol
!=============================================================================== 
! 3.3 Do the interpolation and binning for the given location
!===============================================================================
              call extraction(lon_val,lat_val,alt_val,sol,&
                              GCMlonlen,GCMlatlen,GCMaltlen,GCMtimelen,&
                              GCMlon,GCMlat,GCMalt,GCMtime,&
                              GCM_var,GCMmiss_val,GCMalttype,GCMvarname,extr_value)

              if (extr_value.eq.GCMmiss_val) then
!                write(*,*) "Missing value in gcmfile at lon=",lon_val,"; lat=",lat_val,"; alt=",alt_val,"; sol=",sol 
                solerrcount=solerrcount+1  
                CYCLE sol_bin ! go directly to the next GCM sol of the bin
              endif
              solbinned_value=solbinned_value+extr_value

            enddo sol_bin ! end loop on all GCM sols of the bin
            if ((solcount*LTcount-solerrcount).ne.0) then
              solbinned_value=solbinned_value/(solcount*LTcount-solerrcount)
            else
!              write(*,*)"No GCM value in this sol bin"
              solbinned_value=OBSmiss_val
            endif
            ! Write value to output
            outvar(i,j,k,l)=solbinned_value

            errcount=errcount+solerrcount

            deallocate(sol_list)

          enddo altitude ! end loop on observed altitudes
        enddo longitude ! end loop on observed longitudes
      enddo latitude ! end loop on observed latitudes

      deallocate(int_sol_list)

    enddo Ls ! end loop on observed Solar longitudes

!    write(*,*)"Nb of GCM missing values :",errcount

!===============================================================================  
! 3.4 Write the data in the netcdf output file
!===============================================================================
    status = nf90_put_var(outfid, outvarid, outvar)
    error_text="Error: could not write "//trim(outvarname)//" data in the outfile"
    call status_check(status,error_text)
  
    ! Deallocations before going to the next GCM variable
    deallocate(GCM_var)
    deallocate(OBS_var)
    deallocate(numbin)
    deallocate(outvar)
  ENDDO VAR ! end loop on variables
  
!===================================================================================  
! 4. END OF THE DAY/NIGHT LOOP
!===================================================================================
  IF (dayornight.EQ."dayside") THEN
    ! this is the end of the first loop (on daytime values)
    ! and we still have to loop on nighttime values
    dayornight="nightside"
    deallocate(OBSLT)
    deallocate(OBSLTmax)
    deallocate(OBSLTmin)
    write(*,*)""
    write(*,*) "Beginning the 2nd loop, on nighttime values"; write(*,*)""
    CYCLE DAY_OR_NIGHT
  ELSE ! i.e. dayornight="nightside"
    ! this is the end of the second loop (on nighttime values) 
    ! and thus the end of the program
    EXIT DAY_OR_NIGHT
  ENDIF
ENDDO DAY_OR_NIGHT ! end of the day/night loop that begins in section 1.3

!=================================================================================== 
! 5. CLOSE THE FILES
!===================================================================================
status = nf90_close(gcmfid)
error_text="Error: could not close file "//trim(gcmfile)
call status_check(status,error_text)
status = nf90_close(obsfid)
error_text="Error: could not close file "//trim(obsfile)
call status_check(status,error_text)
status = nf90_close(outfid)
error_text="Error: could not close file "//trim(outfile)
call status_check(status,error_text)

write(*,*)"";write(*,*)"-> Program simu_MCS completed!"

end program simu_MCS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===================================================================================================
! Subroutines
!===================================================================================================

subroutine extraction(lon,lat,alt,sol,&
                  lonlen,latlen,altlen,timelen,&
                  longitude,latitude,altitude,time,&
                  field,missing_value,alttype,varname,value)

implicit none
!================================================================
! Arguments:
!================================================================
real,intent(in) :: lon  ! sought longitude (deg, in [-180:180])
real,intent(in) :: lat  ! sought latitude (deg, in [-90:90])
real,intent(in) :: alt  ! sought altitude (m or Pa)
real,intent(in) :: sol  ! sought date (sols)
integer,intent(in) :: lonlen
integer,intent(in) :: latlen
integer,intent(in) :: altlen
integer,intent(in) :: timelen
real,intent(in) :: longitude(lonlen)
real,intent(in) :: latitude(latlen)
real,intent(in) :: altitude(altlen)
real,intent(in) :: time(timelen)
real,intent(in) :: field(lonlen,latlen,altlen,timelen)
real,intent(in) :: missing_value ! default value for "no data"
character,intent(in) :: alttype ! altitude coord. type:'z' (altitude, m) 
                                !                      'p' (pressure, Pa)
character(len=*),intent(in) :: varname ! variable name (in GCM file)
real,intent(out) :: value

!================================================================
! Local variables:
!================================================================
real,save :: prev_lon=-99 ! previous value of 'lon' routine was called with
real,save :: prev_lat=-99 ! previous value of 'lat' routine was called with
real,save :: prev_alt=-9.e20 ! ! previous value of 'alt'
real,save :: prev_sol=-99 ! previous value of 'sol' routine was called with

! encompasing indexes:
integer,save :: ilon_inf=-1,ilon_sup=-1 ! along longitude
integer,save :: ilat_inf=-1,ilat_sup=-1 ! along latitude
integer,save :: ialt_inf=-1,ialt_sup=-1 ! along altitude
integer,save :: itim_inf=-1,itim_sup=-1 ! along time

! intermediate interpolated values
real :: t_interp(2,2,2) ! after time interpolation
real :: zt_interp(2,2) ! after altitude interpolation
real :: yzt_interp(2) ! after latitude interpolation
real :: coeff ! interpolation coefficient

integer :: i

! By default, set value to missing_value
value=missing_value

!================================================================
! 1. Find encompassing indexes
!================================================================
if (lon.ne.prev_lon) then
  do i=1,lonlen-1
    if (longitude(i).le.lon) then
      ilon_inf=i
    else
      exit
    endif
  enddo
  ilon_sup=ilon_inf+1
endif ! of if (lon.ne.prev_lon)
!write(*,*) 'ilon_inf=',ilon_inf," longitude(ilon_inf)=",longitude(ilon_inf)

if (lat.ne.prev_lat) then
  if (latitude(1).gt.latitude(2)) then
  ! decreasing latitudes, from 90N to (-)90S (LMDZ-like)
    do i=1,latlen-1
      if (latitude(i).ge.lat) then
        ilat_inf=i
      else
        exit
      endif
    enddo
  else
  ! increasing latitudes, from (-)90S to 90N (DYNAMICO-like)
    do i=1,latlen-1
      if (latitude(i).le.lat) then
        ilat_inf=i
      else
        exit
      endif
    enddo
  endif
  ilat_sup=ilat_inf+1
endif ! of if (lat.ne.prev_lat)
!write(*,*) 'ilat_inf=',ilat_inf," latitude(ilat_inf)=",latitude(ilat_inf)

if (alt.ne.prev_alt) then
  if (alttype.eq.'p') then ! pressures are ordered from max to min
    !handle special case for alt not in the altitude(1:altlen) interval
    if ((alt.gt.altitude(1)).or.(alt.lt.altitude(altlen))) then
      ialt_inf=-1
      ialt_sup=-1
      ! return to main program (with value=missing_value)
!      write(*,*)"Problem in extraction : GCM alt p"
      return
    else ! general case
      do i=1,altlen-1
        if (altitude(i).ge.alt) then
          ialt_inf=i
        else
          exit
        endif
      enddo
      ialt_sup=ialt_inf+1
    endif ! of if ((alt.gt.altitude(1)).or.(alt.lt.altitude(altlen)))
  else ! altitudes (m) are ordered from min to max
    !handle special case for alt not in the altitude(1:altlen) interval
    if ((alt.lt.altitude(1)).or.(alt.gt.altitude(altlen))) then
      ialt_inf=-1
      ialt_sup=-1
      ! return to main program (with value=missing_value)
!      write(*,*)"Problem in extraction : GCM alt z"
      return
    else ! general case
      do i=1,altlen-1
        if (altitude(i).le.alt) then
          ialt_inf=i
        else
          exit
        endif
      enddo
      ialt_sup=ialt_inf+1
    endif ! of if ((alt.lt.altitude(1)).or.(alt.gt.altitude(altlen)))
  endif ! of if (alttype.eq.'p') 
endif ! of if (alt.ne.prev_alt)
!write(*,*) 'ialt_inf=',ialt_inf," altitude(ialt_inf)=",altitude(ialt_inf)

if (sol.ne.prev_sol) then
  !handle special case for sol not in the time(1):time(timenlen) interval
  if ((sol.lt.time(1)).or.(sol.gt.time(timelen))) then
    itim_inf=-1
    itim_sup=-1
    ! return to main program (with value=missing_value)
!    write(*,*)"Problem in extraction : GCM sol"
    return
  else ! general case
    do i=1,timelen-1
      if (time(i).le.sol) then
        itim_inf=i
      else
        exit
      endif
    enddo
    itim_sup=itim_inf+1
  endif ! of if ((sol.lt.time(1)).or.(sol.gt.time(timenlen)))
endif ! of if (sol.ne.prev_sol)
!write(*,*) 'itim_inf=',itim_inf," time(itim_inf)=",time(itim_inf)
!write(*,*) 'itim_sup=',itim_sup," time(itim_sup)=",time(itim_sup)

!================================================================
! 2. Interpolate
!================================================================
! check that there are no "missing_value" in the field() elements we need
! otherwise return to main program (with value=missing_value)
if (field(ilon_inf,ilat_inf,ialt_inf,itim_inf).eq.missing_value) then
!  write(*,*)"Error 1 in extraction"
  return
 endif
if (field(ilon_inf,ilat_inf,ialt_inf,itim_sup).eq.missing_value) then
!  write(*,*)"Error 2 in extraction"
  return
 endif
if (field(ilon_inf,ilat_inf,ialt_sup,itim_inf).eq.missing_value) then
!  write(*,*)"Error 3 in extraction"
  return
 endif
if (field(ilon_inf,ilat_inf,ialt_sup,itim_sup).eq.missing_value) then
!  write(*,*)"Error 4 in extraction"
  return
 endif
if (field(ilon_inf,ilat_sup,ialt_inf,itim_inf).eq.missing_value) then
!  write(*,*)"Error 5 in extraction"
  return
 endif
if (field(ilon_inf,ilat_sup,ialt_inf,itim_sup).eq.missing_value) then
!  write(*,*)"Error 6 in extraction"
  return
 endif
if (field(ilon_inf,ilat_sup,ialt_sup,itim_inf).eq.missing_value) then
!  write(*,*)"Error 7 in extraction"
  return
 endif
if (field(ilon_inf,ilat_sup,ialt_sup,itim_sup).eq.missing_value) then
!  write(*,*)"Error 8 in extraction"
  return
 endif
if (field(ilon_sup,ilat_inf,ialt_inf,itim_inf).eq.missing_value) then
!  write(*,*)"Error 9 in extraction"
  return
 endif
if (field(ilon_sup,ilat_inf,ialt_inf,itim_sup).eq.missing_value) then
!  write(*,*)"Error 10 in extraction"
  return
 endif
if (field(ilon_sup,ilat_inf,ialt_sup,itim_inf).eq.missing_value) then
!  write(*,*)"Error 11 in extraction"
  return
 endif
if (field(ilon_sup,ilat_inf,ialt_sup,itim_sup).eq.missing_value) then
!  write(*,*)"Error 12 in extraction"
  return
 endif
if (field(ilon_sup,ilat_sup,ialt_inf,itim_inf).eq.missing_value) then
!  write(*,*)"Error 13 in extraction"
  return
 endif
if (field(ilon_sup,ilat_sup,ialt_inf,itim_sup).eq.missing_value) then
!  write(*,*)"Error 14 in extraction"
  return
 endif
if (field(ilon_sup,ilat_sup,ialt_sup,itim_inf).eq.missing_value) then
!  write(*,*)"Error 15 in extraction"
  return
 endif
if (field(ilon_sup,ilat_sup,ialt_sup,itim_sup).eq.missing_value) then
!  write(*,*)"Error 16 in extraction"
  return
 endif

!================================================================
! 2.1 Linear interpolation in time
!================================================================
coeff=(sol-time(itim_inf))/(time(itim_sup)-time(itim_inf))
t_interp(1,1,1)=field(ilon_inf,ilat_inf,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_inf,ialt_inf,itim_sup)- &
                       field(ilon_inf,ilat_inf,ialt_inf,itim_inf))
t_interp(1,1,2)=field(ilon_inf,ilat_inf,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_inf,ialt_sup,itim_sup)- &
                       field(ilon_inf,ilat_inf,ialt_sup,itim_inf))
t_interp(1,2,1)=field(ilon_inf,ilat_sup,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_sup,ialt_inf,itim_sup)- &
                       field(ilon_inf,ilat_sup,ialt_inf,itim_inf))
t_interp(1,2,2)=field(ilon_inf,ilat_sup,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_sup,ialt_sup,itim_sup)- &
                       field(ilon_inf,ilat_sup,ialt_sup,itim_inf))
t_interp(2,1,1)=field(ilon_sup,ilat_inf,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_inf,ialt_inf,itim_sup)- &
                       field(ilon_sup,ilat_inf,ialt_inf,itim_inf))
t_interp(2,1,2)=field(ilon_sup,ilat_inf,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_inf,ialt_sup,itim_sup)- &
                       field(ilon_sup,ilat_inf,ialt_sup,itim_inf))
t_interp(2,2,1)=field(ilon_sup,ilat_sup,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_sup,ialt_inf,itim_sup)- &
                       field(ilon_sup,ilat_sup,ialt_inf,itim_inf))
t_interp(2,2,2)=field(ilon_sup,ilat_sup,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_sup,ialt_sup,itim_sup)- &
                       field(ilon_sup,ilat_sup,ialt_sup,itim_inf))

!================================================================
! 2.2 Vertical interpolation
!================================================================
if (((varname=='rho').or.(varname=='pressure')).and.(alttype=='z')) then
  ! do the interpolation on the log of the quantity
  coeff=(alt-altitude(ialt_inf))/(altitude(ialt_sup)-altitude(ialt_inf))
  zt_interp(1,1)=log(t_interp(1,1,1))+coeff* &
                             (log(t_interp(1,1,2))-log(t_interp(1,1,1)))
  zt_interp(1,2)=log(t_interp(1,2,1))+coeff* &
                             (log(t_interp(1,2,2))-log(t_interp(1,2,1)))
  zt_interp(2,1)=log(t_interp(2,1,1))+coeff* &
                             (log(t_interp(2,1,2))-log(t_interp(2,1,1)))
  zt_interp(2,2)=log(t_interp(2,2,1))+coeff* &
                             (log(t_interp(2,2,2))-log(t_interp(2,2,1)))
  zt_interp(1:2,1:2)=exp(zt_interp(1:2,1:2))
else ! general case
  coeff=(alt-altitude(ialt_inf))/(altitude(ialt_sup)-altitude(ialt_inf))
  zt_interp(1,1)=t_interp(1,1,1)+coeff*(t_interp(1,1,2)-t_interp(1,1,1))
  zt_interp(1,2)=t_interp(1,2,1)+coeff*(t_interp(1,2,2)-t_interp(1,2,1))
  zt_interp(2,1)=t_interp(2,1,1)+coeff*(t_interp(2,1,2)-t_interp(2,1,1))
  zt_interp(2,2)=t_interp(2,2,1)+coeff*(t_interp(2,2,2)-t_interp(2,2,1))
endif

!================================================================
! 2.3 Latitudinal interpolation
!================================================================
coeff=(lat-latitude(ilat_inf))/(latitude(ilat_sup)-latitude(ilat_inf))
yzt_interp(1)=zt_interp(1,1)+coeff*(zt_interp(1,2)-zt_interp(1,1))
yzt_interp(2)=zt_interp(2,1)+coeff*(zt_interp(2,2)-zt_interp(2,1))

!================================================================
! 2.4 longitudinal interpolation
!================================================================
coeff=(lon-longitude(ilon_inf))/(longitude(ilon_sup)-longitude(ilon_inf))
value=yzt_interp(1)+coeff*(yzt_interp(2)-yzt_interp(1))

end subroutine extraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inidim(outfid,lonlen,latlen,altlen,timelen,lon,lat,alt,time,units_alt,&
                     londimid,latdimid,altdimid,timedimid)
!================================================================
! Purpose:
! Initialize a data file (NetCDF format)
!================================================================
! Remarks:
! The NetCDF file remains open
!================================================================
use netcdf ! NetCDF definitions
implicit none
!================================================================
! Arguments:
!================================================================
integer, intent(in):: outfid
! outfid: [netcdf] file ID
integer, intent(in):: lonlen
! lonlen: longitude length (# of longitude values)
integer, intent(in):: latlen
! latlen: latitude length (# of latitude values)
integer, intent(in):: altlen
! altlen: altitude length (# of altitude values)
integer, intent(in):: timelen
! timelen: time length (# of time values)
real, intent(in):: lon(lonlen)
! lon(): longitude
real, intent(in):: lat(latlen)
! lat(): latitude
real, intent(in):: alt(altlen)
! alt(): altitude
real, intent(in):: time(timelen)
! time(): time (Ls)
character(len=1), intent(in) :: units_alt
! units_alt: altitude coord. type:'z' (altitude, m) 'p' (pressure, Pa)
integer,intent(inout):: londimid
! londimid: [netcdf] lon() (i.e.: longitude) ID in MCS and output files (they are the same)
integer,intent(inout) :: latdimid
! latdimid: [netcdf] lat() ID in MCS and output files (they are the same)
integer,intent(inout):: altdimid
! altdimid: [netcdf] alt() ID in MCS and output files (they are the same)
integer,intent(inout):: timedimid
! timedimid: [netcdf] time() ID in MCS and output files (they are the same)

!================================================================
! Local variables:
!================================================================
integer :: lonvarid,latvarid,altvarid,timevarid,status
! *varid: [netcdf] ID of a variable
! status: [netcdf]  return error code (from called subroutines)
character(len=256) :: error_text
integer :: d ! loop index on dimensions ID

!===============================================================
! 1. Define/write "dimensions" in the same order than MCS file
!    and get their IDs
!================================================================
do d=1,4
  if (altdimid.eq.d) then
    status=nf90_def_dim(outfid, "altitude", altlen, altdimid)
    error_text="Error: could not define the altitude dimension in the outfile"
    call status_check(status,error_text)
  else if (timedimid.eq.d) then
    status=nf90_def_dim(outfid, "time", timelen, timedimid)
    error_text="Error: could not define the time dimension in the outfile"
    call status_check(status,error_text)
  else if (latdimid.eq.d) then
    status=nf90_def_dim(outfid, "latitude", latlen, latdimid)
    error_text="Error: could not define the latitude dimension in the outfile"
    call status_check(status,error_text)
  else if (londimid.eq.d) then
    status=nf90_def_dim(outfid, "longitude", lonlen, londimid)
    error_text="Error: could not define the longitude dimension in the outfile"
    call status_check(status,error_text)
  endif
enddo

!================================================================
! 2. Define "variables" and their attributes
!================================================================
!================================================================
! 2.1 Write "longitude" (data and attributes)
!================================================================

! Insert the definition of the variable
status=nf90_def_var(outfid,"longitude",nf90_float,(/londimid/),lonvarid)
error_text="Error: could not define the longitude variable in the outfile"
call status_check(status,error_text)

! Write the attributes
status=nf90_put_att(outfid,lonvarid,"long_name","longitude")
error_text="Error: could not put the long_name attribute for the longitude variable in the outfile"
call status_check(status,error_text)

status=nf90_put_att(outfid,lonvarid,"units","degrees_east")
error_text="Error: could not put the units attribute for the longitude variable in the outfile"
call status_check(status,error_text)

!================================================================
! 2.2 "latitude"
!================================================================

! Insert the definition of the variable
status=nf90_def_var(outfid,"latitude",nf90_float,(/latdimid/),latvarid)
error_text="Error: could not define the latitude variable in the outfile"
call status_check(status,error_text)

! Write the attributes
status=nf90_put_att(outfid,latvarid,"long_name","latitude")
error_text="Error: could not put the long_name attribute for the latitude variable in the outfile"
call status_check(status,error_text)

status=nf90_put_att(outfid,latvarid,"units","degrees_north")
error_text="Error: could not put the units attribute for the latitude variable in the outfile"
call status_check(status,error_text)

!================================================================
! 2.3 Write "altitude" (data and attributes)
!================================================================

! Insert the definition of the variable
status=nf90_def_var(outfid,"altitude",nf90_float,(/altdimid/),altvarid)
error_text="Error: could not define the altitude variable in the outfile"
call status_check(status,error_text)

! Write the attributes
if (units_alt.eq.'p') then ! pressure coordinate
  status=nf90_put_att(outfid,altvarid,"long_name","pressure")
  error_text="Error: could not put the long_name attribute for the altitude variable in the outfile"
  call status_check(status,error_text)

  status=nf90_put_att(outfid,altvarid,'units',"Pa")
  error_text="Error: could not put the units attribute for the altitude variable in the outfile"
  call status_check(status,error_text)

  status=nf90_put_att(outfid,altvarid,'positive',"down")
  error_text="Error: could not put the positive attribute for the altitude variable in the outfile"
  call status_check(status,error_text)

  status=nf90_put_att(outfid,altvarid,'comment',&
  "The vertical variable is in fact pressure, not altitude. We just call it altitude for easy reading in Ferret and Grads.")
  error_text="Error: could not put the comment attribute for the altitude variable in the outfile"
  call status_check(status,error_text)

else if (units_alt.eq.'z') then ! altitude coordinate
  status=nf90_put_att(outfid,altvarid,"long_name","altitude")
  error_text="Error: could not put the long_name attribute for the altitude variable in the outfile"
  call status_check(status,error_text)

  status=nf90_put_att(outfid,altvarid,'units',"m")
  error_text="Error: could not put the units attribute for the altitude variable in the outfile"
  call status_check(status,error_text)

  status=nf90_put_att(outfid,altvarid,'positive',"up")
  error_text="Error: could not put the positive attribute for the altitude variable in the outfile"
  call status_check(status,error_text)

else
  write(*,*)"I do not understand this unit type ",units_alt," for outfile altitude!"
  stop
end if

!================================================================
! 2.4 "time"
!================================================================

! Insert the definition of the variable
status=nf90_def_var(outfid,"time",nf90_float,(/timedimid/),timevarid)
error_text="Error: could not define the time variable in the outfile"
call status_check(status,error_text)

! Write the attributes
status=nf90_put_att(outfid,timevarid,"long_name","Solar longitude")
error_text="Error: could not put the long_name attribute for the time variable in the outfile"
call status_check(status,error_text)

status=nf90_put_att(outfid,timevarid,"units","days since 0000-01-1 00:00:00")
error_text="Error: could not put the units attribute for the time variable in the outfile"
call status_check(status,error_text)

status=nf90_put_att(outfid,timevarid,"comment",&
"Units is in fact degrees, but set to a dummy value of days for compatibility with Ferret and Grads.")
error_text="Error: could not put the comment attribute for the time variable in the outfile"
call status_check(status,error_text)

!================================================================
! 2.5 End netcdf define mode
!================================================================
status=nf90_enddef(outfid)
error_text="Error: could not end the define mode of the outfile"
call status_check(status,error_text)

!================================================================
! 3. Write "variables" (data of the dimension variables)
!================================================================
! "time"
status=nf90_put_var(outfid,timevarid,time)
error_text="Error: could not write the time variable's data in the outfile"
call status_check(status,error_text)

! "latitude"
status=nf90_put_var(outfid,latvarid,lat)
error_text="Error: could not write the latitude variable's data in the outfile"
call status_check(status,error_text)

! "longitude"
status=nf90_put_var(outfid,lonvarid,lon)
error_text="Error: could not write the longitude variable's data in the outfile"
call status_check(status,error_text)

! "altitude"
status=nf90_put_var(outfid,altvarid,alt)
error_text="Error: could not write the altitude variable's data in the outfile"
call status_check(status,error_text)

end subroutine inidim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ls2sol(ls,sol)

implicit none
!================================================================
!  Arguments:
!================================================================
real,intent(in) :: ls
real,intent(out) :: sol

!================================================================
!  Local:
!================================================================
double precision xref,zx0,zteta,zz
!xref: mean anomaly, zteta: true anomaly, zx0: eccentric anomaly
double precision year_day 
double precision peri_day,timeperi,e_elips
double precision pi,degrad 
parameter (year_day=668.6d0) ! number of sols in a martian year
parameter (peri_day=485.35d0) ! date (in sols) of perihelion
!timeperi: 2*pi*( 1 - Ls(perihelion)/ 360 ); Ls(perihelion)=250.99
parameter (timeperi=1.90258341759902d0)
parameter (e_elips=0.0934d0)  ! eccentricity of orbit
parameter (pi=3.14159265358979d0)
parameter (degrad=57.2957795130823d0)

      if (abs(ls).lt.1.0e-5) then
         if (ls.ge.0.0) then
            sol = 0.0
         else
            sol = real(year_day)
         end if
         return
      end if

      zteta = ls/degrad + timeperi
      zx0 = 2.0*datan(dtan(0.5*zteta)/dsqrt((1.+e_elips)/(1.-e_elips)))
      xref = zx0-e_elips*dsin(zx0)
      zz = xref/(2.*pi)
      sol = real(zz*year_day + peri_day)
      if (sol.lt.0.0) sol = sol + real(year_day)
      if (sol.ge.year_day) sol = sol - real(year_day)


end subroutine ls2sol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_sol_list(solcount,int_sol_list,LTcount,LTave,LTmax,LTmin,lon,f_LT,dayornight,&
                        sol_list)
!================================================================
! Purpose:
! Generate a list that is the combination of two lists :
! - int_sol_list, which is a list of integer values of sol
! - LT_list, which contains a number LTcount of LT values in the
!   interval [LTmin;LTmax] which are evenly distributed around
!   LTave (so that their average is LTave)
!================================================================

implicit none
!================================================================
! Arguments:
!================================================================
integer, intent(in) :: solcount ! nb of integer values of sol
real, intent(in) :: int_sol_list(solcount) ! list of the integer values of sol
integer, intent(in) :: LTcount ! nb of LT per sol to be interpolated
real, intent(in) :: LTave ! average of LT
real, intent(in) :: LTmax, LTmin ! bounds of the LT interval for the interpolation
real, intent(in) :: lon ! longitude value for the transformation into a sol value at lon=0°
external f_LT ! function that changes LT interval for night LT
real f_LT
character (len=10), intent(in) :: dayornight ! to know if we have day or night values

real, intent(out) :: sol_list(solcount*LTcount) ! all the sol values at lon=0° to interpolate

!================================================================
! Local variables:
!================================================================
integer :: N
integer :: a,b,c ! loop iteration indexes
real :: LT_list(LTcount)

N = floor(LTcount/2.)

!================================================================
! 1. Filling of LT_list
!================================================================
SELECT CASE (dayornight)
CASE ("dayside")
  if (abs(LTave-LTmax).le.abs(LTave-LTmin)) then ! LTave is closer to LTmax than to LTmin
  do a=1,N
    LT_list(a)=LTave+a*abs(LTave-LTmax)/(N+1)
    LT_list(N+a)=LTave-a*abs(LTave-LTmax)/(N+1)
  enddo
else ! LTave is closer to LTmin than to LTmax
  do a=1,N
    LT_list(a)=LTave+a*abs(LTave-LTmin)/(N+1)
    LT_list(N+a)=LTave-a*abs(LTave-LTmin)/(N+1)
  enddo
endif
CASE ("nightside")
  if (abs(f_LT(LTave)-f_LT(LTmax)).le.abs(f_LT(LTave)-f_LT(LTmin))) then ! LTave is closer to LTmax than to LTmin
    do a=1,N
      LT_list(a)=LTave+a*abs(f_LT(LTave)-f_LT(LTmax))/(N+1)
      LT_list(a)=mod(LT_list(a),24.) ! reput the LT in a [0;24[ interval if needed
      LT_list(N+a)=LTave-a*abs(f_LT(LTave)-f_LT(LTmax))/(N+1)
      LT_list(N+a)=mod(LT_list(N+a),24.)
    enddo
  else ! LTave is closer to LTmin than to LTmax
    do a=1,N
      LT_list(a)=LTave+a*abs(f_LT(LTave)-f_LT(LTmin))/(N+1)
      LT_list(a)=mod(LT_list(a),24.)
      LT_list(N+a)=LTave-a*abs(f_LT(LTave)-f_LT(LTmin))/(N+1)
      LT_list(N+a)=mod(LT_list(N+a),24.)
    enddo
  endif
END SELECT

if (mod(LTcount,2).eq.1) then ! if LTcount is an odd number
  LT_list(LTcount)=LTave ! add LTave to the list
endif

!================================================================
! 2. Combination of int_sol_list & LT_list into sol_list
!================================================================
c=1
do a=1,solcount
  do b=1,LTcount
    sol_list(c)=int_sol_list(a)+(LT_list(b)-lon/15.)/24.
    c=c+1
  enddo
enddo

end subroutine gen_sol_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine status_check(status,error_text)

use netcdf
implicit none
!================================================================
!  Arguments:
!================================================================
integer,intent(in) :: status
character(len=256),intent(in) :: error_text

if (status.ne.nf90_noerr) then
  write(*,*)trim(error_text)
  write(*,*)trim(nf90_strerror(status))
  stop
endif

end subroutine status_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function LTmod(LT)
!================================================================
! Purpose:
! For night local times management, replace hours  
! from a [0;24[ interval to a [-12;12[ interval in which
! midnight = 0 (in order to ensure continuity when comparing
! night local times)
!================================================================
implicit none
real,intent(in) :: LT

LTmod = 2*mod(LT,12.)-LT
return
end function LTmod
