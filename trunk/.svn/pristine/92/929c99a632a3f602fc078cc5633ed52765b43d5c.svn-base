program aeroptical
! program for computation of aerosol opacities
! author : Antoine Bierjon, April-May 2020
!
!===============================================================================
!     PREFACE
!===============================================================================
! This program takes as inputs a GCM output file (diagfi,stats,concat) that
! contains:
!      - the Mass Mixing Ratios of dust (dustq), water ice (h2o_ice),
!                               stormdust (rdsdustq) & topdust (topdustq)
!      - their effective radius (reffdust, reffice(*),
!                           reffstormdust, refftopdust)
!      - the atmospheric density (rho)
! as well as a wavelength to calibrate the opacities, and the directory
! containing the ASCII files of the aerosols' optical properties.
!
! It outputs a netcdf file containing the opacities [1/km] of the dust,
! water ice and stormdust aerosols calibrated at the prescribed wavelength.
! The type of opacity (extinction or absorption) is given by the user.
!
! (*) due to a high uncertainty of reffice in the GCM, the value is asked
! directly to the user (if the user returns 0, then the program reads the GCM
! file to get reffice)
!
! NOTA BENE:
! 1) if one wanted to add another aerosol to compute, one should look for
! the comments + NEW AER? that are disseminated all along the code to indicate
! the parts of the code that should be modified.
! 2) another enhancement of this program could be the possibility to read its
! own product files and recalibrate the opacities at another wavelength
!===============================================================================


use netcdf

!===============================================================================
! Variable declarations
!===============================================================================

implicit none ! for no implicitly typed variables


! GCM file
character(len=128) :: gcmfile                 ! name of the netcdf GCM input file
integer :: gcmfid                             ! [netcdf] GCM input file ID number
integer :: ierr                               ! [netcdf] subroutine returned error code
character(len=256) :: error_text              ! text to display in case of error
integer :: lonlen,latlen,altlen,timelen       ! nb of grid points along longitude,latitude,altitude,time
integer :: GCM_layers                         ! number of GCM layers
integer :: layerdimout,interlayerdimout       ! "GCM_layers" and "GCM_layers+1" IDs

logical,dimension(:),allocatable :: aerok     ! to know if the needed fields are in the GCM file
integer,dimension(:),allocatable :: mmrvarid  ! stores a MMR variable ID number
integer,dimension(:),allocatable :: reffvarid ! stores a reff variable ID number
integer :: tmpvarid                           ! temporarily stores a variable ID number
real,dimension(:,:,:,:),allocatable :: mmr    ! aerosols mass mixing ratios [kg/kg]
real,dimension(:,:,:,:),allocatable :: reff   ! aerosols effective radii [m]
real :: reffwice_val                          ! value given by the user for reffwice (0 if read in the GCM file) [m]
real,dimension(:,:,:,:),allocatable :: rho    ! atmospheric density [kg.m-3]
integer :: naerkind                           ! nb of aerosols to consider
integer :: iaer                               ! aerosol kind index (1=dust,2=water ice,3=stormdust,4=topdust) ! + NEW AER?
integer :: ilon,ilat,ialt,it                  ! Loop indices for coordinates


! Output file
character(len=128) :: outfile                  ! name of the netcdf output file
integer :: outfid                              ! [netcdf] file ID number
integer :: latdimout,londimout,altdimout,timedimout
                                               ! latdimout: stores latitude dimension ID number
                                               ! londimout: stores longitude dimension ID number
                                               ! altdimout: stores altitude dimension ID number
                                               ! timedimout: stores time dimension ID number
integer :: tmpvaridout                         ! temporarily stores a variable ID number

real :: wvl_val                                ! reference wavelength of the output opacities (given by the user) [m]
character(len=3) :: opatype                    ! opacity type : extinction (ext) or absorption (abs) (given by the user)
integer :: varshape(4)                         ! stores a variable shape (of dimensions' IDs)
character(len=16) :: tmpvarname                ! temporarily stores a variable name
real,dimension(:,:,:,:),allocatable :: opa_aer ! Aerosols opacities [1/km]
real,dimension(:),allocatable :: missval       ! Value to put in outfile when the reff is out of bounds
                                               ! or when there is a mising value in input file


! Optical properties (read in external ASCII files)
character(len=256) :: datadir              ! directory containing the ASCII files
character(len=128) :: optpropfile          ! name of the ASCII optical properties file
logical :: file_ok                         ! to know if the file can be opened
integer :: file_unit = 60                  ! ASCII file ID

integer :: jfile                           ! ASCII file scan index
logical :: endwhile                        ! Reading loop boolean
character(len=132) :: scanline             ! ASCII file scanning line
integer :: read_ok                         ! to know if the line is well read
                   
integer :: nwvl                            ! Number of wavelengths in the domain (VIS or IR)
integer :: nsize                           ! Number of particle sizes stored in the file
integer :: iwvl,isize                      ! Wavelength and Particle size loop indices
real,dimension(:),allocatable :: wvl       ! Wavelength axis [m]
real,dimension(:),allocatable :: radiusdyn ! Particle size axis [m]
real,dimension(:,:),allocatable :: Qext    ! Extinction efficiency coefficient [/]
real,dimension(:,:),allocatable :: omeg    ! Single scattering albedo [/]


! Opacity computation
integer :: iwvl1,iwvl2,isize1,isize2     ! Wavelength and Particle size indices for the interpolation
real :: coeff                            ! Interpolation coefficient
real,dimension(2) :: reffQext            ! Qext after reff interpolation
real :: Qext_val                         ! Qext after both interpolations
real,dimension(2) :: reffomeg            ! omeg after reff interpolation
real :: omeg_val                         ! omeg after both interpolations
real,dimension(:),allocatable :: rho_aer ! Density of the aerosols [kg.m-3]

!===============================================================================
! 0. USER INPUTS
!===============================================================================
write(*,*) "Welcome in the aerosol opacities' computation program !"
write(*,*)

! Ask the GCM file name
WRITE(*,*) "Enter an input file name (GCM diagfi/stats/concat) :"
READ(*,*) gcmfile

! Ask the reffwice to the user
write(*,*)""
write(*,*) "The water ice effective radius computed by the GCM is known to be a bit inaccurate."
write(*,*) "Which value do you want to use for it (in meters) ?"
write(*,*) "(put 0 if you still want the program to read the value in "//trim(gcmfile)//")"
READ(*,*) reffwice_val

! Ask the wavelength to the user
write(*,*)""
WRITE(*,*) "Value of the reference wavelength for the opacities' computation (in meters) ?"
READ(*,*) wvl_val

! Ask the directory containing the optical properties files
write(*,*)""
WRITE(*,*) "In which directory should we look for the optical properties' files ?"
READ(*,'(a)') datadir

! Ask the type of opacity that has to be computed
write(*,*)""
WRITE(*,*) "Do you want to compute extinction or absorption opacities ? (ext/abs)"
READ(*,*) opatype
if ((trim(opatype).ne."ext").and.(trim(opatype).ne."abs")) then
  write(*,*) "opatype = "//opatype
  write(*,*) "Error: the opacity type must be ext or abs"
  stop
endif

!===============================================================================
! 1. GCM FILE & OUTPUT FILE INITIALIZATIONS
!===============================================================================

!==========================================================================
! 1.1 GCM file opening & dimensions - Output file initializations
!==========================================================================

!==========================================================================
! --> 1.1.1 Open the netcdf GCM file given by the user

ierr=nf90_open(gcmfile,nf90_nowrite,gcmfid) ! nowrite mode=the program can only read the opened file
error_text="Error: could not open file "//trim(gcmfile)
call status_check(ierr,error_text)

!==========================================================================
! --> 1.1.2 Creation of the outfile
if (opatype.eq."ext") then
  outfile=gcmfile(1:index(gcmfile, ".nc")-1)//"_OPAext.nc"
else ! opatype.eq."abs"
  outfile=gcmfile(1:index(gcmfile, ".nc")-1)//"_OPAabs.nc"
endif

 
ierr=NF90_CREATE(outfile,IOR(nf90_clobber,nf90_64bit_offset),outfid)
! NB: clobber mode=overwrite any dataset with the same file name !
! 64bit_offset enables the creation of large netcdf files with variables>2GB
error_text="Error: could not create file "//trim(outfile)
call status_check(ierr,error_text)
write(*,*)"";WRITE(*,*)"Output file is: ",trim(outfile);write(*,*)""

!==========================================================================
! --> 1.1.3 Get the dimensions and create them in the output file

! Initialize output file's lat,lon,alt,time and controle dimensions
call inidims(gcmfile,gcmfid,outfile,outfid,&
             lonlen,latlen,altlen,timelen,&
             latdimout,londimout,altdimout,timedimout,&
             GCM_layers,layerdimout,interlayerdimout)

! Initialize output file's aps,bps,ap,bp and phisinit variables
call init2(gcmfid,lonlen,latlen,altlen,GCM_layers,&
           outfid,londimout,latdimout,altdimout,&
           layerdimout,interlayerdimout)

!==========================================================================
! 1.2 GCM aerosols
!==========================================================================

! Number of aerosols
naerkind = 4 ! dust, water ice, stormdust, topdust ! + NEW AER?

! Variables length allocation
allocate(mmrvarid(naerkind))
allocate(reffvarid(naerkind))

! Initialize missing value
allocate(missval(naerkind))
missval(:)=1.e+20

! Initialize aerok to .true. for all aerosols
allocate(aerok(naerkind))
aerok(:)=.true.

!==========================================================================
! --> 1.2.1 DUST

! DUST MASS MIXING RATIO
! Check that the GCM file contains that variable
ierr=nf90_inq_varid(gcmfid,"dustq",mmrvarid(1))
error_text="Failed to find variable dustq in "//trim(gcmfile)&
            //" We'll skip the dust aerosol."
if (ierr.ne.nf90_noerr) then
  write(*,*)trim(error_text)
  aerok(1)=.false. 
endif

! DUST EFFECTIVE RADIUS
if (aerok(1)) then
  ! Check that the GCM file contains that variable
  ierr=nf90_inq_varid(gcmfid,"reffdust",reffvarid(1))
  error_text="Failed to find variable reffdust in "//trim(gcmfile)&
              //" We'll skip the dust aerosol."
  if (ierr.ne.nf90_noerr) then
    write(*,*)trim(error_text)
    aerok(1)=.false.
  endif
endif

!==========================================================================
! --> 1.2.2 WATER ICE 

! WATER ICE MASS MIXING RATIO
! Check that the GCM file contains that variable
ierr=nf90_inq_varid(gcmfid,"h2o_ice",mmrvarid(2))
error_text="Failed to find variable h2o_ice in "//trim(gcmfile)&
            //" We'll skip the water ice aerosol."
if (ierr.ne.nf90_noerr) then
  write(*,*)trim(error_text)
  aerok(2)=.false.
endif

! WATER ICE EFFECTIVE RADIUS
if (aerok(2)) then  
  IF (reffwice_val.eq.0) THEN
    ! Check that the GCM file contains that variable
    ierr=nf90_inq_varid(gcmfid,"reffice",reffvarid(2))
    error_text="Failed to find variable reffice in "//trim(gcmfile)&
                //" We'll skip the water ice aerosol."
    if (ierr.ne.nf90_noerr) then
      write(*,*)trim(error_text)
      aerok(2)=.false.
    endif
  ENDIF
endif

!==========================================================================
! --> 1.2.3 STORMDUST

! STORMDUST MASS MIXING RATIO
! Check that the GCM file contains that variable
ierr=nf90_inq_varid(gcmfid,"rdsdustq",mmrvarid(3))
error_text="Failed to find variable rdsdustq in "//trim(gcmfile)&
            //" We'll skip the stormdust aerosol."
if (ierr.ne.nf90_noerr) then
  write(*,*)trim(error_text)
  aerok(3)=.false.
endif

! STORMDUST EFFECTIVE RADIUS
if (aerok(3)) then
  ! Check that the GCM file contains that variable
  ierr=nf90_inq_varid(gcmfid,"reffstormdust",reffvarid(3))
  error_text="Failed to find variable reffstormdust in "//trim(gcmfile)&
              //" We'll skip the stormdust aerosol."
  if (ierr.ne.nf90_noerr) then
    write(*,*)trim(error_text)
    aerok(3)=.false.
  endif
endif 

!==========================================================================
! --> 1.2.4 TOPDUST

! TOPDUST MASS MIXING RATIO
! Check that the GCM file contains that variable
ierr=nf90_inq_varid(gcmfid,"topdustq",mmrvarid(4))
error_text="Failed to find variable topdustq in "//trim(gcmfile)&
            //" We'll skip the topdust aerosol."
if (ierr.ne.nf90_noerr) then
  write(*,*)trim(error_text)
  aerok(4)=.false.
endif

! TOPDUST EFFECTIVE RADIUS
if (aerok(4)) then
  ! Check that the GCM file contains that variable
  ierr=nf90_inq_varid(gcmfid,"refftopdust",reffvarid(4))
  error_text="Failed to find variable refftopdust in "//trim(gcmfile)&
              //" We'll skip the topdust aerosol."
  if (ierr.ne.nf90_noerr) then
    write(*,*)trim(error_text)
    aerok(4)=.false.
  endif
endif

!==========================================================================
! --> 1.2.5 + NEW AER? 


! Check if there is still at least one aerosol to compute
IF (.NOT.ANY(aerok)) THEN ! + NEW AER?
  write(*,*) "No aerosol among [dust/water ice/stormdust/topdust] was found in the file. Better stop now..."
  stop
ENDIF

!==========================================================================
! 1.3 GCM ATMOSPHERIC DENSITY
!==========================================================================

! Check that the GCM file contains that variable
ierr=nf90_inq_varid(gcmfid,"rho",tmpvarid)
error_text="Failed to find variable rho in "//trim(gcmfile)&
            //" We need it to compute the opacity [1/km]."
call status_check(ierr,error_text)
! Length allocation for each dimension of the 4D variable
allocate(rho(lonlen,latlen,altlen,timelen))
! Load dataset
ierr=nf90_get_var(gcmfid,tmpvarid,rho)
error_text="Failed to load atmospheric density"
call status_check(ierr,error_text)
write(*,*) "Atmospheric density rho loaded from "//trim(gcmfile)


!==========================================================================
! 1.4 Output variable's initializations and datasets loading
!==========================================================================
! Define the ID shape of the output variables
varshape(1)=londimout
varshape(2)=latdimout
varshape(3)=altdimout
varshape(4)=timedimout

! Fill the aerosol density array
allocate(rho_aer(naerkind))
rho_aer(1)=2500. ! dust
rho_aer(2)=920.  ! water ice
rho_aer(3)=2500. ! stormdust
rho_aer(4)=2500. ! topdust
! + NEW AER?

DO iaer = 1, naerkind ! Loop on aerosol kind
  if (aerok(iaer)) then ! check if this aerosol can be computed
    write(*,*)""
    SELECT CASE(iaer)
    CASE(1) ! DUST
      write(*,*)"Computing dust opacity..."
    CASE(2) ! WATER ICE
      write(*,*)"Computing water ice opacity..."
    CASE(3) ! STORMDUST
      write(*,*)"Computing stormdust opacity..."
    CASE(4) ! TOPDUST
      write(*,*)"Computing topdust opacity..."
! + NEW AER?
    END SELECT ! iaer
    
    ! Length allocation of the input variables
    ALLOCATE(mmr(lonlen,latlen,altlen,timelen))
    ALLOCATE(reff(lonlen,latlen,altlen,timelen))
    
    ! Datasets loading
     ! MMR
    ierr=NF90_GET_VAR(gcmfid,mmrvarid(iaer),mmr(:,:,:,:))
    error_text="Failed to load mass mixing ratio"
    call status_check(ierr,error_text)
    write(*,*) "Mass mixing ratio loaded from "//trim(gcmfile)
    ! Get missing value
    ierr=nf90_get_att(gcmfid,mmrvarid(1),"missing_value",missval(1))
    if (ierr.ne.nf90_noerr) then
      missval(1) = 1.e+20
    endif
    
     ! REFF
    if (iaer.ne.2) then
      ! Load dataset
      ierr=NF90_GET_VAR(gcmfid,reffvarid(iaer),reff(:,:,:,:))
      error_text="Failed to load effective radius"
      call status_check(ierr,error_text)
      write(*,*) "Effective radius loaded from "//trim(gcmfile)
    
    else ! water ice special case
      IF (reffwice_val.eq.0) THEN
        ! Load dataset
        ierr=NF90_GET_VAR(gcmfid,reffvarid(iaer),reff(:,:,1,:)) ! reffice is actually a GCM 
                                                              ! surface (column-averaged) variable
        error_text="Failed to load effective radius"
        call status_check(ierr,error_text)
        do ialt=2,altlen
          reff(:,:,ialt,:) = reff(:,:,1,:) ! build up along altitude axis
        enddo
        write(*,*) "Effective radius loaded from "//trim(gcmfile)
      ELSE ! if reffwice_val/=0
        reff(:,:,:,:) = reffwice_val
        write(*,*) "Effective radius loaded from the user input"
      ENDIF
      
    endif ! iaer.ne.2
    
  
    ! Length allocation of the output variable
    ALLOCATE(opa_aer(lonlen,latlen,altlen,timelen))
    
  
!===============================================================================
! 2. OPTICAL PROPERTIES
!===============================================================================
!==========================================================================
! 2.1 Open the ASCII file
!==========================================================================
    IF (wvl_val.lt.5.e-6) THEN
      ! "VISIBLE" DOMAIN
      
      SELECT CASE(iaer)
      CASE(1) ! DUST
        ! Dust file
        optpropfile = "optprop_dustvis_TM_n50.dat"
      CASE(2) ! WATER ICE
        ! Water ice file
        optpropfile = "optprop_icevis_n30.dat"
      CASE(3) ! STORMDUST
        ! Dust file
        optpropfile = "optprop_dustvis_TM_n50.dat"
      CASE(4) ! TOPDUST
        ! Dust file
        optpropfile = "optprop_dustvis_TM_n50.dat"
! + NEW AER?
      END SELECT ! iaer
      
    ELSE ! wvl_val.ge.5.e-6
      ! "INFRARED" DOMAIN
      
      SELECT CASE(iaer)
      CASE(1) ! DUST
        ! Dust file
        optpropfile = "optprop_dustir_n50.dat"
      CASE(2) ! WATER ICE
        ! Water ice file
        optpropfile = "optprop_iceir_n30.dat"
      CASE(3) ! STORMDUST
        ! Dust file
        optpropfile = "optprop_dustir_n50.dat"
      CASE(4) ! TOPDUST
        ! Dust file
        optpropfile = "optprop_dustir_n50.dat"
! + NEW AER?
      END SELECT ! iaer
      
    ENDIF ! wvl_val

    INQUIRE(FILE=trim(datadir)//'/'//trim(optpropfile),EXIST=file_ok)
    if(.not.file_ok) then
      write(*,*)'Problem opening ',trim(optpropfile)
      write(*,*)'It should be in: ',trim(datadir)
      stop 
    endif
    OPEN(UNIT=file_unit,FILE=trim(datadir)//'/'//trim(optpropfile),FORM='formatted')

!==========================================================================
! 2.2 Allocate the optical property table
!==========================================================================
    jfile = 1
    endwhile = .false.
    DO WHILE (.NOT.endwhile)
      READ(file_unit,*,iostat=read_ok) scanline
      if (read_ok.ne.0) then
        write(*,*)'Error reading file ',&
                  trim(datadir)//'/'//trim(optpropfile)
        stop
      endif
      IF ((scanline(1:1).ne.'#').and.(scanline(1:1).ne.' ')) THEN
        BACKSPACE(file_unit)
reading1_seq: SELECT CASE (jfile) ! FIRST READING SEQUENCE
        CASE(1) reading1_seq ! nwvl ----------------------------
            read(file_unit,*,iostat=read_ok) nwvl
            if (read_ok.ne.0) then
              write(*,*)'reading1_seq: Error while reading line: ',&
                        trim(scanline)
              write(*,*)'   of file ',&
                        trim(datadir)//'/'//trim(optpropfile)
              stop
            endif
            jfile = jfile+1
        CASE(2) reading1_seq ! nsize ---------------------------
            read(file_unit,*,iostat=read_ok) nsize
            if (read_ok.ne.0) then
              write(*,*)'reading1_seq: Error while reading line: ',&
                        trim(scanline)
              write(*,*)'   of file ',&
                        trim(datadir)//'/'//trim(optpropfile)
              stop
            endif
            endwhile = .true.
        CASE DEFAULT reading1_seq ! default --------------------
            write(*,*) 'reading1_seq: ',& 
                       'Error while loading optical properties.' 
            stop 
        END SELECT reading1_seq ! ==============================
      ENDIF
    ENDDO ! DO WHILE(.NOT.endwhile)

    ALLOCATE(wvl(nwvl))             ! Wavelength axis
    ALLOCATE(radiusdyn(nsize))      ! Aerosol radius axis
    ALLOCATE(Qext(nwvl,nsize))      ! Extinction efficiency coeff
    ALLOCATE(omeg(nwvl,nsize))      ! Single scattering albedo

!==========================================================================
! 2.3 Read the data
!==========================================================================
    jfile = 1
    endwhile = .false.
    DO WHILE (.NOT.endwhile)
      READ(file_unit,*) scanline
      IF ((scanline(1:1).ne.'#').and.(scanline(1:1).ne.' ')) THEN
        BACKSPACE(file_unit)
reading2_seq: SELECT CASE (jfile) ! SECOND READING SEQUENCE
        CASE(1) reading2_seq ! wvl -----------------------------
            read(file_unit,*,iostat=read_ok) wvl
            if (read_ok.ne.0) then
              write(*,*)'reading2_seq: Error while reading line: ',&
                        trim(scanline)
              write(*,*)'   of file ',&
                        trim(datadir)//'/'//trim(optpropfile)
              stop
            endif
            jfile = jfile+1
        CASE(2) reading2_seq ! radiusdyn -----------------------
            read(file_unit,*,iostat=read_ok) radiusdyn
            if (read_ok.ne.0) then
              write(*,*)'reading2_seq: Error while reading line: ',&
                        trim(scanline)
              write(*,*)'   of file ',&
                        trim(datadir)//'/'//trim(optpropfile)
              stop
            endif
            jfile = jfile+1
        CASE(3) reading2_seq ! Qext ----------------------------
            isize = 1
            DO WHILE (isize .le. nsize)
              READ(file_unit,*,iostat=read_ok) scanline
              if (read_ok.ne.0) then
                write(*,*)'reading2_seq: Error while reading line: ',&
                          trim(scanline)
                write(*,*)'   of file ',&
                          trim(datadir)//'/'//trim(optpropfile)
                stop
              endif
              IF ((scanline(1:1).ne.'#').and.(scanline(1:1).ne.' ')) THEN
                BACKSPACE(file_unit)
                read(file_unit,*) Qext(:,isize)
                isize = isize + 1
              ENDIF
            ENDDO
            jfile = jfile+1
        CASE(4) reading2_seq ! omeg ----------------------------
            isize = 1
            DO WHILE (isize .le. nsize)
              READ(file_unit,*,iostat=read_ok) scanline
              if (read_ok.ne.0) then
                write(*,*)'reading2_seq: Error while reading line: ',&
                          trim(scanline)
                write(*,*)'   of file ',&
                          trim(datadir)//'/'//trim(optpropfile)
                stop
              endif
              IF ((scanline(1:1).ne.'#').and.(scanline(1:1).ne.' ')) THEN
                BACKSPACE(file_unit)
                read(file_unit,*) omeg(:,isize)
                isize = isize + 1
              ENDIF
            ENDDO
            endwhile = .true.
        CASE DEFAULT reading2_seq ! default --------------------
            write(*,*) 'reading2_seq: ',&
                       'Error while loading optical properties.'
            stop
        END SELECT reading2_seq ! ==============================
      ENDIF
    ENDDO

    ! Close the file
    CLOSE(file_unit)
    
    write(*,*)""
    write(*,*) "Wavelength array loaded from file ",trim(datadir)//'/'//trim(optpropfile)
    write(*,*) "ranging from ",wvl(1)," to ",wvl(nwvl)," meters"
    write(*,*) "Effective radius array loaded from file ",trim(datadir)//'/'//trim(optpropfile)
    write(*,*) "ranging from ",radiusdyn(1)," to ",radiusdyn(nsize)," meters"

! + NEW AER? : one may adapt this part to handle the properties' file of the new aerosol

!==========================================================================
! 2.4 Get the optpropfile wavelength values encompassing the ref wavelength
!==========================================================================
    iwvl=1
    DO WHILE ((iwvl.le.nwvl).and.(wvl(iwvl).lt.wvl_val))
      iwvl=iwvl+1
    ENDDO
    if ((iwvl.gt.nwvl).or.(wvl(1).gt.wvl_val)) then
      write(*,*) "ERROR: the reference wavelength is out of the bounds"
      write(*,*) "of the file ",trim(datadir)//'/'//trim(optpropfile)
      write(*,*) "(wvl_inf=",wvl(1),"m ; wvl_sup=",wvl(nwvl),"m)"
      write(*,*) "Ensure you wrote it well (in meters),"
      write(*,*) "or supply a new optical properties' file (change in aeroptical.F90 directly)"
      stop
    endif
    if (wvl(iwvl).eq.wvl_val) then
      ! if the ref wavelength is in the optpropfile
      iwvl1 = iwvl
      iwvl2 = iwvl
    else ! wvl(iwvl)>wvl_val and wvl(iwvl-1)<wvl_val
      iwvl1 = iwvl-1
      iwvl2 = iwvl
    endif


!===============================================================================
! 3. OUTPUT FILE VARIABLES
!===============================================================================
!==========================================================================
! 3.1 Creation of the output individual aerosol opacities
!==========================================================================
    ! Switch to netcdf define mode
    ierr=nf90_redef(outfid)
    error_text="ERROR: Couldn't start netcdf define mode"
    call status_check(ierr,error_text)
    write(*,*)""
    SELECT CASE (iaer)
    CASE(1) ! DUST
      ! Insert the definition of the variable
      tmpvarname="opadust"
      ierr=NF90_DEF_VAR(outfid,trim(tmpvarname),nf90_float,varshape,tmpvaridout)
      error_text="ERROR: Couldn't create "//trim(tmpvarname)//" in the outfile"
      call status_check(ierr,error_text)
      write(*,*) trim(tmpvarname)//" has been created in the outfile"
      
      ! Write the attributes
      if (opatype.eq."ext") then
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Dust extinction opacity at reference wavelength")
      else ! opatype.eq."abs"
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Dust absorption opacity at reference wavelength")
      endif
      
    CASE(2) ! WATER ICE
      ! Insert the definition of the variable
      tmpvarname="opawice"
      ierr=NF90_DEF_VAR(outfid,trim(tmpvarname),nf90_float,varshape,tmpvaridout)
      error_text="ERROR: Couldn't create "//trim(tmpvarname)//" in the outfile"
      call status_check(ierr,error_text)
      write(*,*) trim(tmpvarname)//" has been created in the outfile"
      
      ! Write the attributes
      if (opatype.eq."ext") then
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Water ice extinction opacity at reference wavelength")
      else ! opatype.eq."abs"
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Water ice absorption opacity at reference wavelength")
      endif
      
    CASE(3) ! STORMDUST
      ! Insert the definition of the variable
      tmpvarname="opastormdust"
      ierr=NF90_DEF_VAR(outfid,trim(tmpvarname),nf90_float,varshape,tmpvaridout)
      error_text="ERROR: Couldn't create "//trim(tmpvarname)//" in the outfile"
      call status_check(ierr,error_text)
      write(*,*) trim(tmpvarname)//" has been created in the outfile"
      
      ! Write the attributes
      if (opatype.eq."ext") then
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Stormdust extinction opacity at reference wavelength")
      else ! opatype.eq."abs"
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Stormdust absorption opacity at reference wavelength")
      endif

    CASE(4) ! TOPDUST
      ! Insert the definition of the variable
      tmpvarname="opatopdust"
      ierr=NF90_DEF_VAR(outfid,trim(tmpvarname),nf90_float,varshape,tmpvaridout)
      error_text="ERROR: Couldn't create "//trim(tmpvarname)//" in the outfile"
      call status_check(ierr,error_text)
      write(*,*) trim(tmpvarname)//" has been created in the outfile"
      
      ! Write the attributes
      if (opatype.eq."ext") then
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Topdust extinction opacity at reference wavelength")
      else ! opatype.eq."abs"
        ierr=nf90_put_att(outfid,tmpvaridout,"long_name","Topdust absorption opacity at reference wavelength")
      endif

! + NEW AER?    
    END SELECT ! iaer
    
    ierr=nf90_put_att(outfid,tmpvaridout,"units","opacity/km")
    ierr=nf90_put_att(outfid,tmpvaridout,"refwavelength",wvl_val)
    ierr=nf90_put_att(outfid,tmpvaridout,"missing_value",missval(iaer))
    write(*,*)"with missing value = ",missval(iaer)
    
    ! End netcdf define mode
    ierr=nf90_enddef(outfid)
    error_text="ERROR: Couldn't end netcdf define mode"
    call status_check(ierr,error_text)
      
!==========================================================================
! 3.2 Computation of the opacities
!==========================================================================
    do ilon=1,lonlen
     do ilat=1,latlen
      do ialt=1,altlen
       do it=1,timelen
         ! Get the optpropfile reff values encompassing the GCM reff
         isize=1
         do while((isize.le.nsize).and.(radiusdyn(isize).lt.reff(ilon,ilat,ialt,it)))
           isize=isize+1
         enddo
         if ((isize.gt.nsize).or.(radiusdyn(1).gt.reff(ilon,ilat,ialt,it))) then
!           write(*,*) "WARNING: the GCM reff (",reff(ilon,ilat,ialt,it),"m) is out of the bounds"
!           write(*,*) "of the file ",trim(datadir)//'/'//trim(optpropfile)
!           write(*,*) "(reff_inf=",radiusdyn(1),"m ; reff_sup=",radiusdyn(nsize),"m)"
!           write(*,*) "A missing value will be written for opacity"

           ! NB: this should especially handle cases when reff=0
           opa_aer(ilon,ilat,ialt,it)=missval(iaer)
           
         else if (mmr(ilon,ilat,ialt,it).eq.missval(iaer)) then
           ! if there is a missing value in input file
           opa_aer(ilon,ilat,ialt,it)=missval(iaer)
           
         else
           if (radiusdyn(isize).eq.reff(ilon,ilat,ialt,it)) then
             ! if the GCM reff is exactly in the optpropfile
             isize1 = isize
             isize2 = isize
           else ! radius(isize)>reff and radiusdyn(isize-1)<reff
             isize1 = isize-1
             isize2 = isize
           endif
         
           ! Qext COMPUTATION
           ! Linear interpolation in effective radius
           if (isize2.ne.isize1) then
             coeff = (reff(ilon,ilat,ialt,it)-radiusdyn(isize1))/(radiusdyn(isize2)-radiusdyn(isize1))
           else
             coeff = 0.
           endif
           reffQext(1) = Qext(iwvl1,isize1)+coeff*(Qext(iwvl1,isize2)-Qext(iwvl1,isize1))
           reffQext(2) = Qext(iwvl2,isize1)+coeff*(Qext(iwvl2,isize2)-Qext(iwvl2,isize1))
           ! Linear interpolation in wavelength
           if (iwvl2.ne.iwvl1) then
             coeff = (wvl_val-wvl(iwvl1))/(wvl(iwvl2)-wvl(iwvl1))
           else
             coeff = 0.
           endif
           Qext_val = reffQext(1)+coeff*(reffQext(2)-reffQext(1))
           
           ! COMPUTATION OF THE EXTINCTION OPACITY [1/km]
           opa_aer(ilon,ilat,ialt,it)= 750.*Qext_val*mmr(ilon,ilat,ialt,it)*rho(ilon,ilat,ialt,it)&
                                               / ( rho_aer(iaer) * reff(ilon,ilat,ialt,it) )
             
             
           if (opatype.eq."abs") then
             ! omeg COMPUTATION
             ! Linear interpolation in effective radius
             if (isize2.ne.isize1) then
               coeff = (reff(ilon,ilat,ialt,it)-radiusdyn(isize1))/(radiusdyn(isize2)-radiusdyn(isize1))
             else
               coeff = 0.
             endif
             reffomeg(1) = omeg(iwvl1,isize1)+coeff*(omeg(iwvl1,isize2)-omeg(iwvl1,isize1))
             reffomeg(2) = omeg(iwvl2,isize1)+coeff*(omeg(iwvl2,isize2)-omeg(iwvl2,isize1))
             ! Linear interpolation in wavelength
             if (iwvl2.ne.iwvl1) then
               coeff = (wvl_val-wvl(iwvl1))/(wvl(iwvl2)-wvl(iwvl1))
             else
               coeff = 0.
             endif
             omeg_val = reffomeg(1)+coeff*(reffomeg(2)-reffomeg(1))

             ! COMPUTATION OF THE ABSORPTION OPACITY [1/km]
             opa_aer(ilon,ilat,ialt,it)= (1 - omeg_val) * opa_aer(ilon,ilat,ialt,it)
           endif ! opatype.eq."abs"
                  
         endif
       enddo ! it
      enddo ! ialt
     enddo ! ilat
    enddo ! ilon

!==========================================================================
! 3.3 Writing in the output file
!==========================================================================
    ierr = NF90_PUT_VAR(outfid, tmpvaridout, opa_aer(:,:,:,:))
    error_text="Error: could not write "//trim(tmpvarname)//" data in the outfile"
    call status_check(ierr,error_text)

!==========================================================================
! 3.4 Prepare next loop
!==========================================================================
    DEALLOCATE(mmr)
    DEALLOCATE(reff)
    DEALLOCATE(opa_aer)
    DEALLOCATE(wvl)
    DEALLOCATE(radiusdyn)
    DEALLOCATE(Qext)
    DEALLOCATE(omeg)
    
  endif ! if aerok(iaer)

ENDDO ! iaer

!=============================================================================== 
! 4. Close the files and end the program
!===============================================================================
ierr = nf90_close(gcmfid)
error_text="Error: could not close file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr = nf90_close(outfid)
error_text="Error: could not close file "//trim(outfile)
call status_check(ierr,error_text)

write(*,*)"";write(*,*)"Program aeroptical completed!"

end program aeroptical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine status_check(ierr,error_text)

use netcdf
implicit none
!================================================================
!  Arguments:
!================================================================
integer,intent(in) :: ierr ! NetCDF status number
character(len=256),intent(in) :: error_text

if (ierr.ne.nf90_noerr) then
  write(*,*)trim(error_text)
  write(*,*)trim(nf90_strerror(ierr))
  stop
endif

end subroutine status_check


!*******************************************************************************

subroutine inidims(gcmfile,gcmfid,outfile,outfid,lonlen,latlen,altlen,timelen,&
                   latdimout,londimout,altdimout,timedimout,&
                   GCM_layers,layerdimout,interlayerdimout)

!==============================================================================
! Purpose:
! Read the dimensions of the input file
! and write them in the output file
!==============================================================================
! Remarks:
! The NetCDF files must be open
!==============================================================================
use netcdf

implicit none

!==============================================================================
! Arguments:
!==============================================================================
character(len=128),intent(in) :: gcmfile ! name of the netcdf GCM input file
integer,intent(in) :: gcmfid ! [netcdf] GCM input file ID number
character(len=128),intent(in) :: outfile ! name of the netcdf output file
integer,intent(in) :: outfid ! [netcdf] file ID number

integer,intent(out) ::  lonlen,latlen,altlen,timelen
! nb of grid points along longitude,latitude,altitude,time
integer,intent(out) :: latdimout,londimout,altdimout,timedimout
! latdimout: stores latitude dimension ID number
! londimout: stores longitude dimension ID number
! altdimout: stores altitude dimension ID number
! timedimout: stores time dimension ID number
integer,intent(out) :: GCM_layers ! number of GCM layers
integer,intent(out) :: layerdimout ! "GCM_layers" ID
integer,intent(out) :: interlayerdimout ! "GCM_layers+1" ID

!==============================================================================
! Local variables:
!==============================================================================
integer :: ierr ! [netcdf] subroutine returned error code
character(len=256) :: error_text ! text to display in case of error

integer :: tmpdimid,tmpvarid ! temporary store a dimension/variable ID number
character (len=64) :: long_name,units,positive
! long_name(): [netcdf] long_name attribute
! units(): [netcdf] units attribute
! positive(): [netcdf] positive direction attribute (for altitude)
real, dimension(:), allocatable:: lat,lon,alt,time,ctl
! lat(): array, stores latitude coordinates
! lon(): array, stores longitude coordinates
! alt(): array, stores altitude coordinates
! time(): array, stores time coordinates
! ctl(): array, stores controle variable
integer :: ctllen ! nb of elements in the controle array
integer :: tmpdimidout,tmpvaridout ! temporary stores a dimension/variable ID number

!==============================================================================
! LONGITUDE
!==============================================================================
! Get the dimension in GCM file
ierr=nf90_inq_dimid(gcmfid,"longitude",tmpdimid)
error_text="Error: Dimension <longitude> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=nf90_inquire_dimension(gcmfid,tmpdimid,len=lonlen)
allocate(lon(lonlen))

! Create the dimension in output file
ierr=NF90_DEF_DIM(outfid,"longitude",lonlen,londimout)
error_text="Error: could not define the longitude dimension in the outfile"
call status_check(ierr,error_text)

! Get the field in GCM file
ierr=nf90_inq_varid(gcmfid,"longitude",tmpvarid)
error_text="Error: Field <longitude> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=NF90_GET_VAR(gcmfid,tmpvarid,lon)
error_text="Failed to load longitude"
call status_check(ierr,error_text)

! Create the field in the output file
ierr=NF90_DEF_VAR(outfid,"longitude",nf90_float,(/londimout/),tmpvaridout)
error_text="Error: could not define the longitude variable in the outfile"
call status_check(ierr,error_text)

! Get the field attributes in the GCM file
ierr=nf90_get_att(gcmfid,tmpvarid,"long_name",long_name)
if (ierr.ne.nf90_noerr) then
! if no attribute "long_name", try "title"
  ierr=nf90_get_att(gcmfid,tmpvarid,"title",long_name)
endif
ierr=nf90_get_att(gcmfid,tmpvarid,"units",units)

! Put the field attributes in the output file
ierr=nf90_put_att(outfid,tmpvaridout,"long_name",long_name)
ierr=nf90_put_att(outfid,tmpvaridout,"units",units)

! Write the field values in the output file
ierr=nf90_enddef(outfid) ! end netcdf define mode
error_text="Error: could not end the define mode of the outfile"
call status_check(ierr,error_text)

ierr=NF90_PUT_VAR(outfid,tmpvaridout,lon)
error_text="Error: could not write the longitude field in the outfile"
call status_check(ierr,error_text)

!==============================================================================
! LATITUDE
!==============================================================================
! Switch to netcdf define mode
ierr=nf90_redef(outfid)
error_text="Error: could not switch to define mode in the outfile"
call status_check(ierr,error_text)

! Get the dimension in GCM file
ierr=nf90_inq_dimid(gcmfid,"latitude",tmpdimid)
error_text="Error: Dimension <latitude> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=nf90_inquire_dimension(gcmfid,tmpdimid,len=latlen)
allocate(lat(latlen))

! Create the dimension in output file
ierr=NF90_DEF_DIM(outfid,"latitude",latlen,latdimout)
error_text="Error: could not define the latitude dimension in the outfile"
call status_check(ierr,error_text)

! Get the field in GCM file
ierr=nf90_inq_varid(gcmfid,"latitude",tmpvarid)
error_text="Error: Field <latitude> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=NF90_GET_VAR(gcmfid,tmpvarid,lat)
error_text="Failed to load latitude"
call status_check(ierr,error_text)

! Create the field in the output file
ierr=NF90_DEF_VAR(outfid,"latitude",nf90_float,(/latdimout/),tmpvaridout)
error_text="Error: could not define the latitude variable in the outfile"
call status_check(ierr,error_text)

! Get the field attributes in the GCM file
ierr=nf90_get_att(gcmfid,tmpvarid,"long_name",long_name)
if (ierr.ne.nf90_noerr) then
! if no attribute "long_name", try "title"
  ierr=nf90_get_att(gcmfid,tmpvarid,"title",long_name)
endif
ierr=nf90_get_att(gcmfid,tmpvarid,"units",units)

! Put the field attributes in the output file
ierr=nf90_put_att(outfid,tmpvaridout,"long_name",long_name)
ierr=nf90_put_att(outfid,tmpvaridout,"units",units)

! Write the field values in the output file
ierr=nf90_enddef(outfid) ! end netcdf define mode
error_text="Error: could not end the define mode of the outfile"
call status_check(ierr,error_text)

ierr=NF90_PUT_VAR(outfid,tmpvaridout,lat)
error_text="Error: could not write the latitude field in the outfile"
call status_check(ierr,error_text)

!==============================================================================
! ALTITUDE
!==============================================================================
! Switch to netcdf define mode
ierr=nf90_redef(outfid)
error_text="Error: could not switch to define mode in the outfile"
call status_check(ierr,error_text)

! Get the dimension in GCM file
ierr=nf90_inq_dimid(gcmfid,"altitude",tmpdimid)
error_text="Error: Dimension <altitude> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=nf90_inquire_dimension(gcmfid,tmpdimid,len=altlen)
allocate(alt(altlen))

! Create the dimension in output file
ierr=NF90_DEF_DIM(outfid,"altitude",altlen,altdimout)
error_text="Error: could not define the altitude dimension in the outfile"
call status_check(ierr,error_text)

! Get the field in GCM file
ierr=nf90_inq_varid(gcmfid,"altitude",tmpvarid)
error_text="Error: Field <altitude> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=NF90_GET_VAR(gcmfid,tmpvarid,alt)
error_text="Failed to load altitude"
call status_check(ierr,error_text)

! Create the field in the output file
ierr=NF90_DEF_VAR(outfid,"altitude",nf90_float,(/altdimout/),tmpvaridout)
error_text="Error: could not define the altitude variable in the outfile"
call status_check(ierr,error_text)

! Get the field attributes in the GCM file
ierr=nf90_get_att(gcmfid,tmpvarid,"long_name",long_name)
if (ierr.ne.nf90_noerr) then
! if no attribute "long_name", try "title"
  ierr=nf90_get_att(gcmfid,tmpvarid,"title",long_name)
endif
ierr=nf90_get_att(gcmfid,tmpvarid,"units",units)
ierr=nf90_get_att(gcmfid,tmpvarid,"positive",positive)

! Put the field attributes in the output file
ierr=nf90_put_att(outfid,tmpvaridout,"long_name",long_name)
ierr=nf90_put_att(outfid,tmpvaridout,"units",units)
ierr=nf90_put_att(outfid,tmpvaridout,"positive",positive)

! Write the field values in the output file
ierr=nf90_enddef(outfid) ! end netcdf define mode
error_text="Error: could not end the define mode of the outfile"
call status_check(ierr,error_text)

ierr=NF90_PUT_VAR(outfid,tmpvaridout,alt)
error_text="Error: could not write the altitude field in the outfile"
call status_check(ierr,error_text)

!==============================================================================
! TIME
!==============================================================================
! Switch to netcdf define mode
ierr=nf90_redef(outfid)
error_text="Error: could not switch to define mode in the outfile"
call status_check(ierr,error_text)

! Get the dimension in GCM file
ierr=nf90_inq_dimid(gcmfid,"Time",tmpdimid)
error_text="Error: Dimension <Time> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=nf90_inquire_dimension(gcmfid,tmpdimid,len=timelen)
allocate(time(timelen))

! Create the dimension in output file
ierr=NF90_DEF_DIM(outfid,"Time",timelen,timedimout)
error_text="Error: could not define the time dimension in the outfile"
call status_check(ierr,error_text)

! Get the field in GCM file
ierr=nf90_inq_varid(gcmfid,"Time",tmpvarid)
error_text="Error: Field <Time> is missing in file "//trim(gcmfile)
call status_check(ierr,error_text)

ierr=NF90_GET_VAR(gcmfid,tmpvarid,time)
error_text="Failed to load Time"
call status_check(ierr,error_text)

! Create the field in the output file
ierr=NF90_DEF_VAR(outfid,"Time",nf90_float,(/timedimout/),tmpvaridout)
error_text="Error: could not define the Time variable in the outfile"
call status_check(ierr,error_text)

! Get the field attributes in the GCM file
ierr=nf90_get_att(gcmfid,tmpvarid,"long_name",long_name)
if (ierr.ne.nf90_noerr) then
! if no attribute "long_name", try "title"
  ierr=nf90_get_att(gcmfid,tmpvarid,"title",long_name)
endif
ierr=nf90_get_att(gcmfid,tmpvarid,"units",units)

! Put the field attributes in the output file
ierr=nf90_put_att(outfid,tmpvaridout,"long_name",long_name)
ierr=nf90_put_att(outfid,tmpvaridout,"units",units)

! Write the field values in the output file
ierr=nf90_enddef(outfid) ! end netcdf define mode
error_text="Error: could not end the define mode of the outfile"
call status_check(ierr,error_text)

ierr=NF90_PUT_VAR(outfid,tmpvaridout,time)
error_text="Error: could not write the Time field in the outfile"
call status_check(ierr,error_text)

!==============================================================================
! CONTROLE
!==============================================================================
! Switch to netcdf define mode
ierr=nf90_redef(outfid)
error_text="Error: could not switch to define mode in the outfile"
call status_check(ierr,error_text)

! Get the dimension in GCM file
ierr=nf90_inq_dimid(gcmfid,"index",tmpdimid)
error_text="Dimension <index> is missing in file "//trim(gcmfile)&
           //". We'll skip that one."
if (ierr.ne.nf90_noerr) then
  write(*,*)trim(error_text)
  ierr=nf90_enddef(outfid) ! end netcdf define mode
  error_text="Error: could not end the define mode of the outfile"
  call status_check(ierr,error_text)
else
  ierr=nf90_inquire_dimension(gcmfid,tmpdimid,len=ctllen)
  allocate(ctl(ctllen))

  ! Create the dimension in output file
  ierr=NF90_DEF_DIM(outfid,"index",ctllen,tmpdimidout)
  error_text="Error: could not define the index dimension in the outfile"
  call status_check(ierr,error_text)

  ! Get the field in GCM file
  ierr=nf90_inq_varid(gcmfid,"controle",tmpvarid)
  error_text="Error: Field <controle> is missing in file "//trim(gcmfile)
  call status_check(ierr,error_text)

  ierr=NF90_GET_VAR(gcmfid,tmpvarid,ctl)
  error_text="Failed to load ctl"
  call status_check(ierr,error_text)

  ! Create the field in the output file
  ierr=NF90_DEF_VAR(outfid,"controle",nf90_float,(/tmpdimidout/),tmpvaridout)
  error_text="Error: could not define the controle variable in the outfile"
  call status_check(ierr,error_text)

  ! Get the field attributes in the GCM file
  ierr=nf90_get_att(gcmfid,tmpvarid,"long_name",long_name)
  if (ierr.ne.nf90_noerr) then
  ! if no attribute "long_name", try "title"
    ierr=nf90_get_att(gcmfid,tmpvarid,"title",long_name)
  endif

  ! Put the field attributes in the output file
  ierr=nf90_put_att(outfid,tmpvaridout,"long_name",long_name)

  ! Write the field values in the output file
  ierr=nf90_enddef(outfid) ! end netcdf define mode
  error_text="Error: could not end the define mode of the outfile"
  call status_check(ierr,error_text)

  ierr=NF90_PUT_VAR(outfid,tmpvaridout,ctl)
  error_text="Error: could not write the controle field in the outfile"
  call status_check(ierr,error_text)
endif


!==============================================================================
! Load size of aps() or sigma() (in case it is not altlen)
!==============================================================================
! Switch to netcdf define mode
ierr=nf90_redef(outfid)
error_text="Error: could not switch to define mode in the outfile"
call status_check(ierr,error_text)

! Default is that GCM_layers=altlen
! but for outputs of zrecast, it may be a different value
ierr=nf90_inq_dimid(gcmfid,"GCM_layers",tmpdimid)
if (ierr.ne.nf90_noerr) then
  ! didn't find a GCM_layers dimension; therefore we have:
  GCM_layers=altlen  
else
  ! load value of GCM_layers
  ierr=nf90_inquire_dimension(gcmfid,tmpdimid,len=GCM_layers)
endif

! Create the dimensions in output file
ierr = NF90_DEF_DIM(outfid,"GCM_layers",GCM_layers,layerdimout)
error_text="Error: could not define the GCM_layers dimension in the outfile"
call status_check(ierr,error_text)
ierr = NF90_DEF_DIM(outfid,"GCM_interlayers",GCM_layers+1,interlayerdimout)
error_text="Error: could not define the GCM_interlayers dimension in the outfile"
call status_check(ierr,error_text)

! End netcdf define mode
ierr=nf90_enddef(outfid)
error_text="Error: could not end the define mode of the outfile"
call status_check(ierr,error_text)

end subroutine inidims


!*******************************************************************************

subroutine init2(gcmfid,lonlen,latlen,altlen,GCM_layers, &
                 outfid,londimout,latdimout,altdimout, &
                 layerdimout,interlayerdimout)
!==============================================================================
! Purpose:
! Copy ap() , bp(), aps(), bps(), aire() and phisinit()
! from input file to outpout file
!==============================================================================
! Remarks:
! The NetCDF files must be open
!==============================================================================

use netcdf

implicit none

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: gcmfid  ! NetCDF output file ID
integer, intent(in) :: lonlen ! # of grid points along longitude
integer, intent(in) :: latlen ! # of grid points along latitude
integer, intent(in) :: altlen ! # of grid points along altitude
integer, intent(in) :: GCM_layers ! # of GCM atmospheric layers
integer, intent(in) :: outfid ! NetCDF output file ID
integer, intent(in) :: londimout ! longitude dimension ID
integer, intent(in) :: latdimout ! latitude dimension ID
integer, intent(in) :: altdimout ! altitude dimension ID
integer, intent(in) :: layerdimout ! GCM_layers dimension ID
integer, intent(in) :: interlayerdimout ! GCM_layers+1 dimension ID
!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:),allocatable :: ap,bp ! hybrid vertical coordinates
real,dimension(:),allocatable :: sigma ! sigma levels
real,dimension(:,:),allocatable :: aire ! mesh areas
real,dimension(:,:),allocatable :: phisinit ! Ground geopotential
integer :: ierr
integer :: tmpvarid ! temporary variable ID
logical :: area ! is "aire" available ?
logical :: phis ! is "phisinit" available ?
logical :: hybrid ! are "aps" and "bps" available ?
logical :: apbp ! are "ap" and "bp" available ?

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
ierr=nf90_inq_varid(gcmfid,"aps",tmpvarid)
if (ierr.ne.nf90_noerr) then
  write(*,*) "Ooops. Failed to get aps ID. OK, will look for sigma coord."
  hybrid=.false.
else
  ierr=NF90_GET_VAR(gcmfid,tmpvarid,aps)
  hybrid=.true.
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed reading aps"
  endif

  ! hybrid coordinate bps
!  write(*,*) "bps: altlen=",altlen," GCM_layers=",GCM_layers
  allocate(bps(GCM_layers),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "init2: failed to allocate bps!"
    stop
  endif
  ierr=nf90_inq_varid(gcmfid,"bps",tmpvarid)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed to get bps ID."
  endif
  ierr=NF90_GET_VAR(gcmfid,tmpvarid,bps)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed reading bps"
  endif
endif

! hybrid coordinate ap
allocate(ap(GCM_layers+1),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate ap!"
  stop
else
  ierr=nf90_inq_varid(gcmfid,"ap",tmpvarid)
  if (ierr.ne.nf90_noerr) then
    write(*,*) "Ooops. Failed to get ap ID. OK."
    apbp=.false.
  else
    ierr=NF90_GET_VAR(gcmfid,tmpvarid,ap)
    apbp=.true.
    if (ierr.ne.nf90_noerr) then
      stop "Error: Failed reading ap"
    endif
  endif
endif

! hybrid coordinate bp
allocate(bp(GCM_layers+1),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate bp!"
  stop
else
  ierr=nf90_inq_varid(gcmfid,"bp",tmpvarid)
  if (ierr.ne.nf90_noerr) then
    write(*,*) "Ooops. Failed to get bp ID. OK."
    apbp=.false.
  else
    ierr=NF90_GET_VAR(gcmfid,tmpvarid,bp)
    apbp=.true.
    if (ierr.ne.nf90_noerr) then
      stop "Error: Failed reading bp"
    endif
  endif
endif

! sigma levels (if any)
if (.not.hybrid) then
  allocate(sigma(GCM_layers),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "init2: failed to allocate sigma"
    stop
  endif
  ierr=nf90_inq_varid(gcmfid,"sigma",tmpvarid)
  ierr=NF90_GET_VAR(gcmfid,tmpvarid,sigma)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed reading sigma"
  endif
endif ! of if (.not.hybrid)

! mesh area
allocate(aire(lonlen,latlen),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate aire!"
  stop
endif
ierr=nf90_inq_varid(gcmfid,"aire",tmpvarid)
if (ierr.ne.nf90_noerr) then
  write(*,*)"init2 warning: Failed to get aire ID."
  area = .false.
else
  ierr=NF90_GET_VAR(gcmfid,tmpvarid,aire)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed reading aire"
  endif
  area = .true.
endif

! ground geopotential phisinit
allocate(phisinit(lonlen,latlen),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate phisinit!"
  stop
endif
ierr=nf90_inq_varid(gcmfid,"phisinit",tmpvarid)
if (ierr.ne.nf90_noerr) then
  write(*,*)"init2 warning: Failed to get phisinit ID."
  phis = .false.
else
  ierr=NF90_GET_VAR(gcmfid,tmpvarid,phisinit)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed reading phisinit"
  endif
  phis = .true.
endif

!==============================================================================
! 2. Write
!==============================================================================

!==============================================================================
! 2.1 Hybrid coordinates ap() , bp(), aps() and bps()
!==============================================================================
if(hybrid) then 
! define aps
  ! Switch to netcdf define mode
  ierr=nf90_redef(outfid)
  ! Insert the definition of the variable
  ierr=NF90_DEF_VAR(outfid,"aps",nf90_float,(/layerdimout/),tmpvarid)
  if (ierr.ne.nf90_noerr) then
     stop "init2 Error: Failed to define the variable aps"
  endif
  ! Write the attributes
  ierr=nf90_put_att(outfid,tmpvarid,"long_name","hybrid pressure at midlayers")
  ierr=nf90_put_att(outfid,tmpvarid,"units"," ")
  ! End netcdf define mode
  ierr=nf90_enddef(outfid)

! write aps
  ierr=NF90_PUT_VAR(outfid,tmpvarid,aps)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed to write aps"
  endif

! define bps
  ! Switch to netcdf define mode
  ierr=nf90_redef(outfid)
  ! Insert the definition of the variable
  ierr=NF90_DEF_VAR(outfid,"bps",nf90_float,(/layerdimout/),tmpvarid)
  if (ierr.ne.nf90_noerr) then
     stop "init2 Error: Failed to define the variable bps"
  endif
  ! Write the attributes
  ierr=nf90_put_att(outfid,tmpvarid,"long_name","hybrid sigma at midlayers")
  ierr=nf90_put_att(outfid,tmpvarid,"units"," ")
  ! End netcdf define mode
  ierr=nf90_enddef(outfid)
    
! write bps
  ierr=NF90_PUT_VAR(outfid,tmpvarid,bps)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed to write bps"
  endif

  if (apbp) then
! define ap

    ! Switch to netcdf define mode
    ierr=nf90_redef(outfid)
    ! Insert the definition of the variable
    ierr=NF90_DEF_VAR(outfid,"ap",nf90_float,(/interlayerdimout/),tmpvarid)
    if (ierr.ne.nf90_noerr) then
       stop "init2 Error: Failed to define the variable ap"
    endif
    ! Write the attributes
    ierr=nf90_put_att(outfid,tmpvarid,"long_name","hybrid sigma at interlayers")
    ierr=nf90_put_att(outfid,tmpvarid,"units"," ")
    ! End netcdf define mode
    ierr=nf90_enddef(outfid)

! write ap
    ierr=NF90_PUT_VAR(outfid,tmpvarid,ap)
    if (ierr.ne.nf90_noerr) then
      stop "Error: Failed to write ap"
    endif

! define bp

    ! Switch to netcdf define mode
    ierr=nf90_redef(outfid)
    ! Insert the definition of the variable
    ierr=NF90_DEF_VAR(outfid,"bp",nf90_float,(/interlayerdimout/),tmpvarid)
    if (ierr.ne.nf90_noerr) then
       stop "init2 Error: Failed to define the variable bp"
    endif
    ! Write the attributes
    ierr=nf90_put_att(outfid,tmpvarid,"long_name","hybrid sigma at interlayers")
    ierr=nf90_put_att(outfid,tmpvarid,"units"," ")
    ! End netcdf define mode
    ierr=nf90_enddef(outfid)
  
! write bp
    ierr=NF90_PUT_VAR(outfid,tmpvarid,bp)
    if (ierr.ne.nf90_noerr) then
      stop "Error: Failed to write bp"
    endif
  endif ! of if (apbp)

else

  ! Switch to netcdf define mode
  ierr=nf90_redef(outfid)
  ! Insert the definition of the variable
  ierr=NF90_DEF_VAR(outfid,"sigma",nf90_float,(/layerdimout/),tmpvarid)
  if (ierr.ne.nf90_noerr) then
     stop "init2 Error: Failed to define the variable sigma"
  endif
  ! Write the attributes
  ierr=nf90_put_att(outfid,tmpvarid,"long_name","sigma at midlayers")
  ierr=nf90_put_att(outfid,tmpvarid,"units"," ")
  ! End netcdf define mode
  ierr=nf90_enddef(outfid)
  
! write sigma
  ierr=NF90_PUT_VAR(outfid,tmpvarid,sigma)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed to write sigma"
  endif
endif ! of if (hybrid)

!==============================================================================
! 2.2 aire() and phisinit()
!==============================================================================

if (area) then

  ! Switch to netcdf define mode
  ierr=nf90_redef(outfid)
  ! Insert the definition of the variable
  ierr=NF90_DEF_VAR(outfid,"aire",nf90_float,(/londimout,latdimout/),tmpvarid)
  if (ierr.ne.nf90_noerr) then
     stop "init2 Error: Failed to define the variable aire"
  endif
  ! Write the attributes
  ierr=nf90_put_att(outfid,tmpvarid,"long_name","Mesh area")
  ierr=nf90_put_att(outfid,tmpvarid,"units","m2")
  ! End netcdf define mode
  ierr=nf90_enddef(outfid)
  
  ! write aire
  ierr=NF90_PUT_VAR(outfid,tmpvarid,aire)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed to write aire"
  endif
endif ! of if (area)

if (phis) then

  ! Switch to netcdf define mode
  ierr=nf90_redef(outfid)
  ! Insert the definition of the variable
  ierr=NF90_DEF_VAR(outfid,"phisinit",nf90_float,(/londimout,latdimout/),tmpvarid)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed to define the variable phisinit"
  endif
  ! Write the attributes
  ierr=nf90_put_att(outfid,tmpvarid,"long_name","Ground level geopotential")
  ierr=nf90_put_att(outfid,tmpvarid,"units"," ")
  ! End netcdf define mode
  ierr=nf90_enddef(outfid)

  ! write phisinit
  ierr=NF90_PUT_VAR(outfid,tmpvarid,phisinit)
  if (ierr.ne.nf90_noerr) then
    stop "init2 Error: Failed to write phisinit"
  endif

endif ! of if (phis)


! Cleanup
if (allocated(aps)) deallocate(aps)
if (allocated(bps)) deallocate(bps)
if (allocated(ap)) deallocate(ap)
if (allocated(bp)) deallocate(bp)
if (allocated(sigma)) deallocate(sigma)
if (allocated(phisinit)) deallocate(phisinit)
if (allocated(aire)) deallocate(aire)

end subroutine init2
