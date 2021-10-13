PROGRAM rearrange_startphy

USE netcdf

IMPLICIT NONE

INTEGER :: ncid
INTEGER :: dimids(4)
INTEGER :: varid
INTEGER :: ierr
INTEGER :: nvar
INTEGER :: ndim
CHARACTER(LEN=100) :: varname
INTEGER :: varname_dimids(5)

INTEGER :: i,j,k,nb_cells,nlev,nsoil,nvertex,nlev_p1
INTEGER,ALLOCATABLE :: cell_index(:)
REAL,ALLOCATABLE :: ref_lon(:)
REAL,ALLOCATABLE :: ref_lat(:)
REAL,ALLOCATABLE :: lon(:)
REAL,ALLOCATABLE :: lat(:)
REAL,ALLOCATABLE :: ref_field(:),field(:)
REAL,ALLOCATABLE :: ref_field_3D(:,:),field_3D(:,:)
REAL,ALLOCATABLE :: ref_field_3D_p1(:,:),field_3D_p1(:,:)
REAL,ALLOCATABLE :: ref_field_soil(:,:),field_soil(:,:)
REAL,ALLOCATABLE :: ref_field_vertex(:,:),field_vertex(:,:)
REAL :: lati,loni
REAL :: diff_lat,diff_lon
CHARACTER(LEN=*),PARAMETER :: ref_file="startphy_icosa_ref.nc"
CHARACTER(LEN=*),PARAMETER :: file="startfi.nc"

DOUBLE PRECISION,PARAMETER :: pi=acos(-1.d0)
REAL :: reflatj,reflonj

! load coordinates from files
  ierr=NF90_OPEN(ref_file, NF90_NOWRITE, ncid)
  ierr=NF90_INQ_DIMID(ncid,"physical_points",dimids(1))
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(1), len=nb_cells)
  write(*,*) "nb_cells=",nb_cells,"dimids(1)=",dimids(1)
  allocate(ref_lat(nb_cells),ref_lon(nb_cells))
  
  ierr=NF90_INQ_VARID(ncid,"latitude",varid)
  ierr=NF90_GET_VAR(ncid,varid,ref_lat)
  if (ierr /= nf90_noerr) then
    write(*,*) "cannot load latitude from ",trim(ref_file)
  else
    write(*,*) "ref_lat(1:5)=",ref_lat(1:5)
  endif
  ierr=NF90_INQ_VARID(ncid,"longitude",varid)
  ierr=NF90_GET_VAR(ncid,varid,ref_lon)
  if (ierr /= nf90_noerr) then
    write(*,*) "cannot load longitude from ",trim(ref_file)
  else
    write(*,*) "ref_lon(1:5)=",ref_lon(1:5)
  endif
  ierr=NF90_CLOSE(ncid)
  
  ierr=NF90_OPEN(file, NF90_WRITE, ncid)
  allocate(lat(nb_cells),lon(nb_cells))
  ierr=NF90_INQ_VARID(ncid,"latitude",varid)
  ierr=NF90_GET_VAR(ncid,varid,lat)
  if (ierr /= nf90_noerr) then
    write(*,*) "cannot load lat from ",trim(file)
  else
    write(*,*) "lat(1:5)=",lat(1:5)
  endif
  ierr=NF90_INQ_VARID(ncid,"longitude",varid)
  ierr=NF90_GET_VAR(ncid,varid,lon)
  if (ierr /= nf90_noerr) then
    write(*,*) "cannot load lat from ",trim(file)
  else
    write(*,*) "lon(1:5)=",lon(1:5)
  endif
  
  ! find correspondances between lon/ref_lon & lat/ref_lat
  allocate(cell_index(nb_cells))
  cell_index(:)=0
  do i=1,nb_cells !in case lon and lat are expressed in rad in startphy_icosa_ref /180*pi
    lat(i)=lat(i) /(180.d0/pi) ; lon(i)=lon(i)/(180.d0/pi)
    do j=1,nb_cells
     reflatj=ref_lat(j)
     reflonj=ref_lon(j)

     if (lat(i) ==0.) then
       diff_lat=abs((lat(i)-ref_lat(j))/1.e-8)
     else
       diff_lat=abs((lat(i)-ref_lat(j))/lat(i))
     endif

     if (lon(i) ==0.) then
      diff_lon=abs((lon(i)-ref_lon(j))/1.e-8)
     else
      diff_lon=abs((lon(i)-ref_lon(j))/lon(i))
     endif

     if ((diff_lat <= 1.e-5).and.(diff_lon <= 1.e-5)) then
      cell_index(i)=j
      write(*,*)"j=",j," lat(i)=",lat(i)," ref_lat(j)=",reflatj
      write(*,*)"j=",j," lon(i)=",lon(i)," ref_lon(j)=",reflonj
     endif
    enddo ! of do j=1,nb_cells
    write(*,*) ">> i=",i,"cell_index(i)=",cell_index(i)
    ! sanity check:
    if (cell_index(i)==0) then
      write(*,*) "Error, could not find lon-lat match for i=",i , lat(i), lon(i)
      stop
    endif
    ! writing lat and lon in rad rather than deg
    ierr=NF90_INQ_VARID(ncid,"longitude",varid)
    ierr=NF90_PUT_VAR(ncid,varid,lon)

    ierr=NF90_INQ_VARID(ncid,"latitude",varid)
    ierr=NF90_PUT_VAR(ncid,varid,lat)
  enddo ! of do i=1,nb_cells
  
!  do i=1,nb_cells
!    write(*,*) "i=",i," cell_index(i)=",cell_index(i)
!    write(*,*) "  lat(i)=",lat(i),"ref_lat(cell_index(i))=",ref_lat(cell_index(i))
!  enddo
  
  ! load, rearrange and write variables
  ierr=NF90_INQ_DIMID(ncid,"physical_points",dimids(1))
  ierr=NF90_INQ_DIMID(ncid,"lev",dimids(2))
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(2), len=nlev)
  ierr=NF90_INQ_DIMID(ncid,"subsurface_layers",dimids(3))
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(3), len=nsoil)
  ierr=NF90_INQ_DIMID(ncid,"nvertex",dimids(4))
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(4), len=nvertex)
  ierr=NF90_INQ_DIMID(ncid,"lev_p1",dimids(5))
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(2), len=nlev_p1)
  ierr=NF90_INQUIRE(ncid,nVariables=nvar)
  write(*,*) "nvar=",nvar
  allocate(ref_field(nb_cells),field(nb_cells))
  allocate(ref_field_3D(nb_cells,nlev),field_3D(nb_cells,nlev))
  allocate(ref_field_3D_p1(nb_cells,nlev_p1),field_3D_p1(nb_cells,nlev_p1))
  allocate(ref_field_soil(nb_cells,nsoil),field_soil(nb_cells,nsoil))
  allocate(ref_field_vertex(nvertex,nb_cells),field_vertex(nvertex,nb_cells))
  write(*,*) "dimids :", dimids
  ! loop over variables:
  do k=1,nvar
    varname_dimids(:)=0
    ierr=NF90_INQUIRE_VARIABLE(ncid,k,name=varname,ndims=ndim,dimids=varname_dimids)
      write(*,*) "name: ",trim(varname),"; dimensions:",varname_dimids
    if (ierr /= nf90_noerr) then
      write(*,*) "error for variable k=",k
      write(*,*) nf90_strerror(ierr)
    endif
    ! Processing 2D variables
    if ((ndim==1).and.(varname_dimids(1)==dimids(1))) then
      write(*,*) "processing ",trim(varname)
      ! load field
      ierr=NF90_INQ_VARID(ncid,varname,varid)
      ierr=NF90_GET_VAR(ncid,varid,ref_field)
      ! rearrange field
      do i=1,nb_cells
        field(cell_index(i))=ref_field(i)
      enddo
      ! write field
      ierr=NF90_PUT_VAR(ncid,varid,field)

    ! Processing 3D variables : bounds_lon and bounds_lat
    else if ((ndim==2).and.(varname_dimids(1)==dimids(4))) then
      write(*,*) "processing ",trim(varname)
      ! load field_3D
      ierr=NF90_INQ_VARID(ncid,varname,varid)
      if (ierr /= nf90_noerr) print*,"error inqvarid ",trim(varname)
      ierr=NF90_GET_VAR(ncid,varid,ref_field_vertex)
      if (ierr /= nf90_noerr) print*,"error get_var ",trim(varname)
      ! rearrange field_3D
      do i=1,nb_cells
        field_vertex(:,cell_index(i))=ref_field_vertex(:,i)
      enddo
      ! write field_3D
      ierr=NF90_PUT_VAR(ncid,varid,field_vertex)
      if (ierr /= nf90_noerr) print*,"error putvar ",trim(varname)
 
    ! Processing 3D variables : altitude
    else if ((ndim==2).and.(varname_dimids(2)==dimids(2))) then
      write(*,*) "processing ",trim(varname)
      ! load field_3D
      ierr=NF90_INQ_VARID(ncid,varname,varid)
      if (ierr /= nf90_noerr) print*,"error inqvarid ",trim(varname)
      ierr=NF90_GET_VAR(ncid,varid,ref_field_3D)
      if (ierr /= nf90_noerr) print*,"error get_var ",trim(varname)
      ! rearrange field_3D
      do i=1,nb_cells
        field_3D(cell_index(i),:)=ref_field_3D(i,:)
      enddo
      ! write field_3D
      ierr=NF90_PUT_VAR(ncid,varid,field_3D)
      if (ierr /= nf90_noerr) print*,"error putvar ",trim(varname)

    ! Processing 3D variables : soil
    else if ((ndim==2).and.(varname_dimids(2)==dimids(3))) then
      write(*,*) "processing ",trim(varname)
      ! load field_soil
      ierr=NF90_INQ_VARID(ncid,varname,varid)
      if (ierr /= nf90_noerr) print*,"error inqvarid ",trim(varname)
      ierr=NF90_GET_VAR(ncid,varid,ref_field_soil)
      if (ierr /= nf90_noerr) print*,"error get_var ",trim(varname)
      ! rearrange field_soil
      do i=1,nb_cells
        field_soil(cell_index(i),:)=ref_field_soil(i,:)
      enddo
      ! write field_soil
      ierr=NF90_PUT_VAR(ncid,varid,field_soil)
      if (ierr /= nf90_noerr) print*,"error putvar ",trim(varname)

    ! Processing 3D variables : q2 with lev_p1 as vertical axis
    else if ((ndim==2).and.(varname_dimids(2)==dimids(5))) then
      write(*,*) "processing ",trim(varname)
      ! load field_3D
      ierr=NF90_INQ_VARID(ncid,varname,varid)
      if (ierr /= nf90_noerr) print*,"error inqvarid ",trim(varname)
      ierr=NF90_GET_VAR(ncid,varid,ref_field_3D_p1)
      if (ierr /= nf90_noerr) print*,"error get_var ",trim(varname)
      ! rearrange field_3D
      do i=1,nb_cells
        field_3D_p1(cell_index(i),:)=ref_field_3D_p1(i,:)
      enddo
      ! write field_3D
      ierr=NF90_PUT_VAR(ncid,varid,field_3D_p1)
      if (ierr /= nf90_noerr) print*,"error putvar ",trim(varname)
    endif
  enddo
  
END PROGRAM rearrange_startphy
