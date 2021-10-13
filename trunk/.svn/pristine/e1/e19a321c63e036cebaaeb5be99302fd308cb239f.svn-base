PROGRAM start_archive2icosa

  USE xios
  USE mod_wait
  USE netcdf
  
  IMPLICIT NONE
  INCLUDE "mpif.h"
  INTEGER :: rank
  INTEGER :: size
  INTEGER :: ierr

  CHARACTER(len=*),PARAMETER :: id="client"
  INTEGER :: comm
  TYPE(xios_duration) :: dtime
  CHARACTER(len=15) :: calendar_type
  TYPE(xios_context) :: ctx_hdl

  INTEGER :: n,l
  INTEGER :: src_ibegin, src_iend, src_topo_ibegin, src_topo_iend
  INTEGER :: src_jbegin, src_jend, src_topo_jbegin, src_topo_jend
  INTEGER :: src_ni, src_ni_glo, src_topo_ni, src_topo_ni_glo
  INTEGER :: src_nj, src_nj_glo, src_topo_nj, src_topo_nj_glo
  INTEGER :: src_nlev ! number of vertical layers
  INTEGER :: src_nq=7 ! number of tracers
  INTEGER :: src_nt=1 ! number of time steps
  INTEGER :: src_nsoil=18 ! number of soil layers
  DOUBLE PRECISION,ALLOCATABLE :: lev_values(:) ! vertical axis 
  DOUBLE PRECISION,ALLOCATABLE :: lev_p1_values(:) ! vertical axis 
  DOUBLE PRECISION,ALLOCATABLE :: nq_values(:) ! tracer # axis 
  DOUBLE PRECISION,ALLOCATABLE :: soil_layers_values(:) ! soil axis
  DOUBLE PRECISION,ALLOCATABLE :: src_lon(:) ! mesh center coordinate
  DOUBLE PRECISION,ALLOCATABLE :: src_lat(:)
  DOUBLE PRECISION,ALLOCATABLE :: src_ap(:)
  DOUBLE PRECISION,ALLOCATABLE :: src_bp(:)
  DOUBLE PRECISION,ALLOCATABLE :: src_controle(:)
  DOUBLE PRECISION,ALLOCATABLE :: src_field_2D(:,:)
  DOUBLE PRECISION,ALLOCATABLE :: src_pk(:,:)
  DOUBLE PRECISION,ALLOCATABLE :: src_field_3D(:,:,:)
  DOUBLE PRECISION,ALLOCATABLE :: src_pressure(:,:,:)
  DOUBLE PRECISION,ALLOCATABLE :: src_theta_rhodz(:,:,:)
  DOUBLE PRECISION,ALLOCATABLE :: src_topo_lon(:) ! mesh center coordinate
  DOUBLE PRECISION,ALLOCATABLE :: src_topo_lat(:)
  DOUBLE PRECISION,ALLOCATABLE :: src_topo(:,:)
  DOUBLE PRECISION,ALLOCATABLE :: src_field_3D_soil(:,:,:) ! 3D grid in soil, for tsoil et inertia
  DOUBLE PRECISION,ALLOCATABLE :: src_soil_layers(:) ! soil_depth
  DOUBLE PRECISION,ALLOCATABLE :: src_field_3D_p1(:,:,:) !q2 field
  
  CHARACTER(LEN=*),PARAMETER :: src_file="start_archive_nc4.nc"
  CHARACTER(LEN=*),PARAMETER :: src_topo_file="surface_nc4.nc"
  CHARACTER(LEN=*),PARAMETER :: dst_coord_file="start_icosa_ref.nc"
  CHARACTER(LEN=*),PARAMETER :: src_controle_file="startphy_icosa_ref.nc"
  CHARACTER(LEN=*),PARAMETER :: output_start_file="start_icosa_prefinalize.nc"
  CHARACTER(LEN=*),PARAMETER :: output_startfi_file="startfi_prefinalize.nc"
  DOUBLE PRECISION,ALLOCATABLE :: dst_lon(:),dst_lat(:)
  DOUBLE PRECISION,ALLOCATABLE :: dst_boundslon(:,:) ! mesh corner coordinates
  DOUBLE PRECISION,ALLOCATABLE :: dst_boundslat(:,:)
  INTEGER :: dst_ibegin !, dst_iend
  INTEGER :: dst_ni, dst_ni_glo
  INTEGER :: dst_nvertex
  INTEGER :: ncid
  INTEGER :: dimids(4)
  INTEGER :: varid
  
  INTEGER :: div, remain
  INTEGER :: ts ! time step #
  DOUBLE PRECISION,PARAMETER :: pi=acos(-1.d0)
  DOUBLE PRECISION :: gravity,kappa,preff
! Tracers
  CHARACTER(LEN=11)          :: i_trac,format_string
  INTEGER                    :: i

!!! MPI Initialization
  CALL MPI_INIT(ierr)
  CALL init_wait

!!! XIOS Initialization (get the local communicator)
  CALL xios_initialize(id,return_comm=comm)
! get local rank of MPI process
  CALL MPI_COMM_RANK(comm,rank,ierr)
! get total number of MPI processes
  CALL MPI_COMM_SIZE(comm,size,ierr)

!!! Open files and load sizes and coordinates
  ierr=NF90_OPEN(src_topo_file, NF90_NOWRITE, ncid)
  ierr=NF90_INQ_VARID(ncid,"zMOL",varid)
  ierr=NF90_INQUIRE_VARIABLE(ncid, varid,dimids=dimids)
  write(*,*) "rank=",rank,"dimids=",dimids
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(1), len=src_topo_ni_glo)
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(2), len=src_topo_nj_glo)
  write(*,*) "rank=",rank," src_topo_ni_glo=",src_topo_ni_glo ! longitude
  write(*,*) "rank=",rank," src_topo_nj_glo=",src_topo_nj_glo ! latitude

! assume domain splitup with MPI only along latitudes
  src_topo_ni=src_topo_ni_glo
  src_topo_ibegin=0
  src_topo_iend=src_topo_ibegin+src_ni-1
  write(*,*) "rank=",rank," src_topo_ni=",src_topo_ni
  
  src_topo_jbegin=0
  DO n=0,size-1
    src_topo_nj=src_topo_nj_glo/size
    IF (n<MOD(src_topo_nj_glo,size)) src_topo_nj=src_topo_nj+1
    IF (n==rank) exit
    src_topo_jbegin=src_topo_jbegin+src_topo_nj
  ENDDO
  src_topo_jend=src_topo_jbegin+src_topo_nj-1
  write(*,*) "rank=",rank," src_topo_nj=",src_topo_nj, &
             " src_topo_jbegin=",src_topo_jbegin

  ALLOCATE(src_topo_lon(src_topo_ni))
  ALLOCATE(src_topo_lat(src_topo_nj))
  ALLOCATE(src_topo(src_topo_ni,src_topo_nj))

! load src_topo_lon and src_topo_lat
  ierr=NF90_INQ_VARID(ncid,"longitude",varid)
  ierr=NF90_GET_VAR(ncid,varid, src_topo_lon, &
                    start=(/src_topo_ibegin+1/),count=(/src_topo_ni/))
  WRITE(*,*) rank,":src_topo_lon(1:2)=",src_topo_lon(1:2)
  ierr=NF90_INQ_VARID(ncid,"latitude",varid)
  ierr=NF90_GET_VAR(ncid,varid, src_topo_lat, &
                    start=(/src_topo_jbegin+1/),count=(/src_topo_nj/))
  WRITE(*,*) rank,":src_topo_lat(1:2)=",src_topo_lat(1:2)

! from start_archive.nc file
  ierr=NF90_OPEN(src_file, NF90_NOWRITE, ncid)
  ierr=NF90_INQ_VARID(ncid,"temp",varid)
  ierr=NF90_INQUIRE_VARIABLE(ncid, varid,dimids=dimids)
  write(*,*) "rank=",rank,"dimids=",dimids
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(1), len=src_ni_glo)
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(2), len=src_nj_glo)
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(3), len=src_nlev)
  write(*,*) "rank=",rank," src_ni_glo=",src_ni_glo ! longitude
  write(*,*) "rank=",rank," src_nj_glo=",src_nj_glo ! latitude
  write(*,*) "rank=",rank," src_nlev=",src_nlev ! number of vertical layers
  write(*,*) "rank=",rank," src_nq=",src_nq ! number of tracers
! soil_depth with tsoil variable
  ierr=NF90_INQ_VARID(ncid,"tsoil",varid)
  ierr=NF90_INQUIRE_VARIABLE(ncid, varid,dimids=dimids)
  write(*,*) "rank=",rank,"dimids tsoil =",dimids
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(3), len=src_nsoil)
  write(*,*) "rank=",rank," src_nsoil=",src_nsoil ! number of vertical layers

! assume domain splitup with MPI only along latitudes
  src_ni=src_ni_glo
  src_ibegin=0
  src_iend=src_ibegin+src_ni-1
  write(*,*) "rank=",rank," src_ni=",src_ni
  
  src_jbegin=0
  DO n=0,size-1
    src_nj=src_nj_glo/size
    IF (n<MOD(src_nj_glo,size)) src_nj=src_nj+1
    IF (n==rank) exit
    src_jbegin=src_jbegin+src_nj
  ENDDO
  src_jend=src_jbegin+src_nj-1
  write(*,*) "rank=",rank," src_nj=",src_nj," src_jbegin=",src_jbegin

  ALLOCATE(src_lon(src_ni))
  ALLOCATE(src_lat(src_nj))
  ALLOCATE(src_field_2D(src_ni,src_nj))
  ALLOCATE(src_pk(src_ni,src_nj))
  ALLOCATE(src_field_3D(src_ni,src_nj,src_nlev))
  ALLOCATE(src_field_3D_p1(src_ni,src_nj,src_nlev+1))
  ALLOCATE(src_field_3D_soil(src_ni,src_nj,src_nsoil))
!  ALLOCATE(src_field_soil_layers(src_nsoil))
  ALLOCATE(src_pressure(src_ni,src_nj,src_nlev+1))
  ALLOCATE(src_theta_rhodz(src_ni,src_nj,src_nlev))

! load src_lon and src_lat
  ierr=NF90_INQ_VARID(ncid,"rlonv",varid)
  ierr=NF90_GET_VAR(ncid,varid, src_lon, &
                    start=(/src_ibegin+1/),count=(/src_ni/))
! convert rad to deg
  src_lon(1:src_ni)=src_lon(1:src_ni)*(180.d0/pi)
  WRITE(*,*) rank,":src_lon=",src_lon
  ierr=NF90_INQ_VARID(ncid,"rlatu",varid)
  ierr=NF90_GET_VAR(ncid,varid, src_lat, &
                    start=(/src_jbegin+1/),count=(/src_nj/))
! convert rad to deg
  src_lat(1:src_nj)=src_lat(1:src_nj)*(180.d0/pi)
  WRITE(*,*) rank,":src_lat=",src_lat

! load ap, bp and controle
  ALLOCATE(src_ap(src_nlev+1),src_bp(src_nlev+1))
  ierr=NF90_INQ_VARID(ncid,"ap",varid)
  ierr=NF90_GET_VAR(ncid,varid,src_ap)
  WRITE(*,*) rank,":src_ap(1:5)=",src_ap(1:5)
  ierr=NF90_INQ_VARID(ncid,"bp",varid)
  ierr=NF90_GET_VAR(ncid,varid,src_bp)
  WRITE(*,*) rank,":src_bp(1:5)=",src_bp(1:5)
  
! controle is taken in startphy_icosa_ref as start_archive_nc4 is too old
  ierr=NF90_OPEN(src_controle_file, NF90_NOWRITE, ncid)
  ALLOCATE(src_controle(100))
  ierr=NF90_INQ_VARID(ncid,"controle",varid)
  ierr=NF90_GET_VAR(ncid,varid,src_controle)
  ! day_ini set to 0 as startphy_icosa_ref is a restart
  src_controle(3)=0
  WRITE(*,*) rank,":src_controle(1:5)=",src_controle(1:5)
  gravity=src_controle(7)
  WRITE(*,*) rank,":gravity=",gravity
  kappa=src_controle(9)
  WRITE(*,*) rank,":kappa=",kappa
!  preff=src_controle(18)
! AD: Warning preff is hardcoded because not in controle of restartfi
  preff=610.
  WRITE(*,*) rank,":preff=",preff

! destination coordinates
  ierr=NF90_OPEN(dst_coord_file, NF90_NOWRITE, ncid)
  ierr=NF90_INQ_VARID(ncid,"bounds_lon_mesh",varid)
  ierr=NF90_INQUIRE_VARIABLE(ncid, varid,dimids=dimids)
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(1), len=dst_nvertex)
  ierr=NF90_INQUIRE_DIMENSION(ncid, dimids(2), len=dst_ni_glo)
  write(*,*) "rank=",rank," dst_nvertex=",dst_nvertex ! vertex
  write(*,*) "rank=",rank," dst_ni_glo=",dst_ni_glo ! vertex boundaries

! evenly split into MPI domains
  div    = dst_ni_glo/size
  remain = MOD( dst_ni_glo, size )
  IF (rank < remain) THEN
    dst_ni=div+1 ;
    dst_ibegin=rank*(div+1) ;
  ELSE
    dst_ni=div ;
    dst_ibegin= remain * (div+1) + (rank-remain) * div ;
  ENDIF
  write(*,*) "rank=",rank," dst_ni=",dst_ni

  ALLOCATE(dst_lon(dst_ni))
  ALLOCATE(dst_lat(dst_ni))
  ALLOCATE(dst_boundslon(dst_nvertex,dst_ni))
  ALLOCATE(dst_boundslat(dst_nvertex,dst_ni))

! load dst_lon, dst_lat, dst_boundslon and dst_boundslat
  ierr=NF90_INQ_VARID(ncid,"lon_mesh",varid)
  ierr=NF90_GET_VAR(ncid,varid, dst_lon, &
                    start=(/dst_ibegin+1/), count=(/dst_ni/))
  WRITE(*,*) rank,":dst_lon(1:5)=",dst_lon(1:5)
  ierr=NF90_INQ_VARID(ncid,"lat_mesh",varid)
  ierr=NF90_GET_VAR(ncid,varid, dst_lat, &
                    start=(/dst_ibegin+1/), count=(/dst_ni/))
  WRITE(*,*) rank,":dst_lat(1:5)=",dst_lat(1:5)
  ierr=NF90_INQ_VARID(ncid,"bounds_lon_mesh",varid)
  ierr=NF90_GET_VAR(ncid,varid,dst_boundslon, &
                    start=(/1,dst_ibegin+1/), count=(/dst_nvertex,dst_ni/))
  WRITE(*,*) rank,":dst_boundslon(:,1:2)=",dst_boundslon(:,1:2)
  ierr=NF90_INQ_VARID(ncid,"bounds_lat_mesh",varid)
  ierr=NF90_GET_VAR(ncid,varid, dst_boundslat, &
                    start=(/1,dst_ibegin+1/), count=(/dst_nvertex,dst_ni/))
  WRITE(*,*) rank,":dst_boundslat(:,1:2)=",dst_boundslat(:,1:2)


! Initialize XIOS context
  WRITE(*,*) rank,":CALL xios_context_initialize()"
  CALL xios_context_initialize("test",comm)
  CALL xios_get_handle("test",ctx_hdl)
  CALL xios_set_current_context(ctx_hdl)

! Set XIOS calendar and timestep
  CALL xios_get_calendar_type(calendar_type)
  WRITE(*,*) rank,":calendar_type = ", calendar_type
  dtime%second = 3600
  CALL xios_set_timestep(dtime)

! Set axes
  ! vertical atm axis
  ALLOCATE(lev_values(src_nlev))
  lev_values=(/ (l,l=1,src_nlev) /)
  CALL xios_set_axis_attr("lev",n_glo=src_nlev,value=lev_values)
  ! lev+1 for q2
  ALLOCATE(lev_p1_values(src_nlev+1))
  lev_p1_values=(/ (l,l=1,src_nlev+1) /)
  CALL xios_set_axis_attr("lev_p1",n_glo=src_nlev+1,value=lev_p1_values)
  ! tracers axis
  ALLOCATE(nq_values(src_nq))
  nq_values=(/(l,l=1,src_nq)/)
  CALL xios_set_axis_attr("nq",n_glo=src_nq,value=nq_values)
  ! soil layers axis
  ALLOCATE(soil_layers_values(src_nsoil))
  soil_layers_values=(/ (l,l=1,src_nsoil) /)
  CALL xios_set_axis_attr("soil_layers",n_glo=src_nsoil,value=soil_layers_values)

! Set domains
  CALL xios_set_domain_attr("src_domain_regular", &
                            ni_glo=src_ni_glo, nj_glo=src_nj_glo, &
                            ibegin=src_ibegin, ni=src_ni, &
                            jbegin=src_jbegin, nj=src_nj, &
                            type='rectilinear')
  CALL xios_set_domain_attr("src_domain_regular", &
                             data_dim=2, &
                             data_ibegin=0, data_ni=src_ni, &
                             data_jbegin=0, data_nj=src_nj)
  CALL xios_set_domain_attr("src_domain_regular", &
                            lonvalue_1D=src_lon, &
                            latvalue_1D=src_lat)

  CALL xios_set_domain_attr("src_topo_domain_regular", &
                            ni_glo=src_topo_ni_glo, nj_glo=src_topo_nj_glo, &
                            ibegin=src_topo_ibegin, ni=src_topo_ni, &
                            jbegin=src_topo_jbegin, nj=src_topo_nj, &
                            type='rectilinear')
  CALL xios_set_domain_attr("src_topo_domain_regular", &
                             data_dim=2, &
                             data_ibegin=0, data_ni=src_topo_ni, &
                             data_jbegin=0, data_nj=src_topo_nj)
  CALL xios_set_domain_attr("src_topo_domain_regular", &
                            lonvalue_1D=src_topo_lon, &
                            latvalue_1D=src_topo_lat)

  CALL xios_set_domain_attr("src_domain_regular_clean", &
                            ni_glo=src_ni_glo-1, nj_glo=src_nj_glo, &
                            ibegin=src_ibegin, ni=src_ni-1, &
                            jbegin=src_jbegin, nj=src_nj, &
                            type='rectilinear')
  CALL xios_set_domain_attr("src_domain_regular_clean", &
                             data_dim=2, &
                             data_ibegin=0, data_ni=src_ni-1, &
                             data_jbegin=0, data_nj=src_nj)
  CALL xios_set_domain_attr("src_domain_regular_clean", &
                            lonvalue_1D=src_lon(1:src_ni-1), &
                            latvalue_1D=src_lat)

  CALL xios_set_domain_attr("dst_domain_unstructured", &
                            ni_glo=dst_ni_glo, &
                            ibegin=dst_ibegin, &
                            ni=dst_ni, &
                            type="unstructured")
  CALL xios_set_domain_attr("dst_domain_unstructured", &
                            lonvalue_1D=dst_lon, &
                            latvalue_1D=dst_lat, &
                            bounds_lon_1D=dst_boundslon, &
                            bounds_lat_1D=dst_boundslat, &
                            nvertex=dst_nvertex)

! Finalize XIOS context definition
  WRITE(*,*) rank,":CALL xios_close_context_definition()"
  CALL xios_close_context_definition()
  CALL xios_get_handle("test",ctx_hdl)
  CALL xios_set_current_context(ctx_hdl)

! Temporal loop
  DO ts=1,src_nt
    WRITE(*,*) rank,":ts=",ts
    ! Update calendar
    CALL xios_update_calendar(ts)

    ! Topography
    CALL xios_recv_field("zMOL",src_topo)
    WRITE(*,*) rank,":topo(1:2,1:3)=",src_topo(1:2,1:3)
    ! Send surface geopotential
    CALL xios_send_field("topo",src_topo(:,:)*gravity*1000.)

    ! Surface pressure:
    !! get data using XIOS:
    CALL xios_recv_field("src_ps",src_field_2D)
    WRITE(*,*) rank,":src_ps(1:2,1:3)=",src_field_2D(1:2,1:3)
    !! write data using XIOS
    CALL xios_send_field("ps_clean",src_field_2D(1:src_ni-1,1:src_nj))

    ! compute inter-layer pressures
    DO l=1,src_nlev+1
      src_pressure(:,:,l) = src_ap(l)+src_bp(l)*src_field_2D(:,:)
    ENDDO
    
    ! surface temperature:
    CALL xios_recv_field("src_tsurf",src_field_2D)
    WRITE(*,*) rank,":src_tsurf(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("tsurf_clean",src_field_2D(1:src_ni-1,1:src_nj))

    ! co2 ice coverage:
    CALL xios_recv_field("src_co2ice",src_field_2D)
    WRITE(*,*) rank,":src_co2ice(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("co2ice_clean",src_field_2D(1:src_ni-1,1:src_nj))

    ! emissivity:
    CALL xios_recv_field("src_emis",src_field_2D)
    WRITE(*,*) rank,":src_emis(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("emis_clean",src_field_2D(1:src_ni-1,1:src_nj))

    ! q2 : q2surf for first layer of q2, the rest is q2atm
    CALL xios_recv_field("src_q2surf",src_field_2D)
    WRITE(*,*) rank,":src_q2(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_recv_field("src_q2atm",src_field_3D)
    src_field_3D_p1(:,:,1)=src_field_2D(:,:)
    src_field_3D_p1(:,:,2:src_nlev+1)=src_field_3D(:,:,:)

    CALL xios_send_field("q2_clean",src_field_3D_p1(1:src_ni-1,1:src_nj,1:src_nlev+1))
    
    ! Temperature:
    CALL xios_recv_field("src_temp",src_field_3D)
    ! compute theta_rhodz
    DO l=1,src_nlev
      src_pk(:,:)=((.5/preff)*(src_pressure(:,:,l)+src_pressure(:,:,l+1)))**kappa
      src_theta_rhodz(:,:,l) = src_field_3D(:,:,l) * &
      ((src_pressure(:,:,l)-src_pressure(:,:,l+1))/gravity)/src_pk(:,:)
    ENDDO
    CALL xios_send_field("theta_rhodz_clean", &
                         src_theta_rhodz(1:src_ni-1,1:src_nj,1:src_nlev))

    ! zonal wind
    CALL xios_recv_field("src_u",src_field_3D)
    CALL xios_send_field("u_clean", &
                          src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
    
    ! meridional wind
    CALL xios_recv_field("src_v",src_field_3D)
    CALL xios_send_field("v_clean", &
                          src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
    ! tracers
    DO i=1,src_nq
       if (i < 10) then
           format_string = "(A8,I1)"
       else
           format_string = "(A8,I2)"
       endif
       write (i_trac, format_string) "src_trac", i
       CALL xios_recv_field(i_trac,src_field_3D(:,:,:))
       if (i==1) THEN 
          CALL xios_send_field("co2_clean", &
                        src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
          CALL xios_recv_field("src_co2_surf",src_field_2D)
          CALL xios_send_field("co2_surf_clean", &
                        src_field_2D(1:src_ni-1,1:src_nj))
       elseif (i==2) THEN 
          CALL xios_send_field("dust_number_clean", &
                        src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
          CALL xios_recv_field("src_dust_number_surf",src_field_2D)
          CALL xios_send_field("dust_number_surf_clean", &
                        src_field_2D(1:src_ni-1,1:src_nj))
       elseif (i==3) THEN 
          CALL xios_send_field("dust_mass_clean", &
                        src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
          CALL xios_recv_field("src_dust_mass_surf",src_field_2D)
          CALL xios_send_field("dust_mass_surf_clean", &
                        src_field_2D(1:src_ni-1,1:src_nj))
       elseif (i==4) THEN 
          CALL xios_send_field("ccn_number_clean", &
                        src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
          CALL xios_recv_field("src_ccn_number_surf",src_field_2D)
          CALL xios_send_field("ccn_number_surf_clean", &
                        src_field_2D(1:src_ni-1,1:src_nj))
       elseif (i==5) THEN 
          CALL xios_send_field("ccn_mass_clean", &
                        src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
          CALL xios_recv_field("src_ccn_mass_surf",src_field_2D)
          CALL xios_send_field("ccn_mass_surf_clean", &
                        src_field_2D(1:src_ni-1,1:src_nj))
       elseif (i==6) THEN 
          CALL xios_send_field("h2o_ice_clean", &
                        src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
          CALL xios_recv_field("src_h2o_ice_surf",src_field_2D)
          CALL xios_send_field("h2o_ice_surf_clean", &
                        src_field_2D(1:src_ni-1,1:src_nj))
       elseif (i==7) THEN 
          CALL xios_send_field("h2o_vap_clean", &
                        src_field_3D(1:src_ni-1,1:src_nj,1:src_nlev))
          CALL xios_recv_field("src_h2o_vap_surf",src_field_2D)
          CALL xios_send_field("h2o_vap_surf_clean", &
                        src_field_2D(1:src_ni-1,1:src_nj))
       ENDIF
    ENDDO

    ! soil temperature
    CALL xios_recv_field("src_tsoil",src_field_3D_soil)
    CALL xios_send_field("tsoil_clean", &
                          src_field_3D_soil(1:src_ni-1,1:src_nj,1:src_nsoil))
    ! soil thermal intertia
    CALL xios_recv_field("src_inertiedat",src_field_3D_soil)
    CALL xios_send_field("inertiedat_clean", &
                          src_field_3D_soil(1:src_ni-1,1:src_nj,1:src_nsoil))

    ! ZMEA "zmea Orographie sous-maille"
    CALL xios_recv_field("src_ZMEA",src_field_2D)
    WRITE(*,*) rank,":src_ZMEA(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("ZMEA_clean",src_field_2D(1:src_ni-1,1:src_nj))
    
    ! ZSTD "zstd Orographie sous-maille"
    CALL xios_recv_field("src_ZSTD",src_field_2D)
    WRITE(*,*) rank,":src_ZSTD(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("ZSTD_clean",src_field_2D(1:src_ni-1,1:src_nj))
    
    ! ZSIG "zsig Orographie sous-maille"
    CALL xios_recv_field("src_ZSIG",src_field_2D)
    WRITE(*,*) rank,":src_ZSIG(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("ZSIG_clean",src_field_2D(1:src_ni-1,1:src_nj))
    
    ! ZGAM "zgam Orographie sous-maille"
    CALL xios_recv_field("src_ZGAM",src_field_2D)
    WRITE(*,*) rank,":src_ZGAM(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("ZGAM_clean",src_field_2D(1:src_ni-1,1:src_nj))
    
    ! ZTHE "zthe Orographie sous-maille"
    CALL xios_recv_field("src_ZTHE",src_field_2D)
    WRITE(*,*) rank,":src_ZTHE(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("ZTHE_clean",src_field_2D(1:src_ni-1,1:src_nj))

    ! albedodat "albedo"
    CALL xios_recv_field("src_albedodat",src_field_2D)
    WRITE(*,*) rank,":src_albedodat(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("albedodat_clean",src_field_2D(1:src_ni-1,1:src_nj))

    ! z0 "surface roughness"
    CALL xios_recv_field("src_z0",src_field_2D)
    WRITE(*,*) rank,":src_z0(1:2,1:3)=",src_field_2D(1:2,1:3)
    CALL xios_send_field("z0_clean",src_field_2D(1:src_ni-1,1:src_nj))


  ENDDO ! of DO ts=1,src_nt
  
!! Finalize
  write(*,*) rank,":Finalize: call xios_context_finalize"
  CALL xios_context_finalize()

  write(*,*) rank,":Finalize: call MPI_COMM_FREE"
  CALL MPI_COMM_FREE(comm, ierr)
  write(*,*) rank,":Finalize: call xios_finalize"
  CALL xios_finalize()

  if (rank==0) then
    ! add a couple of things in the "startphy_icosa.nc" file
    write(*,*) rank,"Write controle() to startphy_icosa.nc"
    ierr=NF90_OPEN(output_startfi_file,NF90_WRITE,ncid)
    ierr=NF90_REDEF(ncid) ! switch to define mode
    ierr=NF90_DEF_DIM(ncid,"index",100,dimids(1))
    ierr=NF90_DEF_VAR(ncid,"controle",NF90_DOUBLE,dimids(1),varid)
    ierr=NF90_ENDDEF(ncid) ! switch out of define mode
    ierr=NF90_PUT_VAR(ncid,varid,src_controle(1:100))
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "NetCDF Error:",NF90_STRERROR(ierr)
    endif
    ierr=NF90_CLOSE(ncid)
    ! add a couple of things in the "start_icosa.nc" file
    ierr=NF90_OPEN(output_start_file,NF90_WRITE,ncid)
    ierr=NF90_REDEF(ncid)
    ierr=NF90_DEF_DIM(ncid,"nvertex_u",2,dimids(1))
    ierr=NF90_DEF_VAR(ncid,"iteration",NF90_FLOAT,varid)
    ierr=NF90_ENDDEF(ncid)
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "NetCDF Error:",NF90_STRERROR(ierr)
    endif
    ierr=NF90_PUT_VAR(ncid,varid,0) ! set "iteration" value to 0
    ierr=NF90_CLOSE(ncid)
    
  endif ! of if (rank==0)

  write(*,*) rank,":Finalize: call MPI_FINALIZE"
  CALL MPI_FINALIZE(ierr)

  write(*,*) rank,":my_remap: all is well that ends well!"

END PROGRAM start_archive2icosa
