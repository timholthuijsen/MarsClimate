










      SUBROUTINE surfini(ngrid,piceco2,qsurf)

      USE ioipsl_getin_p_mod, ONLY : getin_p
      use netcdf
      use tracer_mod, only: nqmx, noms
      use geometry_mod, only: longitude, latitude, ! in radians
     &                     cell_area ! for watercaptag diagnosis
      use surfdat_h, only: watercaptag, frost_albedo_threshold,
     &                     albedo_h2o_cap, inert_h2o_ice, albedodat,
     &                     albedice, dryness
      use mod_grid_phy_lmdz, only : klon_glo ! # of physics point on full grid
      use mod_phys_lmdz_para, only : is_master, gather, scatter
      USE comcstfi_h, ONLY: pi
      use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
      use datafile_mod, only: datadir
      IMPLICIT NONE
c=======================================================================
c
c   creation des calottes pour l'etat initial
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------
      include "callkeys.h"

      integer,intent(in) :: ngrid ! number of atmospheric columns
      real,intent(in) :: piceco2(ngrid) ! CO2 ice thickness
      real,intent(inout) :: qsurf(ngrid,nqmx) ! tracer on surface (kg/m2)

      INTEGER ig,icap,iq,alternate
      REAL icedryness ! ice dryness
      
      ! longwatercaptag is watercaptag. Trick for some compilers
      LOGICAL, DIMENSION(100000) :: longwatercaptag
      
! There are 4 different modes for ice distribution:
! icelocationmode = 1 ---> based on data from surface.nc
! icelocationmode = 2 ---> directly predefined for GCM resolutions 32x24 or 64x48
! icelocationmode = 3 ---> based on logical relations for latitude and longitude
! icelocationmode = 4 ---> predefined 64x48 but usable with every
! resolution, and easily adaptable for dynamico 
! For visualisation : > /u/tnalmd/bin/watercaps gcm_txt_output_file
      INTEGER,SAVE :: icelocationmode = 4
       
       
      !in case icelocationmode == 1
      INTEGER i,j
      INTEGER     imd,jmd
      PARAMETER   (imd=360,jmd=180)
      REAL        zdata(imd,jmd)
      REAL        zelat,zelon 

      INTEGER nb_ice(klon_glo,2)   ! number of counts | detected ice for GCM grid
      INTEGER latice(nbp_lat-1,2),lonice (nbp_lon,2) ! number of counts | detected ice along lat & lon axis

      REAL step,count,ratiolat

      INTEGER   ierr,nid,nvarid
      
      REAL,SAVE :: min_icevalue = 500.
      character(len=50) :: string = 'thermal'
      
      character (len=100) :: zedatafile

! problem with nested precompiling flags

      ! to handle parallel cases
      logical watercaptag_glo(ngrid)
      real dryness_glo(ngrid)
      real lati_glo(ngrid)
      real long_glo(ngrid)


c
c=======================================================================
! Initialize watercaptag (default is false)
      watercaptag_glo(:)=.false.

c     water ice outliers
c     ------------------------------------------

      IF ((water) .and. (caps)) THEN
     
c Perennial H20 north cap defined by watercaptag=true (allows surface to be
c hollowed by sublimation in vdifc).

c We might not want albedodat to be modified because it is used to write 
c restart files. Instead, albedo is directly modified when needed (i.e. 
c if we have watercaptag and no co2 ice), below and in albedocaps.F90

c       "Dryness coefficient" controlling the evaporation and
c        sublimation from the ground water ice (close to 1)
c        HERE, the goal is to correct for the fact
c        that the simulated permanent water ice polar caps
c        is larger than the actual cap and the atmospheric
c        opacity not always realistic.

         alternate = 0
         
         if (ngrid .ne. 1) then
           watercaptag(:) = .false.
           longwatercaptag(:) = .false.
         endif
         
         write(*,*) "surfini: Ice dryness ?"
         icedryness=1. ! default value
         call getin_p("icedryness",icedryness)
         write(*,*) "surfini: icedryness = ",icedryness
         dryness (:) = icedryness
         
      ! To be able to run in parallel, we work on the full grid
      ! and dispatch results afterwards

      ! start by geting latitudes and logitudes on full grid
      ! (in serial mode, this is just a copy)
      call gather(latitude,lati_glo)
      call gather(longitude,long_glo)

      if (is_master) then

        IF (ngrid .eq. 1) THEN ! special case for 1d --> do nothing
      
         print*, 'ngrid = 1, do no put ice caps in surfini.F'

        ELSE IF (icelocationmode .eq. 1) THEN
      
         print*,'Surfini: ice caps defined from surface.nc'
            
! This method detects ice as gridded value above min_icevalue in the field "string" from surface.nc
! Typically, it is for thermal inertia above 500 tiu.
! Two conditions are verified:
! 1. GCM ice caps are defined such as area is conserved for a given latitude
! (the approximation is that all points within the GCM latitude resolution have the same area).
! 2. caps are placed to fill the GCM points with the most detected ice first.
      

           
         zedatafile = trim(datadir)
 
        
         ierr=nf90_open(trim(zedatafile)//'/surface.nc',
     &   NF90_NOWRITE,nid)
     
         IF (ierr.NE.nf90_noerr) THEN
       write(*,*)'Error : cannot open file surface.nc '
       write(*,*)'(in phymars/surfini.F)'
       write(*,*)'It should be in :',trim(zedatafile),'/'
       write(*,*)'1) You can set this path in the callphys.def file:'
       write(*,*)'   datadir=/path/to/the/datafiles'
       write(*,*)'2) If necessary, surface.nc (and other datafiles)'
       write(*,*)'   can be obtained online on:'
       write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
       call abort_physic("surfini","missing surface.nc file",1)
         ENDIF
      
      
         ierr=nf90_inq_varid(nid, string, nvarid)
         if (ierr.ne.nf90_noerr) then
          write(*,*) 'surfini error, cannot find ',trim(string)
          write(*,*) ' in file ',trim(zedatafile),'/surface.nc'
          write(*,*)trim(nf90_strerror(ierr))
          call abort_physic("surfini","missing "//trim(string),1)
         endif

         ierr=nf90_get_var(nid, nvarid, zdata)

         if (ierr.ne.nf90_noerr) then
          write(*,*) 'surfini: error failed loading ',trim(string)
          write(*,*)trim(nf90_strerror(ierr))
          call abort_physic("surfini","failed loading "//trim(string),1)
         endif
 
                     
         ierr=nf90_close(nid)
 

         nb_ice(:,1) = 1 ! default: there is no ice
         latice(:,1) = 1
         lonice(:,1) = 1
         nb_ice(:,2) = 0
         latice(:,2) = 0
         lonice(:,2) = 0
         !print*,'jjm,iim',jjm,iim ! jjm =  nb lati , iim = nb longi

         ! loop over the GCM grid - except for poles (ig=1 and ngrid)
         do ig=2,klon_glo-1
      
          ! loop over the surface file grid      
          do i=1,imd
           do j=1,jmd
             zelon = i - 180.
             zelat = 90. - j 
            if ((abs(lati_glo(ig)*180./pi-zelat).le.
     &           90./real(nbp_lat-1)) .and.
     &          (abs(long_glo(ig)*180./pi-zelon).le.
     &           180./real(nbp_lon))) then
              ! count all points in that GCM grid point
              nb_ice(ig,1) = nb_ice(ig,1) + 1
              if (zdata(i,j) > min_icevalue)
                 ! count all detected points in that GCM grid point
     &           nb_ice(ig,2) = nb_ice(ig,2) + 1
             endif
           enddo
          enddo  

        ! projection of nb_ice on GCM lat and lon axes
          latice(1+(ig-2)/nbp_lon,:) =
     &     latice(1+(ig-2)/nbp_lon,:) + nb_ice(ig,:)
          lonice(1+mod(ig-2,nbp_lon),:) = 
     &     lonice(1+mod(ig-2,nbp_lon),:) + nb_ice(ig,:) ! lonice is USELESS ...

         enddo ! of do ig=2,klon_glo-1
     

     
         ! special case for poles
         nb_ice(1,2)   = 1  ! ice prescribed on north pole
         latice(1,:)   = nb_ice(1,:)
         lonice(1,:)   = nb_ice(1,:)
         latice(nbp_lat-1,:) = nb_ice(ngrid,:)
         lonice(nbp_lon,:) = nb_ice(ngrid,:)
      
     
!      print*, 'latice TOT', latice(:,1)
!      print*, 'latice FOUND', latice(:,2)
!      print*, 'lonice TOT', lonice(:,1)
!      print*, 'lonice FOUND', lonice(:,2)
      
!      print*, 'lat ratio', int(real(latice(:,2))/real(latice(:,1))*iim)
!      print*, 'lon ratio', int(real(lonice(:,2))/real(lonice(:,1))*jjm)
      
!      print*,''
!      print*,'sum lat', sum(latice(:,1)), sum(lonice(:,1))
!      print*,'sum lon', sum(latice(:,2)), sum(lonice(:,2))
      
    
         ! loop over GCM latitudes. CONSIDER ONLY NORTHERN HEMISPHERE
         do i=1,(nbp_lat-1)/2
          step  = 1. ! threshold to add ice cap
          count = 0. ! number of ice GCM caps at this latitude
          ! ratiolat is the ratio of area covered by ice within this GCM latitude range
          ratiolat  = real(latice(i,2))/real(latice(i,1))
          !print*,'i',i,(i-1)*iim+2,i*iim+1
     
          ! put ice caps while there is not enough ice,
          ! as long as the threshold is above 20%
          do while ((count.le.ratiolat*nbp_lon).and.(step.ge.0.2))
           count = 0.
           ! loop over GCM longitudes
           do j=1,nbp_lon
            ! if the detected ice ratio in the GCM grid point 
            ! is more than 'step', then add ice
            if (real(nb_ice((i-1)*nbp_lon+1+j,2)) 
     &        / real(nb_ice((i-1)*nbp_lon+1+j,1)) .ge. step) then
                  watercaptag_glo((i-1)*nbp_lon+1+j) = .true.
                  count = count + 1
            endif
           enddo ! of do j=1,nbp_lon
           !print*, 'step',step,count,ratiolat*nbp_lon
           step = step - 0.01
          enddo ! of do while
          !print*, 'step',step,count,ratiolat*nbp_lon

         enddo ! of do i=1,jjm/2
            

        ELSE IF (icelocationmode .eq. 2) THEN
      
         print*,'Surfini: predefined ice caps'
      
         if ((nbp_lon.eq.32).and.((nbp_lat-1).eq.24)) then ! 32x24
           
          print*,'water ice caps distribution for 32x24 resolution'
          longwatercaptag(1:9)    = .true. ! central cap - core
          longwatercaptag(26:33)  = .true. ! central cap
          longwatercaptag(1:33)  = .true. ! central cap
          longwatercaptag(56)  = .true. ! central cap
          longwatercaptag(58)  = .true. ! central cap
          longwatercaptag(60)  = .true. ! central cap
          longwatercaptag(62)  = .true. ! central cap
          longwatercaptag(64)  = .true. ! central cap
!---------------------   OUTLIERS  ----------------------------

         else if ((nbp_lon.eq.64).and.((nbp_lat-1).eq.48)) then ! 64x48

          print*,'water ice caps distribution for 64x48 resolution'
          longwatercaptag(1:65)   = .true. ! central cap - core
          longwatercaptag(75:85)  = .true. ! central cap 
          longwatercaptag(93:114) = .true. ! central cap
!---------------------   OUTLIERS  ----------------------------
          if (.true.) then
          longwatercaptag(136)    = .true. ! outlier, lat = 78.75
          longwatercaptag(138)    = .true. ! outlier, lat = 78.75
          longwatercaptag(140)    = .true. ! outlier, lat = 78.75
          longwatercaptag(142)    = .true. ! outlier, lat = 78.75
          longwatercaptag(161)    = .true. ! outlier, lat = 78.75
          longwatercaptag(163)    = .true. ! outlier, lat = 78.75
          longwatercaptag(165)    = .true. ! outlier, lat = 78.75
          longwatercaptag(183)    = .true. ! outlier, lat = 78.75
          longwatercaptag(185)    = .true. ! outlier, lat = 78.75
          longwatercaptag(187)    = .true. ! outlier, lat = 78.75
          longwatercaptag(189)    = .true. ! outlier, lat = 78.75
          longwatercaptag(191)    = .true. ! outlier, lat = 78.75
          longwatercaptag(193)    = .true. ! outlier, lat = 78.75
          longwatercaptag(194)    = .true. ! outlier, lat = 75
          longwatercaptag(203)    = .true. ! outlier, lat = 75
          longwatercaptag(207)    = .true. ! outlier, lat = 75
          longwatercaptag(244)    = .true. ! outlier, lat = 75
          longwatercaptag(246)    = .true. ! outlier, lat = 75
          longwatercaptag(250)    = .true. ! outlier, lat = 75
          longwatercaptag(252)    = .true. ! outlier, lat = 75
          longwatercaptag(254)    = .true. ! outlier, lat = 75
          longwatercaptag(256)    = .true. ! outlier, lat = 75
          endif
!--------------------------------------------------------------       


            
         else if (klon_glo .ne. 1) then
        
          print*,'No predefined ice location for this resolution :',
     &           nbp_lon,nbp_lat-1
          print*,'Please change icelocationmode in surfini.F'
          print*,'Or add some new definitions ...'
          call abort_physic("surfini",
     &         "no pre-definitions for this resolution",1)
          
         endif

         do ig=1,klon_glo
          if (longwatercaptag(ig)) watercaptag_glo(ig) = .true.
         enddo


        ELSE IF (icelocationmode .eq. 3) THEN
      
         print*,'Surfini: ice caps defined by lat and lon values'

         do ig=1,klon_glo
         
c-------- Towards olympia planitia water caps -----------
c-------------------------------------------------------- 

          if ( ( ( lati_glo(ig)*180./pi .ge. 77.  ) .and. ! cap #2
     .           ( lati_glo(ig)*180./pi .le. 80.  ) .and.
     .           ( long_glo(ig)*180./pi .ge. 110. ) .and.
     .           ( long_glo(ig)*180./pi .le. 181. ) )
     .         .or.

     .         ( ( lati_glo(ig)*180./pi .ge. 75.  ) .and. ! cap #4 (Korolev crater)
     .           ( lati_glo(ig)*180./pi .le. 76.  ) .and.
     .           ( long_glo(ig)*180./pi .ge. 150. ) .and.
     .           ( long_glo(ig)*180./pi .le. 168. ) )
     .         .or.
     .         ( ( lati_glo(ig)*180./pi .ge. 77 ) .and. ! cap #5
     .           ( lati_glo(ig)*180./pi .le. 80.  ) .and.
     .           ( long_glo(ig)*180./pi .ge. -150.) .and.
     .           ( long_glo(ig)*180./pi .le. -110.) ) )
     .         then
             
               if ((alternate .eq. 0)) then  ! 1/2 en 64x48 sinon trop large en lat
              !    watercaptag(ig)=.true.
                  alternate = 1
               else
                  alternate = 0
               endif !end if alternate = 0
               
          endif

c----------- Opposite olympia planitia water cap --------
c-------------------------------------------------------- 

          if ( ( ( lati_glo(ig)*180./pi     .ge.  80 ) .and.
     .         ( lati_glo(ig)*180./pi     .le.  84 ) )
     .         .and.
     .       ( ( long_glo(ig)*180./pi .lt. -95. ) .or.       !!! 32x24
     .         ( long_glo(ig)*180./pi .gt.  85. ) ) ) then   !!! 32x24
!     .     ( ( ( long_glo(ig)*180./pi .ge. -29. ) .and.       !!! 64x48
!     .         ( long_glo(ig)*180./pi .le.  90. ) ) .or.      !!! 64x48
!     .       ( ( long_glo(ig)*180./pi .ge. -77. ) .and.       !!! 64x48
!     .         ( long_glo(ig)*180./pi .le. -70. ) ) ) ) then  !!! 64x48
        !   watercaptag_glo(ig)=.true.
          endif


c -------------------- Central cap ----------------------
c-------------------------------------------------------- 

          if (abs(lati_glo(ig)*180./pi).gt.80)
     .          watercaptag_glo(ig)=.true.
           
c--------------------------------------------------------
c--------------------------------------------------------
         end do ! of (klon_glo)

        ELSE IF (icelocationmode .eq. 4) THEN
      
         print*,'icelocationmode = 4'
         print*,'Surfini: ice caps defined using manual 64x48 settings'
         print*,'(although, it should work with any resolution)'
         call locate_watercaptag(klon_glo,lati_glo,
     &            long_glo,watercaptag_glo)

!         print*,'watercaptag_glo(:), ',watercaptag_glo(:)

        ELSE
      
         print*, 'In surfini.F, icelocationmode is ', icelocationmode
         print*, 'It should be 1, 2, 3 or 4 (default is 4)'
         call abort_physic("surfini","wrong icelocationmode",1)

        ENDIF ! of if (icelocation)
       
       
        ! print caps locations - useful for plots too 
        print*,'surfini: latitude | longitude | ig'
        do ig=1,klon_glo
          dryness_glo(ig) = icedryness

          if (watercaptag_glo(ig)) then
             print*,'surfini: ice water cap', lati_glo(ig)*180./pi,
     &              long_glo(ig)*180./pi, ig
!             write(1,*),ig, lati_glo(ig)*180./pi,
!     &              cell_area(ig)
!             write(2,*), lati_glo(ig)*180./pi,
!     &              long_glo(ig)*180./pi,cell_area(ig)
!             write(3,*), ig, lati_glo(ig)*180./pi,
!     &              long_glo(ig)*180./pi,cell_area(ig)
          endif
        enddo
       
       endif !of if (is_master)
       
       if (ngrid.gt.1) then
        ! Now scatter fields watercaptag and dryness from master to all
        ! (is just a plain copy in serial mode)
        call scatter(dryness_glo,dryness)
        call scatter(watercaptag_glo,watercaptag)
       endif
       ELSE
         watercaptag(:) = .false.
       ENDIF ! (caps & water)
! end of #else of #ifndef MESOSCALE
!      END SUBROUTINE surfini(ngrid,piceco2,qsurf)
      END !SUBROUTINE surfini(ngrid,piceco2,qsurf)

      SUBROUTINE locate_watercaptag(klon_glo,lati_glo,
     &            long_glo,watercaptag_glo)

      USE comcstfi_h, ONLY: pi

      integer, intent(in) :: klon_glo
      real, intent(in) :: lati_glo(klon_glo)
      real, intent(in) :: long_glo(klon_glo)
      logical, intent(out) :: watercaptag_glo(klon_glo)
      integer :: ig,i
!      real, dimension(klon_glo,120) :: wcap
      real, dimension(120,2) :: latedge
      real, dimension(120,2) :: lonedge

! In icelocationmode=2 there are 120 manually predefined grid points where
! watercaptag is true (for the 64x48 resolution). The grid cells corners
! coordinates in latitude and longitude are written below. With this
! routine, we check if the grid cell center is in between any of those
! points. If so, watercaptag = true. 




      latedge(:,1)=(/ 
     & 88.125, 84.375, 84.375, 84.375, 84.375, 84.375,84.375, 84.375, 
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625,
     & 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625,
     & 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625,
     & 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625,
     & 80.625, 80.625, 76.875, 76.875, 76.875, 76.875, 76.875, 76.875,
     & 76.875, 76.875, 76.875, 76.875, 76.875, 76.875, 76.875, 73.125,
     & 73.125, 73.125, 73.125, 73.125, 73.125, 73.125, 73.125, 73.125/)


      latedge(:,2)=(/ 
     & 90. , 88.125, 88.125, 88.125, 88.125, 88.125,88.125, 88.125,   
     & 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125,
     & 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125,
     & 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125,
     & 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125,
     & 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125,
     & 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125,
     & 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125, 88.125,
     & 88.125, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375, 84.375,
     & 84.375, 84.375, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625,
     & 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 80.625, 76.875,
     & 76.875, 76.875, 76.875, 76.875, 76.875, 76.875, 76.875, 76.875/)


      lonedge(:,1)=(/
     &-180.    , -180.    , -177.1875, -171.5625,-165.9375, -160.3125,
     &-154.6875, -149.0625, -143.4375, -137.8125, -132.1875,-126.5625,
     &-120.9375, -115.3125, -109.6875, -104.0625,  -98.4375, -92.8125,
     & -87.1875,  -81.5625,  -75.9375,  -70.3125,  -64.6875, -59.0625,
     & -53.4375,  -47.8125,  -42.1875,  -36.5625,  -30.9375, -25.3125,
     & -19.6875,  -14.0625,   -8.4375,   -2.8125,    2.8125,   8.4375,
     &  14.0625,   19.6875,   25.3125,   30.9375,   36.5625,  42.1875,
     &  47.8125,   53.4375,   59.0625,   64.6875,   70.3125,  75.9375,
     &  81.5625,   87.1875,   92.8125,   98.4375,  104.0625, 109.6875,
     & 115.3125,  120.9375,  126.5625,  132.1875,  137.8125, 143.4375,
     & 149.0625,  154.6875,  160.3125,  165.9375,  171.5625,-132.1875,
     &-126.5625, -120.9375, -115.3125, -109.6875, -104.0625, -98.4375,
     & -92.8125,  -87.1875,  -81.5625,  -75.9375,  -30.9375, -25.3125,
     & -19.6875,  -14.0625,   -8.4375,   -2.8125,    2.8125,   8.4375,
     &  14.0625,   19.6875,   25.3125,   30.9375,   36.5625,  42.1875,
     &  47.8125,   53.4375,   59.0625,   64.6875,   70.3125,  75.9375,
     &  81.5625,   87.1875, -149.0625, -137.8125, -126.5625,-115.3125,
     &  -8.4375,    2.8125,   14.0625,  115.3125,  126.5625, 137.8125,
     & 149.0625,  160.3125,  171.5625, -180.    , -132.1875,-109.6875,
     &  98.4375,  109.6875,  132.1875,  143.4375,  154.6875,165.9375/) 

      lonedge(:,2)=(/ 
     & 180.    , -180.    , -171.5625, -165.9375,-160.3125, -154.6875,
     &-149.0625,-143.4375, -137.8125, -132.1875, -126.5625, -120.9375,
     &-115.3125,-109.6875, -104.0625,  -98.4375,  -92.8125,  -87.1875,
     & -81.5625, -75.9375,  -70.3125,  -64.6875,  -59.0625,  -53.4375,
     & -47.8125, -42.1875,  -36.5625,  -30.9375,  -25.3125,  -19.6875,
     & -14.0625,  -8.4375,   -2.8125,    2.8125,    8.4375,   14.0625,
     &  19.6875,  25.3125,   30.9375,   36.5625,   42.1875,   47.8125,
     &  53.4375,  59.0625,   64.6875,   70.3125,   75.9375,   81.5625,
     &  87.1875,  92.8125,   98.4375,  104.0625,  109.6875,  115.3125,
     & 120.9375, 126.5625,  132.1875,  137.8125,  143.4375,  149.0625,
     & 154.6875, 160.3125,  165.9375,  171.5625,  177.1875, -126.5625,
     &-120.9375,-115.3125, -109.6875, -104.0625,  -98.4375,  -92.8125,
     & -87.1875, -81.5625,  -75.9375,  -70.3125,  -25.3125,  -19.6875,
     & -14.0625,  -8.4375,   -2.8125,    2.8125,    8.4375,   14.0625,
     &  19.6875,  25.3125,   30.9375,   36.5625,   42.1875,   47.8125,
     &  53.4375,  59.0625,   64.6875,   70.3125,   75.9375,   81.5625,
     &  87.1875,  92.8125, -143.4375, -132.1875, -120.9375, -109.6875,
     &  -2.8125,   8.4375,   19.6875,  120.9375,  132.1875,  143.4375,
     & 154.6875, 165.9375,  177.1875, -177.1875, -126.5625, -104.0625,
     & 104.0625, 115.3125,  137.8125,  149.0625,  160.3125,171.5625/)


      watercaptag_glo(:) = .false.
      DO ig=1, klon_glo
        DO i=1, 120
           if ((long_glo(ig)*180./pi.ge.lonedge(i,1))
     &         .and.(long_glo(ig)*180./pi.le.lonedge(i,2))
     &         .and.(lati_glo(ig)*180./pi.ge.latedge(i,1))
     &         .and.(lati_glo(ig)*180./pi.le.latedge(i,2))) then
             watercaptag_glo(ig) = .true.
           endif
        ENDDO !i=1, 120
      ENDDO ! ig=1, klon_glo

      END SUBROUTINE locate_watercaptag
