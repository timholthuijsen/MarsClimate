










      SUBROUTINE concentrations(ngrid,nlayer,nq,
     &                          pplay,pt,pdt,pq,pdq,ptimestep)
                                             
      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,
     &                      igcm_o2, igcm_o3, igcm_h, igcm_h2,
     &                      igcm_oh, igcm_ho2, igcm_n2, igcm_ar,
     &                      igcm_h2o_vap, igcm_n, igcm_no, igcm_no2,
     &                      igcm_n2d, igcm_co2plus, igcm_oplus,
     &                      igcm_o2plus, igcm_coplus, igcm_cplus,
     &                      igcm_nplus, igcm_noplus, igcm_n2plus,
     &                      igcm_hplus, igcm_hco2plus, mmol,
     &                      igcm_he, igcm_elec
      use conc_mod, only: mmean, Akknew, rnew, cpnew
      implicit none

!=======================================================================
! CALCULATION OF MEAN MOLECULAR MASS, Cp, Akk and R
!
! mmean(ngrid,nlayer)	amu
! cpnew(ngrid,nlayer)	J/kg/K
! rnew(ngrid,nlayer)	J/kg/K
! akknew(ngrid,nlayer)	coefficient of thermal concduction
!
! version: April 2012 - Franck Lefevre
!=======================================================================

!     declarations
 
!
! For Fortran 77/Fortran 90 compliance always use line continuation
! symbols '&' in columns 73 and 6
!
! NB: to keep commons aligned, it is better to split them in groups
!     of given types (logical, integer, real, ...)

      COMMON/callkeys_l/callrad,calldifv,calladj,callcond,callsoil      &
     &   ,season,diurnal,lwrite,calllott,callstats,calleofdump          &
     &   ,callnirco2,callnlte,callthermos,callconduct,calleuv           &
     &   ,callmolvis,callmoldiff,thermochem,thermoswater,callemis       &
     &   ,callg2d,linear,rayleigh,tracer                                &
     &   ,scavenging,sedimentation                                      &
     &   ,activice,water,tifeedback,microphys,supersat,caps,photochem   &
     &   ,calltherm,callrichsl,callslope,tituscap,callyamada4,co2clouds &
     &   ,co2useh2o,meteo_flux,activeco2ice,CLFvaryingCO2,spantCO2      &
     &   ,CLFvarying,satindexco2,rdstorm,slpwind,calllott_nonoro        &
     &   ,latentheat_surfwater,gwd_convective_source,startphy_file      &
     &   ,hdo,hdofrac,cap_albedo,temp_dependant_m
     
      COMMON/callkeys_i/iradia,iaervar,ilwd,ilwb,ilwn,ncouche           &
     &   ,nltemodel,nircorr,solvarmod,solvaryear,dustinjection
     
      COMMON/callkeys_r/semi,alphan,euveff,                             &
     &   tke_heat_flux,dustrefir,fixed_euv_value,CLFfixval,             &
     &   coeff_injection,ti_injection,tf_injection,coeff_detrainment
     
      LOGICAL callrad,calldifv,calladj,callcond,callsoil,               &
     &   season,diurnal,lwrite,calllott,calllott_nonoro                 &
     &   ,callstats,calleofdump                                         &
     &   ,callnirco2,callnlte,callthermos,callconduct,                  &
     &    calleuv,callmolvis,callmoldiff,thermochem,thermoswater        &
     &   ,calltherm,callrichsl,callslope,tituscap,callyamada4

      COMMON/aeroutput/dustiropacity

      logical startphy_file

      logical callemis
      logical callg2d
      logical linear
      logical gwd_convective_source

      real semi
      real alphan
      real fixed_euv_value
      real euveff
      real tke_heat_flux
      real coeff_injection ! dust injection scheme coefficient
      real ti_injection ! local time of beginning injection
      real tf_injection ! local time of end injection
      real coeff_detrainment ! rocket dust detrainment coefficient
      real CLFfixval

      integer iaervar
      integer iradia
      integer ilwd
      integer ilwb
      integer ilwn
      integer ncouche
      integer solvarmod   ! model for solar EUV variation
      integer solvaryear  ! mars year for realisticly varying solar EUV 
      integer dustinjection ! dust injection scheme number 

      logical rayleigh
      logical tracer
      logical scavenging
      logical rdstorm ! rocket dust storm parametrization
      logical slpwind ! entrainment by slope wind parametrization
      logical latentheat_surfwater ! latent heat release from ground water ice sublimation/condensation
      logical cap_albedo ! polar cap albedo remains unchanged by water frost deposition
      logical temp_dependant_m ! temperature-dependant water contact parameter
      logical sedimentation
      logical activice,tifeedback,supersat,caps
      logical co2clouds,co2useh2o,meteo_flux,CLFvaryingCO2,satindexco2
      logical activeco2ice
      integer spantCO2
      logical CLFvarying
      logical water
      logical hdo
      logical hdofrac
      logical microphys
      logical photochem
      integer nltemodel
      integer nircorr

      character(len=100) dustiropacity
      real               dustrefir 
  
      integer swrtype ! type of short wave (solar wavelength) radiative
      ! transfer to use 1: Fouquart 2: Toon.
      parameter (swrtype=2)
!      parameter (swrtype=2)

!     input/output

      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of tracers
      real,intent(in) :: pplay(ngrid,nlayer)
      real,intent(in) :: pt(ngrid,nlayer)
      real,intent(in) :: pdt(ngrid,nlayer)
      real,intent(in) :: pq(ngrid,nlayer,nq)
      real,intent(in) :: pdq(ngrid,nlayer,nq)
      real,intent(in) :: ptimestep

!     local variables

      integer       :: i, l, ig, iq
      integer, save :: nbq
      integer,allocatable,save :: niq(:)
      real          :: ni(nq), ntot
      real          :: zq(ngrid, nlayer, nq)
      real          :: zt(ngrid, nlayer)
      real,allocatable,save    :: aki(:)
      real,allocatable,save    :: cpi(:)

      logical, save :: firstcall = .true.

      if (firstcall) then

         ! allocate local saved arrays:
         allocate(aki(nq))
         allocate(cpi(nq))
         allocate(niq(nq))
!        find index of chemical tracers to use
!        initialize thermal conductivity and specific heat coefficients
!        !? values are estimated

         nbq = 0 ! to count number of tracers used in this subroutine

         if (igcm_co2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co2
            aki(nbq) = 3.072e-4
            cpi(nbq) = 0.834e3
         end if
         if (igcm_co /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co
            aki(nbq) = 4.87e-4
            cpi(nbq) = 1.034e3
         end if
         if (igcm_o /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o
            aki(nbq) = 7.59e-4
            cpi(nbq) = 1.3e3
         end if
         if (igcm_o1d /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o1d
            aki(nbq) = 7.59e-4  !?
            cpi(nbq) = 1.3e3    !?
         end if
         if (igcm_o2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o2
            aki(nbq) = 5.68e-4
            cpi(nbq) = 0.9194e3
         end if
         if (igcm_o3 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o3
            aki(nbq) = 3.00e-4  !?
            cpi(nbq) = 0.800e3  !?
         end if
         if (igcm_h /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h
            aki(nbq) = 0.0
            cpi(nbq) = 20.780e3
         end if
         if (igcm_h2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h2
            aki(nbq) = 36.314e-4
            cpi(nbq) = 14.266e3
         end if
         if (igcm_oh /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_oh
            aki(nbq)  = 7.00e-4 !?
            cpi(nbq)  = 1.045e3
         end if
         if (igcm_ho2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ho2
            aki(nbq) = 0.0
            cpi(nbq) = 1.065e3  !?
         end if
         if (igcm_n2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2
            aki(nbq) = 5.6e-4
            cpi(nbq) = 1.034e3
         end if
         if (igcm_ar /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_ar
            aki(nbq) = 3.4e-4
            cpi(nbq) = 5.2e2
            ! aki/cpi from Leslie A. Guildner,
            ! JOURNAL OF RESEARCH of the National Bureau of Standards- 
            ! A. Physics and Chemistry
            ! Vol. 79A, No.2, March-April 1975
         end if
         if (igcm_h2o_vap /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_h2o_vap
            aki(nbq) = 0.0
            cpi(nbq) = 1.870e3
         end if
         if (igcm_n /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_no /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_no
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_no2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_no2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_n2d /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2d
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_he /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_he
            aki(nbq) = 30.e-4
            cpi(nbq) = 5.2e3
            ! aki/cpi from Leslie A. Guildner,
            ! JOURNAL OF RESEARCH of the National Bureau of Standards- 
            ! A. Physics and Chemistry
            ! Vol. 79A, No.2, March-April 1975
         endif
         if(igcm_co2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_co2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_oplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_oplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_o2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_o2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_coplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_coplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_cplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_cplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_nplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_nplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_noplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_noplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_n2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_n2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_hplus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_hplus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_hco2plus /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_hco2plus
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(igcm_elec /= 0) then
            nbq = nbq + 1
            niq(nbq) = igcm_elec
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         
         ! tell the world about it:
         write(*,*) "concentrations: firstcall, nbq=",nbq
         write(*,*) "  niq(1:nbq)=",niq(1:nbq)
         write(*,*) "  aki(1:nbq)=",aki(1:nbq)
         write(*,*) "  cpi(1:nbq)=",cpi(1:nbq)

         firstcall = .false.

      end if ! if (firstcall)

!     update temperature

      do l = 1,nlayer
         do ig = 1,ngrid
            zt(ig,l) = pt(ig,l) + pdt(ig,l)*ptimestep
         end do
      end do

!     update tracers

      do l = 1,nlayer
         do ig = 1,ngrid
            do i = 1,nbq
               iq = niq(i) 
               zq(ig,l,iq) = max(1.e-30, pq(ig,l,iq)
     $                                 + pdq(ig,l,iq)*ptimestep)
            end do
         end do
      end do

!     mmean : mean molecular mass
!     rnew  : specific gas constant

      mmean(:,:)  = 0.

      do l = 1,nlayer
         do ig = 1,ngrid
            do i = 1,nbq
               iq = niq(i) 
               mmean(ig,l) = mmean(ig,l) + zq(ig,l,iq)/mmol(iq)
            end do
            mmean(ig,l) = 1./mmean(ig,l)
            rnew(ig,l) = 8.314/mmean(ig,l)*1.e3     ! J/kg/K		
         end do
      end do

!     cpnew  : specicic heat
!     akknew : thermal conductivity cofficient
      
      cpnew(:,:)  = 0.
      akknew(:,:) = 0.

      do l = 1,nlayer
         do ig = 1,ngrid
            ntot = pplay(ig,l)/(1.381e-23*zt(ig,l))*1.e-6  ! in #/cm3
            do i = 1,nbq
               iq = niq(i) 
               ni(iq) = ntot*zq(ig,l,iq)*mmean(ig,l)/mmol(iq)
               cpnew(ig,l) = cpnew(ig,l) + ni(iq)*cpi(i)
               akknew(ig,l) = akknew(ig,l) + ni(iq)*aki(i)
            end do 
            cpnew(ig,l) = cpnew(ig,l)/ntot
            akknew(ig,l) = akknew(ig,l)/ntot
         end do
!        print*, l, mmean(1,l), cpnew(1,l), rnew(1,l)
      end do

      return
      end 
