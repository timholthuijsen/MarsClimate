module read_profile_mod

implicit none

contains
!=====================================================================================================================!
! SUBROUTINE: read_profile ===========================================================================================!
!=====================================================================================================================!
! Author: Christophe Mathe
! Date: 09/06/2020
!---------------------------------------------------------------------------------------------------------------------!
! Subject:
!---------
!   Read input profile of tracers listed in traceur.def
!---------------------------------------------------------------------------------------------------------------------!
! Comments:
!----------
!     - If no input profile found, q(traceur) and qsurf(traceur) set to 0, except:
!          - q(co2) = 0.95
!          - q(hdo) = q(h2o)*2*155.76e-6*5
!
!     - igcm is not available at this part of the code
!
!     - To ensure that major isotopologue input profile is read before compute minor isotopologue profile, the
!       calculation is performed at the end of the subroutine
!---------------------------------------------------------------------------------------------------------------------!
! Algorithm:
!-----------
!   1. Initialization of q and qsurf to 0
!   2. Get indices of some tracers
!   3. Main
!     3.1. Try to open the input profile
!       3.1.1. Succeed
!       3.1.2. Fail
!         3.1.2.a. Some cases require a special initialization
!   4. Traitment for minor isotopologue if the major isotopologue input profile does not exist
!=====================================================================================================================!
  subroutine read_profile(nb_tracer, nb_layer, qsurf, q)

  use infotrac, only: tname

  implicit none
!---------------------------------------------------------------------------------------------------------------------!
! VARIABLE DECLARATION
!---------------------------------------------------------------------------------------------------------------------!
!  Input arguments:
!------------------
  integer, intent(in) :: &
     nb_layer, & ! number of layer
     nb_tracer  ! number of traceur read from traceur.def 
!---------------------------------------------------------------------------------------------------------------------!
!  Output arguments:
!-------------------
  real, intent(out) :: &
     qsurf(nb_tracer), &    ! kg/m2
     q(nb_layer, nb_tracer) ! kg/kg of atmosphere
!---------------------------------------------------------------------------------------------------------------------!
!  Local:
!--------
  integer :: &
     iq,             & ! loop over nb_tracer
     ilayer,         & ! loop over nb_layer
     ierr,           & ! open file iostat
     indice_h2o_vap, & ! indice of h2o_vap tracer 
     indice_h2o_ice, & ! indice of h2o_ice tracer
     indice_hdo_vap, & ! indice of hdo_vap tracer
     indice_hdo_ice    ! indice of hdo_ice tracer

  character(len=80), dimension(nb_tracer) :: &
     name_tracer ! array of all tracers already read in traceur.def

  logical :: &
     hdo_vap = .false., & ! used to compute hdo_vap profile if its input profile is missing (= .true)
     hdo_ice = .false.    ! used to compute hdo_ice profile if its input profile is missing (= .true)
!=====================================================================================================================!
!=== BEGIN                                                                                                            !
!=====================================================================================================================!
! 1. Initialization of q and qsurf to 0
!---------------------------------------------------------------------------------------------------------------------!
  q(1:nb_layer,1:nb_tracer) = 0. 
  qsurf(1:nb_tracer) = 0. 
  name_tracer(1:nb_tracer) = ""
!---------------------------------------------------------------------------------------------------------------------!
! 2. Get indices of some tracers: igcm is not available at this part of the code
!---------------------------------------------------------------------------------------------------------------------!
  do iq = 1, nb_tracer
    write(name_tracer(iq),"(a)") tname(iq)
    if (trim(name_tracer(iq)) == 'h2o_vap') then
      indice_h2o_vap = iq
    else if (trim(name_tracer(iq)) == 'h2o_ice') then
      indice_h2o_ice = iq
    else if (trim(name_tracer(iq)) == 'hdo_vap') then
      indice_hdo_vap = iq
    else if (trim(name_tracer(iq)) == 'hdo_ice') then
      indice_hdo_ice = iq
    end if
  end do
!---------------------------------------------------------------------------------------------------------------------!
! 3. Main
!---------------------------------------------------------------------------------------------------------------------!
  do iq = 1, nb_tracer
    write(*,*)"  tracer:",trim(name_tracer(iq))
!---------------------------------------------------------------------------------------------------------------------!
! 3.1. Try to open the input profile
!---------------------------------------------------------------------------------------------------------------------!
    open(91,file='profile_'//trim(name_tracer(iq)), status='old', form='formatted', iostat=ierr)
!---------------------------------------------------------------------------------------------------------------------!
! 3.1.1. Succeed: the input file exists
!---------------------------------------------------------------------------------------------------------------------!
    if (ierr.eq.0) then
      read(91,*)qsurf(iq)
      do ilayer = 1, nb_layer
        read(91,*)q(ilayer,iq)
      end do
!---------------------------------------------------------------------------------------------------------------------!
! 3.1.2. Fail: the input file does not exist
!---------------------------------------------------------------------------------------------------------------------!
    else
      write(*,*)"File profile_"//trim(name_tracer(iq))//" not found!"
!---------------------------------------------------------------------------------------------------------------------!
! 3.1.2.a. Some cases require a special initialization
!---------------------------------------------------------------------------------------------------------------------!
      select case(trim(name_tracer(iq)))
        case("co2")
          q(1:nb_layer, iq) = 0.95

        case("hdo_vap")
          hdo_vap = .true.

        case("hdo_ice")
          hdo_ice = .true.

       end select
    end if
    close(91)
  end do ! of do iq = 1, nq
!---------------------------------------------------------------------------------------------------------------------!
! 4. Traitment for minor isotopologue if the major isotopologue input profile does not exist. At this part, ensure that
!    the main isotopologue profile is already read before minors isotopologues
!---------------------------------------------------------------------------------------------------------------------!
  if (hdo_vap.eqv..true. .and. indice_h2o_vap.ne.0) then
    do ilayer = 1, nb_layer
      q(ilayer, indice_hdo_vap) = q(ilayer, indice_h2o_vap)*2*155.76e-6*5
    end do
  end if

  if (hdo_ice.eqv..true. .and. indice_h2o_ice.ne.0) then
    qsurf(indice_hdo_ice) = qsurf(indice_h2o_ice) * 2*155.76e-6*5
    do ilayer = 1, nb_layer
      q(ilayer, indice_hdo_ice) = q(ilayer, indice_h2o_ice) * 2*155.76e-6*5
    end do
  end if
!=====================================================================================================================!
!=== END                                                                                                              !
!=====================================================================================================================!
  end subroutine read_profile
end module read_profile_mod

