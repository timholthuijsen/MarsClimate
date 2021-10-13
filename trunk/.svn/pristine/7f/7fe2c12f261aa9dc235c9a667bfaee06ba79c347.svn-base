       MODULE compute_dtau_mod

        IMPLICIT NONE

        REAL,SAVE :: ti_injection_sol ! time of beginning injection
        REAL,SAVE :: tf_injection_sol ! time of end injection

        REAL,SAVE,ALLOCATABLE :: dtau(:) ! Dust opacity difference (at 610Pa)
                                         ! between GCM and dust scenario

       CONTAINS

        SUBROUTINE compute_dtau(ngrid,nlayer,                           &
                                 zday,pplev,tau_pref_gcm,               &
                                 ptimestep,dustliftday,local_time)

        USE geometry_mod, only: longitude_deg
        USE time_phylmdz_mod, only: dtphys, daysec
        USE comcstfi_h, only: g
        USE tracer_mod, only: alpha_lift,igcm_dust_mass,igcm_dust_number
        USE dimradmars_mod, only: tauvis
        USE dust_param_mod, only: odpref, t_scenario_sol
        
        IMPLICIT NONE
        
        include "callkeys.h"
        
        INTEGER, INTENT(in) :: ngrid
        INTEGER, INTENT(in) :: nlayer
        REAL, INTENT(in) :: zday ! date at lon=0, in fraction of sols
        REAL, INTENT(in) :: pplev(ngrid,nlayer+1) ! pressure (Pa)
        REAL, INTENT(in) :: tau_pref_gcm(ngrid) ! Visible dust opacity column
                            ! at 610Pa as computed in the GCM
        REAL, INTENT(in) :: ptimestep 
        REAL, INTENT(in) :: local_time(ngrid)
        REAL, INTENT(out) :: dustliftday(ngrid) ! Dust injection rate (s-1)
        
        INTEGER :: ig, l
        INTEGER, SAVE :: nb_daystep ! nomber of step a day
        REAL :: tau_pref_target(ngrid) ! dust opacity column at odpref=610 Pa
                ! as extracted from dust scenario
        REAL :: zday_scenario
        REAL,ALLOCATABLE,SAVE :: local_time_prev(:)
        
        LOGICAL, SAVE :: firstcall=.TRUE. ! signals first call to physics
        
        
        IF(firstcall)THEN
                ALLOCATE(local_time_prev(ngrid))
                DO ig=1,ngrid
                   local_time_prev(ig)=modulo(1.+(zday-ptimestep/daysec)&
                                      -INT(zday-ptimestep/daysec)       &
                                      +(longitude_deg(ig)/15)/24,1.)
                ENDDO
                nb_daystep=(daysec/dtphys)
                ! Local time in sol fraction
                ti_injection_sol=ti_injection/24.
                tf_injection_sol=tf_injection/24.
                firstcall=.FALSE.
        ENDIF
        
        ! 1. Obtain tau_pref_target from dust scenario at zday+1
        if (iaervar.eq.1) then
          tau_pref_target = tauvis
        else
          zday_scenario=zday-modulo(zday,1.) ! integer value of the day: the scenario opacity is measured at 14:00
          zday_scenario=zday_scenario+1      ! opacity of the dust scenario is read the day after
          call read_dust_scenario(ngrid,nlayer,zday_scenario,pplev,     &
                                         tau_pref_target)
        endif
       ! for diagnostics
        call WRITEDIAGFI(ngrid,"tau_pref_target", &
                          "target visible dust opacity column at 610Pa", &
                          "",2,tau_pref_target)

        ! 2. Compute dtau() and dustliftday()
        DO ig=1,ngrid
         IF ((local_time(ig).ge.t_scenario_sol).and.                    &
                 (local_time_prev(ig).lt.(t_scenario_sol)))THEN
                 dtau(ig)=tau_pref_target(ig)-tau_pref_gcm(ig)
         ENDIF

        ! Use dtau (when positive) to compute dustliftday
         IF (dtau(ig).LT.0) THEN
             dustliftday(ig)=0.
         ELSE
             dustliftday(ig)=coeff_injection*                           &
                        (dtau(ig)*pplev(ig,1)/odpref)                   &
                        /(daysec*(tf_injection_sol-ti_injection_sol))
         ENDIF
        ENDDO ! of DO ig=1,ngrid

       ! for diagnostics
        call WRITEDIAGFI(ngrid,"dtau","opacity difference wrt scenario",&
                          "",2,dtau)
        call WRITEDIAGFI(ngrid,"dustliftday","dust injection rate",     &
                          "s-1",2,dustliftday)
         
        ! 4. Save local time 
        local_time_prev(1:ngrid)=local_time(1:ngrid)
        
        end subroutine compute_dtau

!=======================================================================
! Initialization of the module variables

        subroutine ini_compute_dtau_mod(ngrid)
       
          implicit none
       
          integer, intent(in) :: ngrid
       
          allocate(dtau(ngrid))
       
        end subroutine ini_compute_dtau_mod
       
        subroutine end_compute_dtau_mod
       
          implicit none
       
          if (allocated(dtau)) deallocate(dtau)

          end subroutine end_compute_dtau_mod       

       END MODULE compute_dtau_mod
