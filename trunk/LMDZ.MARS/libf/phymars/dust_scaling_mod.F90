module dust_scaling_mod

implicit none

contains

  subroutine compute_dustscaling(ngrid,nlayer,naerkind,naerdust, &
                                 zday,pplev, &
                                 tau_pref_scenario,tauscaling, &
                                 dust_rad_adjust,aerosol)
    
    use dust_param_mod, only: dustscaling_mode, odpref
    use dust_rad_adjust_mod, only: compute_dust_rad_adjust
    use dimradmars_mod, only: iaerdust ! dust aerosol indexes 
    
    implicit none
    
    integer,intent(in) :: ngrid ! number of atmospheric columns
    integer,intent(in) :: nlayer ! number of atmospheric layers
    integer,intent(in) :: naerkind ! total number of aerosols
    integer,intent(in) :: naerdust ! number of dust aerosols
    real,intent(in) :: zday
    real,intent(in) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
    real,intent(in) :: tau_pref_scenario(ngrid) ! prescribed visible dust 
                       ! opacity column at odpref reference pressure
    real,intent(out) :: tauscaling(ngrid) ! dust scaling factor
    real,intent(out) :: dust_rad_adjust(ngrid) ! Radiative adjustment 
                          ! factor for dust
    real,intent(inout) :: aerosol(ngrid,nlayer,naerkind) ! opacities
    
    integer :: ig, l , iaer
    real :: taudust(ngrid)

  ! 1. compute/set tauscaling

    if (dustscaling_mode /= 1) then
      ! simple "freedust" case, no effective rescaling using tauscaling, ever
      tauscaling(:) = 1
    endif
    
    if (dustscaling_mode == 1) then
      ! Compute dust column opacity using aerosol() dusts
      taudust(:) = 0
      do iaer=1,naerdust ! loop on all dust aerosols
        do l=1,nlayer
          do ig=1,ngrid
              taudust(ig)=taudust(ig)+aerosol(ig,l,iaerdust(iaer))
          enddo
        enddo
      enddo

    elseif (dustscaling_mode == 2) then
      ! Compute dust column opacity using only background dust
      taudust(:) = 0
      do l=1,nlayer
        do ig=1,ngrid
            taudust(ig)=taudust(ig)+aerosol(ig,l,iaerdust(1))
        enddo
      enddo

    endif ! of if (dustscaling_mode == 1) elseif (dustscaling_mode == 2)
      
  ! 2. compute the scaling factors (tauscaling or dust_rad_adjust)
    if (dustscaling_mode==1) then
      ! GCM v5.3 style: tauscaling is computed so that
      ! aerosol() opacities correspond to the prescribed tau_pref_scenario()
      tauscaling(:)=tau_pref_scenario(:)*pplev(:,1)/odpref/taudust(:)
    elseif (dustscaling_mode==2) then
      ! GCM v6 style, compute dust_rad_adjust
      call compute_dust_rad_adjust(ngrid,nlayer,zday,pplev, &
                                   taudust,dust_rad_adjust)
    endif

  ! 3. Apply dust aerosol opacities rescaling
    if (dustscaling_mode <=1) then
      do iaer=1,naerdust
        do l=1,nlayer
          do ig=1,ngrid
            aerosol(ig,l,iaerdust(iaer)) = max(1E-20, &
                      aerosol(ig,l,iaerdust(iaer))* tauscaling(ig))
          enddo
        enddo
      enddo
    else ! duscaling_mode==2, use dust_rad_adjust
      do iaer=1,naerdust
        do l=1,nlayer
          do ig=1,ngrid
            aerosol(ig,l,iaerdust(iaer)) = max(1E-20, &
                      aerosol(ig,l,iaerdust(iaer))*dust_rad_adjust(ig))
          enddo
        enddo
      enddo
    endif
  end subroutine compute_dustscaling

end module dust_scaling_mod
