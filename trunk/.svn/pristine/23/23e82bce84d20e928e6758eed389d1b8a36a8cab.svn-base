!****************************************************************
!
!     Photochemical routine 
!
!     Authors: Franck Lefevre, Francisco Gonzalez-Galindo
!     -------
!
!     Version: 14/11/2020
!
!     ASIS scheme : for details on the method see
!     Cariolle et al., Geosci. Model Dev., 10, 1467-1485, 2017.
!
!*****************************************************************

subroutine photochemistry(nlayer, nq, nesp, ionchem, deutchem,                 &
                          nb_reaction_3_max, nb_reaction_4_max, nphot,         &
                          nb_phot_max, nphotion,                               &
                          jonline, ig, lswitch, zycol, sza, ptimestep, press,  &
                          alt, temp, temp_elect, dens, zmmean,                 &
                          dist_sol, zday,                                      &
                          surfdust1d, surfice1d, jo3, jh2o,em_no,em_o2,        &
                          tau, iter)

use param_v4_h, only: jion

implicit none

include "callkeys.h"

!===================================================================
!     inputs:
!===================================================================

integer, intent(in) :: nlayer ! number of atmospheric layers
integer, intent(in) :: nq     ! number of tracers in traceur.def
integer, intent(in) :: nesp   ! number of traceurs in chemistry
integer, intent(in) :: nb_reaction_3_max   
                              ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max
                              ! number of bimolecular reactions
integer, intent(in) :: nphot
                              ! number of photodissociations
integer, intent(in) :: nb_phot_max
                              ! number of reactions treated numerically as photodissociations
integer, intent(in) :: nphotion
                              ! number of photoionizations
logical, intent(in) :: ionchem! switch for ion chemistry
logical, intent(in) :: deutchem
                              ! switch for deuterium chemistry
logical, intent(in) :: jonline! switch for on-line calculation of photolysis rates
integer :: ig                 ! grid point index
      
real :: sza                   ! solar zenith angle (deg)
real :: ptimestep             ! physics timestep (s)
real :: press(nlayer)         ! pressure (hpa)
real :: alt(nlayer)           ! altitude (km)
real :: temp(nlayer)          ! temperature (k)
real :: temp_elect(nlayer)    ! electronic temperature (K)
real :: dens(nlayer)          ! density (cm-3)
real :: zmmean(nlayer)        ! mean molar mass (g/mole)
real :: dist_sol              ! sun distance (au) 
real :: zday                  ! date (time since Ls=0, in martian days)
real :: surfdust1d(nlayer)    ! dust surface area (cm2/cm3)
real :: surfice1d(nlayer)     ! ice surface area (cm2/cm3)
real :: tau                   ! dust optical depth

!===================================================================
!     input/output:
!===================================================================
      
real :: zycol(nlayer,nq)       ! chemical species volume mixing ratio

!===================================================================
!     output:
!===================================================================
      
integer :: iter(nlayer)        ! iteration counter
real    :: jo3(nlayer)         ! photodissociation rate o3 -> o1d
real    :: jh2o(nlayer)        ! photodissociation rate h2o -> h + oh
real    :: em_no(nlayer)       ! NO nightglow volume emission rate
real    :: em_o2(nlayer)       ! O2 nightglow volume emission rate

!===================================================================
!     local:
!===================================================================

integer :: phychemrat            ! (physical timestep)/(nominal chemical timestep)
integer :: j_o3_o1d, j_h2o, ilev, iesp
integer :: lswitch
logical, save :: firstcall = .true.
logical :: jionos                ! switch for J parameterization

! tracer indexes in the chemistry:
! Species always considered in the chemistry
integer,parameter :: i_co2  =  1
integer,parameter :: i_co   =  2
integer,parameter :: i_o    =  3
integer,parameter :: i_o1d  =  4
integer,parameter :: i_o2   =  5
integer,parameter :: i_o3   =  6
integer,parameter :: i_h    =  7
integer,parameter :: i_h2   =  8
integer,parameter :: i_oh   =  9
integer,parameter :: i_ho2  = 10
integer,parameter :: i_h2o2 = 11
integer,parameter :: i_h2o  = 12
integer,parameter :: i_n    = 13
integer,parameter :: i_n2d  = 14
integer,parameter :: i_no   = 15
integer,parameter :: i_no2  = 16
integer,parameter :: i_n2   = 17
!Species that may or not be included
!Their indices will be defined, if species are to be considered, in firstcall
integer,save      :: i_co2plus  = 0
integer,save      :: i_oplus    = 0
integer,save      :: i_o2plus   = 0
integer,save      :: i_noplus   = 0
integer,save      :: i_coplus   = 0
integer,save      :: i_cplus    = 0
integer,save      :: i_n2plus   = 0
integer,save      :: i_nplus    = 0
integer,save      :: i_hplus    = 0
integer,save      :: i_hco2plus = 0
integer,save      :: i_hcoplus  = 0
integer,save      :: i_h2oplus  = 0
integer,save      :: i_h3oplus  = 0
integer,save      :: i_ohplus   = 0
integer,save      :: i_elec     = 0
integer,save      :: i_hdo      = 0
integer,save      :: i_od       = 0
integer,save      :: i_d        = 0
integer,save      :: i_hd       = 0
integer,save      :: i_do2      = 0
integer,save      :: i_hdo2     = 0

!integer,parameter :: i_co2plus = 18
!integer,parameter :: i_oplus   = 19
!integer,parameter :: i_o2plus  = 20
!integer,parameter :: i_noplus  = 21
!integer,parameter :: i_coplus  = 22
!integer,parameter :: i_cplus   = 23
!integer,parameter :: i_n2plus  = 24
!integer,parameter :: i_nplus   = 25
!integer,parameter :: i_hplus   = 26
!integer,parameter :: i_hco2plus= 27
!integer,parameter :: i_hcoplus = 28
!integer,parameter :: i_h2oplus = 29
!integer,parameter :: i_h3oplus = 30
!integer,parameter :: i_ohplus  = 31
!integer,parameter :: i_elec    = 32
!integer,parameter :: i_hdo     = 33
!integer,parameter :: i_od      = 34
!integer,parameter :: i_d       = 35
!integer,parameter :: i_hd      = 36
!integer,parameter :: i_do2     = 37
!integer,parameter :: i_hdo2    = 38

integer :: i_last

!Tracer indexes for photionization coeffs

integer,parameter :: induv_co2 = 1
integer,parameter :: induv_o2  = 2
integer,parameter :: induv_o   = 3
integer,parameter :: induv_h2o = 4
integer,parameter :: induv_h2  = 5
integer,parameter :: induv_h2o2= 6
integer,parameter :: induv_o3  = 7
integer,parameter :: induv_n2  = 8
integer,parameter :: induv_n   = 9
integer,parameter :: induv_no  = 10
integer,parameter :: induv_co  = 11
integer,parameter :: induv_h   = 12
integer,parameter :: induv_no2 = 13

integer :: ilay

real :: ctimestep           ! standard timestep for the chemistry (s) 
real :: dt_guess            ! first-guess timestep (s) 
real :: dt_corrected        ! corrected timestep (s) 
real :: time                ! internal time (between 0 and ptimestep, in s)

real, dimension(nlayer,nesp)            :: rm   ! mixing ratios 
real (kind = 8), dimension(nesp)        :: cold ! number densities at previous timestep (molecule.cm-3) 
real (kind = 8), dimension(nlayer,nesp) :: c    ! number densities at current timestep (molecule.cm-3) 
real (kind = 8), dimension(nesp)        :: cnew ! number densities at next timestep (molecule.cm-3) 
 
! reaction rates
 
real (kind = 8), dimension(nlayer,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nlayer,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nlayer,nb_reaction_4_max) :: v_4
logical :: hetero_dust, hetero_ice

! matrix

real (kind = 8), dimension(nesp,nesp) :: mat, mat1
integer, dimension(nesp)              :: indx
integer                               :: code

! production and loss terms (for first-guess solution only)

real (kind = 8), dimension(nesp) :: prod, loss

!===================================================================
!     initialisation of the reaction indexes
!===================================================================

if (firstcall) then
   !Indices for species that may or not be considered
   i_last = i_n2
   if(ionchem) then
      i_co2plus  = i_last + 1
      i_oplus    = i_last + 2
      i_o2plus   = i_last + 3
      i_noplus   = i_last + 4
      i_coplus   = i_last + 5
      i_cplus    = i_last + 6
      i_n2plus   = i_last + 7
      i_nplus    = i_last + 8
      i_hplus    = i_last + 9
      i_hco2plus = i_last + 10
      i_hcoplus  = i_last + 11
      i_h2oplus  = i_last + 12
      i_h3oplus  = i_last + 13
      i_ohplus   = i_last + 14
      i_elec     = i_last + 15
      !Update last used index
      i_last = i_elec
   endif
   if(deutchem) then
      i_hdo  = i_last + 1
      i_od   = i_last + 2
      i_d    = i_last + 3
      i_hd   = i_last + 4
      i_do2  = i_last + 5
      i_hdo2 = i_last + 6
      i_last = i_hdo2
   endif
   print*,'photochemistry: check number of indices',i_last,nesp
   print*,'photochemistry: initialize indexes'
   call indice(nb_reaction_3_max,nb_reaction_4_max,                 &
               nb_phot_max, ionchem, deutchem, i_co2, i_co, i_o,    &
               i_o1d, i_o2, i_o3, i_h,i_h2, i_oh, i_ho2, i_h2o2,    &
               i_h2o, i_n, i_n2d, i_no, i_no2, i_n2, i_co2plus,     &
               i_oplus, i_o2plus, i_noplus, i_coplus, i_cplus,      &
               i_n2plus, i_nplus, i_hplus, i_hco2plus, i_hcoplus,   &
               i_h2oplus, i_h3oplus, i_ohplus, i_elec,              &
               i_hdo, i_od, i_d, i_hd, i_do2, i_hdo2)
   firstcall = .false.
end if

!===================================================================
!     initialisation of mixing ratios and densities       
!===================================================================

call gcmtochim(nlayer, ionchem, deutchem, nq, zycol, lswitch,  &
               nesp, i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h, &
               i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
               i_n, i_n2d, i_no, i_no2, i_n2, i_co2plus,       &
               i_oplus, i_o2plus, i_noplus, i_coplus, i_cplus, &
               i_n2plus, i_nplus, i_hplus, i_hco2plus,         &
               i_hcoplus, i_h2oplus, i_h3oplus, i_ohplus,      &
               i_elec, i_hdo, i_od, i_d, i_hd, i_do2, i_hdo2,  &
               dens, rm, c) 

!===================================================================
!     photolysis rates
!===================================================================

jionos = .true.

if (jonline) then
   if (sza <= 113.) then ! day at 300 km
      call photolysis_online(nlayer, deutchem, nb_phot_max, alt,        &
                             press, temp, zmmean,                       &
                             i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,  &
                             i_h2, i_oh, i_ho2, i_h2o2, i_h2o,          &
                             i_n, i_n2d, i_no, i_no2, i_n2, nesp, rm,   &
                             tau, sza, dist_sol, v_phot)

      !Calculation of photoionization rates, if needed
      if (jionos .and. ionchem) then
         call jthermcalc_e107(ig,nlayer,2,c,nesp,temp,alt,sza,zday)
         do ilay=1,lswitch-1
            call phdisrate(ig,nlayer,2,sza,ilay)
         end do
         !CO2 photoionization
         v_phot(:,nphot+ 1) = jion(induv_co2,:,1)
         v_phot(:,nphot+ 2) = jion(induv_co2,:,2)
         v_phot(:,nphot+ 3) = jion(induv_co2,:,2)
         v_phot(:,nphot+ 4) = jion(induv_co2,:,3)
         v_phot(:,nphot+ 5) = jion(induv_co2,:,3)
         v_phot(:,nphot+ 6) = jion(induv_co2,:,4)
         v_phot(:,nphot+ 7) = jion(induv_co2,:,4)
         !O2 photoionization
         v_phot(:,nphot+ 8) = jion(induv_o2,:,1)
         !O photoionization
         v_phot(:,nphot+ 9) = jion(induv_o,:,1)
         !NO photoionization
         v_phot(:,nphot+10) = jion(induv_no,:,1)
         !CO photoionization
         v_phot(:,nphot+11) = jion(induv_co,:,1)
         v_phot(:,nphot+12) = jion(induv_co,:,2)
         v_phot(:,nphot+13) = jion(induv_co,:,2)
         !N2 photoionization
         v_phot(:,nphot+14) = jion(induv_n2,:,1)
         v_phot(:,nphot+15) = jion(induv_n2,:,2)
         v_phot(:,nphot+16) = jion(induv_n2,:,2)
         !N photoionization
         v_phot(:,nphot+17) = jion(induv_n,:,1)
         !H photoionization
         v_phot(:,nphot+18) = jion(induv_h,:,1)
      end if
      
   else ! night
      v_phot(:,:) = 0.
   end if
!else if(jparam) then
!   call jthermcalc_e107(ig,nlayer,2,c,nesp,temp,alt,sza,zday)
!   do ilay=1,lswitch-1
!      call phdisrate(ig,nlayer,2,sza,ilay)
!   enddo
!   v_phot(:,1)=jdistot(2,:)
!   v_phot(:,2)=jdistot_b(2,:)
!   v_phot(:,3)=jdistot(1,:)
!   v_phot(:,4)=jdistot_b(1,:)
!   v_phot(:,5)=jdistot(7,:)
!   v_phot(:,6)=jdistot_b(7,:)
!   v_phot(:,7)=jdistot(4,:)
!   v_phot(:,8)=jdistot(6,:)
!   v_phot(:,10)=jdistot(5,:)
!   v_phot(:,11)=jdistot(10,:)
!   v_phot(:,12)=jdistot(13,:)
!   v_phot(:,13)=jdistot(8,:)
!   v_phot(:,14)=jion(1,:,1)
!   v_phot(:,15)=jion(1,:,2)
!   v_phot(:,16)=jion(1,:,2)
!   v_phot(:,17)=jion(1,:,3)
!   v_phot(:,18)=jion(1,:,3)
!   v_phot(:,19)=jion(1,:,4)
!   v_phot(:,20)=jion(1,:,4)
!   v_phot(:,21)=jion(2,:,1)
!   v_phot(:,22)=jion(3,:,1)
!   v_phot(:,23)=jion(10,:,1)
!   v_phot(:,24)=jion(11,:,1)
!   v_phot(:,25)=jion(11,:,2)
!   v_phot(:,26)=jion(11,:,2)
!   v_phot(:,27)=jion(8,:,1)
!   v_phot(:,28)=jion(8,:,2)
!   v_phot(:,29)=jion(8,:,2)
!   v_phot(:,30)=jion(9,:,1)
!   v_phot(:,31)=jion(12,:,1)
else
   tau = tau*7./press(1) ! dust in the lookup table is at 7 hpa
   call photolysis(nlayer, nb_phot_max, lswitch, press, temp, sza, tau,  &
                   zmmean, dist_sol,rm(:,i_co2), rm(:,i_o3), v_phot)
end if
! save o3 and h2o photolysis for output

j_o3_o1d = 5
jo3(:) = v_phot(:,j_o3_o1d)
j_h2o = 7
jh2o(:) = v_phot(:,j_h2o)

!===================================================================
!     reaction rates                                     
!===================================================================
!     switches for heterogeneous chemistry
!     hetero_ice  : reactions on ice clouds
!     hetero_dust : reactions on dust    
!===================================================================

hetero_dust = .false.
hetero_ice  = .false.

call reactionrates(nlayer, ionchem, deutchem,                             &
                   nb_reaction_3_max, nb_reaction_4_max,                  &
                   nb_phot_max, nphotion, lswitch, dens, c(:,i_co2),      &
                   c(:,i_o2), c(:,i_o), c(:,i_n2), press, temp,           &
                   temp_elect, hetero_dust, hetero_ice,                   &
                   surfdust1d, surfice1d, v_phot, v_3, v_4)

!===================================================================
!     ctimestep : standard chemical timestep (s), defined as 
!                 the fraction phychemrat of the physical timestep                           
!===================================================================

phychemrat = 1

ctimestep = ptimestep/real(phychemrat)

!print*, "ptimestep  = ", ptimestep
!print*, "phychemrat = ", phychemrat
!print*, "ctimestep  = ", ctimestep
!stop

!===================================================================
!     loop over levels         
!===================================================================

do ilev = 1,lswitch - 1

!  initializations

   time = 0.
   iter(ilev) = 0
   dt_guess = ctimestep
   cold(:) = c(ilev,:)

!  internal loop for the chemistry

   do while (time < ptimestep)

   iter(ilev) = iter(ilev) + 1    ! iteration counter
  
!  first-guess: fill matrix

   call fill_matrix(ilev, mat1, prod, loss, c, nesp, nlayer,              &
                    nb_reaction_3_max, nb_reaction_4_max, nb_phot_max,    &
                    v_phot, v_3, v_4)

!  adaptative evaluation of the sub time step

   call define_dt(nesp, dt_corrected, dt_guess, ctimestep, cold(:), c(ilev,:),  &
                  mat1, prod, loss, dens(ilev))

   if (time + dt_corrected > ptimestep) then
      dt_corrected = ptimestep - time
   end if

!  if (dt_corrected /= dt_guess) then  ! the timestep has been modified

!  form the matrix identity + mat*dt_corrected

   mat(:,:) = mat1(:,:)*dt_corrected
   do iesp = 1,nesp
      mat(iesp,iesp) = 1. + mat(iesp,iesp)
   end do

!  solve the linear system  M*Cn+1 = Cn (RHS in cnew, then replaced by solution)
   cnew(:) = c(ilev,:)

#ifdef LAPACK
   call dgesv(nesp,1,mat,nesp,indx,cnew,nesp,code)
#else
   write(*,*) "photochemistry error, missing LAPACK routine dgesv"
   call abort_physic("photochemistry","missing LAPACK routine",1)
#endif

!  end if

!  eliminate small values

   where (cnew(:)/dens(ilev) < 1.e-30)
      cnew(:) = 0.
   end where

!  update concentrations

   cold(:)   = c(ilev,:)
   c(ilev,:) = cnew(:)

!  force charge neutrality (mod fgg, july 2019)

   if (ionchem) then
      if(c(ilev,i_elec).ne.c(ilev,i_co2plus)+c(ilev,i_oplus)+c(ilev,i_o2plus)+&
           c(ilev,i_noplus)+c(ilev,i_coplus)+c(ilev,i_cplus)+c(ilev,i_n2plus)+&
           c(ilev,i_nplus)+c(ilev,i_hplus)+c(ilev,i_hco2plus)+                &
           c(ilev,i_hcoplus)+c(ilev,i_h2oplus)+c(ilev,i_h3oplus)+             &
           c(ilev,i_ohplus)) then
         c(ilev,i_elec) = c(ilev,i_co2plus)+c(ilev,i_oplus)+c(ilev,i_o2plus)+ &
              c(ilev,i_noplus)+c(ilev,i_coplus)+c(ilev,i_cplus)+              &
              c(ilev,i_n2plus)+c(ilev,i_nplus)+c(ilev,i_hplus)+               &
              c(ilev,i_hco2plus)+c(ilev,i_hcoplus)+c(ilev,i_h2oplus)+         &
              c(ilev,i_h3oplus)+c(ilev,i_ohplus)
         !      write(*,*)'photochemistry/359'
         !      write(*,*)'Forcing charge neutrality at ilev,',ilev,' ig=',ig
      end if
   end if

   cnew(:)   = 0.

!  increment internal time

   time = time + dt_corrected
   dt_guess = dt_corrected     ! first-guess timestep for next iteration

   end do ! while (time < ptimestep)

end do ! ilev

!add calculation of NO and O2 nightglows
em_no(:)=c(:,i_o)*c(:,i_n)*v_4(:,26)   !2.8e-17*(300./temp(:)))**0.5
em_o2(:)=0.75*c(:,i_o)*c(:,i_o)*c(:,i_co2)*v_3(:,1)   !2.5*9.46e-34*exp(485./temp(:))*dens(:)

!===================================================================
!     save chemical species for the gcm       
!===================================================================

call chimtogcm(nlayer, ionchem, deutchem, nq, zycol, lswitch,  &
               nesp, i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h, &
               i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
               i_n, i_n2d, i_no, i_no2, i_n2, i_co2plus,       &
               i_oplus, i_o2plus, i_noplus, i_coplus, i_cplus, &
               i_n2plus, i_nplus, i_hplus, i_hco2plus,         &
               i_hcoplus, i_h2oplus, i_h3oplus, i_ohplus,      &
               i_elec, i_hdo, i_od, i_d, i_hd, i_do2, i_hdo2,  &
               dens, c) 
contains

!================================================================

 subroutine define_dt(nesp, dtnew, dtold, ctimestep, cold, ccur, mat1, &
                      prod, loss, dens)

!================================================================
! iterative evaluation of the appropriate time step dtnew
! according to curvature criterion based on
! e = 2 Rtol [r Cn+1 -(1-r) Cn + Cn-1 ]/[(1+r) Cn]
! with r = (tn - tn-1)/(tn+1 - tn)
!================================================================

implicit none

! input

integer :: nesp  ! number of species in the chemistry

real :: dtold, ctimestep
real (kind = 8), dimension(nesp)      :: cold, ccur
real (kind = 8), dimension(nesp,nesp) :: mat1
real (kind = 8), dimension(nesp)      :: prod, loss
real                                  :: dens

! output

real :: dtnew

! local

real (kind = 8), dimension(nesp)      :: cnew
real (kind = 8), dimension(nesp,nesp) :: mat
real (kind = 8) :: atol, ratio, e, es, coef

integer                  :: code, iesp, iter
integer, dimension(nesp) :: indx

real :: dttest

! parameters

real (kind = 8), parameter :: dtmin   = 10.      ! minimum time step (s)
real (kind = 8), parameter :: vmrtol  = 1.e-11   ! absolute tolerance on vmr
real (kind = 8), parameter :: rtol    = 0.05     ! rtol recommended value : 0.1-0.02
integer,         parameter :: niter   = 3        ! number of iterations
real (kind = 8), parameter :: coefmax = 2.
real (kind = 8), parameter :: coefmin = 0.1 
logical                    :: fast_guess = .true.

dttest = dtold   ! dttest = dtold = dt_guess

atol = vmrtol*dens ! absolute tolerance in molecule.cm-3

do iter = 1,niter

if (fast_guess) then

! first guess : fast semi-implicit method

   do iesp = 1, nesp
      cnew(iesp) = (ccur(iesp) + prod(iesp)*dttest)/(1. + loss(iesp)*dttest)
   end do

else

! first guess : form the matrix identity + mat*dt_guess

   mat(:,:) = mat1(:,:)*dttest
   do iesp = 1,nesp
      mat(iesp,iesp) = 1. + mat(iesp,iesp)
   end do

! form right-hand side (RHS) of the system

   cnew(:) = ccur(:)

! solve the linear system  M*Cn+1 = Cn (RHS in cnew, then replaced by solution)

#ifdef LAPACK
   call dgesv(nesp,1,mat,nesp,indx,cnew,nesp,code)
#else
   write(*,*) "photochemistry error, missing LAPACK routine dgesv"
   call abort_physic("photochemistry","missing LAPACK routine",1)
#endif

end if

! ratio old/new subtimestep

ratio = dtold/dttest

! e : local error indicatocitr

e = 0.

do iesp = 1,nesp
   es = 2.*abs((ratio*cnew(iesp) - (1. + ratio)*ccur(iesp) + cold(iesp))   &
         /(1. + ratio)/max(ccur(iesp)*rtol,atol))

   if (es > e) then
      e = es
   end if
end do

! timestep correction

coef = max(coefmin, min(coefmax,0.8/sqrt(e)))

dttest = max(dtmin,dttest*coef)
dttest = min(ctimestep,dttest)

end do ! iter

! new timestep

dtnew = dttest

end subroutine define_dt

!======================================================================

 subroutine reactionrates(nlayer, ionchem, deutchem ,                         &
                          nb_reaction_3_max, nb_reaction_4_max, &
                          nb_phot_max, nphotion, lswitch, dens, co2, o2, o,      &
                          n2, press, t, t_elect, hetero_dust, hetero_ice,        &
                          surfdust1d, surfice1d, v_phot, v_3, v_4)
 
!================================================================
! compute reaction rates                                        !
!----------------------------------------------------------------
! reaction               type                array              !
!----------------------------------------------------------------
! A + B    --> C + D     bimolecular         v_4                !
! A + A    --> B + C     quadratic           v_3                !
! A + C    --> B + C     quenching           v_phot             !
! A + ice  --> B + C     heterogeneous       v_phot             !
!================================================================

use comcstfi_h
use photolysis_mod, only : nphot

implicit none

!----------------------------------------------------------------------
!     input
!----------------------------------------------------------------------

integer, intent(in)     :: nlayer            ! number of atmospheric layers
integer, intent(in)     :: nb_reaction_3_max ! number of quadratic reactions
integer, intent(in)     :: nb_reaction_4_max ! number of bimolecular reactions
integer, intent(in)     :: nb_phot_max       ! number of reactions treated numerically as photodissociations
integer, intent(in)     :: nphotion          ! number of photoionizations
logical, intent(in)     :: ionchem
logical, intent(in)     :: deutchem
integer                 :: lswitch           ! interface level between lower
                                             ! atmosphere and thermosphere chemistries
real, dimension(nlayer) :: dens              ! total number density (molecule.cm-3)
real, dimension(nlayer) :: press             ! pressure (hPa)
real, dimension(nlayer) :: t                 ! temperature (K)
real, dimension(nlayer) :: t_elect           ! electronic temperature (K)
real, dimension(nlayer) :: surfdust1d        ! dust surface area (cm2.cm-3)
real, dimension(nlayer) :: surfice1d         ! ice surface area (cm2.cm-3)
real (kind = 8), dimension(nlayer) :: co2    ! co2 number density (molecule.cm-3)
real (kind = 8), dimension(nlayer) :: o2     ! o2 number density (molecule.cm-3)
real (kind = 8), dimension(nlayer) :: o      ! o number density (molecule.cm-3)
real (kind = 8), dimension(nlayer) :: n2     ! n2 number density (molecule.cm-3)
logical :: hetero_dust, hetero_ice           ! switches for heterogeneous chemistry

!----------------------------------------------------------------------
!     output
!----------------------------------------------------------------------

real (kind = 8), dimension(nlayer,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nlayer,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nlayer,nb_reaction_4_max) :: v_4

!----------------------------------------------------------------------
!     local
!----------------------------------------------------------------------

integer          :: ilev
integer          :: nb_phot, nb_reaction_3, nb_reaction_4
real :: ak0, ak1, xpo, rate, rate1, rate2
real :: k1a0, k1b0, k1ainf, k1a, k1b, fc, fx, x, y, gam
real, dimension(nlayer) :: deq
real, dimension(nlayer) :: a001, a002, a003,                           &
                           b001, b002, b003, b004, b005, b006, b007,   &
                           b008, b009, b010, b011, b012,               &
                           c001, c002, c003, c004, c005, c006, c007,   &
                           c008, c009, c010, c011, c012, c013, c014,   &
                           c015, c016, c017, c018,                     &
                           d001, d002, d003, d004, d005, d006, d007,   &
                           d008, d009, d010, d011, d012,               &
                           e001, e002,                                 &
                           f001, f002, f003, f004, f005, f006, f007,   &
                           f008, f009, f010, f011, f012, f013, f014,   &
                           f015, f016, f017, f018, f019, f020, f021,   &
                           f022, f023, f024, f025, f026, f027, f028,   &
                           f029, f030, f031, f032,                     &
                           i001, i002, i003, i004, i005, i006,         &
                           i007, i008, i009, i010, i011, i012,         &
                           i013, i014, i015, i016, i017, i018, i019,   &
                           i020, i021, i022, i023, i024, i025, i026,   &
                           i027, i028, i029, i030, i031, i032, i033,   &
                           i034, i035, i036, i037, i038, i039, i040,   &
                           i041, i042, i043, i044, i045, i046, i047,   &
                           i048, i049, i050, i051, i052, i053, i054,   &
                           i055, i056, i057, i058, i059, i060, i061,   &
                           i062, i063,                                 &
                           h001, h002, h003, h004, h005

!----------------------------------------------------------------------
!     initialisation
!----------------------------------------------------------------------

      nb_phot       = nphot + nphotion ! initialised to the number of photolysis + number of photoionization rates
      nb_reaction_3 = 0
      nb_reaction_4 = 0

!----------------------------------------------------------------------
!     oxygen reactions 
!----------------------------------------------------------------------

!---  a001: o + o2 + co2 -> o3 + co2

!     jpl 2003
!
!     co2/n2 efficiency as a third body = 2.075
!     from sehested et al., j. geophys. res., 100, 1995.

      a001(:) = 2.075*6.0e-34*(t(:)/300.)**(-2.4)*dens(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a001(:)

!---  a002: o + o + co2 -> o2 + co2

!     Tsang and Hampson, J. Chem. Phys. Ref. Data, 15, 1087, 1986

!     a002(:) = 2.5*5.2e-35*exp(900./t(:))*dens(:)

!     Campbell and Gray, Chem. Phys. Lett., 18, 607, 1973

!     a002(:) = 1.2e-32*(300./t(:))**(2.0)*dens(:)  ! yung expression
      a002(:) = 2.5*9.46e-34*exp(485./t(:))*dens(:) ! nist expression

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = a002(:)

!---  a003: o + o3 -> o2 + o2

!     jpl 2003

      a003(:) = 8.0e-12*exp(-2060./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a003(:)

!----------------------------------------------------------------------
!     o(1d) reactions
!----------------------------------------------------------------------

!---  b001: o(1d) + co2  -> o + co2

!     jpl 2006

      b001(:) = 7.5e-11*exp(115./t(:))
   
      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b001(:)*co2(:)

!---  b002: o(1d) + h2o  -> oh + oh

!     jpl 2006
 
      b002(:) = 1.63e-10*exp(60./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b002(:)

!---  b003: o(1d) + h2  -> oh + h

!     jpl 2011

      b003(:) = 1.2e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b003(:)

!---  b004: o(1d) + o2  -> o + o2

!     jpl 2006

      b004(:) = 3.3e-11*exp(55./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b004(:)*o2(:)
    
!---  b005: o(1d) + o3  -> o2 + o2

!     jpl 2003

      b005(:) = 1.2e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b005(:)
    
!---  b006: o(1d) + o3  -> o2 + o + o

!     jpl 2003

      b006(:) = 1.2e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b006(:)
    
!---  b007: o(1d) + ch4 -> ch3 + oh

!     jpl 2003

      b007(:) = 1.5e-10*0.75

!---  b008: o(1d) + ch4 -> ch3o + h

!     jpl 2003

      b008(:) = 1.5e-10*0.20
!
!---  b009: o(1d) + ch4 -> ch2o + h2

!     jpl 2003

      b009(:) = 1.5e-10*0.05

      if(deutchem) then
!
!---     b010: o(1d) + hdo -> od + oh
         b010(:) = b002(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = b010(:)

!
!---     b011: o(1d) + hd -> h + od

         !Laurent et al., 1995

         b011(:) = 1.3e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = b011(:)

!
!---     b012: o(1d) + hd -> d + oh

         !Laurent et al., 1995

         b012(:) = 1.0e-10
         
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = b012(:)

      endif  !deutchem
!----------------------------------------------------------------------
!     hydrogen reactions
!----------------------------------------------------------------------

!---  c001: o + ho2 -> oh + o2

!     jpl 2003

      c001(:) = 3.0e-11*exp(200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c001(:)

!---  c002: o + oh -> o2 + h

!     jpl 2011

      c002(:) = 1.8e-11*exp(180./t(:))

!     c002(:) = c002(:)*1.12 ! FL li et al. 2007

!     robertson and smith, j. chem. phys. a 110, 6673, 2006

!     c002(:) = 11.2e-11*t(:)**(-0.32)*exp(177./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c002(:)

!---  c003: h + o3 -> oh + o2

!     jpl 2003

      c003(:) = 1.4e-10*exp(-470./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c003(:)

!---  c004: h + ho2 -> oh + oh

!     jpl 2006

      c004(:) = 7.2e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c004(:)

!---  c005: h + ho2 -> h2 + o2

!     jpl 2006

      c005(:) = 6.9e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c005(:)

!---  c006: h + ho2 -> h2o + o

!     jpl 2006

      c006(:) = 1.6e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c006(:)

!---  c007: oh + ho2 -> h2o + o2

!     jpl 2003

!     canty et al., grl, 2006 suggest to increase this rate
!     by 20%. not done here.

      c007(:) = 4.8e-11*exp(250./t(:))

!     c007(:) = c007(:)*0.9 ! FL li et al. 2007

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c007(:)

!---  c008: ho2 + ho2 -> h2o2 + o2

!     jpl 2015

      c008(:) = 3.0e-13*exp(460./t(:))

!     christensen et al., grl, 13, 2002

!     c008(:) = 1.5e-12*exp(19./t(:))

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c008(:)

!---  c009: oh + h2o2 -> h2o + ho2

!     jpl 2006

      c009(:) = 1.8e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c009(:)

!---  c010: oh + h2 -> h2o + h

!     jpl 2006

      c010(:) = 2.8e-12*exp(-1800./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c010(:)

!---  c011: h + o2 + co2 -> ho2 + co2

!     jpl 2011
!     co2/n2 efficiency as a third body = 2.4
!     from ashman and haynes, 27th symposium on combustion, 1998.

      do ilev = 1,lswitch-1
!        ak0 = 3.1*2.4*4.4e-32*(t(ilev)/300.)**(-1.3) ! FL li et al 2017
         ak0 = 2.4*4.4e-32*(t(ilev)/300.)**(-1.3)
         ak1 = 7.5e-11*(t(ilev)/300.)**(0.2)

         rate = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2)
         c011(ilev) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c011(:)

!---  c012: o + h2o2 -> oh + ho2

!     jpl 2003

      c012(:) = 1.4e-12*exp(-2000./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c012(:)

!---  c013: oh + oh -> h2o + o

!     jpl 2006

      c013(:) = 1.8e-12

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c013(:)

!---  c014: oh + o3 -> ho2 + o2

!     jpl 2003

      c014(:) = 1.7e-12*exp(-940./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c014(:)

!---  c015: ho2 + o3 -> oh + o2 + o2

!     jpl 2003

      c015(:) = 1.0e-14*exp(-490./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c015(:)

!---  c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2

!     jpl 2011

      c016(:) = 2.5*2.1e-33*exp(920./t(:))*dens(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c016(:)

!---  c017: oh + oh + co2 -> h2o2 + co2

!     jpl 2003

      do ilev = 1,lswitch-1
         ak0 = 2.5*6.9e-31*(t(ilev)/300.)**(-1.0)
         ak1 = 2.6e-11*(t(ilev)/300.)**(0.0)

         rate = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2)
         c017(ilev) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c017(:)

!---  c018: h + h + co2 -> h2 + co2

!     baulch et al., 2005

      c018(:) = 2.5*1.8e-30*(t(:)**(-1.0))*dens(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c018(:)


!----------------------------------------------------------------------
!     nitrogen reactions
!----------------------------------------------------------------------

!---  d001: no2 + o -> no + o2

!     jpl 2006

      d001(:) = 5.1e-12*exp(210./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d001(:)

!---  d002: no + o3 -> no2 + o2

!     jpl 2006

      d002(:) = 3.0e-12*exp(-1500./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d002(:)

!---  d003: no + ho2 -> no2 + oh

!     jpl 2011

      d003(:) = 3.3e-12*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d003(:)

!---  d004: n + no -> n2 + o

!     jpl 2011

      d004(:) = 2.1e-11*exp(100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d004(:)

!---  d005: n + o2 -> no + o

!     jpl 2011

      d005(:) = 1.5e-11*exp(-3600./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d005(:)

!---  d006: no2 + h -> no + oh

!     jpl 2011

      d006(:) = 4.0e-10*exp(-340./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d006(:)

!---  d007: n + o -> no
 
      d007(:) = 2.8e-17*(300./t(:))**0.5

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d007(:)

!---  d008: n + ho2 -> no + oh

!     brune et al., j. chem. phys., 87, 1983 

      d008(:) = 2.19e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d008(:)

!---  d009: n + oh -> no + h

!     atkinson et al., j. phys. chem. ref. data, 18, 881, 1989 

      d009(:) = 3.8e-11*exp(85./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d009(:)

!---  d010: n(2d) + o  -> n + o

!     herron, j. phys. chem. ref. data, 1999

      d010(:) = 3.3e-12*exp(-260./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = d010(:)*o(:)

!---  d011: n(2d) + n2  -> n + n2

!     herron, j. phys. chem. ref. data, 1999

      d011(:) = 1.7e-14

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = d011(:)*n2(:)

!---  d012: n(2d) + co2  -> no + co

!     herron, j. phys. chem. ref. data, 1999

      d012(:) = 3.6e-13

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d012(:)

!----------------------------------------------------------------------
!     carbon reactions
!----------------------------------------------------------------------

!---  e001: oh + co -> co2 + h

!     jpl 2003

!     e001(:) = 1.5e-13*(1 + 0.6*press(:)/1013.)

!     mccabe et al., grl, 28, 3135, 2001

!     e001(:) = 1.57e-13 + 3.54e-33*dens(:)

!     jpl 2015

      do ilev = 1,lswitch-1

!        branch 1 : oh + co -> h + co2

         rate1 = 1.5e-13*(t(ilev)/300.)**(0.0)

!        branch 2 : oh + co + m -> hoco + m

         ak0 = 5.9e-33*(t(ilev)/300.)**(-1.0)
         ak1 = 1.1e-12*(t(ilev)/300.)**(1.3)
         rate2 = (ak0*dens(ilev))/(1. + ak0*dens(ilev)/ak1)
         xpo = 1./(1. + alog10((ak0*dens(ilev))/ak1)**2)

         e001(ilev) = rate1 + rate2*0.6**xpo
      end do

!     joshi et al., 2006

!     do ilev = 1,lswitch-1
!        k1a0 = 1.34*2.5*dens(ilev)                                  &
!              *1/(1/(3.62e-26*t(ilev)**(-2.739)*exp(-20./t(ilev)))  &
!              + 1/(6.48e-33*t(ilev)**(0.14)*exp(-57./t(ilev))))     ! typo in paper corrected
!        k1b0 = 1.17e-19*t(ilev)**(2.053)*exp(139./t(ilev))          &
!             + 9.56e-12*t(ilev)**(-0.664)*exp(-167./t(ilev))
!        k1ainf = 1.52e-17*t(ilev)**(1.858)*exp(28.8/t(ilev))        &
!               + 4.78e-8*t(ilev)**(-1.851)*exp(-318./t(ilev))
!        x = k1a0/(k1ainf - k1b0)
!        y = k1b0/(k1ainf - k1b0)
!        fc = 0.628*exp(-1223./t(ilev)) + (1. - 0.628)*exp(-39./t(ilev))  &
!           + exp(-t(ilev)/255.)
!        fx = fc**(1./(1. + (alog(x))**2))                           ! typo in paper corrected
!        k1a = k1a0*((1. + y)/(1. + x))*fx
!        k1b = k1b0*(1./(1.+x))*fx
!        e001(ilev) = k1a + k1b
!     end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e001(:)

!---  e002: o + co + m -> co2 + m

!     tsang and hampson, 1986.

      e002(:) = 2.5*6.5e-33*exp(-2184./t(:))*dens(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e002(:)


!----------------------------------------------------------------------
!     deuterium reactions
!     only if deutchem = true
!----------------------------------------------------------------------

      if(deutchem) then

!---     f001: od + oh -> hdo + o

         ! Bedjanian, Y.; Le Bras, G.; and Poulet, G., 
         !Kinetic Study of OH + OH and OD + OD Reactions, 
         !J. Phys. Chem. A, 103, pp. 7017 - 7025, 1999

         f001(:) = 1.5e-12
      
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f001(:)

!---     f002: od + h2 -> HDO + H

         !Talukdar, R.K. et al., 
         !Kinetics of hydroxyl radical reactions with isotopically labeled hydrogen, 
         !J. Phys. Chem., 100, pp.   3037 - 3043, 1996

         f002(:) = 7.41e-15

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f002(:)

!---     f003: od + ho2 -> hdo + o2
      
         f003(:) = c007(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f003(:)

!---     f004: od + h2o2 -> hdo + ho2

         !Vaghjiani, G.L. & Ravishankara, A.R., 
         !Reactions of OH and OD with H2O2 and D2O2, 
         !J. Phys. Chem., 93, pp. 7833 - 7837, 1989;

         f004(:) = 1.79e-12

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f004(:)

!---     f005: o + od -> o2 + d

         !Following Yung+1998, rate equal to that of O + OH -> O2 + H (c002)

         f005(:) = c002(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f005(:)

!---     f006: od + h2 -> h2o + d

         !Following Yung+1998, rate equal to that of OH + H2 -> H2O + H (c010)

         f006(:) = c010(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f006(:)

!---     f007: od + h -> oh + d

         !Rate following Yung+1988 and the rate of the inverse reaction (f012) from Atahan et al. J. Chem. Phys. 2005

         f007(:) = 1.61e-10*((298./t(:))**0.32)*exp(-16./t(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f007(:)

!---     f008: co + od -> co2 + d

         !Following Yung+1988 rate equal to that of reaction CO + OH -> CO2 + H

         f008(:) = e001(:) 

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f008(:)

!---     f009: o3 + d -> o2 + od

         !Rate from NIST, Yu & Varandas, J. Chem. Soc. Faraday Trans., 93, 2651-2656, 1997

         f009(:) = 7.41e-11*exp(-379./t(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f009(:)

!---     f010: HO2 + D -> OH + OD

         !Following Yung+1988, rate equal to 0.71 times the rate for reaction HO2 + H -> OH + OH (c004)

         f010(:) = 0.71*c004(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f010(:)

!---     f011: HO2 + D -> HDO + O

         !Following Yung+1988, rate equal to 0.71 times the rate for reaction HO2 + H -> H2O + O (c006)

         f011(:) = 0.71*c006(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f011(:)

!---     f012: OH + D -> OD + H

         !Rate from NIST, Atahan et al. J. Chem. Phys. 2005

         f012(:) = 1.16e-10*((298./t(:))**0.32)*exp(-16./t(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f012(:)

!---     f013: h + d + co2 -> hd + co2

         !According to Yung et al. 1988, rate equal to that of H + H + CO2 
         !(reaction c018). Source: baulch et al., 2005

         f013(:) = c018(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f013(:)

!---     f014: D + HO2 -> HD + O2
      
         !According to Yung et al., rate equal to 0.71 times the rate of
         ! H + HO2 -> H2 + O2 (reaction c005, source JPL 2019)

         f014(:) = 0.71*c005(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f014(:)
      
!---     f015: OH + HD -> HDO + H

         !Talukdar et al., Kinetics of hydroxyl radical reactions with 
         !isotopically labeled hydrogen, J. Phys. Chem. 100, 3037-3043, 1996

         f015(:) = 8.5e-13*exp(-2130./t(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f015(:)

!---     f016: OH + HD -> H2O + D

         !Talukdar et al., 1996

         f016(:) = 4.15e-12*exp(-2130./t(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f016(:)

!---     f017: D + O2 + CO2 -> DO2 + CO2

         !Breshears et al., Room-temperature rate constants for the reaction 
         !D + O2 + M -> DO2 + M, M=Ar, D2, CO2 and F2. J. Chem. Soc. Faraday 
         !Trans., 87, 2337-2355 (1991)

         f017(:) = 1.59e-31

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f017(:)

!---     f018: OD + O3 -> DO2 + O2

         !According to Yung+1988, rate equal to that of reaccion 
         !OH + O3 -> HO2 + O2, (reaction c014)

         f018(:) = c014(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f018(:)

!---     f019: D + HO2 -> DO2 + H

         !Yung et al., 1988

         f019(:) = 1.0e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f019(:)

!---     f020: O + DO2 -> OD + O2

         !According to Yung+1988, rate equal to that of O + HO2 -> OH + O2
         ! -> reaction c001

         f020(:) = c001(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f020(:)
      
!---     f021: H + DO2 -> OH + OD

         !According to Yung+1988, rate equal to that of H + HO2 -> OH + OH
         ! -> reaction c004

         f021(:) = c004(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f021(:)

!---     f022: H + DO2 -> HD + O2

         !According to Yung+1988, rate equal to that of H + HO2 -> H2 + O2
         ! -> reaction c005

         f022(:) = c005(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f022(:)

!---     f023: H + DO2 -> HDO + O

         !According to Yung+1988, rate equal to that of H + HO2 -> H2O + O
         ! -> reaction c006

         f023(:) = c006(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f023(:)

!---     f024: H + DO2 -> HO2 + D

         !Yung+1988

         f024(:) = 1.85e-10*exp(-890./t(:))
         
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f024(:)

!---     f025: OH + DO2 -> HDO + O2

         !According to Yung+1988, rate equal to that of OH + HO2 -> H2O + O2
         ! -> reaction c007

         f025(:) = c007(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f025(:)

!---     f026: DO2 + O3 -> OD + O2 + O2

         !According to Yung+1988, rate equal to that of the reaction
         ! HO2 + O3 -> OH + O2 + O2 -> reaction c015

         f026(:) = c015(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f026(:)

!---     f027: OD + OH + CO2 -> HDO2 + CO2

         !According to Yung+1988, rate equal to that of the reaction
         ! OH + OH + CO2 -> H2O2 + CO2 (reaction c017)

         f027(:) = c017(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f027(:)

!---     f028: DO2 + HO2 -> HDO2 + O2

         !According to Yung+1988, rate equal to that of HO2 + HO2 -> H2O2 + O2
         ! (reaction c008)

         f028(:) = c008(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f028(:)

!---     f029: O + HDO2 -> OD + HO2

         !According to Yung+1988, rate half that of O + H2O2 -> OH + HO2
         ! (reaction c012)

         f029(:) = 0.5*c012(:)
      
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f029(:)

!---     f030: O + HDO2 -> OH + DO2

         !According to Yung+1988, rate half that of O + H2O2 -> OH + HO2
         ! (reaction c012)

         f030(:) = 0.5*c012(:)
      
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f030(:)

!---     f031: OH + HDO2 -> HDO + HO2

         !According to Yung+1988, rate half that of OH + H2O2 -> H2O + HO2
         ! (reaction c009)

         f031(:) = 0.5*c009(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f031(:)


!---     f032: OH + HDO2 -> H2O + DO2

         !According to Yung+1988, rate half that of OH + H2O2 -> H2O + HO2
         ! (reaction c009)

         f032(:) = 0.5*c009(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = f032(:)

      endif !deutchem

!----------------------------------------------------------------------
!     ionospheric reactions
!     only if ionchem=true
!----------------------------------------------------------------------

      if (ionchem) then

!---     i001: co2+ + o2 -> o2+ + co2

!        aninich, j. phys. chem. ref. data 1993

         i001(:) = 5.5e-11*(300./t_elect(:))**0.5

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i001(:)

!---     i002: co2+ + o -> o+ + co2

!        UMIST database

         i002(:) = 9.6e-11
      
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i002(:)

!---     i003: co2+ + o -> o2+ + co

!        UMIST database

         i003(:) = 1.64e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i003(:)

!---     i004: o2+ + e- -> o + o

!        Alge et al., J. Phys. B, At. Mol. Phys. 1983

         i004(:) = 2.0e-7*(300./t_elect(:))**0.7

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i004(:)

!---     i005: o+ + co2 -> o2+ + co

!        UMIST database

         i005(:) = 9.4e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i005(:)


!---     i006: co2+ + e- -> co + o

!        UMIST database

         i006(:) = 3.8e-7*(300./t_elect(:))**0.5

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i006(:)


!---     i007: co2+ + no -> no+ + co2

!        UMIST database

         i007(:) = 1.2e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i007(:)

!---     i008: o2+ + no -> no+ + o2

!        UMIST database

         i008(:) = 4.6e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i008(:)

!---     i009: o2+ + n2 -> no+ + no
      
!        Fox & Sung 2001

         i009(:) = 1.0e-15
      
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i009(:)

!---     i010: o2+ + n -> no+ + o

!        Fox & Sung 2001

         i010(:) = 1.0e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i010(:)

!---     i011: o+ + n2 -> no+ + n

!        Fox & Sung 2001

         i011(:) = 1.2e-12 * (300./t_elect(:))**0.45

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i011(:)

!---     i012: no+ + e -> n + o

!        UMIST database

         i012(:) = 4.3e-7*(300./t_elect(:))**0.37

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i012(:)


!---     i013: co+ + co2 -> co2+ + co

!        UMIST database

         i013(:) = 1.0e-9

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i013(:)


!---     i014: co+ + o -> o+ + co

!        UMIST database

         i014(:) = 1.4e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i014(:)

!---     i015: c+ + co2 -> co+ + co

!        UMIST database

         i015(:) = 1.1e-9

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i015(:)


!---     i016: N2+ + co2 -> co2+ + N2

!        Fox & Song 2001

         i016(:) = 9.0e-10*(300./t_elect(:))**0.23

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i016(:)


!---     i017: N2+ + o -> no+ + N

!        Fox & Song 2001

         i017(:) = 1.33e-10*(300./t_elect(:))**0.44

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i017(:)

!---     i018: N2+ + co -> co+ + N2

!        UMIST

         i018(:) = 7.4e-11

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i018(:)

!---     i019: N2+ + e -> N + N

!        UMIST

         i019(:) = 7.7e-7*(300./t_elect(:))**0.3

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i016(:)

!---     i020: N2+ + o -> o+ + N2

!        Fox & Song 2001

         i020(:) = 7.0e-12*(300./t_elect(:))**0.23

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i020(:)

!---     i021: N+ + co2 -> co2+ + N

!        UMIST

         i021(:) = 7.5e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i021(:)

!---     i022: CO+ + H -> H+ + CO

!        Fox & Sung 2001

         i022(:) = 4.0e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i022(:)

!---     i023: O+ + H -> H+ + O

!        UMIST

         i023(:) = 5.66e-10*((t_elect(:)/300.)**0.36)*exp(8.6/t_elect(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i023(:)

!---     i024: H+ + O -> O+ + H

!        UMIST

         i024(:) = 6.86e-10*((t_elect(:)/300.)**0.26)*exp(-224.3/t_elect(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i024(:)

!---     i025: CO+ + H2 -> HCO2+ + H

!        UMIST

         i025(:) = 9.5e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i025(:)

!---     i026: HCO2+ + e -> H + CO2

!        UMIST

         i026(:) = 1.75e-8*((300./t_elect(:))**0.5)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i026(:)

!---     i027+i028: HCO2+ + e -> H + O + CO

!        UMIST
         !Reaction splitted in 2: i027: 0.5 (HCO2+ + e-) -> H
         !i028: 0.5 (HCO2+ + e-) -> O + CO

         i027(:) = 8.1e-7*((300./t_elect(:))**0.64)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i027(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i027(:)

!---     i029: HCO2+ + e -> OH + CO

!        UMIST

         i029(:) = 3.2e-7*((300./t_elect(:))**0.64)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i029(:)

!---     i030: HCO2+ + e -> H + CO2

         i030(:) = 6.0e-8*((300./t_elect(:))**0.64)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i030(:)

!---     i031: HCO2+ + O -> HCO+ + O2

!        UMIST

         i031(:) = 1.e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i031(:)

!---     i032: HCO2+ + CO -> HCO+ + CO2

!        UMIST, from Prassad & Huntress 1980

         i032(:) = 7.8e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i032(:)

!---     i033: H+ + CO2 -> HCO+ + O

!        UMIST, from Smith et al., Int. J. Mass Spectrom. Ion Proc., 117, 457-473(1992) 

         i033(:) = 3.5e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i033(:)


!---     i034: CO2+ + H -> HCO+ + O

!        Seen in Fox 2015, from Borodi et al., Int. J. Mass Spectrom. 280, 218-225, 2009

         i034(:) = 4.5e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i034(:)

!---     i035: CO+ + H2 -> HCO+ + H

         !UMIST, from Scott et al., J. Chem. Phys., 106, 3982-3987(1997)

         i035(:) = 7.5e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i035(:)

!---     i036: HCO+ + e- -> CO + H

         !UMIST, from Mitchell, Phys. Rep., 186, 215 (1990)

         i036(:) = 2.4e-7 *((300./t_elect(:))**0.69)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i036(:)

!---     i037: CO2+ + H2O -> H2O+ + CO2

         !UMIST, from Karpas, Z., Anicich, V.G., and Huntress, W.T., Chem. Phys. Lett., 59, 84 (1978)

         i037(:) = 2.04e-9 *((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i037(:)

!---     i038: CO+ + H2O -> H2O+ + CO

         !UMIST, from Huntress, W.T., McEwan, M.J., Karpas, Z., and Anicich, V.G., Astrophys. J. Supp. Series, 44, 481 (1980)

         i038(:) = 1.72e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i038(:)

!---     i039: O+ + H2O -> H2O+ + O

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i039(:) = 3.2e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i039(:)

!---     i040: N2+ + H2O -> H2O+ + N2

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i040(:) = 2.3e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i040(:)

!---     i041: N+ + H2O -> H2O+ + N

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i041(:) = 2.8e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i041(:)


!---     i042: H+ + H2O -> H2O+ + H

         !UMIST, from D. Smith, P. Spanel and C. A. Mayhew, Int. J. Mass Spectrom. Ion Proc., 117, 457-473(1992)

         i042(:) = 6.9e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i042(:)

!---     i043: H2O+ + O2 -> O2+ + H2O

         !UMIST, from A. B. Raksit and P. Warneck, J. Chem. Soc. Faraday Trans., 76, 1084-1092(1980)

         i043(:) = 4.6e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i043(:)

!---     i044: H2O+ + CO -> HCO+ + OH

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i044(:) = 5.0e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i044(:)

!---     i045: H2O+ + O -> O2+ + H2

         !UMIST, from Viggiano, A.A, Howarka, F., Albritton, D.L., Fehsenfeld, F.C., Adams, N.G., and Smith, D., Astrophys. J., 236, 492 (1980)

         i045(:) = 4.0e-11
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i045(:)

!---     i046: H2O+ + NO -> NO+ + H2O

         !UMIST, from A. B. Raksit and P. Warneck, J. Chem. Soc. Faraday Trans., 76, 1084-1092(1980)

         i046(:) = 2.7e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i046(:)

!---     i047: H2O+ + e- -> H + H + O

         !UMIST, from Rosen, S., Derkatch, A., Semaniak, J., et al., 2000, Far. Disc., 115, 295
         
         i047(:) = 3.05e-7*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i047(:)

!---     i048: H2O+ + e- -> H + OH
         
         !UMIST, from Rosen, S., Derkatch, A., Semaniak, J., et al., 2000, Far. Disc., 115, 295

         i048(:) = 8.6e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i048(:)

!---     i049: H2O+ + e- -> O + H2

         !UMIST, from Rosen, S., Derkatch, A., Semaniak, J., et al., 2000, Far. Disc., 115, 295

         i049(:) = 3.9e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i049(:)

!---     i050: H2O+ + H2O -> H3O+ + OH

         !UMIST, from Huntress, W.T. and Pinizzotto, R.F., J. Chem. Phys., 59, 4742 (1973)

         i050(:) = 2.1e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i050(:)


!---     i051: H2O+ + H2 -> H3O+ + H

         !UMIST, from A. B. Raksit and P. Warneck, J. Chem. Soc. Faraday Trans., 76, 1084-1092(1980)

         i051(:) = 6.4e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i051(:)

!---     i052: HCO+ + H2O -> H3O+ + CO

         !UMIST, from Adams, N.G., Smith, D., and Grief, D., Int. J. Mass Spectrom. Ion Phys., 26, 405 (1978)

         i052(:) = 2.5e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i052(:)

!---     i053: H3O+ + e -> OH + H + H

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874

         i053(:) = 3.05e-7*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i053(:)

!---     i054: H3O + + e -> H2O + H

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874
         
         i054(:) = 7.09e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i054(:)

!---     i055: H3O+ + e -> OH + H2

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874

         i055(:) = 5.37e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i055(:)

!---     i056: H3O+ + e -> O + H2 + H

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874

         i056(:) = 5.6e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i056(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i056(:)

!---     i057: O+ + H2 -> OH+ + H 

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i057(:) = 1.7e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i057(:)

!---     i058: OH+ + O -> O2+ + H

         !UMIST, from Prasad & Huntress, 1980, ApJS, 43, 1

         i058(:) = 7.1e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i058(:)

!---     i059: OH+ + CO2 -> HCO2+ + O

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i059(:) = 1.44e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i059(:)

!---     i060: OH+ + CO -> HCO+ + O

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i060(:) = 1.05e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i060(:)

!---     i061: OH+ + NO -> NO+ + OH (tasa de reacciÃ³n UMIST 3.59e-10)

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i061(:) = 3.59e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i061(:)

!---     i062: OH+ + H2 -> H2O+ + H (tasa de reacciÃ³n UMIST 1.01e-9,

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i062(:) = 1.01e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i062(:)

!---     i063: OH+ + O2 -> O2+ + OH (tasa de reacciÃ³n UMIST 5.9e-10

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i063(:) = 5.9e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i063(:)

      end if   !ionchem

!----------------------------------------------------------------------
!     heterogeneous chemistry 
!----------------------------------------------------------------------

      if (hetero_ice) then

!        k = (surface*v*gamma)/4 (s-1)
!        v = 100*sqrt(8rt/(pi*m))  (cm s-1)
 
!---     h001: ho2 + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.025
      
         gam = 0.025
         h001(:) = surfice1d(:)       &
                   *100.*sqrt(8.*8.31*t(:)/(33.e-3*pi))*gam/4.
 
!        h002: oh + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.03
 
         gam = 0.03
         h002(:) = surfice1d(:)       &
                   *100.*sqrt(8.*8.31*t(:)/(17.e-3*pi))*gam/4.

!---     h003: h2o2 + ice -> products
 
!        gamma = 0.    test value
 
         gam = 0.
         h003(:) = surfice1d(:)       &
                   *100.*sqrt(8.*8.31*t(:)/(34.e-3*pi))*gam/4.
      else
         h001(:) = 0.
         h002(:) = 0.
         h003(:) = 0.
      end if

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h001(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h002(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h003(:)

      if (hetero_dust) then
 
!---     h004: ho2 + dust -> products
 
!        jacob, 2000: gamma = 0.2
!        see dereus et al., atm. chem. phys., 2005
 
         gam = 0.2
         h004(:) = surfdust1d(:)  &
                   *100.*sqrt(8.*8.31*t(:)/(33.e-3*pi))*gam/4.
 
!---     h005: h2o2 + dust -> products
 
!        gamma = 5.e-4
!        see dereus et al., atm. chem. phys., 2005
 
         gam = 5.e-4
         h005(:) = surfdust1d(:)  &
                   *100.*sqrt(8.*8.31*t(:)/(34.e-3*pi))*gam/4.
      else
         h004(:) = 0.
         h005(:) = 0.
      end if
 
      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h004(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h005(:)

end subroutine reactionrates

!======================================================================

 subroutine fill_matrix(ilev, mat, prod, loss, c, nesp, nlayer,            &
                        nb_reaction_3_max, nb_reaction_4_max, nb_phot_max, &
                        v_phot, v_3, v_4)

!======================================================================
! filling of the jacobian matrix
!======================================================================

use types_asis

implicit none

! input

integer             :: ilev    ! level index
integer             :: nesp    ! number of species in the chemistry
integer, intent(in) :: nlayer  ! number of atmospheric layers
integer, intent(in) :: nb_reaction_3_max 
                               ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max
                               ! number of bimolecular reactions
integer, intent(in) :: nb_phot_max
                               ! number of processes treated numerically as photodissociations

real (kind = 8), dimension(nlayer,nesp)              :: c    ! number densities
real (kind = 8), dimension(nlayer,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nlayer,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nlayer,nb_reaction_4_max) :: v_4

! output

real (kind = 8), dimension(nesp,nesp), intent(out) :: mat  ! matrix
real (kind = 8), dimension(nesp), intent(out)      :: prod, loss

! local

integer :: iesp
integer :: ind_phot_2,ind_phot_4,ind_phot_6
integer :: ind_3_2,ind_3_4,ind_3_6
integer :: ind_4_2,ind_4_4,ind_4_6,ind_4_8
integer :: iphot,i3,i4

real(kind = 8) :: eps, eps_4  ! implicit/explicit coefficient

! initialisations 

mat(:,:) = 0.
prod(:)  = 0.
loss(:)  = 0.

! photodissociations
! or reactions a + c -> b + c
! or reactions a + ice -> b + c

do iphot = 1,nb_phot_max

  ind_phot_2 = indice_phot(iphot)%z2
  ind_phot_4 = indice_phot(iphot)%z4
  ind_phot_6 = indice_phot(iphot)%z6

  mat(ind_phot_2,ind_phot_2) = mat(ind_phot_2,ind_phot_2) + indice_phot(iphot)%z1*v_phot(ilev,iphot)
  mat(ind_phot_4,ind_phot_2) = mat(ind_phot_4,ind_phot_2) - indice_phot(iphot)%z3*v_phot(ilev,iphot)
  mat(ind_phot_6,ind_phot_2) = mat(ind_phot_6,ind_phot_2) - indice_phot(iphot)%z5*v_phot(ilev,iphot)

  loss(ind_phot_2) = loss(ind_phot_2) + indice_phot(iphot)%z1*v_phot(ilev,iphot)
  prod(ind_phot_4) = prod(ind_phot_4) + indice_phot(iphot)%z3*v_phot(ilev,iphot)*c(ilev,ind_phot_2)
  prod(ind_phot_6) = prod(ind_phot_6) + indice_phot(iphot)%z5*v_phot(ilev,iphot)*c(ilev,ind_phot_2)

end do

! reactions a + a -> b + c 

do i3 = 1,nb_reaction_3_max

  ind_3_2 = indice_3(i3)%z2
  ind_3_4 = indice_3(i3)%z4
  ind_3_6 = indice_3(i3)%z6

  mat(ind_3_2,ind_3_2) = mat(ind_3_2,ind_3_2) + indice_3(i3)%z1*v_3(ilev,i3)*c(ilev,ind_3_2)
  mat(ind_3_4,ind_3_2) = mat(ind_3_4,ind_3_2) - indice_3(i3)%z3*v_3(ilev,i3)*c(ilev,ind_3_2)
  mat(ind_3_6,ind_3_2) = mat(ind_3_6,ind_3_2) - indice_3(i3)%z5*v_3(ilev,i3)*c(ilev,ind_3_2)

  loss(ind_3_2) = loss(ind_3_2) + indice_3(i3)%z1*v_3(ilev,i3)*c(ilev,ind_3_2)
  prod(ind_3_4) = prod(ind_3_4) + indice_3(i3)%z3*v_3(ilev,i3)*c(ilev,ind_3_2)*c(ilev,ind_3_2)
  prod(ind_3_6) = prod(ind_3_6) + indice_3(i3)%z5*v_3(ilev,i3)*c(ilev,ind_3_2)*c(ilev,ind_3_2)

end do

! reactions a + b -> c + d 

eps = 1.d-10

do i4 = 1,nb_reaction_4_max

  ind_4_2 = indice_4(i4)%z2
  ind_4_4 = indice_4(i4)%z4
  ind_4_6 = indice_4(i4)%z6
  ind_4_8 = indice_4(i4)%z8

  eps_4 = abs(c(ilev,ind_4_2))/(abs(c(ilev,ind_4_2)) + abs(c(ilev,ind_4_4)) + eps)
  eps_4 = min(eps_4,1.0)

  mat(ind_4_2,ind_4_2) = mat(ind_4_2,ind_4_2) + indice_4(i4)%z1*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4) 
  mat(ind_4_2,ind_4_4) = mat(ind_4_2,ind_4_4) + indice_4(i4)%z1*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)
  mat(ind_4_4,ind_4_2) = mat(ind_4_4,ind_4_2) + indice_4(i4)%z3*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_4,ind_4_4) = mat(ind_4_4,ind_4_4) + indice_4(i4)%z3*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)   
  mat(ind_4_6,ind_4_2) = mat(ind_4_6,ind_4_2) - indice_4(i4)%z5*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_6,ind_4_4) = mat(ind_4_6,ind_4_4) - indice_4(i4)%z5*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)
  mat(ind_4_8,ind_4_2) = mat(ind_4_8,ind_4_2) - indice_4(i4)%z7*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_8,ind_4_4) = mat(ind_4_8,ind_4_4) - indice_4(i4)%z7*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)


  loss(ind_4_2) = loss(ind_4_2) + indice_4(i4)%z1*v_4(ilev,i4)*c(ilev,ind_4_4)
  loss(ind_4_4) = loss(ind_4_4) + indice_4(i4)%z3*v_4(ilev,i4)*c(ilev,ind_4_2)
  prod(ind_4_6) = prod(ind_4_6) + indice_4(i4)%z5*v_4(ilev,i4)*c(ilev,ind_4_2)*c(ilev,ind_4_4)
  prod(ind_4_8) = prod(ind_4_8) + indice_4(i4)%z7*v_4(ilev,i4)*c(ilev,ind_4_2)*c(ilev,ind_4_4)

end do

end subroutine fill_matrix

!================================================================

 subroutine indice(nb_reaction_3_max, nb_reaction_4_max,                &
                   nb_phot_max, ionchem, deutchem, i_co2, i_co, i_o,    &
                   i_o1d, i_o2, i_o3, i_h,i_h2, i_oh, i_ho2, i_h2o2,    &
                   i_h2o, i_n, i_n2d, i_no, i_no2, i_n2, i_co2plus,     &
                   i_oplus, i_o2plus, i_noplus, i_coplus, i_cplus,      &
                   i_n2plus, i_nplus, i_hplus, i_hco2plus, i_hcoplus,   &
                   i_h2oplus, i_h3oplus, i_ohplus, i_elec,              &
                   i_hdo, i_od, i_d, i_hd, i_do2, i_hdo2)

!================================================================
! set the "indice" arrays used to fill the jacobian matrix      !
!----------------------------------------------------------------
! reaction                                   type               !
!----------------------------------------------------------------
! A + hv   --> B + C     photolysis          indice_phot        ! 
! A + B    --> C + D     bimolecular         indice_4           !
! A + A    --> B + C     quadratic           indice_3           !
! A + C    --> B + C     quenching           indice_phot        !
! A + ice  --> B + C     heterogeneous       indice_phot        !
!================================================================

use types_asis

implicit none

! input

integer :: i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
           i_n, i_n2d, i_no, i_no2, i_n2,                   &
           i_co2plus, i_oplus, i_o2plus, i_noplus, i_coplus, &
           i_cplus, i_n2plus, i_nplus, i_hplus, i_hco2plus, &
           i_hcoplus, i_h2oplus, i_h3oplus, i_ohplus, i_elec, &
           i_hdo, i_od, i_d, i_hd, i_do2, i_hdo2
integer, intent(in) :: nb_reaction_3_max
                       ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max
                       ! number of bimolecular reactions
integer, intent(in) :: nb_phot_max
                       ! number of processes treated numerically as photodissociations
logical, intent(in) :: ionchem
logical, intent(in) :: deutchem

! local

integer :: nb_phot, nb_reaction_3, nb_reaction_4
integer :: i_dummy

allocate (indice_phot(nb_phot_max))
allocate (indice_3(nb_reaction_3_max))
allocate (indice_4(nb_reaction_4_max))

i_dummy = 1

nb_phot       = 0
nb_reaction_3 = 0
nb_reaction_4 = 0

!===========================================================
!      O2 + hv -> O + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2, 2.0, i_o, 0.0, i_dummy) 

!===========================================================
!      O2 + hv -> O + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2, 1.0, i_o, 1.0, i_o1d) 

!===========================================================
!      CO2 + hv -> CO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co, 1.0, i_o) 

!===========================================================
!      CO2 + hv -> CO + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co, 1.0, i_o1d) 

!===========================================================
!      O3 + hv -> O2 + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o3, 1.0, i_o2, 1.0, i_o1d) 

!===========================================================
!      O3 + hv -> O2 + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o3, 1.0, i_o2, 1.0, i_o) 

!===========================================================
!      H2O + hv -> H + OH
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o, 1.0, i_h, 1.0, i_oh) 

!===========================================================
!      H2O2 + hv -> OH + OH
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 2.0, i_oh, 0.0, i_dummy) 

!===========================================================
!      HO2 + hv -> OH + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 1.0, i_oh, 1.0, i_o) 

!===========================================================
!      H2 + hv -> H + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2, 1.0, i_h, 1.0, i_h) 

!===========================================================
!      NO + hv -> N + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_no, 1.0, i_n, 1.0, i_o) 

!===========================================================
!      NO2 + hv -> NO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_no2, 1.0, i_no, 1.0, i_o) 

!===========================================================
!      N2 + hv -> N + N
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_n2, 1.0, i_n2d, 1.0, i_n) 

!Only if deuterium chemistry included
if (deutchem) then
!===========================================================
!      HDO + hv -> OD + H
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_hdo, 1.0, i_h, 1.0, i_od) 

!===========================================================
!      HDO + hv -> D + OH
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_hdo, 1.0, i_d, 1.0, i_oh) 
   
endif !deutchem

!Only if ion chemistry included
if (ionchem) then

!===========================================================
!      CO2 + hv -> CO2+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co2plus, 1.0, i_elec)

!===========================================================
!      CO2 + hv -> O+ + CO + e-
!===========================================================
!We divide this reaction in two

!0.5 CO2 + hv -> CO
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_co, 0.0, i_dummy)

!0.5 CO2 + hv -> O+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_oplus, 1.0, i_elec)

!===========================================================
!      CO2 + hv -> CO+ + O + e-
!===========================================================
!We divide this reaction in two

!0.5 CO2 + hv -> O
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_o, 0.0, i_dummy)

!0.5 CO2 + hv -> CO+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_coplus, 1.0, i_elec)

!===========================================================
!      CO2 + hv -> C+ + O2 + e-
!===========================================================
!We divide this reaction in two

!0.5 CO2 + hv -> O2
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_o2, 0.0, i_dummy)

!0.5 CO2 + hv -> C+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_cplus, 1.0, i_elec)

!===========================================================
!      O2 + hv -> O2+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_o2, 1.0, i_o2plus, 1.0, i_elec)

!===========================================================
!      O + hv -> O+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_o, 1.0, i_oplus, 1.0, i_elec)

!===========================================================
!      NO + hv -> NO+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_no, 1.0, i_noplus, 1.0, i_elec)

!===========================================================
!      CO + hv -> CO+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_co, 1.0, i_coplus, 1.0, i_elec)

!===========================================================
!      CO + hv -> C+ + O + e-
!===========================================================
!We divide this reaction in two

!0.5 CO + hv -> O
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co, 1.0, i_o, 0.0, i_dummy)

!0.5 CO + hv -> C+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co, 1.0, i_cplus, 1.0, i_elec)

!===========================================================
!      N2 + hv -> N2+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_n2, 1.0, i_n2plus, 1.0, i_elec)

!===========================================================
!      N2 + hv -> N+ + N + e-
!===========================================================
!We divide this reaction in two

!0.5 N2 + hv -> N
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_n2, 1.0, i_n, 0.0, i_dummy)

!0.5 N2 + hv -> N+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_n2, 1.0, i_nplus, 1.0, i_elec)

!===========================================================
!      N + hv -> N+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_n, 1.0, i_nplus, 1.0, i_elec)

!===========================================================
!      H + hv -> H+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_h, 1.0, i_hplus, 1.0, i_elec)

end if   !ionchem

!===========================================================
!      a001 : O + O2 + CO2 -> O3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_o2, 1.0, i_o3, 0.0, i_dummy) 

!===========================================================
!      a002 : O + O + CO2 -> O2 + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_o, 1.0, i_o2, 0.0, i_dummy) 

!===========================================================
!      a003 : O + O3 -> O2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_o3, 2.0, i_o2, 0.0, i_dummy) 

!===========================================================
!      b001 : O(1D) + CO2 -> O + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o1d, 1.0, i_o, 0.0, i_dummy) 

!===========================================================
!      b002 : O(1D) + H2O -> OH + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_h2o, 2.0, i_oh, 0.0, i_dummy) 

!===========================================================
!      b003 : O(1D) + H2 -> OH + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_h2, 1.0, i_oh, 1.0, i_h) 

!===========================================================
!      b004 : O(1D) + O2 -> O + O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o1d, 1.0, i_o, 0.0, i_dummy) 

!===========================================================
!      b005 : O(1D) + O3 -> O2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_o3, 2.0, i_o2, 0.0, i_dummy) 

!===========================================================
!      b006 : O(1D) + O3 -> O2 + O + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_o3, 1.0, i_o2, 2.0, i_o) 


if(deutchem) then
!===========================================================
!      b010 : O(1D) + HDO -> OD + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_hdo, 1.0, i_od, 1.0, i_oh)


!===========================================================
!      b011 : O(1D) + HD -> H + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_hd, 1.0, i_h, 1.0, i_od)


!===========================================================
!      b012 : O(1D) + HD -> D + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_hd, 1.0, i_d, 1.0, i_oh)

endif !deutchem

!===========================================================
!      c001 : O + HO2 -> OH + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_ho2, 1.0, i_oh, 1.0, i_o2) 

!===========================================================
!      c002 : O + OH -> O2 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_oh, 1.0, i_o2, 1.0, i_h) 

!===========================================================
!      c003 : H + O3 -> OH + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_o3, 1.0, i_oh, 1.0, i_o2) 

!===========================================================
!      c004 : H + HO2 -> OH + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 2.0, i_oh, 0.0, i_dummy) 

!===========================================================
!      c005 : H + HO2 -> H2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 1.0, i_h2, 1.0, i_o2) 

!===========================================================
!      c006 : H + HO2 -> H2O + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 1.0, i_h2o, 1.0, i_o) 

!===========================================================
!      c007 : OH + HO2 -> H2O + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_ho2, 1.0, i_h2o, 1.0, i_o2) 

!===========================================================
!      c008 : HO2 + HO2 -> H2O2 + O2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_ho2, 1.0, i_h2o2, 1.0, i_o2) 

!===========================================================
!      c009 : OH + H2O2 -> H2O + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_h2o2, 1.0, i_h2o, 1.0, i_ho2) 

!===========================================================
!      c010 : OH + H2 -> H2O + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_h2, 1.0, i_h2o, 1.0, i_h) 

!===========================================================
!      c011 : H + O2 + CO2 -> HO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_o2, 1.0, i_ho2, 0.0, i_dummy) 

!===========================================================
!      c012 : O + H2O2 -> OH + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_h2o2, 1.0, i_oh, 1.0, i_ho2) 

!===========================================================
!      c013 : OH + OH -> H2O + O
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_oh, 1.0, i_h2o, 1.0, i_o) 

!===========================================================
!      c014 : OH + O3 -> HO2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_o3, 1.0, i_ho2, 1.0, i_o2) 

!===========================================================
!      c015 : HO2 + O3 -> OH + O2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ho2, 1.0, i_o3, 1.0, i_oh, 2.0, i_o2) 

!===========================================================
!      c016 : HO2 + HO2 + CO2 -> H2O2 + O2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_ho2, 1.0, i_h2o2, 1.0, i_o2) 

!===========================================================
!      c017 : OH + OH + CO2 -> H2O2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_oh, 1.0, i_h2o2, 0.0, i_dummy) 

!===========================================================
!      c018 : H + H + CO2 -> H2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_h, 1.0, i_h2, 0.0, i_dummy) 


!===========================================================
!      d001 : NO2 + O -> NO + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no2, 1.0, i_o, 1.0, i_no, 1.0, i_o2) 

!===========================================================
!      d002 : NO + O3 -> NO2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no, 1.0, i_o3, 1.0, i_no2, 1.0, i_o2) 

!===========================================================
!      d003 : NO + HO2 -> NO2 + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no, 1.0, i_ho2, 1.0, i_no2, 1.0, i_oh) 

!===========================================================
!      d004 : N + NO -> N2 + O 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_no, 1.0, i_n2, 1.0, i_o) 

!===========================================================
!      d005 : N + O2 -> NO + O 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_o2, 1.0, i_no, 1.0, i_o) 

!===========================================================
!      d006 : NO2 + H -> NO + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no2, 1.0, i_h, 1.0, i_no, 1.0, i_oh) 

!===========================================================
!      d007 : N + O -> NO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_o, 1.0, i_no, 0.0, i_dummy) 

!===========================================================
!      d008 : N + HO2 -> NO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_ho2, 1.0, i_no, 1.0, i_oh) 

!===========================================================
!      d009 : N + OH -> NO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_oh, 1.0, i_no, 1.0, i_h) 

!===========================================================
!      d010 : N(2D) + O -> N + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_n2d, 1.0, i_n, 0.0, i_dummy) 

!===========================================================
!      d011 : N(2D) + N2 -> N + N2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_n2d, 1.0, i_n, 0.0, i_dummy) 

!===========================================================
!      d012 : N(2D) + CO2 -> NO + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n2d, 1.0, i_co2, 1.0, i_no, 1.0, i_co) 

!===========================================================
!      e001 : CO + OH -> CO2 + H 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_oh, 1.0, i_co2, 1.0, i_h) 

!===========================================================
!      e002 : CO + O + M -> CO2 + M 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_o, 1.0, i_co2, 0.0, i_dummy) 
if(deutchem) then
!===========================================================
!      f001 : OD + OH -> HDO + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_oh, 1.0, i_hdo, 1.0, i_o) 

!===========================================================
!      f002 : OD + H2 -> HDO + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_h2, 1.0, i_hdo, 1.0, i_h)

!===========================================================
!      f003 : OD + HO2 -> HDO + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_ho2, 1.0, i_hdo, 1.0, i_o2)

!===========================================================
!      f004 : OD + H2O2 -> HDO + HO2
!===========================================================
   
   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_h2o2, 1.0, i_hdo, 1.0, i_ho2)

!===========================================================
!      f005 : O + OD -> O2 + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_od, 1.0, i_o2, 1.0, i_d)

!===========================================================
!      f006 : OD + H2 -> H2O + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h2, 1.0, i_od, 1.0, i_h2o, 1.0, i_d)

!===========================================================
!      f007 : OD + H -> OH + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_od, 1.0, i_oh, 1.0, i_d)

!===========================================================
!      f008 : CO + OD -> CO2 + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_od, 1.0, i_co2, 1.0, i_d)

!===========================================================
!      f009 : O3 + D -> O2 + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o3, 1.0, i_d, 1.0, i_o2, 1.0, i_od)

!===========================================================
!      f010 : HO2 + D -> OH + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_ho2, 1.0, i_d, 1.0, i_oh, 1.0, i_od)

!===========================================================
!      f011 : HO2 + D -> HDO + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_ho2, 1.0, i_d, 1.0, i_hdo, 1.0, i_o)

!===========================================================
!      f012 : OH + D -> H + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_d, 1.0, i_h, 1.0, i_od)

!===========================================================
!      f013 : H + D + CO2 -> HD + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_d, 1.0, i_hd, 0.0, i_dummy)

!===========================================================
!      f014 : D + HO2 -> HD + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_d, 1.0, i_ho2, 1.0, i_hd, 1.0, i_o2)


!===========================================================
!      f015 : OH + HD -> HDO + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hd, 1.0, i_hdo, 1.0, i_h)


!===========================================================
!      f016 : OH + HD -> H2O + D 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hd, 1.0, i_h2o, 1.0, i_d)


!===========================================================
!      f017 : D + O2 + CO2 -> DO2 + CO2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_d, 1.0, i_o2, 1.0, i_do2, 0.0, i_dummy)


!===========================================================
!      f018 : OD + O3 -> DO2 + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_o3, 1.0, i_do2, 1.0, i_o2)


!===========================================================
!      f019 : D + HO2 -> DO2 + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_d, 1.0, i_ho2, 1.0, i_do2, 1.0, i_h)


!===========================================================
!      f020 : O + DO2 -> OD + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_do2, 1.0, i_od, 1.0, i_o2)


!===========================================================
!      f021 : H + DO2 -> OH + OD
!===========================================================
   
   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_od, 1.0, i_oh)


!===========================================================
!      f022 : H + DO2 -> HD + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_hd, 1.0, i_o2)


!===========================================================
!      f023 : H + DO2 -> HDO + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_hdo, 1.0, i_o)


!===========================================================
!      f024 : H + DO2 -> HO2 + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_ho2, 1.0, i_d)


!===========================================================
!      f025 : OH + DO2 -> HDO + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_do2, 1.0, i_hdo, 1.0, i_o2)


!===========================================================
!      f026 : DO2 + O3 -> OD + O2 + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_do2, 1.0, i_o3, 1.0, i_od, 2.0, i_o2)


!===========================================================
!      f027 : OD + OH + CO2 -> HDO2 + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_oh, 1.0, i_hdo2, 0.0, i_dummy)


!===========================================================
!      f028 : DO2 + HO2 -> HDO2 + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_do2, 1.0, i_ho2, 1.0, i_hdo2, 1.0, i_o2)

!===========================================================
!      f029 : O + HDO2 -> OD + HO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_hdo2, 1.0, i_od, 1.0, i_ho2)


!===========================================================
!      f030 : O + HDO2 -> OH + DO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_hdo2, 1.0, i_oh, 1.0, i_do2)


!===========================================================
!      f031 : OH + HDO2 -> HDO + HO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hdo2, 1.0, i_hdo, 1.0, i_ho2)


!===========================================================
!      f032 : OH + HDO2 -> H2O + DO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hdo2, 1.0, i_h2o, 1.0, i_do2)


endif !deutchem

!Only if ion chemistry
if (ionchem) then

!===========================================================
!      i001 : CO2+ + O2 -> O2+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_o2, 1.0, i_o2plus, 1.0, i_co2)

!===========================================================
!      i002 : CO2+ + O -> O+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_o, 1.0, i_oplus, 1.0, i_co2)

!===========================================================
!      i003 : CO2+ + O -> O2+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_co)

!===========================================================
!      i004 : O2+ + e- -> O + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_elec, 2.0, i_o, 0.0, i_dummy)

!===========================================================
!      i005 : O+ + CO2 -> O2+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_co2, 1.0, i_o2plus, 1.0, i_co)

!===========================================================
!      i006 : CO2+ + e -> CO + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_elec, 1.0, i_co, 1.0, i_o)

!===========================================================
!      i007 : CO2+ + NO -> NO+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_no, 1.0, i_noplus, 1.0, i_co2)

!===========================================================
!      i008 : O2+ + NO -> NO+ + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_no, 1.0, i_noplus, 1.0, i_o2)

!===========================================================
!      i009 : O2+ + N2 -> NO+ + NO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_n2, 1.0, i_noplus, 1.0, i_no)

!===========================================================
!      i010 : O2+ + N -> NO+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_n, 1.0, i_noplus, 1.0, i_o)

!===========================================================
!      i011 : O+ + N2 -> NO+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_n2, 1.0, i_noplus, 1.0, i_n)

!===========================================================
!      i012 : NO+ + e -> N + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_noplus, 1.0, i_elec, 1.0, i_n, 1.0, i_o)

!===========================================================
!      i013 : CO+ + CO2 -> CO2+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_co2, 1.0, i_co2plus, 1.0, i_co)

!===========================================================
!      i014 : CO+ + O -> O+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_o, 1.0, i_oplus, 1.0, i_co)

!===========================================================
!      i015 : C+ + CO2 -> CO+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_cplus, 1.0, i_co2, 1.0, i_coplus, 1.0, i_co)

!===========================================================
!      i016 : N2+ + CO2 -> CO2+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_co2, 1.0, i_co2plus, 1.0, i_n2)

!===========================================================
!      i017 : N2+ + O -> NO+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_o, 1.0, i_noplus, 1.0, i_n)

!===========================================================
!      i018 : N2+ + CO -> CO+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_co, 1.0, i_coplus, 1.0, i_n2)

!===========================================================
!      i019 : N2+ + e -> N + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_elec, 2.0, i_n, 0.0, i_dummy)

!===========================================================
!      i020 : N2+ + O -> O+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_o, 1.0, i_oplus, 1.0, i_n2)

!===========================================================
!      i021 : N+ + CO2 -> CO2+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_nplus, 1.0, i_co2, 1.0, i_co2plus, 1.0, i_n)

!===========================================================
!      i022 : CO+ + H -> H+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_h, 1.0, i_hplus, 1.0, i_co)

!===========================================================
!      i023 : O+ + H -> H+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_h, 1.0, i_hplus, 1.0, i_o)

!===========================================================
!      i024 : H+ + O -> O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hplus, 1.0, i_o, 1.0, i_oplus, 1.0, i_h)

!===========================================================
!      i025 : CO2+ + H2 -> HCO2+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_h2, 1.0, i_hco2plus, 1.0, i_h)

!===========================================================
!      i026 : HCO2+ + e -> H + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_elec, 1.0, i_h, 1.0, i_co2)

!===========================================================
!      i027 : HCO2+ + e -> H + O + CO 
!===========================================================
!We divide this reaction in two

!0.5HCO2+ + 0.5e -> H

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(.5, i_hco2plus, 0.5, i_elec, 1.0, i_h, 0.0, i_dummy)

!0.5 HCO2+ + 0.5 e -> O + CO

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(0.5, i_hco2plus, 0.5, i_elec, 1.0, i_o, 1.0, i_co)

!===========================================================
!      i029 : HCO2+ + e -> OH + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_elec, 1.0, i_oh, 1.0, i_co)


!===========================================================
!      i030 : HCO2+ + e -> H + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_elec, 1.0, i_h, 1.0, i_co2)


!===========================================================
!      i031 : HCO2+ + O -> HCO+ + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_o, 1.0, i_hcoplus, 1.0, i_o2)


!===========================================================
!      i032 : HCO2+ + CO -> HCO+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_co, 1.0, i_hcoplus, 1.0, i_co2)


!===========================================================
!      i033 : H+ + CO2 -> HCO+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hplus, 1.0, i_co2, 1.0, i_hcoplus, 1.0, i_o)


!===========================================================
!      i034 : CO2+ + H -> HCO+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_h, 1.0, i_hcoplus, 1.0, i_o)


!===========================================================
!      i035 : CO+ + H2 -> HCO+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_h2, 1.0, i_hcoplus, 1.0, i_h)


!===========================================================
!      i036 : HCO+ + e- -> CO + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hcoplus, 1.0, i_elec, 1.0, i_co, 1.0, i_h)

!===========================================================
!      i037 : CO2+ + H2O -> H2O+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_co2)

!===========================================================
!      i038 : CO+ + H2O -> H2O+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_co)

!===========================================================
!      i039 : O+ + H2O -> H2O+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_o)

!===========================================================
!      i040 : N2+ + H2O -> H2O+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_n2)

!===========================================================
!      i041 : N+ + H2O -> H2O+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_nplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_n)

!===========================================================
!      i042 : H+ + H2O -> H2O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_h)

!===========================================================
!      i043 : H2O+ + O2 -> O2+ + H2O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_o2, 1.0, i_o2plus, 1.0, i_h2o)

!===========================================================
!      i044 : H2O+ + CO -> HCO+ + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_co, 1.0, i_hcoplus, 1.0, i_oh) 

!===========================================================
!      i045 : H2O+ + O -> O2+ + H2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_h2)

!===========================================================
!      i046 : H2O+ + NO -> NO+ + H2O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_no, 1.0, i_noplus, 1.0, i_h2o)

!===========================================================
!      i047 : H2O+ + e- -> H + H + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_elec, 2.0, i_h, 1.0, i_o)

!===========================================================
!      i048 : H2O+ + e- -> H + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_elec, 1.0, i_h, 1.0, i_oh)

!===========================================================
!      i049 : H2O+ + e- -> H2 + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_elec, 1.0, i_h2, 1.0, i_o)

!===========================================================
!      i050 : H2O+ + H2O -> H3O+ + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_h2o, 1.0, i_h3oplus, 1.0, i_oh)

!===========================================================
!      i051 : H2O+ + H2 -> H3O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_h2, 1.0, i_h3oplus, 1.0, i_h)

!===========================================================
!      i052 : HCO+ + H2O -> H3O+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hcoplus, 1.0, i_h2o, 1.0, i_h3oplus, 1.0, i_co)

!===========================================================
!      i053: H3O+ + e -> OH + H + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h3oplus, 1.0, i_elec, 1.0, i_oh, 2.0, i_h)

!===========================================================
!      i054: H3O+ + e -> H2O + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h3oplus, 1.0, i_elec, 1.0, i_h2o, 1.0, i_h)

!===========================================================
!      i055: H3O+ + e -> HO + H2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h3oplus, 1.0, i_elec, 1.0, i_oh, 1.0, i_h2)

!===========================================================
!      i056: H3O+ + e -> O + H2 + H
!===========================================================
!We divide this reaction in two

!0.5H3O+ + 0.5e -> O

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(0.5, i_h3oplus, 0.5, i_elec, 1.0, i_o, 0.0, i_dummy)

!0.5H3O+ + 0.5e -> H2 + H

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(0.5, i_h3oplus, 0.5, i_elec, 1.0, i_h2, 1.0, i_h)

!===========================================================
!      i057: O+ + H2 -> OH+ + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_h2, 1.0, i_ohplus, 1.0, i_h)

!===========================================================
!      i058: OH+ + O -> O2+ + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_h)

!===========================================================
!      i059: OH+ + CO2 -> HCO2+ + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_co2, 1.0, i_hco2plus, 1.0, i_o)

!===========================================================
!      i060: OH+ + CO -> HCO+ + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_co, 1.0, i_hcoplus, 1.0, i_o)

!===========================================================
!      i061: OH+ + NO -> NO+ + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_no, 1.0, i_noplus, 1.0, i_oh)

!===========================================================
!      i062: OH+ + H2 -> H2O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_h2, 1.0, i_h2oplus, 1.0, i_h)

!===========================================================
!      i063: OH+ + O2 -> O2+ + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_o2, 1.0, i_o2plus, 1.0, i_oh)

end if    !ionchem

!===========================================================
!      h001: HO2 + ice -> products
!            treated as
!            HO2 -> 0.5 H2O + 0.75 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 0.5, i_h2o, 0.75, i_o2) 

!===========================================================
!      h002: OH + ice -> products
!            treated as
!            OH -> 0.5 H2O + 0.25 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_oh, 0.5, i_h2o, 0.25, i_o2) 

!===========================================================
!      h003: H2O2 + ice -> products
!            treated as
!            H2O2 -> H2O + 0.5 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 1.0, i_h2o, 0.5, i_o2) 

!===========================================================
!      h004: HO2 + dust -> products
!            treated as
!            HO2 -> 0.5 H2O + 0.75 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 0.5, i_h2o, 0.75, i_o2) 

!===========================================================
!      h005: H2O2 + dust -> products
!            treated as
!            H2O2 -> H2O + 0.5 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 1.0, i_h2o, 0.5, i_o2) 

!===========================================================
!  check dimensions 
!===========================================================

print*, 'nb_phot       = ', nb_phot
print*, 'nb_reaction_4 = ', nb_reaction_4
print*, 'nb_reaction_3 = ', nb_reaction_3

if ((nb_phot /= nb_phot_max)             .or.  &
    (nb_reaction_3 /= nb_reaction_3_max) .or.  &
    (nb_reaction_4 /= nb_reaction_4_max)) then
   print*, 'wrong dimensions in indice' 
   call abort_physic("indice","wrong array dimensions",1)
end if  

end subroutine indice

!*****************************************************************

      subroutine gcmtochim(nlayer, ionchem, deutchem, nq, zycol,     &
                           lswitch, nesp,                            &
                           i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h, &
                           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,         &
                           i_n, i_n2d, i_no, i_no2, i_n2,            &
                           i_co2plus, i_oplus, i_o2plus, i_noplus,   &
                           i_coplus, i_cplus, i_n2plus, i_nplus,     &
                           i_hplus, i_hco2plus, i_hcoplus, i_h2oplus,&
                           i_h3oplus, i_ohplus, i_elec, i_hdo, i_od, &
                           i_d, i_hd, i_do2, i_hdo2, dens, rm, c) 
        
!*****************************************************************

      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,         &
     &                      igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,  &
     &                      igcm_ho2, igcm_h2o2, igcm_h2o_vap,           &
     &                      igcm_n, igcm_n2d, igcm_no, igcm_no2, igcm_n2,&
     &                      igcm_co2plus, igcm_oplus, igcm_o2plus,       &
     &                      igcm_noplus, igcm_coplus, igcm_cplus,        &
     &                      igcm_n2plus, igcm_nplus, igcm_hplus,         &
     &                      igcm_hco2plus, igcm_hcoplus, igcm_h2oplus,   &
     &                      igcm_h3oplus, igcm_ohplus, igcm_elec,        &
     &                      igcm_hdo_vap, igcm_od, igcm_d, igcm_hd,      &
     &                      igcm_do2, igcm_hdo2

      implicit none

      include "callkeys.h"

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     input:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      integer, intent(in) :: nlayer ! number of atmospheric layers
      integer, intent(in) :: nq     ! number of tracers in the gcm
      logical, intent(in) :: ionchem
      logical, intent(in) :: deutchem
      integer :: nesp               ! number of species in the chemistry
      integer :: lswitch            ! interface level between chemistries

      integer :: i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,           &
                 i_h2, i_oh, i_ho2, i_h2o2, i_h2o,                   &
                 i_n, i_n2d, i_no, i_no2, i_n2,                      &
                 i_co2plus, i_oplus, i_o2plus, i_noplus, i_coplus,   &
                 i_cplus, i_n2plus, i_nplus, i_hplus, i_hco2plus,    &
                 i_hcoplus, i_h2oplus, i_h3oplus, i_ohplus, i_elec,  &
                 i_hdo, i_od, i_d, i_hd, i_do2, i_hdo2

      real :: zycol(nlayer,nq)      ! volume mixing ratios in the gcm
      real :: dens(nlayer)          ! total number density (molecule.cm-3) 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     output:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real, dimension(nlayer,nesp)            :: rm ! volume mixing ratios 
      real (kind = 8), dimension(nlayer,nesp) :: c  ! number densities (molecule.cm-3) 
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     local:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer      :: l, iesp
      logical,save :: firstcall = .true.
      
!     first call initializations

      if (firstcall) then

!       identify the indexes of the tracers we need

         if (igcm_co2 == 0) then
            write(*,*) "gcmtochim: Error; no CO2 tracer !!!"
            call abort_physic("gcmtochim","missing co2 tracer",1)
         endif
         if (igcm_co == 0) then
            write(*,*) "gcmtochim: Error; no CO tracer !!!"
            call abort_physic("gcmtochim","missing co tracer",1)
         end if
         if (igcm_o == 0) then
            write(*,*) "gcmtochim: Error; no O tracer !!!"
            call abort_physic("gcmtochim","missing o tracer",1)
         end if
         if (igcm_o1d == 0) then
            write(*,*) "gcmtochim: Error; no O1D tracer !!!"
            call abort_physic("gcmtochim","missing o1d tracer",1)
         end if
         if (igcm_o2 == 0) then
            write(*,*) "gcmtochim: Error; no O2 tracer !!!"
            call abort_physic("gcmtochim","missing o2 tracer",1)
         end if
         if (igcm_o3 == 0) then
            write(*,*) "gcmtochim: Error; no O3 tracer !!!"
            call abort_physic("gcmtochim","missing o3 tracer",1)
         end if
         if (igcm_h == 0) then
            write(*,*) "gcmtochim: Error; no H tracer !!!"
            call abort_physic("gcmtochim","missing h tracer",1)
         end if
         if (igcm_h2 == 0) then
            write(*,*) "gcmtochim: Error; no H2 tracer !!!"
            call abort_physic("gcmtochim","missing h2 tracer",1)
         end if
         if (igcm_oh == 0) then
            write(*,*) "gcmtochim: Error; no OH tracer !!!"
            call abort_physic("gcmtochim","missing oh tracer",1)
         end if
         if (igcm_ho2 == 0) then
            write(*,*) "gcmtochim: Error; no HO2 tracer !!!"
            call abort_physic("gcmtochim","missing ho2 tracer",1)
         end if
         if (igcm_h2o2 == 0) then
            write(*,*) "gcmtochim: Error; no H2O2 tracer !!!"
            call abort_physic("gcmtochim","missing h2o2 tracer",1)
         end if
         if (igcm_n == 0) then
            write(*,*) "gcmtochim: Error; no N tracer !!!"
            call abort_physic("gcmtochim","missing n tracer",1)
         end if
         if (igcm_n2d == 0) then
            write(*,*) "gcmtochim: Error; no N2D tracer !!!"
            call abort_physic("gcmtochim","missing n2d tracer",1)
         end if
         if (igcm_no == 0) then
            write(*,*) "gcmtochim: Error; no NO tracer !!!"
            call abort_physic("gcmtochim","missing no tracer",1)
         end if
         if (igcm_no2 == 0) then
            write(*,*) "gcmtochim: Error; no NO2 tracer !!!"
            call abort_physic("gcmtochim","missing no2 tracer",1)
         end if
         if (igcm_n2 == 0) then
            write(*,*) "gcmtochim: Error; no N2 tracer !!!"
            call abort_physic("gcmtochim","missing n2 tracer",1)
         end if
         if (igcm_h2o_vap == 0) then
            write(*,*) "gcmtochim: Error; no water vapor tracer !!!"
            call abort_physic("gcmtochim","missing h2o_vap tracer",1)
         end if
         if(deutchem) then
            if (igcm_hdo_vap == 0) then
               write(*,*) "gcmtochim: Error; no HDO tracer !!!"
               call abort_physic("gcmtochim","missing hdo_vap tracer",1)
            end if
            if (igcm_od == 0) then
               write(*,*) "gcmtochim: Error; no OD tracer !!!"
               call abort_physic("gcmtochim","missing od tracer",1)
            end if
            if (igcm_d == 0) then
               write(*,*) "gcmtochim: Error; no D tracer !!!"
               call abort_physic("gcmtochim","missing d tracer",1)
            end if
            if (igcm_hd == 0) then
               write(*,*) "gcmtochim: Error; no HD tracer !!!"
               call abort_physic("gcmtochim","missing hd tracer",1)
            end if
            if (igcm_do2 == 0) then
               write(*,*) "gcmtochim: Error; no DO2 tracer !!!"
               call abort_physic("gcmtochim","missing do2 tracer",1)
            end if
            if (igcm_hdo2 == 0) then
               write(*,*) "gcmtochim: Error; no HDO2 tracer !!!"
               call abort_physic("gcmtochim","missing hdo2 tracer",1)
            end if
         endif !deutchem
         if (ionchem) then
            if (igcm_co2plus == 0) then
               write(*,*) "gcmtochim: Error; no CO2+ tracer !!!"
               call abort_physic("gcmtochim","missing co2plus tracer",1)
            end if
            if (igcm_oplus == 0) then
               write(*,*) "gcmtochim: Error; no O+ tracer !!!"
               call abort_physic("gcmtochim","missing oplus tracer",1)
            end if
            if (igcm_o2plus == 0) then
               write(*,*) "gcmtochim: Error; no O2+ tracer !!!"
               call abort_physic("gcmtochim","missing o2plus tracer",1)
            end if
            if (igcm_noplus == 0) then
               write(*,*) "gcmtochim: Error; no NO+ tracer !!!"
               call abort_physic("gcmtochim","missing noplus tracer",1)
            endif
            if (igcm_coplus == 0) then
               write(*,*) "gcmtochim: Error; no CO+ tracer !!!"
               call abort_physic("gcmtochim","missing coplus tracer",1)
            endif
            if (igcm_cplus == 0) then
               write(*,*) "gcmtochim: Error; no C+ tracer !!!"
               call abort_physic("gcmtochim","missing cplus tracer",1)
            endif
            if (igcm_n2plus == 0) then
               write(*,*) "gcmtochim: Error; no N2+ tracer !!!"
               call abort_physic("gcmtochim","missing n2plus tracer",1)
            endif
            if (igcm_nplus == 0) then
               write(*,*) "gcmtochim: Error; no N+ tracer !!!"
               call abort_physic("gcmtochim","missing nplus tracer",1)
            endif
            if (igcm_hplus == 0) then
               write(*,*) "gcmtochim: Error; no H+ tracer !!!"
               call abort_physic("gcmtochim","missing hplus tracer",1)
            endif
            if (igcm_hco2plus == 0) then
               write(*,*) "gcmtochim: Error; no HCO2+ tracer !!!"
               call abort_physic("gcmtochim","missing hco2plus tracer",1)
            endif
            if (igcm_hcoplus == 0) then
               write(*,*) "gcmtochim: Error; no HCO+ tracer !!!"
               call abort_physic("gcmtochim","missing hcoplus tracer",1)
            endif
            if (igcm_h2oplus == 0) then
               write(*,*) "gcmtochim: Error; no H2O+ tracer !!!"
               call abort_physic("gcmtochim","missing h2oplus tracer",1)
            endif
            if (igcm_h3oplus == 0) then
               write(*,*) "gcmtochim: Error; no H3O+ tracer !!!"
               call abort_physic("gcmtochim","missing h3oplus tracer",1)
            endif
            if (igcm_ohplus == 0) then
               write(*,*) "gcmtochim: Error; no OH+ tracer !!!"
               call abort_physic("gcmtochim","missing ohplus tracer",1)
            endif
            if (igcm_elec == 0) then
               write(*,*) "gcmtochim: Error; no e- tracer !!!"
               call abort_physic("gcmtochim","missing elec tracer",1)
            end if
         end if  ! ionchem
         firstcall = .false.
      end if ! of if (firstcall)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     initialise mixing ratios
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l = 1,nlayer
         rm(l,i_co2)  = zycol(l, igcm_co2)
         rm(l,i_co)   = zycol(l, igcm_co)
         rm(l,i_o)    = zycol(l, igcm_o)
         rm(l,i_o1d)  = zycol(l, igcm_o1d)
         rm(l,i_o2)   = zycol(l, igcm_o2)
         rm(l,i_o3)   = zycol(l, igcm_o3)
         rm(l,i_h)    = zycol(l, igcm_h)
         rm(l,i_h2)   = zycol(l, igcm_h2)
         rm(l,i_oh)   = zycol(l, igcm_oh)
         rm(l,i_ho2)  = zycol(l, igcm_ho2)
         rm(l,i_h2o2) = zycol(l, igcm_h2o2)
         rm(l,i_h2o)  = zycol(l, igcm_h2o_vap)
         rm(l,i_n)    = zycol(l, igcm_n)
         rm(l,i_n2d)  = zycol(l, igcm_n2d)
         rm(l,i_no)   = zycol(l, igcm_no)
         rm(l,i_no2)  = zycol(l, igcm_no2)
         rm(l,i_n2)   = zycol(l, igcm_n2)
      enddo

      if (deutchem) then
         do l = 1,nlayer
            rm(l,i_hdo)  = zycol(l, igcm_hdo_vap)
            rm(l,i_od)   = zycol(l, igcm_od)
            rm(l,i_d)    = zycol(l, igcm_d)
            rm(l,i_hd)   = zycol(l, igcm_hd)
            rm(l,i_do2)  = zycol(l, igcm_do2)
            rm(l,i_hdo2)  = zycol(l, igcm_hdo2)
         end do
      endif

      if (ionchem) then
         do l = 1,nlayer
            rm(l,i_co2plus)  = zycol(l, igcm_co2plus)
            rm(l,i_oplus)    = zycol(l, igcm_oplus)
            rm(l,i_o2plus)   = zycol(l, igcm_o2plus)
            rm(l,i_noplus)   = zycol(l, igcm_noplus)
            rm(l,i_coplus)   = zycol(l, igcm_coplus)
            rm(l,i_cplus)    = zycol(l, igcm_cplus)
            rm(l,i_n2plus)   = zycol(l, igcm_n2plus)
            rm(l,i_nplus)    = zycol(l, igcm_nplus)
            rm(l,i_hplus)    = zycol(l, igcm_hplus)
            rm(l,i_hco2plus) = zycol(l, igcm_hco2plus)
            rm(l,i_hcoplus)  = zycol(l, igcm_hcoplus)
            rm(l,i_h2oplus)  = zycol(l, igcm_h2oplus)
            rm(l,i_h3oplus)  = zycol(l, igcm_h3oplus)
            rm(l,i_ohplus)   = zycol(l, igcm_ohplus)
            rm(l,i_elec)     = zycol(l, igcm_elec)
         end do 
      end if

      where (rm(:,:) < 1.e-30)
         rm(:,:) = 0.
      end where

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     initialise number densities
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iesp = 1,nesp
         do l = 1,nlayer
            c(l,iesp) = rm(l,iesp)*dens(l)
         end do
      end do

      end subroutine gcmtochim

!*****************************************************************
 
      subroutine chimtogcm(nlayer, ionchem, deutchem, nq, zycol,      &
                           lswitch, nesp,                             &
                           i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,  &
                           i_h2, i_oh, i_ho2, i_h2o2, i_h2o,          &
                           i_n, i_n2d, i_no, i_no2, i_n2,             &
                           i_co2plus, i_oplus, i_o2plus, i_noplus,    &
                           i_coplus, i_cplus, i_n2plus, i_nplus,      &
                           i_hplus, i_hco2plus, i_hcoplus, i_h2oplus, &
                           i_h3oplus, i_ohplus, i_elec, i_hdo, i_od,  &
                           i_d, i_hd, i_do2, i_hdo2, dens, c) 
 
!*****************************************************************
 
      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,          &
                            igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,   &
                            igcm_ho2, igcm_h2o2, igcm_h2o_vap,            &
                            igcm_n, igcm_n2d, igcm_no, igcm_no2, igcm_n2, &
                            igcm_co2plus, igcm_oplus, igcm_o2plus,        &
                            igcm_noplus, igcm_coplus, igcm_cplus,         &
                            igcm_n2plus, igcm_nplus, igcm_hplus,          &
                            igcm_hco2plus, igcm_hcoplus, igcm_h2oplus,    &
                            igcm_h3oplus, igcm_ohplus, igcm_elec,         &
                            igcm_hdo_vap, igcm_od, igcm_d, igcm_hd,       &
                            igcm_do2, igcm_hdo2

      implicit none

#include "callkeys.h"

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     inputs:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      integer, intent(in) :: nlayer  ! number of atmospheric layers
      integer, intent(in) :: nq      ! number of tracers in the gcm
      logical, intent(in) :: ionchem
      logical, intent(in) :: deutchem
      integer :: nesp                ! number of species in the chemistry
      integer :: lswitch             ! interface level between chemistries
      integer :: i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h,       &
                 i_h2, i_oh, i_ho2, i_h2o2, i_h2o,               &
                 i_n, i_n2d, i_no, i_no2, i_n2,                  &
                 i_co2plus, i_oplus, i_o2plus, i_noplus,         &
                 i_coplus, i_cplus, i_n2plus, i_nplus,           &
                 i_hplus, i_hco2plus, i_hcoplus, i_h2oplus,      &
                 i_h3oplus, i_ohplus, i_elec, i_hdo, i_od, i_d,  &
                 i_hd, i_do2, i_hdo2

      real :: dens(nlayer)     ! total number density (molecule.cm-3) 
      real (kind = 8), dimension(nlayer,nesp) :: c  ! number densities (molecule.cm-3) 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     output:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
      real zycol(nlayer,nq)  ! volume mixing ratios in the gcm

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     local:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer l
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     save mixing ratios for the gcm
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l = 1,lswitch-1
         zycol(l, igcm_co2)     = c(l,i_co2)/dens(l) 
         zycol(l, igcm_co)      = c(l,i_co)/dens(l) 
         zycol(l, igcm_o)       = c(l,i_o)/dens(l) 
         zycol(l, igcm_o1d)     = c(l,i_o1d)/dens(l)
         zycol(l, igcm_o2)      = c(l,i_o2)/dens(l) 
         zycol(l, igcm_o3)      = c(l,i_o3)/dens(l) 
         zycol(l, igcm_h)       = c(l,i_h)/dens(l)  
         zycol(l, igcm_h2)      = c(l,i_h2)/dens(l) 
         zycol(l, igcm_oh)      = c(l,i_oh)/dens(l) 
         zycol(l, igcm_ho2)     = c(l,i_ho2)/dens(l) 
         zycol(l, igcm_h2o2)    = c(l,i_h2o2)/dens(l)
         zycol(l, igcm_h2o_vap) = c(l,i_h2o)/dens(l)
         zycol(l, igcm_n)       = c(l,i_n)/dens(l)
         zycol(l, igcm_n2d)     = c(l,i_n2d)/dens(l)
         zycol(l, igcm_no)      = c(l,i_no)/dens(l)
         zycol(l, igcm_no2)     = c(l,i_no2)/dens(l)
         zycol(l, igcm_n2)      = c(l,i_n2)/dens(l)
      enddo
      
      if(deutchem) then
         do l=1,lswitch-1
            zycol(l, igcm_hdo_vap) = c(l,i_hdo)/dens(l)
            zycol(l, igcm_od)      = c(l,i_od)/dens(l)
            zycol(l, igcm_d)       = c(l,i_d)/dens(l)
            zycol(l, igcm_hd)      = c(l,i_hd)/dens(l)
            zycol(l, igcm_do2)     = c(l,i_do2)/dens(l)
            zycol(l, igcm_hdo2)    = c(l,i_hdo2)/dens(l)
         end do
      endif !deutchem

      if (ionchem) then
         do l = 1,lswitch-1
            zycol(l, igcm_co2plus) = c(l,i_co2plus)/dens(l)
            zycol(l, igcm_oplus)   = c(l,i_oplus)/dens(l)
            zycol(l, igcm_o2plus)  = c(l,i_o2plus)/dens(l)
            zycol(l, igcm_noplus)  = c(l,i_noplus)/dens(l)
            zycol(l, igcm_coplus)  = c(l,i_coplus)/dens(l)
            zycol(l, igcm_cplus)   = c(l,i_cplus)/dens(l)
            zycol(l, igcm_n2plus)  = c(l,i_n2plus)/dens(l)
            zycol(l, igcm_nplus)   = c(l,i_nplus)/dens(l)
            zycol(l, igcm_hplus)   = c(l,i_hplus)/dens(l)
            zycol(l, igcm_hco2plus)= c(l,i_hco2plus)/dens(l)
            zycol(l, igcm_hcoplus) = c(l,i_hcoplus)/dens(l)
            zycol(l, igcm_h2oplus) = c(l,i_h2oplus)/dens(l)
            zycol(l, igcm_h3oplus) = c(l,i_h3oplus)/dens(l)
            zycol(l, igcm_ohplus)  = c(l,i_ohplus)/dens(l)
            zycol(l, igcm_elec)    = c(l,i_elec)/dens(l)
         end do
      end if 

      end subroutine chimtogcm

end subroutine photochemistry
