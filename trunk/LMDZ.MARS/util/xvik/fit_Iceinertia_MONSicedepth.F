       PROGRAM SUPERFIT

c   ------------------------------------------------------------------
c      #### Special version with ice depth and ice inertia #####
c      Program used to compute best-fit cap albedo, ICE DEPTH  and total CO2
c      inventory 
c      Make use of 2 annual GCM  run  performed with 2 different cap albedos. 
c      Based on Hourdin et al. (JGR, 1995) 
c      Modified to account for non-linear sensitivity of polar cap mass
c      To cap albedo  (FF,2000)
c      Modified/cleaned-up to use 4 input files and minimise cap albedos,
c      MONS-derived ice depth and total pressure (EM,2008)
c      Modified to minimise wrt ice thermal inertia coefficients (EM, 2009)
c   ------------------------------------------------------------------
       
       implicit none

       integer ngcm,nvl1,nsol,time_unit

c ***** THE FOLLOWING LINES MUST BE ADAPTED FOR EACH CASES : **************
!       parameter(ngcm=10704)  ! size of the 1 year raw GCM data
!       parameter(ngcm=2672) ! "*4 & *6" files
       parameter(ngcm=8028) ! "*3 & *5" files
       parameter(nvl1=669)  ! size of the Vl1 Ps Observation
c *************************************************************


c      size the smoothed data (1/sol) which are used for the simulation 
       parameter(nsol=669)  ! size the smoothed data (1/sol)
c      parameter(nsol=445)  ! size the smoothed data (1/sol)

c      OBSERVED VL1 pressure
c      ~~~~~~~~~~~~~~~~~~~~~
       real pvl_obs(nvl1), solvl(nvl1), lsvl(nvl1)

c      Raw Data from 4 simulations with 2 different albedo/iced :
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c runs are numbered 1-4 such than runs 1&2 are at same ice depth
c                             and runs 3&4 are at same ice depth
c                             and runs 1&3 are at same ice inertia
c                             and runs 2&4 are at same ice inertia
       real iceigcmN(4) , iceigcmS(4) 
       real icedgcm(4)
       integer refrun ! run number of the 'reference' run 
       ! reference values for functions (e.g. values of run #1)
       real iceigcm_refN, icedgcm_ref
       real iceigcm_refS

       real ptot(4)       
       real patm(ngcm,4), pcapn(ngcm,4),pcaps(ngcm,4)
       real pvl_gcm(ngcm,4)     ! Simulated VL1 pressure
       real solgcm(ngcm)        ! time in runs in sol
       real lsgcm(ngcm)         ! time in runs in ls
       integer year_xvik        ! year of xvik files to use for fit
       integer plot_test        ! =1 if output file is to plot COST(DN,DS), =2 for COST(IN,IS)s
       
c      Smoothed Data from the 2 simulations with 2 different albedo :
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       real box
       integer n ,i
       real sol(nsol)       ! time in runs
       real ls(nsol)
       real solconv
       real pvl_sm(nsol,4)  ! Simulated smooth VL1 pressure
       real patm_sm(nsol,4), pcapn_sm(nsol,4),pcaps_sm(nsol,4)
       real alphavl1(nsol) ! pvl/patm ratio
c      Mass cap  sensitivity to thermal inertia (derivative)     
       real dpdicein(nsol),dpdiceis(nsol) 
       real dpdicein_fltr(nsol),dpdiceis_fltr(nsol) 
c      Mass cap  sensitivity to ice inertia (derivative)
       real dpden(nsol),dpdes(nsol) 
       real dpden_fltr(nsol),dpdes_fltr(nsol) 

       real iceis, icein, ptry ,pvl1
       real iceds, icedn
! ice inertia range to explore, and maximum allowed difference between N&S ice inertias:
       real iceimin,iceimax,maxiceidiff
       real deltaicei ! step in ice thermal inertia
! ice depth range to explore, and maximum allowed difference between N&S ice depths:
       real icedmin,icedmax,maxiceddiff
       real deltaiced ! step in ice depth coefficient 
! total pressure range to explore, and pressure step:
       real pmin,pmax,deltap
       real cost,cost4plot
       real iceis4plot,ps4plot,icein4plot !,iceds4plot
       real iceds4plot,icedn4plot
       real costmin, pfit, iceinfit, iceisfit
       real pcapn_new, pcaps_new
       real icednfit, icedsfit
       logical fit 
       
       real fonc, fonc2n, fonc2s

c      Inputs files (from xvik) variables  :
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
       character(len=33) filename1, filename2
       character(len=99) dset_DminImin,dset_DmaxImin
       character(len=99) dset_DminImax,dset_DmaxImax 
       integer ie


c----------------------------------------------------------------------
c      Initialisation :
c      ~~~~~~~~~~~~~~
c      Implicit function : behavior of Pcap = fonc(alb)
c                          (e.g. fonc(alb) = alb  if linear)
!!       fonc(albn) = exp(3.46*albn) ! for albedo, both north & south
       fonc(icein) = icein
c       fonc2n(icedn) = log(icedn +10)
       fonc2n(icedn) = icedn ! for northern ice depth coefficient
c       fonc2s(iceds) = log(iceds +1)
       fonc2s(iceds) = iceds ! for southern ice depth coefficient

c ***** THE FOLLOWING LINES MUST BE ADAPTED FOR EACH CASES : **************
       
       write(*,*) 'Program written for xvik outputs files ''xpsol'' ',
     &'and ''xprestot'' from xvik.F'
       write(*,*) 'Enter extreme parameters values used in xvik.F'
       write(*,*) 'Minimum ice depth coefficient'
       read(*,*)   icedmin
       write(*,*)  icedmin
       write(*,*) 'Maximum ice depth coefficient'
       read(*,*)   icedmax
       write(*,*)  icedmax      
       write(*,*) 'Minimum ice thermal inertia '
       read(*,*)   iceimin
       write(*,*)  iceimin
       write(*,*) 'Maximum ice thermal inertia '
       read(*,*)   iceimax
       write(*,*)  iceimax
       
    
       maxiceidiff = 3000  ! maximum allowed difference between N and S best fit ice thermal inertia
       deltaicei= 20	 ! step in ice thermal inertia
       maxiceddiff = 30.e-4  ! maximum allowed difference between N and S best fit ice depth
       deltaiced= 0.5e-4 ! step in ice depth coefficient
       pmin = 690  ! minimum allowed  best fit CO2 total pressure (Pa)
       pmax = 720 ! maximum allowed  best fit CO2 total pressure (Pa)
       deltap=1.0 ! step in pressure

c      File to be read : 
c      ~~~~~~~~~~~~~~~
c runs are numbered 1-4 such than runs 1&2 are at same ice depth
c                             and runs 3&4 are at same ice depth
c                             and runs 1&3 are at same ice inertia
c                             and runs 2&4 are at same ice inertia
       ! First run:
       iceigcmN(1)= iceimin     ! Northern Cap ice thermal inertia
       iceigcmS(1)= iceimin     ! Southern Cap ice thermal inertia
       icedgcm(1) = icedmin     ! Ice depth coefficient
       ! Second run:
       iceigcmN(2) = iceimax    ! Northern Cap ice thermal inertia
       iceigcmS(2) = iceimax    ! Southern Cap ice thermal inertia
       icedgcm(2)  = icedmin    ! Ice depth coefficient
       ! Third run:
       iceigcmN(3) = iceimin    ! Northern Cap ice thermal inertia
       iceigcmS(3) = iceimin    ! Southern Cap ice thermal inertia
       icedgcm(3)  = icedmax    ! Ice depth coefficient
       ! Fourth run:
       iceigcmN(4) = iceimax    ! Northern Cap ice thermal inertia
       iceigcmS(4) = iceimax    ! Southern Cap ice thermal inertia
       icedgcm(4)  = icedmax    ! Ice depth coefficient
       
       ! set reference values: those of one of the files (here file 1)
       refrun=1
       ! ice inertia of the reference run used to simulate function
       iceigcm_refN=iceigcmN(refrun) 
       ! ice inertia of the reference run used to simulate function
       iceigcm_refS=iceigcmS(refrun) 
       ! icedepth coefficient of the reference run used to simulate function
       icedgcm_ref=icedgcm(refrun)
       
       write(*,*) 'Path to xvik outputs with minimum ice depth ', 
     &'coefficient and minimum ice thermal inertia ?'
       read (*,'(a)')  dset_DminImin
       write(*,*)  dset_DminImin
       write(*,*) 'Path to xvik outputs with minimum ice depth ', 
     &'coefficient and maximum ice thermal inertia ?'
       read (*,'(a)')  dset_DminImax
       write(*,*)  dset_DminImax
       write(*,*) 'Path to xvik outputs with maximum ice depth ', 
     &'coefficient and minimum ice thermal inertia ?'
       read (*,'(a)')  dset_DmaxImin
       write(*,*)  dset_DmaxImin
       write(*,*) 'Path to xvik outputs with maximum ice depth ', 
     &'coefficient and maximum ice thermal inertia ?'
       read (*,'(a)')  dset_DmaxImax
       write(*,*)  dset_DmaxImax
       
       write(*,*) 'Which year of xvik outputs do you want to use ?'
       read (*,*)  year_xvik
       write(*,*)  year_xvik
      
       write(*,*) 'Xvik files in sol (1), ls (2) or both (3)'
       read(*,*) time_unit       
       write(*,*) time_unit
       

       write(*,*) 'COST being the model/obs difference'
       write(*,*) 'IN ans IS being the ice thermal inertia of ', 
     &'Northern and Southern Cap'
       write(*,*) 'DN ans DS being the ice depth coefficient of ', 
     &'of Northern and Southern Cap'
       write(*,*) 'Do you want to use output file ',
     &'''minimization.txt'' to plot COST(DN,DS) (1) or COST(IN,IS) (2)'   
       read(*,*) plot_test
       write(*,*) plot_test

c      For four simulation, files :  "time (VL1 sol) , Ps Vl1 (Pa)"
c                                    "time (VL1 sol), Patm,PcapN,PcapS (Pa)"



       write(filename1,fmt='(a6,i1)') 'xpsol1',year_xvik
       write(filename2,fmt='(a8,i1)') 'xprestot',year_xvik
       

       ! Dmin Imin files
       
       write(*,*) 'Opening ', trim(dset_DminImin)//'/'//trim(filename1)
       open(21,file=trim(dset_DminImin)//'/'//trim(filename1),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DminImin)//'/'
     &//trim(filename1)
          stop 
       endif
       
       write(*,*) 'Opening ', trim(dset_DminImin)//'/'//trim(filename2)
       open(11,file=trim(dset_DminImin)//'/'//trim(filename2),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DminImin)//'/'
     &//trim(filename2)
          stop
       endif  
       
       ! Dmin Imax files
           
       write(*,*) 'Opening ', trim(dset_DminImax)//'/'//trim(filename1)
       open(22,file=trim(dset_DminImax)//'/'//trim(filename1),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DminImax)//'/'
     &//trim(filename1)
          stop
       endif
       
       write(*,*) 'Opening ', trim(dset_DminImax)//'/'//trim(filename2)
       open(12,file=trim(dset_DminImax)//'/'//trim(filename2),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DminImax)//'/'
     &//trim(filename2)
          stop
       endif        
       
       ! Dmax Imin files
       
       write(*,*) 'Opening ', trim(dset_DmaxImin)//'/'//trim(filename1)
       open(23,file=trim(dset_DmaxImin)//'/'//trim(filename1),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DmaxImin)//'/'
     &//trim(filename1)
          stop
       endif
       
       write(*,*) 'Opening ', trim(dset_DmaxImin)//'/'//trim(filename2)
       open(13,file=trim(dset_DmaxImin)//'/'//trim(filename2),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DmaxImin)//'/'
     &//trim(filename2)
          stop
       endif        
       
      
       ! Dmax Imax files
      
       write(*,*) 'Opening ', trim(dset_DmaxImax)//'/'//trim(filename1)
       open(24,file=trim(dset_DmaxImax)//'/'//trim(filename1),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DmaxImax)//'/'
     &//trim(filename1)
          stop
       endif
       
       write(*,*) 'Opening ', trim(dset_DmaxImax)//'/'//trim(filename2)
       open(14,file=trim(dset_DmaxImax)//'/'//trim(filename2),iostat=ie)
       
       if (ie.ne.0) then 
          write(*,*) 'Error opening file ',trim(dset_DmaxImax)//'/'
     &//trim(filename2)
          stop
       endif  
      
       

c *************************************************************
c      Observed Smooth Viking Curves
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       open(30,file='VL1')
       do n=1,nsol
           read(30,*) solvl(n), lsvl(n), pvl_obs(n)
       end do
       close(30)

c      Opening output file      

       !open(33, file = 'minimization.txt')


c      reading Simulation results 
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~


      do n=1,ngcm

c       Reading Viking Lander 1 simulated pressure for 4 runs
        
        if (time_unit == 1) then
        read(21,*) solgcm(n),  pvl_gcm(n,1)
        read(22,*) solgcm(n),  pvl_gcm(n,2)
        read(23,*) solgcm(n),  pvl_gcm(n,3)
        read(24,*) solgcm(n),  pvl_gcm(n,4)

c       Reading atmospheric and "cap" total pressure for all runs
        read(11,*) solgcm(n), patm(n,1), pcapn(n,1),pcaps(n,1)
        read(12,*) solgcm(n), patm(n,2), pcapn(n,2),pcaps(n,2)
        read(13,*) solgcm(n), patm(n,3), pcapn(n,3),pcaps(n,3)
        read(14,*) solgcm(n), patm(n,4), pcapn(n,4),pcaps(n,4)
	
        elseif (time_unit == 2) then
        read(21,*) lsgcm(n),  pvl_gcm(n,1)
        read(22,*) lsgcm(n),  pvl_gcm(n,2)
        read(23,*) lsgcm(n),  pvl_gcm(n,3)
        read(24,*) lsgcm(n),  pvl_gcm(n,4)

c       Reading atmospheric and "cap" total pressure for all runs
        read(11,*) lsgcm(n), patm(n,1), pcapn(n,1),pcaps(n,1)
        read(12,*) lsgcm(n), patm(n,2), pcapn(n,2),pcaps(n,2)
        read(13,*) lsgcm(n), patm(n,3), pcapn(n,3),pcaps(n,3)
        read(14,*) lsgcm(n), patm(n,4), pcapn(n,4),pcaps(n,4)
        call ls2sol(lsgcm(n),solconv)
        solgcm(n) = solconv

        elseif (time_unit == 3) then
        read(21,*) solgcm(n), lsgcm(n),  pvl_gcm(n,1)
        read(22,*) solgcm(n), lsgcm(n),  pvl_gcm(n,2)
        read(23,*) solgcm(n), lsgcm(n),  pvl_gcm(n,3)
        read(24,*) solgcm(n), lsgcm(n),  pvl_gcm(n,4)

c       Reading atmospheric and "cap" total pressure for all runs
        read(11,*) solgcm(n), lsgcm(n), patm(n,1), pcapn(n,1),pcaps(n,1)
        read(12,*) solgcm(n), lsgcm(n), patm(n,2), pcapn(n,2),pcaps(n,2)
        read(13,*) solgcm(n), lsgcm(n), patm(n,3), pcapn(n,3),pcaps(n,3)
        read(14,*) solgcm(n), lsgcm(n), patm(n,4), pcapn(n,4),pcaps(n,4)
	
	else 
	write(*,*) 'Wrong integer for xvik files format :',
     &' must be 1, 2 or 3'
        stop

        endif
	
c       Checking total CO2 inventory for all runs :
        do i=1,4
          if(n.eq.1) then
            ptot(i) = patm(n,i)+pcapn(n,i)+pcaps(n,i)
            write(*,*) 'For run = ',i,'  Ptot= ',ptot(i)
          else
           if(abs(patm(n,i)+pcapn(n,i)+pcaps(n,i)-ptot(i)).gt.3)then
              write(*,*)'total pressure not constant for run i= ',i
              write(*,*) 'n=',1,' ptot=',ptot(i)
              write(*,*)'n=',n,' ptot=',patm(n,i)+pcapn(n,i)+pcaps(n,i)
           end if 
          end if 
        end do
      end do

	close(11)
	close(12)
	close(13)
	close(14)
	close(21)
	close(22)
	close(23)
	close(24)
	
	
c      Smoothing simulated GCM Viking 1 pressure curves
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do n=1,nsol
        sol(n)=float(n)
      end do
      

c    Running average for both runs (with "boxsize" box in days)
      box=20.
      do i=1,4
         call runave(solgcm,pvl_gcm(1,i),ngcm,669.,box,
     &               sol,pvl_sm(1,i),nsol)
         call runave(solgcm,patm(1,i),ngcm,669.,box,
     &               sol,patm_sm(1,i),nsol)
         call runave(solgcm,pcapn(1,i),ngcm,669.,box,
     &               sol,pcapn_sm(1,i),nsol)
         call runave(solgcm,pcaps(1,i),ngcm,669.,box,
     &               sol,pcaps_sm(1,i),nsol)
      end do 
      

c      Computing mass cap  sensitivity to albedo (derivative) and alphaVl1
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      alphaVL1 = Pvl1/patm at a given time (see Hourdin et al., JGR,  1995)

       do n=1,nsol
         ! Evaluate derivative of pressure wrt Northern inertia, as the
         ! mean of the 2 values obtained using runs 1&2 and 3&4
         dpdicein(n)=(((pcapn_sm(n,1)-pcapn_sm(n,2))/
     &               (fonc(iceigcmN(1))-fonc(iceigcmN(2))))+
     &             ((pcapn_sm(n,3)-pcapn_sm(n,4))/
     &               (fonc(iceigcmN(3))-fonc(iceigcmN(4)))))*(1./2.)
         ! Evaluate derivative of pressure wrt Southern inertia, as the
         ! mean of the 2 values obtained using runs 1&2 and 3&4
         dpdiceis(n)=(((pcaps_sm(n,1)-pcaps_sm(n,2))/
     &               (fonc(iceigcmS(1))-fonc(iceigcmS(2))))+
     &             ((pcaps_sm(n,3)-pcaps_sm(n,4))/
     &               (fonc(iceigcmS(3))-fonc(iceigcmS(4)))))*(1./2.)
         ! Evaluate derivative of pressure wrt Northern ice depth coefficient,
         ! as the mean of the 2 values obtained using runs 1&3 and 2&4
         dpden(n)=(((pcapn_sm(n,1)-pcapn_sm(n,3))/
     &               (fonc2n(icedgcm(1))-fonc2n(icedgcm(3))))+
     &             ((pcapn_sm(n,2)-pcapn_sm(n,4))/
     &               (fonc2n(icedgcm(2))-fonc2n(icedgcm(4)))))*(1./2.)
         ! Evaluate derivative of pressure wrt Southern ice depth coefficient,
         ! as the mean of the 2 values obtained using runs 1&3 and 2&4
         dpdes(n)=(((pcaps_sm(n,1)-pcaps_sm(n,3))/
     &               (fonc2s(icedgcm(1))-fonc2s(icedgcm(3))))+
     &             ((pcaps_sm(n,2)-pcaps_sm(n,4))/
     &               (fonc2s(icedgcm(2))-fonc2s(icedgcm(4)))))*(1./2.)
         ! Evaluate alphaVL1 coefficient, as the mean of alphaVL1 coefficients
         ! of all 4 runs:
         alphavl1(n) = (1./4.)*(pvl_sm(n,1)/patm_sm(n,1)
     &      +pvl_sm(n,2)/patm_sm(n,2)+pvl_sm(n,3)/patm_sm(n,3)
     &      +pvl_sm(n,4)/patm_sm(n,4))
         !write(91,*) sol(n), dpden(n), dpdes(n)
c         write(*,*) 'pcapn_sm(n,1),pcapn_sm(n,2)'
c    &                ,pcapn_sm(n,1),pcapn_sm(n,2)
c         write(*,*)'fonc(albgcmN(1)),fonc(albgcmS(2))',
c    &    fonc(albgcmN(1)),fonc(albgcmS(2))
c         write(*,*)'albgcmN(1)),albgcmS(2)',
c    &    albgcmN(1),albgcmS(2)
         !write(90,*)pvl_sm(n,1)/patm_sm(n,1),pvl_sm(n,2)/patm_sm(n,2),
c    &              pvl_sm(n,3)/patm_sm(n,3),pvl_sm(n,4)/patm_sm(n,4)
       end do ! of do n=1,nsol

c      -------------------------------------------------
c      Smooth the derivative like Frederic did ? :
         call runave(sol,dpdicein,nsol,669.,60.,
     &               sol,dpdicein_fltr,nsol)
         call runave(sol,dpdiceis,nsol,669.,60.,
     &               sol,dpdiceis_fltr,nsol)
         call runave(sol,dpden,nsol,669.,60.,
     &               sol,dpden_fltr,nsol)
         call runave(sol,dpdes,nsol,669.,60.,
     &               sol,dpdes_fltr,nsol)
       do n=1,nsol
          dpdiceis(n)=dpdiceis_fltr(n)
          dpdicein(n)=dpdicein_fltr(n)
          dpdes(n)=dpdes_fltr(n)
          dpden(n)=dpden_fltr(n)
          !write(92,*) sol(n), dpden_fltr(n), dpdes_fltr(n)
       end do
c      -------------------------------------------------


c  Compute best fit parameters by minimizing Cost function 
c  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  STUPID ultra-robust minimization method: waisting CPU to save human time...

      open(33, file =  'minimization.txt')
c 'minimzation.txt' will contain data to plot COST(DN,DS) or COST(IN,IS) depending on plot_test


c to plot COST(DN,DS) (plot_test=1), loops must be in this order : DN->DS->IN->IS

      costmin = 1.e30 ! initialization
                  
      
      if(plot_test.eq.1) then
c     to plot COST(DN,DS) (plot_test=1), loops must be in this order : DN->DS->IN->IS            
        write(*,*)'   DN        DS      cost       IN       IS       Ps'
        icedn=icedmin ! initialization		
        do while (icedn.le.icedmax) ! loop on northern ice depth coefficient
!         albn=albmin ! initialization
!         do while (albn.le.albmax) ! loop on northern cap albedo
         iceds=max(icedn-maxiceddiff,icedmin) ! initialization
         do while (iceds.le.min(icedn+maxiceddiff,icedmax))
          cost4plot = 1.e30 ! initializationqsub 
!          iceds=max(icedn-maxiceddiff,icedmin) ! initialization
!          do while (iceds.le.min(icedn+maxiceddiff,icedmax))
          icein=iceimin ! initialization
          do while (icein.le.iceimax) ! loop on northern ice inertia
           iceis=max(icein-maxiceidiff,iceimin) ! initialization
           do while (iceis.le.min(icein+maxiceidiff,iceimax))
            ptry=pmin ! initialization
            do while (ptry.le.pmax) ! loop on total pressure
              cost =0. !initialization
              do n=1,nsol  
	        !write(*,*) n
               ! Pressure corresponding to Northern cap
               pcapn_new=pcapn_sm(n,refrun)+
     &            (fonc(icein) -fonc(iceigcm_refN)) * dpdicein(n)
     &          + (fonc2n(icedn) -fonc2n(icedgcm_ref)) * dpden(n)
               pcapn_new= max(pcapn_new,0.)
               ! Pressure corresponding to Southern Cap
               pcaps_new=pcaps_sm(n,refrun)+
     &            (fonc(iceis) -fonc(iceigcm_refS)) * dpdiceis(n)
     &          + (fonc2s(iceds) -fonc2s(icedgcm_ref)) * dpdes(n)
               pcaps_new= max(pcaps_new,0.)
               ! Pressure at VL1 site
               pvl1 = alphavl1(n)*
     &                (ptry - pcapn_new - pcaps_new)
               ! cumulative squared differences between predicted
               ! and observed pressures at VL1 site
               cost= cost+ ( pvl_obs(n) - pvl1)**2 
	        
              end do ! of do n=1,nsol
              ! store parameters which lead to minimum cost
              ! (i.e. best fit so far)
              if(cost.lt.costmin) then
                 costmin = cost
                 pfit=ptry
                 iceinfit=icein
                 iceisfit=iceis
                 icednfit=icedn
                 icedsfit=iceds
              end if
	   
              if(cost.lt.cost4plot) then
c                RMS, best pressure and best albedos for these icedsivities value 
                 cost4plot = cost
                 iceis4plot=iceis
                ! iceds4plot=iceds
                 icein4plot=icein
                 ps4plot=ptry
              end if
              ptry=ptry+deltap ! increment ptry
            end do ! of do while (ptry.le.pmax)
            iceis=iceis+deltaicei ! increment iceis
           end do ! of do while (iceis.le.min(icein+maxiceidiff,iceimax))
           icein=icein+deltaicei ! increment icein
          enddo ! of do while (icein.le.iceimax)
!           iceds=iceds+deltaiced ! increment iceds
!          end do ! of do while (iceds.le.min(icedn+maxiceddiff,icedmax))
!          write(*,fmt='(1pe9.2,f5.2,f9.3,1pe9.2,
!     &                  f5.2,f7.2)') 
           write(*,fmt='(6(1pe10.3))') ! output to screen
     &          icedn,iceds,sqrt(cost4plot/float(nsol)),
     &          icein4plot,iceis4plot,ps4plot
           write(33,fmt='(6(1pe10.3))') ! output to file
     &          icedn,iceds,sqrt(cost4plot/float(nsol)),
     &          icein4plot,iceis4plot,ps4plot
          iceds=iceds+deltaiced ! increment iceds
         enddo ! of do while (iceds.le.min(icedn+maxiceddiff,icedmax))
         write(33,*)' ' ! blank line in output file
!          albn=albn+deltaalb ! increment albn
!         end do ! of do while (albn.le.albmax)
         icedn=icedn+deltaiced
        end do ! do while (icedn.le.icedmax)

      elseif(plot_test.eq.2) then    
c     to plot COST(IN,IS) (plot_test=2), loops must be in this order : IN->IS->DN->DS            
        write(*,*)'   DN        DS      cost       IN       IS       Ps'
        write(33,*)'   DN        DS      cost       IN       IS     
     &Ps'	
        icein=iceimin ! initialization		
        do while (icein.le.iceimax) ! loop on northern ice depth coefficient
!         albn=albmin ! initialization
!         do while (albn.le.albmax) ! loop on northern cap albedo
         iceis=max(icein-maxiceidiff,iceimin)! initialization
         do while (iceis.le.min(icein+maxiceidiff,iceimax))
          cost4plot = 1.e30 ! initializationqsub 
!          iceds=max(icedn-maxiceddiff,icedmin) ! initialization
!          do while (iceds.le.min(icedn+maxiceddiff,icedmax))
          icedn=icedmin ! initialization
          do while (icedn.le.icedmax) ! loop on northern ice inertia
           iceds=max(icedn-maxiceddiff,icedmin) ! initialization
           do while (iceds.le.min(icedn+maxiceddiff,icedmax))
            ptry=pmin ! initialization
            do while (ptry.le.pmax) ! loop on total pressure
              cost =0. !initialization
              do n=1,nsol  
	        !write(*,*) n
               ! Pressure corresponding to Northern cap
               pcapn_new=pcapn_sm(n,refrun)+
     &            (fonc(icein) -fonc(iceigcm_refN)) * dpdicein(n)
     &          + (fonc2n(icedn) -fonc2n(icedgcm_ref)) * dpden(n)
               pcapn_new= max(pcapn_new,0.)
               ! Pressure corresponding to Southern Cap
               pcaps_new=pcaps_sm(n,refrun)+
     &            (fonc(iceis) -fonc(iceigcm_refS)) * dpdiceis(n)
     &          + (fonc2s(iceds) -fonc2s(icedgcm_ref)) * dpdes(n)
               pcaps_new= max(pcaps_new,0.)
               ! Pressure at VL1 site
               pvl1 = alphavl1(n)*
     &                (ptry - pcapn_new - pcaps_new)
               ! cumulative squared differences between predicted
               ! and observed pressures at VL1 site
               cost= cost+ ( pvl_obs(n) - pvl1)**2 
	        
              end do ! of do n=1,nsol
              ! store parameters which lead to minimum cost
              ! (i.e. best fit so far)
              if(cost.lt.costmin) then
                 costmin = cost
                 pfit=ptry
                 iceinfit=icein
                 iceisfit=iceis
                 icednfit=icedn
                 icedsfit=iceds
              end if
	   
              if(cost.lt.cost4plot) then
c                RMS, best pressure and best albedos for these icedsivities value 
                 cost4plot = cost
                 iceds4plot=iceds
                ! iceds4plot=iceds
                 icedn4plot=icedn
                 ps4plot=ptry
              end if
              ptry=ptry+deltap ! increment ptry
            end do ! of do while (ptry.le.pmax)
            iceds=iceds+deltaiced ! increment iceis
           end do ! of do while (iceis.le.min(icein+maxiceidiff,iceimax))
           icedn=icedn+deltaiced ! increment icein
          enddo ! of do while (icein.le.iceimax)
!           iceds=iceds+deltaiced ! increment iceds
!          end do ! of do while (iceds.le.min(icedn+maxiceddiff,icedmax))
!          write(*,fmt='(1pe9.2,f5.2,f9.3,1pe9.2,
!     &                  f5.2,f7.2)') 
           write(*,fmt='(6(1pe10.3))') ! output to screen
     &          icedn4plot,iceds4plot,sqrt(cost4plot/float(nsol)),
     &          icein,iceis,ps4plot
           write(33,fmt='(6(1pe10.3))') ! output to file
     &          icedn4plot,iceds4plot,sqrt(cost4plot/float(nsol)),
     &          icein,iceis,ps4plot           
          iceis=iceis+deltaicei ! increment iceds
         enddo ! of do while (iceds.le.min(icedn+maxiceddiff,icedmax))
         write(33,*)' ' ! blank line in output file
!          albn=albn+deltaalb ! increment albn
!         end do ! of do while (albn.le.albmax)
         icein=icein+deltaicei
        end do ! do while (icedn.le.icedmax)        
      else 
       write(*,*) 'Wrong integer for plot (must be 1 or 2)'
       stop
       
      endif
                
      close(33)  
                     
      write(*,*) 'Best fit Ptot=', pfit
      write(*,*) 'Best fit In=', iceinfit
      write(*,*) 'Best fit IS=', iceisfit
      write(*,*) 'Best fit Dn=', icednfit
      write(*,*) 'Best fit DS=', icedsfit
      write(*,*) 'RMS difference model/obs=',sqrt(costmin/float(nsol)) 

c  Synthethic VL1 pressure curves for information 
c  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      open(41,file='pvlfit')
      open(42,file='xprestotfit')
      do n=1,nsol
           pcapn_new=pcapn_sm(n,refrun)+
     &     (fonc(iceinfit) -fonc(iceigcm_refN)) * dpdicein(n)
     &      + (fonc2n(icednfit) -fonc2n(icedgcm_ref)) * dpden(n)
           pcapn_new= max(pcapn_new,0.)
           pcaps_new=pcaps_sm(n,refrun)+
     &     (fonc(iceisfit) -fonc(iceigcm_refS)) * dpdiceis(n)
     &      + (fonc2s(icedsfit) -fonc2s(icedgcm_ref)) * dpdes(n)
           pcaps_new= max(pcaps_new,0.)
	
      call sol2ls(sol(n),solconv)   
      
           write(41,*) sol(n), solconv, alphavl1(n)*
     &        ( pfit - pcapn_new - pcaps_new)
           write(42,*) sol(n),solconv, pfit - pcapn_new - pcaps_new,
     &           pcapn_new, pcaps_new, pfit
      end do
      close(41)
      close(42)      
     
      end 

c *****************************************************************
       SUBROUTINE RUNAVE (x,y,nmax,xperiod,xave,xsmooth,ysmooth,nsm)


       IMPLICIT NONE 

c -----------------------------------------------------------
c      Computing runnig average ysmooth on xsmooth coordinate for a periodic 
c      field y in coordinate x PERIODIC between 0 and xperiod
c      Averaging box size : xave
c      F.Forget 1999
c -----------------------------------------------------------

       integer nmax,nsm
       real y(nmax), ysmooth(nsm)
       real x(nmax), xsmooth(nsm)
       real xperiod,xave
       integer n,i, nave, imin

       integer nbig
       parameter (nbig=99999999)
       real xx(nbig)

c -----------------------------------------------------------


      if (nbig.lt.3*nmax) stop 'Must increase nbig in runave'

c     Reindexation des donnees
      do n=1,nmax
          xx(n) = x(n) - xperiod
          xx(n+nmax) = x(n) 
          xx(n+2*nmax) = x(n) +xperiod
      end do

c     Moyenne glissante
      imin=1
      do n=1,nsm
          ysmooth(n) =0.
          nave =0
          do i=imin,3*nmax
             if (xx(i).ge.(xsmooth(n)-0.5*xave)) then 
                if (xx(i).gt.(xsmooth(n)+0.5*xave)) goto 999
                ysmooth(n) =  ysmooth(n) + y(mod(i-1,nmax)+1)
                nave = nave +1
             end if
          end do
 999      continue
          imin = i -1 -nave
         ysmooth(n) = ysmooth(n)/ float (nave)
      end do
      end 

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
      
      subroutine sol2ls(sol,Ls)
!==============================================================================
! Purpose: 
! Convert a date/time, given in sol (martian day),
! into solar longitude date/time, in Ls (in degrees),
! where sol=0 is (by definition) the northern hemisphere
!  spring equinox (where Ls=0).
!==============================================================================
! Notes:
! Even though "Ls" is cyclic, if "sol" is greater than N (martian) year,
! "Ls" will be increased by N*360
! Won't work as expected if sol is negative (then again,
! why would that ever happen?)
!==============================================================================

      implicit none

!==============================================================================
! Arguments:
!==============================================================================
      real,intent(in) :: sol
      real,intent(out) :: Ls

!==============================================================================
! Local variables:
!==============================================================================
      real year_day,peri_day,timeperi,e_elips,twopi,degrad
      data year_day /669./            ! # of sols in a martian year
      data peri_day /485.0/           
      data timeperi /1.9082314/ 
      data e_elips  /0.093358/
      data twopi       /6.2831853/    ! 2.*pi
      data degrad   /57.2957795/      ! pi/180

      real zanom,xref,zx0,zdx,zteta,zz

      integer count_years
      integer iter

!==============================================================================
! 1. Compute Ls
!==============================================================================

      zz=(sol-peri_day)/year_day
      zanom=twopi*(zz-nint(zz))
      xref=abs(zanom)

!  The equation zx0 - e * sin (zx0) = xref, solved by Newton
      zx0=xref+e_elips*sin(xref)
      do iter=1,20 ! typically, 2 or 3 iterations are enough
         zdx=-(zx0-e_elips*sin(zx0)-xref)/(1.-e_elips*cos(zx0))
         zx0=zx0+zdx
         if(abs(zdx).le.(1.e-7)) then
!            write(*,*)'iter:',iter,'     |zdx|:',abs(zdx)
             exit
         endif 
      enddo

      if(zanom.lt.0.) zx0=-zx0

      zteta=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
      Ls=zteta-timeperi

      if(Ls.lt.0.) then
         Ls=Ls+twopi
      else
         if(Ls.gt.twopi) then
            Ls=Ls-twopi
         endif
      endif

      Ls=degrad*Ls
! Ls is now in degrees

!==============================================================================
! 1. Account for (eventual) years included in input date/time sol
!==============================================================================

      count_years=0 ! initialize
      zz=sol  ! use "zz" to store (and work on) the value of sol
      do while (zz.ge.year_day)
          count_years=count_years+1
          zz=zz-year_day
      enddo

! Add 360 degrees to Ls for every year
      if (count_years.ne.0) then
         Ls=Ls+360.*count_years
      endif


      end subroutine sol2ls

