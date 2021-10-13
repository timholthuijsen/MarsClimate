










      PROGRAM xvik

      USE filtreg_mod, ONLY: inifilr
      USE comconst_mod, ONLY: dtvr,g,r,pi
     
      
      IMPLICIT NONE
      
      
c=======================================================================
c
c  Pression au site Viking
c
c=======================================================================


c-----------------------------------------------------------------------
c   declarations:
c-----------------------------------------------------------------------


      include "dimensions.h"
      include "paramet.h"
      include "comdissip.h"
      include "comgeom2.h"
      include "netcdf.inc"      


      INTEGER itau,nbpas,nbpasmx 
      PARAMETER(nbpasmx=1000000)
      REAL temps(nbpasmx)
      INTEGER unitlec
      INTEGER i,j,l,jj
      REAL constR

c   Declarations NCDF:
c   -----------------
      CHARACTER*100  varname
      INTEGER ierr,nid,nvarid,dimid
      LOGICAL nc
      INTEGER start_ps(3),start_temp(4),start_co2ice(3)
      INTEGER count_ps(3),count_temp(4),count_co2ice(3)

c   declarations pour les points viking:
c   ------------------------------------
      INTEGER ivik(2),jvik(2),ifile(2),iv
      
      REAL, PARAMETER ::  lonvik1 = -47.95
      REAL, PARAMETER ::  latvik1 =  22.27
      REAL, PARAMETER ::  lonvik2 =  134.29
      REAL, PARAMETER ::  latvik2 =  47.67
      
      REAL, PARAMETER :: phivik1 = -3637
      REAL, PARAMETER :: phivik2 = -4505
      
      
      REAL lonvik(2),latvik(2),phivik(2),phisim(2)
      REAL unanj

c   variables meteo:
c   ----------------
      REAL vnat(iip1,jjm,llm),unat(iip1,jjp1,llm)
      REAL t(iip1,jjp1,llm),ps(iip1,jjp1),pstot, phis(iip1,jjp1)
      REAL co2ice(iip1,jjp1), captotN,captotS
      real t7(iip1,jjp1) ! temperature in 7th atmospheric layer

      REAL zp1,zp2,zp2_sm,zu,zv,zw(0:1,0:1,2),zalpha,zbeta

      LOGICAL firstcal
      INTEGER*4 day0

      REAL ziceco2(iip1,jjp1)
      REAL day,zt,sollong,sol,dayw,dayw_ls
      REAL airtot1,gh

      INTEGER ii,iyear,kyear

      CHARACTER*2 chr2

       
c   declarations de l'interface avec mywrite:
c   -----------------------------------------

      CHARACTER file*80
      CHARACTER pathchmp*80,pathsor*80,nomfich*80
      
      INTEGER Time_unit
      

c   externe:
c   --------

      EXTERNAL iniconst,inigeom,covcont,mywrite
      EXTERNAL exner,pbar
      EXTERNAL coordij,moy2
      EXTERNAL SSUM
      REAL SSUM
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c-----------------------------------------------------------------------
c   initialisations:
c-----------------------------------------------------------------------

      chr2="0"
      unanj=669.
      print*,'WARNING!!!',unanj,'Jours/an'
      nc=.true.
      
      phivik(1) = phivik1
      phivik(2) = phivik2
      
      print *, 'COORDVIKIIIN', latvik, lonvik
      print*, 'LES PHIVIK', phivik
      
      



      WRITE(*,*) 'Chemin des fichiers histoires'
      READ (*,'(a)')  pathchmp
      WRITE(*,*) 'Chemin des fichiers sorties'
      READ (*,'(a)')  pathsor
      
      WRITE(*,*) 'Fichiers de sortie en sol (1) 
     &,en ls (2) ,les deux (3)'
      READ (*,*)  Time_unit
      
      
      write (*,*)'>>>>>>>>>>>>>>>>', phivik,g
      DO iv=1,2
         phivik(iv)=phivik(iv)*3.73
      END DO

c-----------------------------------------------------------------------
c   ouverture des fichiers xgraph:
c-----------------------------------------------------------------------
      ifile(1)=12
      ifile(2)=13
      kyear=-1
      unitlec=11
      
      
      print*,'Entrer un fichier NC (sans le .nc)'
      READ(5,'(a)',err=9999) nomfich
      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   grande boucle sur les fichiers histoire:
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      firstcal=.true.
      DO WHILE(len_trim(nomfich).GT.0.AND.len_trim(nomfich).LT.50)
      PRINT *,'>>>  nomfich : ',trim(nomfich)

c----------------------------------------------------------------------
c   Ouverture des fichiers histoire:
c----------------------------------------------------------------------

      file=pathchmp(1:len_trim(pathchmp))//'/'//
     s     nomfich(1:len_trim(nomfich))
      PRINT*,'file.nc: ', file(1:len_trim(file))//'.nc'
      PRINT*,'timestep ',dtvr

      IF(nc) THEN
      ierr= NF_OPEN(file(1:len_trim(file))//'.nc',NF_NOWRITE,nid)        
      ELSE
         PRINT*,'Ouverture binaire ',file
         OPEN(unitlec,file=file,status='old',form='unformatted',
     .   iostat=ierr)
      ENDIF

c----------------------------------------------------------------------
c   initialisation de la physique:
c----------------------------------------------------------------------

      CALL readhead_NC(file(1:len_trim(file))//'.nc',day0,phis,constR)

      WRITE (*,*) 'day0 = ' , day0

      CALL conf_gcm( 99, .TRUE. )
      CALL iniconst
      CALL inigeom

c----------------------------------------------------------------------
c   Lecture temps :
c----------------------------------------------------------------------


      ierr= NF_INQ_DIMID (nid,"Time",dimid)
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'xvik: Le champ <Time> est absent'
          CALL abort
        ENDIF

      ierr= NF_INQ_DIMLEN (nid,dimid,nbpas)

      ierr = NF_INQ_VARID (nid, "Time", nvarid)
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, temps)
        IF (ierr.NE.NF_NOERR) THEN
          PRINT*, 'xvik: Lecture echouee pour <Time>'
          CALL abort
        ENDIF

        PRINT*,'temps(1:10)',(temps(itau),itau=1,10)
        
        
c-----------------------------------------------------------------------
c   coordonnees des point Viking:
c   --------------------------------------------------------------------

      lonvik(1) = lonvik1 * pi/180.
      latvik(1) = latvik1 * pi/180.
      lonvik(2) = lonvik2 * pi/180.
      latvik(2) = latvik2 * pi/180.
      
                    
c----------------------------------------------------------------------   
c   ponderations pour les 4 points autour de Viking
c----------------------------------------------------------------------


      DO iv=1,2
        ! locate index of GCM grid points near VL
         do i=1,iim
           ! we know longitudes are ordered -180...180
           if ((lonvik(iv).ge.rlonu(i)).and.
     &         (lonvik(iv).le.rlonu(i+1))) then
             ivik(iv)=i
             exit
           endif
         enddo
         do j=1,jjm-1
           !we know tha latitudes are ordered 90...-90
           if ((latvik(iv).le.rlatv(j)).and.
     &         (latvik(iv).ge.rlatv(j+1))) then
             jvik(iv)=j
             exit
           endif
         enddo
         zalpha=(lonvik(iv)-rlonu(ivik(iv)))/
     s          (rlonu(ivik(iv)+1)-rlonu(ivik(iv)))
         zbeta=(latvik(iv)-rlatv(jvik(iv)))/
     s          (rlatv(jvik(iv)+1)-rlatv(jvik(iv)))
         zw(0,0,iv)=(1.-zalpha)*(1.-zbeta)
         zw(1,0,iv)=zalpha*(1.-zbeta)
         zw(0,1,iv)=(1.-zalpha)*zbeta
         zw(1,1,iv)=zalpha*zbeta
      ENDDO

c----------------------------------------------------------------------
c   altitude reelle et modele aux points Viking
c----------------------------------------------------------------------


      DO iv=1,2
         phisim(iv)=0.
         DO jj=0,1
            j=jvik(iv)+jj
            DO ii=0,1
               i=ivik(iv)+ii
               phisim(iv)=phisim(iv)+zw(ii,jj,iv)*phis(i,j)
            ENDDO
         ENDDO
      ENDDO
      PRINT*,'relief aux points Viking pour les sorties:',phivik
           

c----------------------------------------------------------------------
c   lectures des etats:
c   -------------------------------------------------------------------

       airtot1=1./(SSUM(ip1jmp1,aire,1)-SSUM(jjp1,aire,iip1))

c======================================================================
c   debut de la boucle sur les etats dans un fichier histoire:
c======================================================================


       count_ps=(/iip1,jjp1,1/)
       count_co2ice=(/iip1,jjp1,1/)
       count_temp=(/iip1,jjp1,llm,1/)
       
       DO itau=1,nbpas

       start_ps=(/1,1,itau/)
       start_co2ice=(/1,1,itau/)
       start_temp=(/1,1,1,itau/)
       
c----------------------------------------------------------------------       
c   lecture drs des champs:
c----------------------------------------------------------------------


ccccccccc  LECTURE Ps ccccccccccccccccccccccccccc


          ierr = NF_INQ_VARID (nid, "ps", nvarid)
          ierr = NF_GET_VARA_DOUBLE(nid, nvarid,start_ps,count_ps, ps)
          IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'xvik: Lecture echouee pour <ps>'
            CALL abort
          ENDIF
          
          PRINT*,'ps',ps(iip1/2,jjp1/2)

ccccccccc  LECTURE Temperature ccccccccccccccccccccccccccc


          ierr = NF_INQ_VARID (nid, "temp", nvarid)
          ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start_temp,count_temp, t)
          IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'xvik: Lecture echouee pour <temp>'
            ! Ehouarn: proceed anyways
            ! CALL abort
            write(*,*)'--> Setting temperature to zero !!!'
            t(1:iip1,1:jjp1,1:llm)=0.0
            write(*,*)'--> looking for temp7 (temp in 7th layer)'
            ierr=NF_INQ_VARID(nid,"temp7", nvarid)
            if (ierr.eq.NF_NOERR) then
            write(*,*) "    OK, found temp7 variable"
            ierr=NF_GET_VARA_DOUBLE(nid,nvarid,start_ps,count_ps,t7)
              if (ierr.ne.NF_NOERR) then
                write(*,*)'xvik: failed loading temp7 !'
                stop
              endif
            else ! no 'temp7' variable
              write(*,*)'  No temp7 variable either !'
              write(*,*)'  Will have to to without ...'
              t7(1:iip1,1:jjp1)=0.0
            endif
          ELSE ! t() was successfully loaded, copy 7th layer to t7()
            t7(1:iip1,1:jjp1)=t(1:iip1,1:jjp1,7)
          ENDIF



ccccccccc  LECTURE co2ice ccccccccccccccccccccccccccc


          ierr = NF_INQ_VARID (nid, "co2ice", nvarid)
          ierr = NF_GET_VARA_DOUBLE(nid,nvarid,start_co2ice,
     &    count_co2ice,  co2ice)
          IF (ierr.NE.NF_NOERR) THEN
            PRINT*, 'xvik: Lecture echouee pour <co2ice>'
            CALL abort
          ENDIF

c----------------------------------------------------------------------
c Gestion du temps
c ---------------------------------------------------------------------

          day=temps(itau)
          PRINT*,'day ',day
          sol=day+day0
          iyear=sol/unanj
          WRITE (*,*) 'iyear',iyear
          sol=sol-iyear*unanj

c----------------------------------------------------------------------
c Ouverture / fermeture des fichiers
c ---------------------------------------------------------------------

          IF (iyear.NE.kyear) THEN
             WRITE(chr2(1:1),'(i1)') iyear+1
             WRITE (*,*) 'iyear bis',iyear
             WRITE (*,*) 'chr2'
             WRITE (*,*)  chr2
             IF(iyear.GE.9) WRITE(chr2,'(i2)') iyear+1
             kyear=iyear
             DO ii=1,2
                CLOSE(10+ifile(ii))
                CLOSE(2+ifile(ii))
                CLOSE(4+ifile(ii))
                CLOSE(6+ifile(ii))
                CLOSE(8+ifile(ii))
                CLOSE(16+ifile(ii))
                CLOSE(12+ifile(ii))
                CLOSE(14+ifile(ii))
                CLOSE(97)
                CLOSE(98)
             ENDDO
             CLOSE(5+ifile(1))
             OPEN(ifile(1)+10,file='xpsol1'//chr2,form='formatted')
             OPEN(ifile(2)+10,file='xpsol2'//chr2,form='formatted')                                  
             OPEN(97,file='xprestot'//chr2,form='formatted')

          ENDIF
 
          dayw = sol
          call sol2ls(sol,sollong)
          dayw_ls = sollong
          
          
          
c----------------------------------------------------------------------
c Calcul de la moyenne de pression planetaire
c ---------------------------------------------------------------------


          pstot=0.
          captotS=0.
          captotN=0.
          DO j=1,jjp1
             DO i=1,iim
                pstot=pstot+aire(i,j)*ps(i,j)
             ENDDO
          ENDDO
 
              DO j=1,jjp1/2
                 DO i=1,iim
                    captotN = captotN  +aire(i,j)*co2ice(i,j)
                 ENDDO
              ENDDO
              DO j=jjp1/2+1, jjp1
                 DO i=1,iim
                    captotS = captotS  +aire(i,j)*co2ice(i,j)
                 ENDDO
              ENDDO


c --------------Ecriture fichier sortie xprestot----------------------- 
c  Sol ou ls ou les deux 
c  Ps_moy_planetaire (Pa)
c  Pequivalente de glace de CO2 au Nord (si entierement sublimee) (Pa)
c  Pequivalente de glace de CO2 au Sud (si entierement sublimee) (Pa) 


         IF(Time_unit == 1) THEN
              WRITE(97,'(4e16.6)') dayw,pstot*airtot1
     &       , captotN*g*airtot1, captotS*g*airtot1       
 
         ELSEIF (Time_unit == 2) THEN    
              WRITE(97,'(4e16.6)') dayw_ls,pstot*airtot1
     &       , captotN*g*airtot1, captotS*g*airtot1
     
         ELSE 
             WRITE(97,'(5e16.6)') dayw,dayw_ls,pstot*airtot1
     &       , captotN*g*airtot1,captotS*g*airtot1
     
                    
         ENDIF           

c----------------------------------------------------------------------
c boucle sur les sites vikings:
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c interpolation de la temperature dans la 7eme couche, de la pression
c de surface et des vents aux points viking.
c----------------------------------------------------------------------

         IF(.NOT.firstcal) THEN
          
          DO iv=1,2

             zp1=0.
             zp2=0.
             zp2_sm=0.
             zt=0.

             DO jj=0,1
             
                j=jvik(iv)+jj
                
                DO ii=0,1
                
                   i=ivik(iv)+ii
                   zt=zt+zw(ii,jj,iv)*t7(i,j)
                   zp1=zp1+zw(ii,jj,iv)*log(ps(i,j)) ! interpolate in log(P)
                   WRITE (*,*) 'ps autour iv',ps(i,j),iv

                ENDDO
             ENDDO
             
             zp1=exp(zp1) ! because of the bilinear interpolation in log(P)
             WRITE (*,*) 'constR ',constR 
             WRITE (*,*) 'zt ',zt
             gh=constR*zt            
             
c---------------------------------------------------------------------- 
c  pression au sol extrapolee a partir de la temp. 7eme couche
c----------------------------------------------------------------------
           
             if (gh.eq.0) then ! if we don't have temperature values
               ! assume a scale height of 10km
               zp2=zp1*exp(-(phivik(iv)-phisim(iv))/(3.73*1.e4))
             else
               zp2=zp1*exp(-(phivik(iv)-phisim(iv))/gh)
             endif
            
          WRITE (*,*) 'iv,pstot,zp2, zp1, phivik(iv),phisim(iv),gh'
          WRITE (*,*) iv,pstot*airtot1,zp2,zp1,phivik(iv),phisim(iv),gh
             

c ------Ecriture 2 fichiers (1 pour Vl1, 1 pour VL2) sortie xpsol ------
c  Sol ou ls ou les deux
c  Ps site VLi (i=1,2) a  l'altitude GCM (Pa)
c  Ps site VLi (i=1,2) a  l'altitude exacte  (interpolee) (Pa)
              
             IF(Time_unit == 1) THEN
             	WRITE(ifile(iv)+10,'(3e15.5)') dayw,zp2,zp1
             ELSEIF (Time_unit == 2) THEN    
                WRITE(ifile(iv)+10,'(3e15.5)') dayw_ls,zp2,zp1
             ELSE   
                WRITE(ifile(iv)+10,'(4e15.5)') dayw,dayw_ls,zp2,zp1
             ENDIF     
              
          ENDDO

         ENDIF
         
         firstcal=.false.


c======================================================================
c   Fin de la boucle sur les etats du fichier histoire:
c======================================================================

      ENDDO

      ierr= NF_CLOSE(nid)

      PRINT*,'Fin du fichier',nomfich
      print*,'Entrer un nouveau fichier NC 
     &(sans le .nc) ou return pour finir'
      READ(5,'(a)',err=9999) nomfich


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c   Fin de la boucle sur les fichiers histoire:
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ENDDO

      PRINT*,'relief du point V1',.001*phis(ivik(1),jvik(1))/g
      PRINT*,'relief du point V2',.001*phis(ivik(2),jvik(2))/g
      DO iv=1,2
         PRINT*,'Viking',iv,'   i=',ivik(iv),'j  =',jvik(iv)
         WRITE(6,7777)
     s   (rlonv(i)*180./pi,i=ivik(iv)-1,ivik(iv)+2)
         print*
         DO j=jvik(iv)-1,jvik(iv)+2
            WRITE(6,'(f8.1,10x,5f7.1)')
     s   rlatu(j)*180./pi,(phis(i,j)/(g*1000.),i=ivik(iv)-1,ivik(iv)+2)
         ENDDO
         print*
         print*,'zw'
         write(6,'(2(2f10.4/))') ((zw(ii,jj,iv),ii=0,1),jj=0,1)
         print*,'altitude interpolee (km) ',phisim(iv)/1000./g
      ENDDO
      PRINT*,'R=',r
 9999  PRINT*,'Fin '

7777  FORMAT ('latitude/longitude',4f7.1)



      END

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
