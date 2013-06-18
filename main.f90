Program Pbl_flux
implicit real (A-H,O-Z)

 real(8) bx2(7), bix2(11)
 integer n
 character(17)  datechar
 real(8) airtemp, airhum, windspeed, soiltemp, pressure
 real(8) hflux1,hflux2,z0t1,z0t2,mflux1,mflux2
 real(8) tet2,tet1, tetz0T1, tetz0T2, tetz0t3,tetz0t4, tempdiff
 real(8) z0t, nut
 real(8) ro, cp
 real(8) deltat1, deltat2
 real(8) Us(1000),GRAD(1000),difftemp(1000,1000)
 integer method	! 1= simple iteration, 2= hords, 3 = no z0t correction , 4 = create function  with z0correction
 open(90,file ='istok2012time.txt',status='old')
 open(91,file ='istok2012temp.txt',status='old')
 open(92,file ='istok2012hum.txt',status='old')
 open(93,file ='istok2012wind.txt',status='old')
 open(94,file ='istok2012soil.txt',status='old')
 open(89,file = 'istok2012pres.txt',status='old')
 open(95,file ='istok2012flux.txt')
 open(96,file = 'istok2012diftemp.txt')
 open(196,file = 'func.STD',access = 'direct', recl = 8)
 open(197,file = 'difflux.STD',access = 'direct',recl =8)
 read(*,*)ifilelngth
 
 anu = 28.98
 method = 4 !
 Rconst = 8.31
 cp = 1006
 errmark = -999.9
 20 format(a17)
 21 format(f6.2)
 call PBLDAT
 bix2 = 0
 bx2(1) = 0
 bx2(2) =  0
 bx2(3) =  0
 bx2(4) =  0
 bx2(5) =  0
 bx2(6) =  2
 bx2(7) =  0.02
 itdrag = 10   
   ro = 1.204
 do i = 1,ifilelngth

  read(90,20)datechar
  read(91,*)airtemp
  read(89,*) pressure
  read(92,*)airhum
  read(93,*)windspeed
  read(94,*)soiltemp
  write(*,*) airtemp, soiltemp, 'wooops'
 ! read(*,*) k
  if(abs((airtemp -errmark)*(airhum-errmark)*(windspeed-errmark)*(soiltemp - errmark)).gt.0.1) then

   ro = pressure*100./(anu*Rconst*0.5*(airtemp+soiltemp +273.15*2))
  TET2     = (airtemp+273.15)*(100000./100000)**0.286
      TET1     = (soiltemp  +273.15)*(100000./100000)**0.286
	  if(abs(tet2-tet1).ge. 50.) then
	  write(*,*) 'ERROR!!!!!!!'
	  read(*,*)k
	  endif
	
 !tet2 = airtemp + 273.15
 !tet1 = soiltemp +273.15
 bx2(1) = windspeed
 bx2(2) = TET2
 tetz0t2 = tet1+0.03
 tetz0t1 = tet1
 iter2 = 1
 deltat1 = 0

if( method == 1) then
 !! SIMPLE ITERATION METHOD
do 	while( abs(tetz0t1-tetz0t2).gt. 0.002)

	  bx2(2) = TET2
	  bx2(3) = tetz0t1
	!  if(abs(tet2-tet1).ge. 50.) then
	!  write(*,*) 'ERROR!!!!!!!'
	!  endif
	 ! write(*,*) tet2,tet1
	 iter2 = iter2+1

	  call dragvl(bx2,bix2,itdrag)
	  c_u = bix2(8)
	  c_t = bix2(9)
	   if(iter2 == 1) then
	   hflux2   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t1)
	   mflux2  = ro*c_u*c_u*windspeed*windspeed
	   z0t2  = bix2(6)
	   endif
       hflux1   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t1)
	   mflux1  = ro*c_u*c_u*windspeed*windspeed
	   z0t1  = bix2(6)
	   tetz0t2= tetz0t1
		nut = 0.000068*(tetz01  - 273.15) + 0.0237

	 tetz0t1 = tet1 - hflux1*z0t1/(nut)
	 deltaT2 = abs(tetz0t1-tetz0t2)
	 deltaFI = abs((hflux1*z0t1 - hflux2*z0t2)/(deltaT2*nut))
	 if(deltaFI.gt.0) then
	 write(*,*) datechar,'new thetaz0t', tetz0t1,'new flux',hflux1,'diff', deltaFI
	 	! read(*,*)k
	 endif
	 if ( abs(deltat2- deltat1).lt.0.002) then
	   tetz0t1 = (tetz0t1	+  tetz0t2)/2.

	!   read(*,*) k 
	   else 
	   tetz0t1 = (tetz0t1	+  tetz0t2)/2.
	  endif
	  hflux2=hflux1
	  z0t2=z0t1

	  deltat1= deltat2
	! write(*,*) datechar,tet2, tetz0t1, (tet2-tetz0t1), hflux, iter2
	 enddo
	 !! END OF SIMLE ITERATION METHOD
	 endif

	 if( method == 2) then 
	 ! Chords method
       tet2 = airtemp + 273.15
 tet1 = soiltemp +273.15
 tetz0t2 = tet1
 tetz0t1 = (tet1+tet2)/2. +1
 iter2 = 1
 nut = 0.000068*(tetz02  - 273.15) + 0.0237
 do while (abs(tetz0t1-tetz0t2).gt.0.002) 
 if(iter2 ==1) then

	bx2(3) = tetz0t1
	call dragvl(bx2,bix2,itdrag)
	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t1 = bix2(6)
	  hflux1   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t1)
	   mflux1  = ro*c_u*c_u*windspeed*windspeed
	 fx1 = Tet1 - hflux1*z0t1/nut -  tetz0t1
	 write(*,*) 'FX1',fx1, hflux1,z0t1,	 hflux1*z0t1/nut
	 
	 endif
	 iter2=iter2+1

	 bx2(3)  = tetz0t2
	  call dragvl(bx2,bix2,itdrag)
	  	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t2 = bix2(6)
	  hflux2   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t2)
	   mflux2  = ro*c_u*c_u*windspeed*windspeed
	   fx2 = tet1 - hflux2*z0t2/nut -  tetz0t2

	   tetz0t3 = tetz0t1 - (tetz0t2 - tetz0t1)*fx1/(fx2-fx1)
	   if(tetz0t3.gt. 350)then
	   write(*,*)datechar,tetz0t3, 'TEMPDIFF',tet1-tetz0t3, TET2- TET1
	   read(*,*) k 
	   endif
	   tetz0t1 = tetz0t2
	   hflux1=hflux2
	   mflux1 = mflux2
	   z0t1=z0t2
	   tetz0t2 = tetz0t3
	enddo
	  bx2(3)  = tetz0t3
	  call dragvl(bx2,bix2,itdrag)
	  	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t2 = bix2(6)
	  hflux2   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t2)
	   mflux2  = ro*c_u*c_u*windspeed*windspeed
	endif
!	if(iter2 ==1) then
!	method =3
!	endif
	  if( method ==3) then
	  bx2(3) = tet1
	  tetz0t3 = tet1
	   call dragvl(bx2,bix2,itdrag)
	  	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t2 = bix2(6)
	  hflux2   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t2)
	   mflux2  = ro*c_u*c_u*windspeed*windspeed
	  endif

	  write(*,*) datechar,'diff',(tetz0t3-tet1), hflux2, iter2
	 	  tempdiff = tetz0t3 - tet1
	  write(95,*) datechar,mflux2, hflux2, tempdiff
!	  if(abs(tempdiff).ge.8.0) then
!	  write(*,*)' correction is too high', datechar
!	  write(*,*) 'u = ',windspeed,'tet2 =', tet2
!	  write(*,*) 'tet0 =',tet1 , 'tetz0t =', tetz0t3
!	  write(*,*) 'Grad1 =', tet2 -tet1,'Grad2 = ',tet2 - tetz0t3
!	  read(*,*) k
!	  endif  
	  tempdiff = tetz0t3 - tet1
!	  write(96,*) datechar,tempdiff
!	  if(datechar == '19.08.2012 03:00')  then
!	  read(*,*) k
!	  endif 
	   

 !  write(*,20)datechar
 ! write(*,*)airtemp
 ! write(*,*)airhum
 ! write(*,*)windspeed
 ! write(*,*)soiltemp
 !read(*,*)j
 endif
 enddo
 close(90)
 close(91)
 close(92)
 close(93)
 close(94)
 close(95)
 close(96)
 if(method == 4) then
 n = 1
 do i = 1,100
  do j = 1, 100
     Us(i) = 1. + 10.*i/100.
	 GRAD(j) = -20. + 40.*j/100.
	  tet2 = 18 + 273.15

 tet1 =  grad(j) + tet2
 tetz0t2 = tet1
 tetz0t1 = (tet1+tet2)/2.
 tetz0t3 = (tet1+tet2)/2.
 tetz0t4  = tet2
 iter2 = 1
 nut = 0.000068*(tetz02  - 273.15) + 0.0237
 bx2(1) = Us(i)
 !write(*,*)'HI'
 goto 111
 do while (abs(tetz0t1-tetz0t2).gt.0.002) 
 if(iter2 ==1) then

	bx2(3) = tetz0t1
	call dragvl(bx2,bix2,itdrag)
	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t1 = bix2(6)
	  hflux1   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t1)
	   mflux1  = ro*c_u*c_u*windspeed*windspeed
	 fx1 = Tet1 - hflux1*z0t1/nut -  tetz0t1
!	 write(*,*) 'FX1',fx1, hflux1,z0t1,	 hflux1*z0t1/nut
	 
	 endif
	 iter2=iter2+1

	 bx2(3)  = tetz0t2
	  call dragvl(bx2,bix2,itdrag)
	  	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t2 = bix2(6)
	  hflux2   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t2)
	   mflux2  = ro*c_u*c_u*windspeed*windspeed
	   fx2 = tet1 - hflux2*z0t2/nut -  tetz0t2

	   tetz0t3 = tetz0t1 - (tetz0t2 - tetz0t1)*fx1/(fx2-fx1)
!	   if(tetz0t3.gt. 350)then
!	   write(*,*)datechar,tetz0t3, 'TEMPDIFF',tet1-tetz0t3, TET2- TET1 !
	 !  read(*,*) k 														!
	  ! endif
	   tetz0t1 = tetz0t2
	   hflux1=hflux2
	   mflux1 = mflux2
	   z0t1=z0t2
	   tetz0t2 = tetz0t3
	   write(*,*) iter2
	enddo
  !  goto 222
111 CONTINUE !Dihotomy method
   do while( abs(tetz0t1 -tetz0t2).gt.0.001)
    bx2(3) = tetz0t2
	  call dragvl(bx2,bix2,itdrag)
	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t2 = bix2(6)
	  hflux2   = -cp*ro*c_u*c_t*US(i)*(TET2-tetz0t2)
	   mflux2  = ro*c_u*c_u*Us(i)*Us(i)
	 fx2 = Tet1 - hflux2*z0t2/nut -  tetz0t2

   bx2(3) = tetz0t1
 !  write(*,*)'HI1'
    call dragvl(bx2,bix2,iterdrag)
	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t1 = bix2(6)
	  hflux1   = -cp*ro*c_u*c_t*Us(i)*(TET2-tetz0t1)
	   mflux1  = ro*c_u*c_u*Us(i)*Us(i)
	 fx1 = Tet1 - hflux1*z0t1/nut -  tetz0t1
	! write(*,*) tetz0t2,tetz0t1
	  if(fx1*fx2.le. 0.) then
	  tetz0t4 = tetz0t1
	   tetz0t1 = (tetz0t4 + tetz0t2)/2.
	   tetz0t2 = tetz0t2
	   !write(*,*) 'FIRSt',fx1,fx2,hflux2, tetz0t2,tetz0t1
	   
	   else
	   tetz0t3 = tetz0t1
	   tetz0t1 = (tetz0t4 + tetz0t1)/2.
	   tetz0t2 = tetz0t3
	   tetz0t4 = tetz0t4
	  ! write(*,*)'Second',fx1,fx2, tetz0t2, tetz0t1
	   endif
	   iter2 =iter2+1   
	enddo
  	   !read(*,*)k
	  bx2(3)  = tetz0t1
	  call dragvl(bx2,bix2,itdrag)
	  	 c_u = bix2(8)
	  c_t = bix2(9)
	  z0t2 = bix2(6)
	  hflux2   = -cp*ro*c_u*c_t*windspeed*(TET2-tetz0t3)
         bx2(3) = TET1
         call dragvl(bx2,bix2,itdrag)
         c_u = bix2(8)
         c_t = bix2(9)
!222 CONTINUE
         hflux1 = -cp*ro*c_u*c_t*windspeed*(TET2-TET1)
	 ! write(*,*) hflux2, tet2-tetz0t2,grad(j),Us(i), i,j 
	   mflux2  = ro*c_u*c_u*windspeed*windspeed
	   difftemp(i,j) = tetz0t3 - tet1
	   write(196,rec =n) difftemp(i,j)/z0t2
          write(197,rec=n) (hflux1-hflux2)
	   n =n+1
	   enddo 
	   enddo
	endif
	  


 end




 SUBROUTINE PBLDAT

!     PBLDAT assigns values to parameters used in surface layer parameterization
!     (see subr. dragvl and subr., called by dragvl)

      IMPLICIT real (A-H,O-Z)
!
!*=====================================================================
!*       INITIALIZATION OF PBL COMMON BLOCKS (FBL,DRAG)               =
!*    FBL(1)     :  FOR LINEAR EXTRAPOLATION OF WIND ONTO SURFACE     =
!*    FBL(2)     :  ---------------------------------------------     =
!*    FBL(3)     :  PRESCRIBED MINIMAL VALUE OF SURFACE WIND MODULE   =
!*    FBL(4)     :  HEIGHT OF CONSTANT FLUX LAYER ( M )               =
!*    FBL(5)     :  CRITICAL VALUE OF RELATIVE HUMIDITY IN            =
!*               :  CALCULATION OF EQUIVALENT POTENTIAL TEMPERATURE   =
!*    FBL(6)     :  SEA WATER TEMPERATURE UNDER ICE ( DEG. K )        =
!*    FBL(7)     :  PRECRIBED MINIMAL VALUE OF SEA ICE SURFACE        =
!*                  TEMPERATURE ( DEG. K )                            =
!*    FBL(8)     :  RESERVED                                          =
!*    FBL(9)     :  PBL WIND TURNING ANGLE OVER OCEAN ( DEG )         =
!*    FBL(10)    :  --------------------------- ICE ----------------- =
!*    FBL(11)    :  --------------------------- SNOW COVERED SOIL --- =
!*    FBL(12)    :  --------------------------- BARED SOIL ---------- =
!*    FBL(13)    :  --------------------------- IN TROPICS ---------- =
!*    FBL(14)    :  BOUNDARIES OF TROPICS ( RADIANS )                 =
!*    FBL(15)    :  RESERVED
!*    FBL(16)    :  (SNOW COVERED SOIL HEAT CONDUCTIVITY)/DEPTH ----- =
!*    FBL(17)    :  (SEA ICE HEAT CONDUCTIVITY)/(SEA ICE DEPTH) ----- =
!*    FBL(18)    :  (SOIL HEAT CONDUCTIVITY)/DEPTH, CAL/(DEG.*M**2*S) =
!*    FBL(19)    :  BULK SOIL HEAT CAPACITY, CAL/(DEG.*M**2)          =
!*    FBL(20)    :  (SOIL MOISTURE CONDUCTIVITY)/DEPTH,SEC**-1        =
!*=====================================================================
      COMMON /FBL/ FBL(20)
  !    COMMON /BL1/ AZ0P,VEG,WL,DZZ
  !    COMMON /HYDR/ SNCR,WLMMX,WSMAX,WSL,WSG,CEFF, &
   !   &              CA,DZL,DZG,BMIN,WSDL,WSDG,DMIN,DMAX,D,HR,WIINF, &
   !   &              ZRM,ZRMM,TRM,FLXMIN,TOMIN,HSNold
 !     COMMON /VEGSW/ CCQ,CK,SWW,TL
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &               alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &               zeta_st,a_m,a_t,b_m,b_t
      real FBL,AZ0P,VEG,WL,DZZ,SNCR,WLMMX,WSMAX,WSL,WSG,CEFF, &
      & CA,DZL,DZG,BMIN,WSDL,WSDG,DMIN,DMAX,D,HR,WIINF,ZRM,ZRMM,TRM, &
      & FLXMIN,TOMIN,CCQ,CK,SWW,TL, &
      & vkc,anu,z0min,chc,p_m,p_t,q_m,q_t,alpha_m,alpha_t,beta_m, &
      & beta_t,gamma_m,gamma_t,zeta_st,a_m,a_t,b_m,b_t

      FBL(1) = 1.
      FBL(2) = 1.
      FBL(3) = 1.
      FBL(4) = 70.
      FBL(5) = 1.
      FBL(6) = 271.5
      FBL(7) = 243.2
      FBL(8) = 0.
      FBL(9) = 20.
      FBL(10) = 10.
      FBL(11) = 30.
      FBL(12) = 30.
      FBL(13) = 0.
      FBL(14) = 20.*3.1415926536/180.
!     FBL(8)...FBL(14) WILL BE NOT USED IN THE MULTILEVEL PBL
      FBL(15) = 0.
      FBL(16) = 0.04/2.
      FBL(17) = 0.5/3.
      FBL(18) = 0.6/2.
      FBL(19) = 4.0E+04
      FBL(20) = 1.0E-6/5.
!------------------------------------------------------------------
!*    subroutine dragvl constants
!------------------------------------------------------------------
      vkc = 0.4
      anu = 0.000015
      z0min = 1.5e-5
      chc = 0.0132
      zeta_st = -2.
      p_m = - 0.25
      p_t = - 0.5
      q_m = - 1./3.
      q_t = - 1./3.
      alpha_m = 1.
      alpha_t = 1.
      beta_m = 4.7
      beta_t = 4.7
      gamma_m = 16.
      gamma_t = 16.
      a_m = (1. + gamma_m*(p_m/q_m - 1.)*zeta_st) &
      &      *(1. - gamma_m*zeta_st)**(p_m/q_m - 1.)
      b_m = a_m*gamma_m*(p_m/q_m)/(1. + gamma_m*(p_m/q_m - 1.)*zeta_st)
      a_t = (1. + gamma_t*(p_t/q_t - 1.)*zeta_st) &
      &      *(1. - gamma_t*zeta_st)**(p_t/q_t - 1.)
      b_t = a_t*gamma_t*(p_t/q_t)/(1. + gamma_t*(p_t/q_t - 1.)*zeta_st)

      RETURN
      END subroutine  PBLDAT


      SUBROUTINE dragvl(bx2,bix2,itdrag)
      
!     dragvl calculates exchange coefficients in aerodynamic formulas for surface
!     sensible heat, latent heat and momentum fluxes at the surface
!     following Businger-Dayer interpolation formulas, Beljaars parameterization and others

!     Input variables:
!     bx2(1) - the module of wind speed at the height Z, m/s
!     bx2(2) - the potential temperature at the height Z, K
!     bx2(3) - the potential temperature at the surface (water, snow, soil, etc.), K
!     bx2(4) - the specific humidity at the height Z, kg/kg
!     bx2(5) - the specific humidity at the surface (water, snow, soil, etc.), kg/kg
!     bx2(6) - the height Z, m
!     bx2(7) - the surface roughness parameter, m
!     itdrag - the number of iterations in dragvl, numbers from 5 to 10 are usually acceptable

!     Output variables:
!     bix2(1) - zeta = z/(monin - obukhov length)
!     bix2(2) - ...
!     bix2(3) - ...
!     bix2(4) - ...
!     bix2(5) - ...
!     bix2(6) - ...
!     bix2(7) - ...
!     bix2(8) - c_u, the non-dimensional exchange coefficient for momentum
!     bix2(9) - c_t, the non-dimensional exchange coefficient for temperature and humidity

!     Sensible_heat_flux = 
!     -specific_heat_of_air_at_constant_pressure*air_density*c_u*c_t*wind_speed_module*
!     (potent_temp_air_at_Z-potent_temp_air_at_surface)

!     For unfrozen surfaces
!     Latent_heat_flux = -heat_water_vapor_transition*air_density*c_u*c_t*wind_speed_module*
!     (spec_hum_air_at_Z-spec_hum_air_at_surf) 

!     or for frozen surfaces

!     Latent_heat_flux = -heat_ice_vapor_transition*air_density*c_u*c_t*wind_speed_module*
!     (spec_hum_air_at_Z-spec_hum_air_at_surf)

      real bx(7),bix(11)

      real(8), intent(in) :: bx2(7)
      real(8), intent(out):: bix2(11)
	  COMMON /FBL/ FBL(20)
          
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t,grav
      real u_a,tpot_a,tpot_s,q_a,q_s,z,z0m,z0t,tvir_a, &
      & tvir_s,dtvir,alam,c_u,c_t,ust,hfl,re,xx,zeta,zeta0m, &
      & zeta0t,smo,psi_m,psi_t,phi_m,phi_t
      INTEGER iter,itdrag
     
         
      grav = 9.81
      
      bx=sngl(bx2)  

      u_a = amax1(bx(1),0.1) ! speed 
      tpot_a = bx(2) ! potential temperature - air
      tpot_s = bx(3) ! potential temperature - surface
      q_a = bx(4)
      q_s = bx(5)
      z = bx(6)
      if(bx(7) .gt. 0.) then
          z0m = bx(7)
      else
          z0m = z0min
      end if
      tvir_a = tpot_a * (1. + 0.61 * q_a) ! virtual temp - air
     	tvir_s = tpot_s * (1. + 0.61 * q_s) ! virtual temp - surf
      dtvir = tvir_a - tvir_s  ! for smo
      alam = grav/(0.5*(tpot_a + tpot_s))  ! for smo
      c_u = vkc/log(z/z0m) ! for smo
      c_t = c_u
	  !open(110,file = 'z0t.txt')
      do iter = 1,itdrag
          ust = c_u * u_a  ! for smo
          hfl =  - c_u * c_t * u_a * dtvir ! for smo Hflux
          if(bx(7) .lt. 0.) z0m = amax1(chc*ust**2/grav, z0min)
          re = ust*z0m/anu
		  
          if(re .le. 0.111) then
              xx = - 2.43	  
          else
              if(0.111 .lt. re .and. re .le. 16.3) then
                  xx = 0.83*log(re) - 0.6
              else
                  xx = 0.49 * re**0.45
              end if
          end if
          z0t = amax1(z0m*exp(-xx), z0min)
!         z0t = z0m	 
        if(iter.eq. itdrag) then
		  !write(*,*) re,z0t, 'CAUTION!!'
		 
          endif
          if(abs(dtvir) .lt. 1.e-5) then
              zeta = 0
              zeta0m = 0
              zeta0t = 0
          else
              smo = - ust**3/(alam*vkc*hfl) ! vkc = 0.4
              smo = sign(1.,smo)*amax1(abs(smo),1.e-3) !
              zeta = z/smo
              zeta0m = z0m/smo
              zeta0t = z0t/smo
          end if
          c_u = vkc/amax1( (log(z/z0m) - psi_m(zeta,zeta0m)), 0.1)
          c_t = vkc/amax1( (log(z/z0t) - psi_t(zeta,zeta0t)), 0.1)
           
      end do
!         c_t = c_u
      bix(1) = zeta
      bix(2) = alam*z*dtvir/u_a**2
      bix(3) = re
      bix(4) = log(z0m/z0t)
      bix(5) = z0m																																																																																																									
      bix(6) = z0t
      bix(7) = 0.
      bix(8) = c_u
      bix(9) = c_t
      bix(10) = vkc*ust*z/phi_m(zeta)
      bix(11) = phi_m(zeta) / phi_t(zeta)

      bix2=dble(bix)

      return
      END     
                        

      real function phi_m(zeta)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,a
      a=1.0

!*   universal function for momentum: businger et al. (1971) formulation
!*   for zeta .ge. zeta_st and large et al. (1994) formulation otherwise
!*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
!         y = alpha_m + beta_m*zeta     !  businger et al.
          y = 1. + zeta*(a+0.667*(6.-0.35*zeta)*exp(-0.35*zeta)) ! b&h
      else
          if(zeta .ge. zeta_st) then
              y = (1. - gamma_m*zeta)**(p_m)
          else
              y = (a_m - b_m*zeta)**(q_m)
          end if
      end if
      phi_m = y

      return
      end

      real function phi_t(zeta)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,a
      a=1.0

!*   universal function for heat: businger et al. (1971) formulation
!*   for zeta .ge. zeta_st and large et al. (1994) formulation otherwise
!*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
!         y = alpha_t + beta_t*zeta
          y = 1. + zeta * &
          & (a*sqrt(1.+2.*zeta*a/3.)+0.667*(6.-0.35*zeta)*exp(-0.35*zeta)) ! b&h
      else
          if(zeta .ge. zeta_st) then
              y = (1. - gamma_t*zeta)**(p_t)
          else
              y = (a_t - b_t*zeta)**(q_t)
          end if
      end if
      phi_t = y

      return
      end

      real function psi_m1(x)

      real x,pi
      pi = 3.1415926
!*   integrated universal function for momentum based on businger et al.
!*   (1971) formulation for zeta_st .le. zeta .le. 0.
!*    p_m = - 1./4. (paulson, 1970)
      psi_m1 = log(0.125*(1. + x**2)*(1. + x)**2) -2.*atan(x) + pi/2.
      return
      end

      real function psi_m2(a,x)
      real x,a
!*   integrated universal function for momentum based on large et al.
!*   (1994) formulation for zeta .le. zeta_st
!*    only for q_m = -1./3.
      psi_m2 = (0.5*log(x**2 + x + 1.) &
      &         + sqrt(3.)*atan((2.*x + 1.)/sqrt(3.))) &
      &         / sign(1.,a)*abs(a)**(1./3.)
      return
      end

      real function psi_m(zeta,zeta0)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &                alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &                zeta_st,a_m,a_t,b_m,b_t
      real zeta,y,zeta0,x,x_st,psi_m2,psi_m1,x0,x0_st,a
      
      a=1.0
!*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
!         y = (1.-alpha_m)*log(zeta/zeta0) - beta_m*(zeta-zeta0)
          y = - (a*zeta + &
          &  0.667*(zeta-(5./0.35))*exp(-0.35*zeta)+(0.667*5./0.35))  ! b&h
      else
          if(zeta .ge. zeta_st) then
              x = (1. - gamma_m*zeta)**(- p_m)
              y = psi_m1(x)
          else
              x0 = 1. - (b_m/a_m)*zeta
              x = sign(1.,x0)*abs(x0)**(- q_m)
              x0_st = 1. - (b_m/a_m)*zeta_st
              x_st = sign(1.,x0_st)*abs(x0_st)**(- q_m)
!              y = psi_m1(zeta_st) + psi_m2(a_m,x) - psi_m2(a_m,x_st)
              y = psi_m1(x_st) + psi_m2(a_m,x) - psi_m2(a_m,x_st)
          end if
      end if
      psi_m = y
      return
      end

      real function psi_t1(x)
      real x
!*   integrated universal function for heat based on businger et al.
!*   (1971) formulation for zeta_st .le. zeta .le. 0.
!*    p_t = - 1./2. (paulson, 1970)
      psi_t1 = log(0.25*(1. + x)**2) !in (paulson, 1970) log(0.25*(1. + x**2)**2)!
      return
      end

      real function psi_t(zeta,zeta0)
      common /p_drag/ vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &               alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &               zeta_st,a_m,a_t,b_m,b_t
      real vkc,anu,z0min,chc,p_m,p_t,q_m,q_t, &
      &               alpha_m,alpha_t,beta_m,beta_t,gamma_m,gamma_t, &
      &               zeta_st,a_m,a_t,b_m,b_t
      real zeta,zeta0,y,x,x0,psi_t1,x0_st,x_st,psi_m2,a
      a=1.0
!*   for zeta .ge. 0. : beljaars and holtslag (1991)
      if(zeta .ge. 0.) then
!         y = (1.-alpha_t)*log(zeta/zeta0) - beta_t*(zeta-zeta0)
          y = - (a*(1.+(2./3.)*zeta*a)**(1.5) &
          & + 0.667*(zeta-(5./0.35))*exp(-0.35*zeta)+((0.667*5./0.35)-1.))
      else
          if(zeta .ge. zeta_st) then
              x = (1. - gamma_t*zeta)**(- p_t)
              y = psi_t1(x)
          else
              x0 = 1. - (b_t/a_t)*zeta
              x = sign(1.,x0)*abs(x0)**(- q_t)
              x0_st = 1. - (b_t/a_t)*zeta_st
              x_st = sign(1.,x0_st)*abs(x0_st)**(- q_t)
              y = psi_t1(x_st) + psi_m2(a_t,x) - psi_m2(a_t,x_st)
          end if
      end if
      psi_t = y											
      return
      end
