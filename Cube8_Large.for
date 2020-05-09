c     program 'Cube8.for' **********************************************
c     Finds ground state of six electrons inside a bubble. ************
c     Uses simple diffusion. *******************************************
	implicit real*8 (a-h,o-z)
      real*8 m,k_bar
	integer*4 t1,t2,clock_rate,clock_max
	h_bar=10.5459d0
	m=9.10956d-4
      pi=2*dacos(0.0d0)
c     ******************************************************************
c     Open file to read in parameters. ********************************* 
c     ******************************************************************
	open(1,file='Cube8.inp',status='unknown')
c     ******************************************************************
c     Open file for summary. *******************************************
c     ******************************************************************
	open(4,file='Cube8_summ.dat',status='unknown')
	write(4,*)
	close(4)
c     ******************************************************************
c     Read in parameters. **********************************************
c     ******************************************************************
c     Read in range of maximum of bubble radius. ***********************
	read(1,*) r_1_min,r_1_max,r_1_step 
c     Read in range of minimum of bubble radius. ***********************
      read(1,*) r_2_min,r_2_max,r_2_step 
c     Read in space step. **********************************************
      read(1,*) space_step  ! Space step.
c     Read in time step. ***********************************************
	read(1,*) time_step  
c     Read in number of iterations between saves. **********************
      read(1,*) n_save 
	close(1) 
	i1_max=idnint((r_1_max-r_1_min)/r_1_step)
	i2_max=idnint((r_2_max-r_2_min)/r_2_step)		  
      do i1=0,i1_max
        r_0=10.1  ! Starting guess for distance from end of the lobe. 
c       Note that this must not be a multiple of space_step or an error 
c       will occur because the potential energy can become infinite. ***
        r_1=r_1_min+dfloat(i1)*r_1_step
	  do i2=0,i2_max
	    call cpu_time(start_time)
c     ******************************************************************
          r_2=r_2_min+dfloat(i2)*r_2_step
c     Calculate volume and surface. ************************************
c     ******************************************************************
          call volume_surface(r_1,r_2,surface,volume)
c         **************************************************************
c         Calculate the energy. ****************************************
c         **************************************************************
          call system_clock ( t1, clock_rate, clock_max )
          call energy_calc(space_step,time_step,r_1,r_2,r_0,
     &	n_save,surface,e,k_bar,v_bar, r_2_min)
          call system_clock ( t2, clock_rate, clock_max )
            time = dble(t2-t1)/dble(clock_rate)
          write(*,147) r_1,r_2,r_0,e,v_bar,8*e-4*v_bar+surface*0.375,
     &                  volume,surface, time
147       format(' r_1',f7.1,' r_2',f7.1,' r_0',f8.3,
     &    '  e',f9.1,'  v_bar',f9.1,'  total e',f9.1, 'volume',f12.2,
     &          'surface', f12.2, 'time', f10.2)		 
	    open(4,file='Cube8_summ.dat',status='unknown',
     &	position='append')
	    
          write(4,148) r_1,r_2,r_0,e,v_bar,8*e-4*v_bar+surface*0.375,    
     &     volume, surface, time
148       format(3f8.3,f10.1,f10.1,f10.1,f12.2,f12.2,f8.3)
          close(4)
		write(*,*) 'Elapsed time = ',dble(t2-t1)/dble(clock_rate)		 
        end do
c       Insert a blank line into file.  ********************************
        open(4,file='Cube8_summ.dat',status='unknown',position='append')
        write(4,*)
	  close(4)       
      end do
	stop
      end
c     ******************************************************************
c     Set shape of the bubble
c     ******************************************************************
c      subroutine radius(r_1,r_2,theta,phi)
c            r_bubble=-0.5*r_2+1.5*r_1-
c     &	  1.5*(r_1-r_2)*(dsin(theta)**4*(dcos(phi)**4+
c     &      dsin(phi)**4)+dcos(theta)**4)
c      end
c     ******************************************************************
c     ******************************************************************
      subroutine energy_calc(space_step,time_step,r_1,r_2,r_0,
     &n_save,surface,e,k_bar,v_bar, r_2_min)  
c     This subroutine finds the energy and wave function for given values 
c     of r_1,r_2 and r_0. It then returns a new value of r_0_new.  
c     ******************************************************************
c     ******************************************************************
      implicit real*8 (a-h,o-z)
	real*8 psi(-200:200,-200:200,-200:200)
	real*8 psi_p(-200:200,-200:200,-200:200)
	real*8 v(-200:200,-200:200,-200:200)
      integer*4 i(-200:200,-200:200,-200:200)
      real*8 e_array(20000),v_bar_array(20000),r_0_array(20000)
      real*8 m,k_bar
	h_bar=10.5459d0
	m=9.10956d-4
      pi=2*dacos(0.0d0)

      
c     ******************************************************************
c     Set initial wave function. ***************************************
c     If r_2 equals r_2_min use a new initial wave function.  **********
c     Otherwise use wave function from last value of r_2. **************
c     Uninitialized values within the bubble are set to 1. *************
c     ******************************************************************
      n_1=r_1/space_step+1
      write(*,*) r_2, r_2_min
c      if (r_2 .gt. r_2_min) then
c            do n_x=-n_1,n_1
c                  do n_y=-n_1,n_1
c                        do n_z=-n_1,n_1
c                              if (psi(n_x, n_y, n_z) .eq. 0) then    
c                              psi(n_x,n_y,n_z)=1.0
c                              i(n_x,n_y,n_z)=1.0
c                              end if
c                        end do
c                  end do
c            end do 
c      else
            do n_x=-n_1,n_1
                  do n_y=-n_1,n_1
                        do n_z=-n_1,n_1    
                              psi(n_x,n_y,n_z)=1.0
                              i(n_x,n_y,n_z)=1
                        end do
                  end do
            end do 
c      end if

      do n_x=-n_1,n_1
	  do n_y=-n_1,n_1
          do n_z=-n_1,n_1
c           Check to see if outside bubble. ****************************
            eps=1.0d-9  ! To avoid problems with the datan2 function. 
	      x=n_x*space_step+eps
	      y=n_y*space_step+eps
	      z=n_z*space_step+eps
            r=dsqrt(x**2+y**2+z**2)
	      theta=dacos(z/r)
            phi=datan2(y,x)
c           Calculate radius of bubble in direction of walker. *********
c           ************************************************************

c           THIS DETERMINES THE SHAPE OF THE BUBBLE

c           ************************************************************

	      r_bubble=-0.5*r_2+1.5*r_1-
     &	  1.5*(r_1-r_2)*(dsin(theta)**4*(dcos(phi)**4+
     &      dsin(phi)**4)+dcos(theta)**4)
c           Check to see if outside. ***********************************
	      if (r .gt. r_bubble) then 
		    psi(n_x,n_y,n_z)=0.0 
	        i(n_x,n_y,n_z)=0
            end if    
	    end do 
	  end do 
	end do 
c     ******************************************************************
c     Begin iteration. *************************************************
c     ******************************************************************  
      n_iter=0
20    continue
      n_iter=n_iter+1
c     ******************************************************************
c     Calculate potential energy at all positions. *********************
c     ******************************************************************
	do n_x=-n_1,n_1
        do n_y=-n_1,n_1
	    do n_z=-n_1,n_1 
            call potential_energy(r_1,r_0,space_step,n_x,n_y,n_z,
     &      v(n_x,n_y,n_z))
	    end do 
	  end do 
	end do 
c     ******************************************************************
c     Change psi by diffusion. *****************************************
c     This also makes sure psi is zero outside the bubble. *************
c     ******************************************************************
	do n_x=-n_1,n_1
        do n_y=-n_1,n_1
	    do n_z=-n_1,n_1 
            del_sq_psi=(psi(n_x+1,n_y,n_z)+psi(n_x-1,n_y,n_z)+
     &	  psi(n_x,n_y+1,n_z)+psi(n_x,n_y-1,n_z)+psi(n_x,n_y,n_z+1)
     &      +psi(n_x,n_y,n_z-1)-6*psi(n_x,n_y,n_z))/space_step**2
	      psi_p(n_x,n_y,n_z)=i(n_x,n_y,n_z)*(psi(n_x,n_y,n_z)+
     &	  (h_bar/(2*m))*del_sq_psi*time_step) 
	    end do 
	  end do 
	end do 
c     ******************************************************************
c     Change psi due to the potential. *********************************
c     ******************************************************************
	do n_x=-n_1,n_1
        do n_y=-n_1,n_1
	    do n_z=-n_1,n_1 
            psi(n_x,n_y,n_z)=psi_p(n_x,n_y,n_z)*
     &	  dexp(-v(n_x,n_y,n_z)*time_step/h_bar)
	    end do 
	  end do 
	end do 
c     ******************************************************************
c     Calculate energy, potential energy and new r_0. ******************
c     ******************************************************************
      sum_psi=0.0		   
	sum_k=0.0
      sum_v=0.0
      sum_r_0=0.0
      do n_x=-n_1,n_1
        do n_y=-n_1,n_1
          do n_z=-n_1,n_1
            sum_psi=sum_psi+psi(n_x,n_y,n_z)**2
            del_sq_psi=(psi(n_x+1,n_y,n_z)+psi(n_x-1,n_y,n_z)+
     &      psi(n_x,n_y+1,n_z)+psi(n_x,n_y-1,n_z)+
     &      psi(n_x,n_y,n_z+1)+psi(n_x,n_y,n_z-1)-
     &      6*psi(n_x,n_y,n_z))/space_step**2
	      sum_k=sum_k-psi(n_x,n_y,n_z)*(h_bar**2/(2*m))*del_sq_psi
	      sum_v=sum_v+psi(n_x,n_y,n_z)**2*v(n_x,n_y,n_z)  
            r=dsqrt(space_step**2*(n_z**2+n_y**2+n_x**2))
c     ******************************************************************
c     NOTE: used to be r = space_step*r_0
c     return to this, not sure why r_0 is calculated this way
c     Must ask Humphrey
c     ******************************************************************
	      sum_r_0=sum_r_0+psi(n_x,n_y,n_z)**2.0*(r_1-r)
	    end do
        end do
	end do
      k_bar=sum_k/sum_psi
  	v_bar=sum_v/sum_psi
	r_0=sum_r_0/sum_psi
      e=k_bar+v_bar
	e_total=8*e-4*v_bar+surface*0.375
	if (mod(n_iter,10) .eq. 0)
     &write(*,100) n_iter,r_0,k_bar,v_bar,e,e_total
100   format(' n_iter',i6,'  r_0',f7.3,'  k',f8.1,'  v',f8.1,'  e',
     &f8.1,'  e_tot',f8.1)        
c     ******************************************************************
c     Renormalize. ***************************************************
c     ******************************************************************
      factor=1/dsqrt(sum_psi)
      do n_x=-n_1,n_1
        do n_y=-n_1,n_1
          do n_z=-n_1,n_1
	      psi(n_x,n_y,n_z)=psi(n_x,n_y,n_z)*factor
	    end do
        end do
	end do
c     ****************************************************************      
c     Save psi. ******************************************************	
c     ****************************************************************      
      if (mod(n_iter,n_save) .eq. 0) then  
        open(3,file='Cube8_psi.dat',status='unknown')
	  do n_x=-n_1,n_1
	    do n_z=-n_1,n_1
		  write(3,140) n_x*space_step,n_z*space_step,
     &	  psi_p(n_x,0,n_z),psi(n_x,0,n_z)
140         format(2f7.1,2es14.6) 
          end do
		write(3,140)
	  end do	   
	  close(3)
      endif 
c     ******************************************************************
c     Decide whether result has converged. *****************************
c     ******************************************************************
	e_array(n_iter)=e
	v_bar_array(n_iter)=v_bar
	r_0_array(n_iter)=r_0
      if (n_iter .lt. 100) goto 20
      n_iter_p=(0.7*n_iter/n_save)*n_save
	delta_e=dabs(e_array(n_iter)/e_array(n_iter_p)-1)
	if (delta_e .gt. 0.001) goto 20
	delta_v_bar=dabs(v_bar_array(n_iter)/v_bar_array(n_iter_p)-1)
	if (delta_v_bar .gt. 0.001) goto 20
	delta_r_0=dabs(r_0_array(n_iter)/r_0_array(n_iter_p)-1)
	if (delta_r_0 .gt. 0.001) goto 20
	return
	goto 20
      return
	end
c     ******************************************************************
c     ******************************************************************
      subroutine potential_energy(r_1,r_0,space_step,n_x,n_y,
     &n_z,v)
c     ******************************************************************      
c     ******************************************************************
      implicit real*8 (a-h,o-z)
      e=480.325d0
      pi=2*dacos(0.0d0)
      x=n_x*space_step
	y=n_y*space_step
	z=n_z*space_step
      r_bar=(r_1-r_0)
      cartesian = r_bar/(3.0**.5)
	v=0.0d0
      x_1=-cartesian
	y_1=-cartesian
	z_1=-cartesian
	v=v+e**2/dsqrt((x-x_1)**2+(y-y_1)**2+(z-z_1)**2)
      x_2=-cartesian
	y_2=cartesian
	z_2=cartesian
	v=v+e**2/dsqrt((x-x_2)**2+(y-y_2)**2+(z-z_2)**2)
      x_3=cartesian
	y_3=-cartesian
	z_3=cartesian
	v=v+e**2/dsqrt((x-x_3)**2+(y-y_3)**2+(z-z_3)**2)
      x_4=-cartesian
	y_4=-cartesian
	z_4=cartesian
	v=v+e**2/dsqrt((x-x_4)**2+(y-y_4)**2+(z-z_4)**2)
      x_5=cartesian
	y_5=cartesian
	z_5=-cartesian
      v=v+e**2/dsqrt((x-x_5)**2+(y-y_5)**2+(z-z_5)**2)
      x_6=-cartesian
	y_6=cartesian
	z_6=-cartesian
      v=v+e**2/dsqrt((x-x_6)**2+(y-y_6)**2+(z-z_6)**2)
      x_7=cartesian
	y_7=-cartesian
	z_7=-cartesian
      v=v+e**2/dsqrt((x-x_7)**2+(y-y_7)**2+(z-z_7)**2)
      return
	end
c     ******************************************************************
c     ******************************************************************
      subroutine volume_surface(r_1,r_2,surface,volume)
c     ******************************************************************
c     ******************************************************************
	implicit real*8 (a-h,o-z)
      real*8 x(0:200,0:200),y(0:200,0:200),z(0:200,0:200),r(0:200,0:200) 
	real*8 theta(0:100),cos_theta(0:100),sin_theta(0:100)
	real*8 phi(0:100),cos_phi(0:100),sin_phi(0:100)
	pi=2*dacos(0.0d0)
c     Set up array for sines and cosines. ****************************** 
      n_theta_max=100
      d_theta=pi/n_theta_max
      n_phi_max=100
	d_phi=2*pi/n_phi_max
	do n_theta=0,n_theta_max
        theta(n_theta)=n_theta*d_theta
	  cos_theta(n_theta)=dcos(theta(n_theta))
        sin_theta(n_theta)=dsin(theta(n_theta))
	end do
      do n_phi=0,n_phi_max
	  phi(n_phi)=n_phi*d_phi
        cos_phi(n_phi)=dcos(phi(n_phi))
	  sin_phi(n_phi)=dsin(phi(n_phi))
	end do
c     Calculate x,y,z and r. ******************************************** 
	do n_theta=0,n_theta_max
        do n_phi=0,n_phi_max
	    r(n_theta,n_phi)=-0.5*r_2+1.5*r_1-
     &	1.5*(r_1-r_2)*(sin_theta(n_theta)**4*(cos_phi(n_phi)**4+
     &    sin_phi(n_phi)**4)+cos_theta(n_theta)**4)
	    x(n_theta,n_phi)=r(n_theta,n_phi)*
     *	sin_theta(n_theta)*dcos(phi(n_phi))
	    y(n_theta,n_phi)=r(n_theta,n_phi)*
     *	sin_theta(n_theta)*dsin(phi(n_phi))
	    z(n_theta,n_phi)=r(n_theta,n_phi)*cos_theta(n_theta)
        end do
      end do
c     ******************************************************************
c     Find volume of bubble. *******************************************
c     ******************************************************************
      volume=0.0d0
 	do n_theta=0,n_theta_max-1
        do n_phi=0,n_phi_max-1
          r_average=(r(n_theta,n_phi)+r(n_theta+1,n_phi)+
     &	r(n_theta,n_phi+1)+r(n_theta+1,n_phi+1))/4.d0
	    volume=volume+(r_average**3/3.d0)*
     &	sin_theta(n_theta)*d_theta*d_phi
        end do
     	end do
c     ******************************************************************
c     Find surface area of bubble. *************************************
c     ******************************************************************
	surface=0.0d0
 	do n_theta=0,n_theta_max-1
        do n_phi=0,n_phi_max-1
c    		Variation with respect to phi.  ******************************
          f_x=0.5*(x(n_theta,n_phi+1)-x(n_theta,n_phi)+
     &    x(n_theta+1,n_phi+1)-x(n_theta+1,n_phi))
          f_y=0.5*(y(n_theta,n_phi+1)-y(n_theta,n_phi)+
     &    y(n_theta+1,n_phi+1)-y(n_theta+1,n_phi))
          f_z=0.5*(z(n_theta,n_phi+1)-z(n_theta,n_phi)+
     &    z(n_theta+1,n_phi+1)-z(n_theta+1,n_phi))
c    		Variation with respect to theta.  
          g_x=0.5*(x(n_theta+1,n_phi)-x(n_theta,n_phi)+
     &    x(n_theta+1,n_phi+1)-x(n_theta,n_phi+1))
          g_y=0.5*(y(n_theta+1,n_phi)-y(n_theta,n_phi)+
     &    y(n_theta+1,n_phi+1)-y(n_theta,n_phi+1))
          g_z=0.5*(z(n_theta+1,n_phi)-z(n_theta,n_phi)+
     &    z(n_theta+1,n_phi+1)-z(n_theta,n_phi+1))
c         Calculate cross product. *************************************
          a_x=f_y*g_z-f_z*g_y		
          a_y=f_z*g_x-f_x*g_z		
          a_z=f_x*g_y-f_y*g_x
          d_surface=dsqrt(a_x**2+a_y**2+a_z**2)				
          surface=surface+d_surface
        end do
     	end do
	return 
	end
