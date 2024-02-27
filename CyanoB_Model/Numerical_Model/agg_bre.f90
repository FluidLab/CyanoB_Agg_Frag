!Research project: Fluid dynamics of aquatic microorganisms
!Code title: Simulation of Microcystis colony formation under turbulent shear
!Author: Yuri Sinzato - PhD Student - University of Amsterdal
!Date (start): 03-08-2022
!Date (current version): 19-06-2023 
!!!!!!!!!!!!!!!!! 
!	Description:
!!!!!!!!!!!!!!!!! 
! Routine to simulate the dynamics of Microcystis colony formation under:
! - Uniform homogeneous turbulent shear flow carachterized by tubulent dissipation rate epsilon and kinematic viscosity nu
! - No cell division
! - Aggregation is also turbulent driven, with a stickness coeficient alpha. Model described in "Burd and Jackson, 2009, Particle Aggregation"
! - Break-up driven by erosion - Double population model, adapted from Babler 2008
! - Colony size distribution described by frequency as a function of radius n(r) and mass fraction (mass-weighted distribution n_m(r))
! - Solution by classes method (discritized distribution function)
! - Two population kinds: division-formed colonies, with large critical size, and aggregation-formed colonies, with small critical size

! Physical imput paramenters
! m1 : Mass of a single cell
! r1 : Radius of a single cell
! epsilon : turbulent energy dissipation rate (m2/s3)
! nu : kinematic viscosity
! d_f : 3-D Fractal dimension defined by (m'/m1) = (r'/r1)^d_f
! tau_h = sqrt(nu/epsilon) : reference hydrodynamic time
! alpha : Stickiness
! phi_o : Initial equivalent volume fraction = (Total number of cells*Single cell volume/Suspension volume)
! q : Breakup law exponent
! m_cr : Breakup law critical mass for the average shear rate
! t_exp = Minimum simulation time in dimensional units (seconds)

! Variables in dimensional form, denoted by a ' whenever there is a non-dimensional form:
! csd_f' : frequency weighted colony size distribution defined by: csd_f'(r')*dr' = number of colonies with size between [r;r+dr] per volume of suspension
! csd_m' : mass weighted colony size distribution defined by: csd_m'(r')*ds' = m(r')*csd_f'(r')*dr' = mass of colonies with size between [r;r+dr] per volume of suspension
! r' : radius of a colony 
! r : relative radius of a colony (r'/r1)
! m(r)= r^df: relative mass of a colony with radius r
! t_d : dimensional time 
! t : non-dimensional time t = t_d*sqrt(epsilon/nu)
! K_b : Brakage rate function

! Populations
!csd_m(1,1:N_r) : Aggregation-formed colonies. Can be broken and aggregated
!csd_m(2,1:N_r) : Division-formed colonies. Can be broken, but not aggregated. Fragments become kind 1
  
! Numerical paramenters (all in non-dimensional form)
! N_a : Number of size segments that are aggregativelly active
! N_r : Total number of size segments 
! h_r : Size segment length
! N_tmax : Maximum number of time steps
! Dt : Time step
! Dt_i : Initial time step
! Dt_l : Lower limit of time step
! Dt_u : Upper limit of time step
! tol_csd : Tolerance in csd solver
 
!Lapack path C:\Users\YZS10\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc\LocalState\rootfs\usr\lib\x86_64-linux-gnu\lapack
!Blas path C:\Users\YZS10\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc\LocalState\rootfs\usr\lib\x86_64-linux-gnu\blas


!!!!!!!!!!!!!!!!! 
!	Module to read parameters file and define global variables 
!!!!!!!!!!!!!!!!!!
module input_parameters
	implicit none
	save
		integer::N_a,N_tmax,N_r
		double precision::h_r,t_exp,Dt_i,Dt_l,Dt_u,tol_csd,mu,alpha,phi_o,d_f,epsilon,nu
		double precision,dimension(1:2)::ep_s,q	

	contains

		subroutine read_parameters()
			implicit none
			character(len=12)::nome_var
			CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A,F12.6,A,F12.6,A)"   !! Formato do nome dos arquivos com os perfis

		!! Lê o arquivo de parâmetros
		OPEN(UNIT=10,FILE = 'Numerical_Model/Temporary/parameters.txt')
		read(10,*) nome_var,N_a
		read(10,*) nome_var,N_r
		read(10,*) nome_var,h_r
		read(10,*) nome_var,t_exp
		read(10,*) nome_var,N_tmax
		read(10,*) nome_var,Dt_i
		read(10,*) nome_var,Dt_l
		read(10,*) nome_var,Dt_u
		read(10,*) nome_var,tol_csd
		read(10,*) nome_var,mu
		read(10,*) nome_var,alpha
		read(10,*) nome_var,phi_o
		read(10,*) nome_var,ep_s(1)
		read(10,*) nome_var,ep_s(2)
		read(10,*) nome_var,d_f
		read(10,*) nome_var,q(1)
		read(10,*) nome_var,q(2)
		read(10,*) nome_var,epsilon
		read(10,*) nome_var,nu
		close(10)	

		end subroutine


	end module input_parameters

module globals
	use input_parameters
	implicit none
	save	

	double precision,allocatable::r(:),m(:),K_b(:,:),beta(:,:)

	end module globals

program agg_bre
	use input_parameters
	use globals
	implicit none


		!! All variable are non-dimensional
		
		!!!!!!!!!!!!!!!!! 
		!	Declaring variables:
		!!!!!!!!!!!!!!!!! 
		double precision,allocatable::csd_f(:,:),csd_m(:,:)
		double precision,allocatable::csd_f_t(:,:),csd_m_t(:,:),t_csd(:),r_ave_t(:)
		double precision::t,t_d,Dt,Dt_real,cell_con,col_con
		integer::i,j,ifault,debug

		call read_parameters()
		allocate(r(1:N_r),m(1:N_r),csd_f(1:2,1:N_r),csd_m(1:2,1:N_r),K_b(1:2,1:N_r),&
		&beta(1:N_a,1:N_a),csd_f_t(1:N_r,0:100),csd_m_t(1:N_r,0:100),t_csd(0:100),r_ave_t(0:100))
		OPEN(UNIT=20,FILE = 'Numerical_Model/Temporary/r_data.txt')
		OPEN(UNIT=40,FILE = 'Numerical_Model/Temporary/csd_data.txt')
		
		!!!!!!!!!!!!!!!!! 
		!	Coordinates vector
		!!!!!!!!!!!!!!!!! 
		
		do i=1,N_r
		 r(i) = (i-0.5d+0)*h_r
		 m(i) = r(i)**d_f
		end do	

		!!!!!!!!!!!!!!!!! 
		!	Brekage rates - Calculate K_b values at the center of each segments. K_b Function is defined by calcK_b
		!!!!!!!!!!!!!!!!! 

		do i=1,N_r
		 K_b(1,i) = calcK_b(1,i)
		 K_b(2,i) = calcK_b(2,i)
		end do	
	
		!!!!!!!!!!!!!!!!! 
		!	Collision rates - Calculate beta values at the center of each segments. beta Function is defined by calcbeta
		!!!!!!!!!!!!!!!!! 

		do i=1,N_a
			do j=1,N_a
		 	 beta(i,j) = calcbeta(r(i),r(j))
			end do	
		end do	
		
		!!!!!!!!!!!!!!!!! 
		!	Initial distribution
		!!!!!!!!!!!!!!!!! 
		
		t = 0
		t_d = t*(nu/epsilon)**0.5
		Dt = Dt_i
		call read_csd(csd_m)
		! uniform distribution
		do i=1,N_r
			csd_f(1,i) = csd_m(1,i)/m(i)
			csd_f(2,i) = csd_m(2,i)/m(i)
		end do	
		cell_con = calcconc(csd_m)
		col_con = calccol_con(csd_f)
		! Write on file
		
		write(20,"(A15,T20,A15,T40,A15,T60,A15,T80,A15,T100,A15,T120,A15,T140,A15,T160,A15,T180,A15,T200,A15,T220,A15&
		&,T240,A15,T260,A15,T280)",advance='yes')'#(1)t(s)','(2)r10',&
		&'(3)rv1v','(4)rm1m','(5)col_con','(6)cell_con','(7)bla-bla','(8)frac_kind1','(9)kind1_r_per50','(10)kind1_r_per25',&
		&'(11)kind1_r_per75','(12)kind2_r_per50','(13)kind2_r_per25','(14)kind2_r_per75'
		write(20,"(F16.6,T20,F16.10,T40,F16.10,T60,F16.10,T80,F16.10,T100,F16.10,T120,F16.10,T140,F16.10,T160,F16.10,T180&
		&,F16.10,T200,F16.10,T220,F16.10,T240,F16.10,T260,F16.10,T280)",advance='yes')&
		&t_d,calcr_ave_f(csd_f),calcr_ave_m(csd_m),calcr_ave_m(csd_m),col_con,cell_con,0.0d+0&
		&,calc_frac(csd_m),calcr_per(csd_m,0.50d+0,1),calcr_per(csd_m,0.25d+0,1),calcr_per(csd_m,0.75d+0,1)&
		&,calcr_per(csd_m,0.50d+0,2),calcr_per(csd_m,0.25d+0,2),calcr_per(csd_m,0.75d+0,2)
		write(40,"(A50)",advance='yes')'#Normalized size distribution per time organized in blocks' 
		write(40,"(A15,T20,A15,T40,A15,T60,A15,T80,A19,T100,A19,T120)",advance='yes')'#(1)Block number','(2)t(s)','(3)r'&
		&,'(4)csd_f(norm)','(5)csd_m(kind1-norm)','(6)csd_m(kind2-norm)' 
		do i=1,N_r
			write(40,"(I15,T20,F16.6,T40,F16.10,T60,F16.10,T80,F16.10,T100,F16.10,T120)",advance='yes')&
			&0,t_d,r(i),csd_f(1,i)/col_con,csd_m(1,i)/cell_con,csd_m(2,i)/cell_con
		end do		
		write(40,"(A2,T4)",advance='yes')' '		
		write(40,"(A2,T4)",advance='yes')' '

		!! Clear terminal
		call execute_command_line ("clear")
		print*,'Running...:'
		
		!!!!!!!!!!!!!!!!! 
		!	Time iteration
		!!!!!!!!!!!!!!!!! 
		
		
		do i=1,N_tmax
			! Call Runge-Kutta-Fehlberg to advance the Colony size distribution by a time step Dt_real
			call rkf45(csd_m,Dt,Dt_real)
			! Update current non-dimensional time
			t = t+Dt_real
			! Calculate dimensional time in seconds
			t_d = t*(nu/epsilon)**0.5
			! Calculate frequency weighted distribution
			do j=1,N_r
				csd_f(1,j) = csd_m(1,j)/m(j)
				csd_f(2,j) = csd_m(2,j)/m(j)
			end do
			! Calculate concentrations
			cell_con = calcconc(csd_m)
			col_con = calccol_con(csd_f)			
			! Write results in txt. file
			write(20,"(F16.6,T20,F16.10,T40,F16.10,T60,F16.10,T80,F16.10,T100,F16.10,T120,F16.10,T140,F16.10,T160,F16.10,T180&
			&,F16.10,T200,F16.10,T220,F16.10,T240,F16.10,T260,F16.10,T280)",advance='yes')&
			&t_d,calcr_ave_f(csd_f),calcr_ave_m(csd_m),calcr_ave_m(csd_m),col_con,cell_con,0.0d+0&
			&,calc_frac(csd_m),calcr_per(csd_m,0.50d+0,1),calcr_per(csd_m,0.25d+0,1),calcr_per(csd_m,0.75d+0,1)&
			&,calcr_per(csd_m,0.50d+0,2),calcr_per(csd_m,0.25d+0,2),calcr_per(csd_m,0.75d+0,2)
			if (floor(t_d/300) /= floor((t_d-Dt_real*(nu/epsilon)**0.5)/300)) then    ! Print every 300 s
				do j=1,N_r
					write(40,"(I15,T20,F16.6,T40,F16.10,T60,F16.10,T80,F16.10,T100,F16.10,T120)",advance='yes')&
					&i,t_d,r(j),csd_f(1,j)/col_con,csd_m(1,j)/cell_con,csd_m(2,j)/cell_con
				end do		
				write(40,"(A2,T4)",advance='yes')' '		
					write(40,"(A2,T4)",advance='yes')' '
			end if
			! Print simulation status in terminal
			call execute_command_line ("clear")
			print*,'Running... :','i= ',100*t_d/t_exp,'%  ;  ','t = ',t_d, 's'
			if (t_d>t_exp) then
				print*,'Reached experimental time'
				exit
			end if 
			
		end do		
			
		!!!!!!!!!!!!!!!!! 
		!	!Finish results
		!!!!!!!!!!!!!!!!!	
		close(20)
		close(40)
		call execute_command_line ("clear")
		print*,'Running...  Finished!'		
		deallocate(r,m,csd_f,csd_m,K_b,beta)


	CONTAINS

		!!!!!!!!!!!!!!!!! 
		!	Breakage rate function - Bäbler 2008 - Gaussian turbulent distribution
		!!!!!!!!!!!!!!!!! 

		function calcK_b(j,i)
		use input_parameters
		use globals
			implicit none
			double precision::calcK_b,l_cr,c1,m_sc
			integer::i,j
			
			c1 = 1
			
			! i is the size bin
			! j is the population kind

			l_cr = max((epsilon/ep_s(j))**(-q(j)/d_f),1.0d+0)
			m_sc = (7.5d+0)**0.5*(r(i)/l_cr)**(-d_f/2/q(j))
			if (r(i) .LE. (1.0d+0+0.5*h_r)) then ! No breakup for single cells
				calcK_b = 0
!			else if (r(i) .LE. 2.0d+0*l_cr) then   ! Breakup function up to particles twice the critical size
!				calcK_b = c1*sqrt(2.0/15/3.14159)*exp(-m_sc**2)/erf(m_sc)
			else    ! Regularization for particles larger than twice the critical size
!				m_sc = (7.5d+0)**0.5*(2.0)**(-d_f/2/q)
				calcK_b = c1*sqrt(2.0/15/3.14159)*exp(-m_sc**2)/erf(m_sc)
			end if
			return
		end function calcK_b

		!!!!!!!!!!!!!!!!! 
		!	Collision rate function for simple turbulent shear
		!!!!!!!!!!!!!!!!! 

		function calcbeta(ra,rb)
		use input_parameters
		use globals
			implicit none
			double precision::calcbeta,ra,rb
			
			calcbeta=0.163*(ra +rb)**3
			return 
		end function calcbeta

		!!!!!!!!!!!!!!!!! 
		!	Total mass concentration
		!!!!!!!!!!!!!!!!! 

		function calcconc(csd_m)
		use input_parameters
		use globals
			implicit none
			double precision::calcconc
			double precision,dimension(1:2,1:N_r)::csd_m
			integer::i
			
			calcconc = 0
			do i=1,N_r
		 	 calcconc = calcconc + csd_m(1,i)*h_r + csd_m(2,i)*h_r
			end do
			return
		end function calcconc

		!!!!!!!!!!!!!!!!! 
		!	Total concentration of colonies
		!!!!!!!!!!!!!!!!! 

		function calccol_con(csd_f)
		use input_parameters
		use globals
			implicit none
			double precision::calccol_con
			double precision,dimension(1:2,1:N_r)::csd_f
			integer::i
			
			calccol_con = 0
			do i=1,N_r
		 	 calccol_con = calccol_con + csd_f(1,i)*h_r + csd_f(2,i)*h_r
			end do
			return
		end function calccol_con


		!!!!!!!!!!!!!!!!! 
		!	Frequence-weighted average colony radius
		!!!!!!!!!!!!!!!!! 

		function calcr_ave_f(csd_f)
		use input_parameters
		use globals
			implicit none
			double precision::calcr_ave_f,mean,norm
			double precision,dimension(1:2,1:N_r)::csd_f
			integer::i
			
			mean = 0
			norm = 0
			do i=1,N_r
		 	  mean = mean + csd_f(1,i)*h_r*r(i) + csd_f(2,i)*h_r*r(i)
		 	  norm = norm + csd_f(1,i)*h_r + csd_f(2,i)*h_r 
			end do
			calcr_ave_f = mean/norm
			return
		end function calcr_ave_f


		!!!!!!!!!!!!!!!!! 
		!	Mass-weighted average colony radius
		!!!!!!!!!!!!!!!!! 

		function calcr_ave_m(csd_m)
		use input_parameters
		use globals
			implicit none
			double precision::calcr_ave_m,mean,norm
			double precision,dimension(1:2,1:N_r)::csd_m
			integer::i
			
			mean = 0
			norm = 0
			do i=1,N_r
		 	  mean = mean + csd_m(1,i)*h_r*r(i) +csd_m(2,i)*h_r*r(i)
		 	  norm = norm + csd_m(1,i)*h_r + csd_m(2,i)*h_r
			end do
			calcr_ave_m = mean/norm
			return
		end function calcr_ave_m

		!!!!!!!!!!!!!!!!! 
		!	Mass-fraction of population kind 1
		!!!!!!!!!!!!!!!!! 

		function calc_frac(csd_m)
		use input_parameters
		use globals
			implicit none
			double precision::calc_frac,total,kind1
			double precision,dimension(1:2,1:N_r)::csd_m
			integer::i
			
			total = 0
			kind1 = 0
			do i=1,N_r
		 	  total = total + csd_m(1,i)*h_r + csd_m(2,i)*h_r
		 	  kind1 = kind1 + csd_m(1,i)*h_r
			end do
			calc_frac = kind1/total
			return
		end function calc_frac

		!!!!!!!!!!!!!!!!! 
		!	Radius coresponding to percentile p of Mass-weighted size distribution for each kind
		!!!!!!!!!!!!!!!!! 

		function calcr_per(csd_m,p,j)
		use input_parameters
		use globals
			implicit none
			double precision,dimension(1:2,1:N_r)::csd_m
			double precision,dimension(1:N_r)::norm_csd
			double precision,dimension(0:N_r)::cum_csd
			double precision::calcr_per,p,sum
			integer::i,j

			sum = 0
			! Normalize
			do i=1,N_r
			  norm_csd(i) = csd_m(j,i) + 1.0d-10
			  sum = sum + norm_csd(i)	
			 end do		
			do i=1,N_r
			  norm_csd(i) = norm_csd(i)/sum
			  end do
			! Cumulative sum
			! Find percentile using linear interpolation
			cum_csd(0) = 0
			do i=1,N_r
			  cum_csd(i) = cum_csd(i-1) + norm_csd(i)
			  if (cum_csd(i) > p) then
			    calcr_per = ((i-1)*h_r*(cum_csd(i)-p)+i*h_r*(p - cum_csd(i-1)))/(cum_csd(i)-cum_csd(i-1))
			    exit
			    end if
			  end do
			return
		end function calcr_per

		!!!!!!!!!!!!!!!!! 
		!	Read initial distribution from txt file
		!!!!!!!!!!!!!!!!! 
		subroutine read_csd(csd_m)
		use input_parameters
		use globals
			implicit none
			double precision,dimension(1:2,1:N_r),intent(out)::csd_m
			CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(I15,T20,F16.6,T40,F16.10,T60,F16.10,T80,F16.10,T100,F16.10,T120)"   !! Formato do nome dos arquivos com os perfis
			character(len=200)::pretext
			real::t_dum,r_dum,csdf_dum
			integer::i,i_blo

	
			OPEN(UNIT=30,FILE = 'Numerical_Model/Temporary/csd_initial.txt')
			read(30,"(A)") pretext
			read(30,"(A)") pretext
			do i=1,N_r
				read(30,*) i_blo,t_dum,r_dum,csdf_dum,csd_m(1,i),csd_m(2,i)
			end do  		  	
			close(30)

		end subroutine read_csd		


		!!!!!!!!!!!!!!!!! 
		!	Partial time derivative of csd_m (Mass-weighted colony size distribution)
		!!!!!!!!!!!!!!!!! 

		subroutine derivative(dndt,csd_m)
		use input_parameters
		use globals
			implicit none
			double precision,dimension(1:2,1:N_r),intent(in)::csd_m
			double precision,dimension(1:2,1:N_r),intent(out)::dndt
			double precision::n_col,m_col,m_bre,r_com,r_des
			integer::i,j,debug,i_des,i_single
	
			!!! Initialize derivative
			do i=1,N_r
				dndt(1,i) = 0
				dndt(2,i) = 0
			end do		
						

			!!! Aggregation - Calculate aggregative terms by performing collision between all segments
			!!! Only Kind 1 can aggregate
			do i=1,N_a
				do j=1,N_a
				   ! Number of colisions
				   n_col = 0.5*phi_o*6.0*alpha*beta(i,j)*csd_m(1,i)*csd_m(1,j)/m(i)/m(j)/3.14159
				   m_col = n_col*(m(i)+m(j))
				   ! Losing colonies
				   dndt(1,i) = dndt(1,i)- n_col*m(i)
				   dndt(1,j) = dndt(1,j)- n_col*m(j)
				   ! Destination colonies
				   r_des = (m(i)+m(j))**(1.0d+0/d_f)	
				   i_des = nint(r_des/h_r)	
				   dndt(1,i_des) = dndt(1,i_des) + (r(i_des+1) - r_des)*m_col/h_r
				   dndt(1,i_des+1) = dndt(1,i_des+1) + (r_des - r(i_des))*m_col/h_r   

				end do		
			end do		
			
			!!! Breakage - Calculate breakage terms by performing break of between all segments
			!!! Fragment distribution assumed to be of a single cell + the remaining colony
			!! Kind 1 - fragments remain in kind 1
			i_single = ceiling(1.0d+0/h_r)
			if (((1.0d+0 - r(i_single)).GE. 0.0d+0).OR.(i_single .eq. 1)) then
				do i=(i_single+2),N_r
					! Losing colonies	
					m_bre = K_b(1,i)*csd_m(1,i)
					dndt(1,i) = dndt(1,i) - m_bre
					!Fragments - Model 1 - Erosion
					r_com = (m(i)-1)**(1.0d+0/d_f)
					!dndt(1,i_single) = dndt(1,i_single) + (r(i_single+1) - 1.0d+0)*m_bre/m(i)/h_r
					!dndt(1,i_single+1) = dndt(1,i_single+1) + (1.0d+0-r(i_single))*m_bre/m(i)/h_r
					!dndt(1,i) = dndt(1,i) + (r_com-r(i-1))*(m(i)-1.0d+0)*m_bre/m(i)/h_r
					!dndt(1,i-1) = dndt(1,i-1) + (r(i)-r_com)*(m(i)-1.0d+0)*m_bre/m(i)/h_r
					!Fragments - Model 2 - Two equal fragments 
					r_com = (m(i)*0.5)**(1.0d+0/d_f)
					i_des = nint(r_com/h_r)
					dndt(1,i_des) = dndt(1,i_des) + (r(i_des+1) - r_com)*m_bre/h_r
				    dndt(1,i_des+1) = dndt(1,i_des+1) + (r_com - r(i_des))*m_bre/h_r 
				end do					
			else
				do i=(i_single+1),N_r
					! Losing colonies	
					m_bre = K_b(1,i)*csd_m(1,i)
					dndt(1,i) = dndt(1,i) - m_bre
					!Fragments - Model 1 - Erosion 
					r_com = (m(i)-1)**(1.0d+0/d_f)
					!dndt(1,i_single-1) = dndt(1,i_single-1) + (r(i_single) - 1.0d+0)*m_bre/m(i)/h_r
					!dndt(1,i_single) = dndt(1,i_single) + (1.0d+0-r(i_single-1))*m_bre/m(i)/h_r
					!dndt(1,i) = dndt(1,i) + (r_com-r(i-1))*(m(i)-1.0d+0)*m_bre/m(i)/h_r
					!dndt(1,i-1) = dndt(1,i-1) + (r(i)-r_com)*(m(i)-1.0d+0)*m_bre/m(i)/h_r
					!Fragments - Model 2 - Two equal fragments 
					r_com = (m(i)*0.5)**(1.0d+0/d_f)
					i_des = nint(r_com/h_r)
					dndt(1,i_des) = dndt(1,i_des) + (r(i_des+1) - r_com)*m_bre/h_r
				    dndt(1,i_des+1) = dndt(1,i_des+1) + (r_com - r(i_des))*m_bre/h_r   
				end do
			end if
			!! Kind 2 - single cell fragments become kind 1'
			i_single = ceiling(1.0d+0/h_r)
			if (((1.0d+0 - r(i_single)).GE. 0.0d+0).OR.(i_single .eq. 1)) then
				do i=(i_single+2),N_r
					! Losing colonies	
					m_bre = K_b(2,i)*csd_m(2,i)
					dndt(2,i) = dndt(2,i) - m_bre
					!Fragments
					r_com = (m(i)-1)**(1.0d+0/d_f)
					dndt(1,i_single) = dndt(1,i_single) + (r(i_single+1) - 1.0d+0)*m_bre/m(i)/h_r
					dndt(1,i_single+1) = dndt(1,i_single+1) + (1.0d+0-r(i_single))*m_bre/m(i)/h_r
					dndt(2,i) = dndt(2,i) + (r_com-r(i-1))*(m(i)-1.0d+0)*m_bre/m(i)/h_r
					dndt(2,i-1) = dndt(2,i-1) + (r(i)-r_com)*(m(i)-1.0d+0)*m_bre/m(i)/h_r
				end do					
			else
				do i=(i_single+1),N_r
					! Losing colonies	
					m_bre = K_b(2,i)*csd_m(2,i)
					dndt(2,i) = dndt(2,i) - m_bre
					!Fragments
					r_com = (m(i)-1)**(1.0d+0/d_f)
					dndt(1,i_single-1) = dndt(1,i_single-1) + (r(i_single) - 1.0d+0)*m_bre/m(i)/h_r
					dndt(1,i_single) = dndt(1,i_single) + (1.0d+0-r(i_single-1))*m_bre/m(i)/h_r
					dndt(2,i) = dndt(2,i) + (r_com-r(i-1))*(m(i)-1.0d+0)*m_bre/m(i)/h_r
					dndt(2,i-1) = dndt(2,i-1) + (r(i)-r_com)*(m(i)-1.0d+0)*m_bre/m(i)/h_r
				end do
			end if
			! !! Kind 1 - fragments remain in kind 1
			! i_single = nint(1.0d+0/h_r)
			! do i=(i_single+2),N_r
			! 	! Losing colonies	
			! 	m_bre = K_b(1,i)*csd_m(1,i)
			! 	dndt(1,i) = dndt(1,i) - m_bre
			! 	!Fragments
			! 	r_com = (m(i)-1)**(1.0d+0/d_f)
			! 	dndt(1,i_single) = dndt(1,i_single) + (r(i_single+1) - 1.0d+0)*m_bre/m(i)/h_r
			! 	dndt(1,i_single+1) = dndt(1,i_single+1) + (1.0d+0-r(i_single))*m_bre/m(i)/h_r
			! 	dndt(1,i) = dndt(1,i) + (r_com-r(i-1))*(m(i)-1.0d+0)*m_bre/m(i)/h_r
			! 	dndt(1,i-1) = dndt(1,i-1) + (r(i)-r_com)*(m(i)-1.0d+0)*m_bre/m(i)/h_r
			! end do	
			! !! Kind 2 - single cell fragments become kind 1'
			! i_single = nint(1.0d+0/h_r)
			! do i=(i_single+2),N_r
			! 	! Losing colonies	
			! 	m_bre = K_b(2,i)*csd_m(2,i)
			! 	dndt(2,i) = dndt(2,i) - m_bre
			! 	!Fragments
			! 	r_com = (m(i)-1)**(1.0d+0/d_f)
			! 	dndt(1,i_single) = dndt(1,i_single) + (r(i_single+1) - 1.0d+0)*m_bre/m(i)/h_r
			! 	dndt(1,i_single+1) = dndt(1,i_single+1) + (1.0d+0-r(i_single))*m_bre/m(i)/h_r
			! 	dndt(2,i) = dndt(2,i) + (r_com-r(i-1))*(m(i)-1.0d+0)*m_bre/m(i)/h_r
			! 	dndt(2,i-1) = dndt(2,i-1) + (r(i)-r_com)*(m(i)-1.0d+0)*m_bre/m(i)/h_r
			! end do					

		end subroutine derivative		

		!!!!!!!!!!!!!!!!! 
		!	Runge-Kutta-Fehlberg method
		!!!!!!!!!!!!!!!!! 
		subroutine rkf45(csd_m,Dt,Dt_real)
		use input_parameters
		use globals
			implicit none
			double precision,dimension(1:2,1:N_r),intent(inout)::csd_m
			double precision,intent(inout)::Dt
			double precision,intent(out)::Dt_real
			double precision,dimension(1:2,1:N_r)::k1,k2,k3,k4,k5,k6,n_y,n_z
			double precision::n_col,m_col,m_bre,Dn_max,s,Dt_no_conv
			integer::i,debug
	
  		  ! Advance in colony size distribution based on previous distribution na for a given time step Dt
    		  ! Runge-Kutta-Fehlberg method
    		  ! Adpative time step
    			s = 0
    			i = 0
			do i=1,10
				call derivative(k1,csd_m)
				call derivative(k2,csd_m+Dt*k1/4)
				call derivative(k3,csd_m+Dt*(3 * k1 / 32 + 9 * k2 / 32))
				call derivative(k4,csd_m+Dt*(1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197))
				call derivative(k5,csd_m+Dt*(439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104))
				call derivative(k6,csd_m+Dt*(- 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40))
				n_y =  Dt*(25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4101 - k5 / 5)
     	  			n_z = Dt*(16 * k1 / 135 + 6656 * k3 / 12825 + 28561 * k4 / 56430 - 9 * k5 / 50 + 2 * k6 / 55)
      	 			Dn_max = maxval(abs(n_z - n_y))
    	  			s = (tol_csd*calcconc(csd_m)/ Dn_max) ** 0.2
				Dt_real = Dt
				Dt_no_conv = Dt
				Dt = max(Dt_l, min(Dt_u, s * Dt)) 
				if (s>0.9) then
					csd_m = csd_m + n_y
					exit
				end if
				if (i==10) then
					Dt_real = Dt_no_conv
					csd_m = csd_m + n_y
					print*,'RKF45 did not reach converge'
					exit
				end if
				if (isnan(Dt)) then !In case nan value for Dt_l, use lower bound 
					Dt = Dt_l
				end if
			end do
		end subroutine rkf45		



	end program agg_bre
