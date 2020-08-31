Program CG_1D

! Modules
USE Weight_Leg_Lobat_mod
USE Lagrange_mod
USE Matrix_inverse_mod

!IMPLICIT NONE
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Np          		: number of global grid points
! N			: polynomial order
! Ne			: number of elements
! u			: velocity
! M, Me, De		: respectively global mass matrix, element mass matrix, differentiation element matrix
! Miv			: inverse of global mass matrix
! w			: weight values
! R, Re			: global residual vector, element residual vector
! Xg			: globalgrid points
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

INTEGER						:: i, j, k1, k2,I1, I2, J1, n1
DOUBLE PRECISION 				:: u = 2.0
INTEGER, PARAMETER 				:: N = 2, Qd = 3, Ne=50, Np=101
INTEGER 					:: e, k, u1, v, TN
DOUBLE PRECISION 				:: p, q1, outp1, outp2
DOUBLE PRECISION, DIMENSION(0:Qd) 		:: w
DOUBLE PRECISION , DIMENSION(0:Np-1,0:Np-1) 	:: M, Miv
DOUBLE PRECISION , DIMENSION(0:N,0:N)		:: Me, De
INTEGER , DIMENSION(0:N) 			:: intm                 
DOUBLE PRECISION, DIMENSION(0:Qd) 		:: X
DOUBLE PRECISION, DIMENSION(0:N) 		:: xroots, X1, w1
DOUBLE PRECISION, DIMENSION(0:Np-1) 		:: R, RG, q, q0, qn, qh
DOUBLE PRECISION, DIMENSION(0:N) 		:: Re, qe0, fo
DOUBLE PRECISION, DIMENSION(0:Np-1)		:: Xg
DOUBLE PRECISION				:: a = 0.1, dx, Tf = 0.5, ax = -1, bx = 1, dtest, dt, el
DOUBLE PRECISION, DIMENSION(:),	ALLOCATABLE	:: t

OPEN(unit=1,FILE='1Dw.txt',ACTION='READWRITE')
OPEN(unit=2,FILE='ExactDw.txt',ACTION='READWRITE')

CALL weight_LobattoPoints(Qd, xroots, X, w)

! mass matrix at each element
do i = 0, N
   
   do j = 0, N
      do k = 0, Qd

	 	! call of the basis function( Lagrange polynomial) 

	 	call Lagrange_deriv(N, i, X(k), xroots,outp1, outp2)
	 
     		p = outp1
	 	call Lagrange_deriv(N, j, X(k), xroots,outp1, outp2)

	 	q1 = outp1
         
	 	Me(i,j) = Me(i,j) + w(k)*p*q1
      enddo
   enddo
enddo

Me = (1.0D0/2.0)*Me

! DSS matrix

do e = 1, Ne
        do i = 0, N

            ! call of array intma
            call intma(e,N,intm)
	    	u1 = intm(i)
            
            do j = 0, N

                ! call of array intma
                call intma(e,N,intm)
	        	v = intm(j)
                
                M(u1,v) = M(u1,v) + ((bx-ax)/Ne)*Me(i,j)
            enddo
        enddo
enddo

CALL weight_LobattoPoints(N, xroots, X1, w1)

! Differentiation element matrix

do i = 0, N
   
   do j = 0, N
      do k = 0, N

	 ! call of the basis function 

	 call Lagrange_deriv(N, i, X1(k), X1,outp1, outp2)
	 
         p = outp1
	 
	 call Lagrange_deriv(N, j, X1(k), X1,outp1, outp2)

	 q1 = outp2
         
	 De(i,j) = De(i,j) - w1(k)*p*q1
      enddo
   enddo
enddo

! Inverse matrix of global matrix

CALL Matrix_inverse(Np, M, Miv)

! computation of the global grid points
! spatial stuff
do i = 0, Np-1
	Xg(i) = -1.0 + 2.0*i/(Np-1)
enddo

dx = (bx-ax)/(Np-1)

! time stuff
dtest = a*dx/abs(u)
TN = int(Tf/dtest+1.0)

dt = Tf/TN

ALLOCATE(t(0:TN))

do i = 0, TN
t(i) = i*dt
CALL qinit(Xg-u*t(i),q0)
write(2,*) q0
enddo


! Initial condition

CALL qinit(Xg,q0)

! writting of the initial values into a file
CALL qinit_el(ax, el)
q0(0) = el
CALL qinit_el(bx, el)
q0(Np-1) = el

write(1,*) q0

! initialisation of the solution to q0

q = q0

! computation of the solution of 1D wave equations using CG

do n1 = 0, TN
   
   do e = 1, Ne
	call intma(e,N,intm)	
	
	qe0 = q(intm)

	! call of the function f(q)

	CALL fe(u,qe0, fo)

	! computation of element residual vector
	do i = 0, N
	    do j = 0, N
	    	Re(i) = Re(i) -De(i,j)*fo(j)	    
	    enddo 
	enddo
   	fo = 0.0D0   ! reinitialisation of the output of the function fe = u*qe to 0
	
	! computation of global residual vector
   	do i = 0, N

	    I1 = intm(i)
	    
	    R(I1) = R(I1) + Re(i)
	    
    	enddo
	Re = 0.0D0  ! reinitialisation of the element residual vector to 0
   enddo

   !reinitialisation of global residual using inverse mass matrix
   
   do I2 = 0, Np-1
       do J1 = 0, Np-1
       RG(I2) = RG(I2) + Miv(I2,J1)*R(J1)	
       
       enddo
   enddo
   R = 0.0D0  ! reinitialisation of the global residual vector to 0

   ! solution of the wave equation at time n+1

   qn = q + dt*RG
   
   RG = 0.0D0
   q = 0.0D0
   
   ! initial conditions
   CALL qinit_el(ax-u*t(n1+1), el)    	
   qn(0) = el
   CALL qinit_el(bx-u*t(n1+1), el)
   qn(Np-1) = el
   write(1,*) qn
   
   RG = 0.0D0
   q = 0.0D0
   q = qn
enddo

CONTAINS

! Array intma

SUBROUTINE intma(e, N, intmm)

!==================================================
! e	: element
! intmm	: array that contains the values of intma
!==================================================

IMPLICIT NONE

INTEGER                                     :: s,t,r, N
INTEGER , INTENT(IN)                        :: e
INTEGER, DIMENSION(0:N), INTENT(OUT)        :: intmm


t = (e-1)*N
r = N*e

do s = t, r
   intmm(s-t) = s
enddo

END SUBROUTINE intma


SUBROUTINE fe(u,q, fo)

IMPLICIT NONE 		

DOUBLE PRECISION			:: u
DOUBLE PRECISION, DIMENSION(0:N)	:: q
DOUBLE PRECISION, DIMENSION(0:N)	:: fo

fo = u*q

END SUBROUTINE fe

! Initial conditions

SUBROUTINE qinit(X,q0)

IMPLICIT NONE 

DOUBLE PRECISION, DIMENSION(0:Np-1)	:: q0, X

q0 = exp(-(4*X)**2)


END SUBROUTINE qinit

SUBROUTINE qinit_el(X,q0)

IMPLICIT NONE 

DOUBLE PRECISION	:: q0, X

q0 = exp(-(4*X)**2)


END SUBROUTINE qinit_el

! Matrix form printing

subroutine print_matrix(b,n,m)
IMPLICIT NONE

INTEGER			::n,m,i
DOUBLE PRECISION	::b(n,m) 

do i=1,n 
    print '(20f6.2)',b(i,1:m) 
enddo

end subroutine print_matrix

end program CG_1D
