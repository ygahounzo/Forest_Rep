!include "Lobatto_mod.f95"
!include "Lagrange_mod.f95"

program Mass_matrix_DSS

USE Weight_Leg_Lobat_mod
USE Lagrange_mod

INTEGER 				:: N = 1, Q = 2, u, v
INTEGER 				:: s,e, i, j, k, Ne=2, Np=3
DOUBLE PRECISION 			:: p, q1, outp1, outp2
DOUBLE PRECISION, DIMENSION(0:2) 	:: w
DOUBLE PRECISION , DIMENSION(0:2,0:2) 	:: M
DOUBLE PRECISION , DIMENSION(0:1,0:1)	:: Me
INTEGER , DIMENSION(0:1) 		:: intm
DOUBLE PRECISION, DIMENSION(0:2) 	:: xi
DOUBLE PRECISION, DIMENSION(0:1) 	:: xroots

CALL weight_LobattoPoints(Q, xroots, xi, w)

! mass matrix at each element
do i = 0, N
   
   do j = 0, N
      do k = 0, Q

	 ! call of the basis function 

	 call Lagrange_deriv(N, i, xi(k), xroots,outp1, outp2)
	 
         p = outp1
	 call Lagrange_deriv(N, j, xi(k), xroots,outp1, outp2)

	 q1 = outp1
         
	 Me(i,j) = Me(i,j) + w(k)*p*q1
      enddo
   enddo
enddo

Me = (1.0D0/2.0)*Me
print*,'Me = '
call print_matrix(Me,Q,Q)

! DSS matrix
do e = 1, Ne
        do i = 0, N

            ! call of array intma
            call intma(e,N,intm)
	    u = intm(i)
            
            do j = 0, N

                ! call of array intma
                call intma(e,N,intm)
	        v = intm(j)
                
                M(u,v) = M(u,v) + Me(i,j)
            enddo
        enddo
enddo

print*, 'M = '
call print_matrix(M,Np,Np)

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


! Matrix form printing

subroutine print_matrix(b,n,m)

IMPLICIT NONE
INTEGER					::n,m,i
DOUBLE PRECISION			::b(n,m) 

do i=1,n 
print '(20f6.2)',b(i,1:m) 
enddo

endsubroutine print_matrix


end program Mass_matrix_DSS
