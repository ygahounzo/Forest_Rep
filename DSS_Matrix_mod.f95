MODULE DSS_Mass_matrix_mod

USE Lagrange_mod

USE Weight_Leg_Lobat_mod


IMPLICIT NONE

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

! mass matrix at each element

SUBROUTINE Mass_matrix_element(N, Q, Me)

!=================================================================
! N	: order of polynomials
! Q 	: number of LGL integration points
! Me	: mass matrix of an element 
!==================================================================

IMPLICIT NONE

INTEGER						:: i, j, k
INTEGER						:: N, Q
DOUBLE PRECISION, DIMENSION(0:Q) 		:: w
DOUBLE PRECISION , DIMENSION(0:N,0:N)		:: Me
DOUBLE PRECISION, DIMENSION(0:Q) 		:: xi
DOUBLE PRECISION, DIMENSION(0:Q-1) 		:: xroots
DOUBLE PRECISION				:: outp1, outp2

CALL weight_LobattoPoints(Q, xroots, xi, w)


do i = 0, N
   
   do j = 0, N
      do k = 0, Q

	 ! call of the basis function 
	 CALL Lagrange_deriv(N, i, xi(k1), xroots, outp1, outp2)
         p = outp1
	 CALL Lagrange_deriv(N, i, xi(k1), xroots, outp1, outp2)
	 q1 = outp1
         
	 Me(i,j) = Me(i,j) + w(k)*p*q1

      enddo
   enddo
enddo

Me = (1.0D0/2.0)*Me

END SUBROUTINE Mass_matrix_element


! DSS matrix


SUBROUTINE DSS_Mass_matrix(N, Ne, Q, Me, M)

!=================================================================
! N	: order of polynomials
! Ne	: number of elements
! Q 	: number of LGL integration points
! Me	: mass matrix of an element 
! M	: global mass matrix
!=================================================================

IMPLICIT NONE

INTEGER						:: i, j, k, N, Q, Ne, u, v

DOUBLE PRECISION , DIMENSION(0:N*Ne,0:N*Ne) 	:: M
DOUBLE PRECISION , DIMENSION(0:N,0:N)		:: Me
INTEGER , DIMENSION(0:N) 			:: intm
DOUBLE PRECISION, DIMENSION(0:Q) 		:: xi

CALL Mass_matrix_element(N, Q, Me)

do e = 1, Ne
        do i = 0, N

            ! call of array intma

            call intma(e, N, intm)
	    u = intm(i)
            
            do j = 0, N

                ! call of array intma

                call intma(e, N, intm)

	        v = intm(j)
                
                M(u,v) = M(u,v) + Me(i,j)

            enddo
        enddo
enddo

END SUBROUTINE DSS_Mass_matrix

MODULE DSS_Mass_matrix_mod
