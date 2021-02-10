Module Weight_Leg_Lobat_mod

IMPLICIT NONE

CONTAINS

! Legendre and Lobatto polynomials

SUBROUTINE Lobatto_deriv_poly(P,x,outp1,outp2)

!=========================================================
! P	: Order of the polynomials
! x	: represent variable
! outp1	: P_th order of Lobatto polynomial
! oupt2	: derivative of P_th order of Lobatto polynomial
! oupt3	: P_th order of Legendre polynomial
!=========================================================

IMPLICIT NONE
INTEGER                           	      :: i, P
DOUBLE PRECISION,INTENT(IN)                   :: x
DOUBLE PRECISION,DIMENSION(0:P+1)             :: Li, dLi, ddLi
DOUBLE PRECISION,DIMENSION(0:P+1)             :: Bi,dBi
DOUBLE PRECISION, INTENT(OUT)                 :: outp1, outp2

Li(0) 	= 1.0    		! lengendre polynomial order 0
dLi(0) 	= 0.0   		! first derivative of lengendre polynomial order 0
ddLi(0) = 0.0  			! second derivative of lengendre polynomial order 0
Li(1) 	= x      		! lengendre polynomial order 1
dLi(1) 	= 1.0   		! first derivative of lengendre polynomial order 1
ddLi(1) = 0.0  			! second derivative of lengendre polynomial order 1

do i = 2,P+1

Li(i) 	= ((2.0*i-1.0)/i)*x*Li(i-1) - ((i-1.0)/i)*Li(i-2)       ! lengendre polynomial order i

dLi(i)	= i*Li(i-1) + x*dLi(i-1)                               	! first derivative of lengendre polynomial order i

ddLi(i) = (i+1.0)*dLi(i-1) + x*ddLi(i-1)                      	! second derivative of lengendre polynomial order i

Bi(i) 	= (1.0-x**2)*dLi(i-1)                                   ! lobatto polynomial order i

dBi(i) 	= -2.0*x*dLi(i-1) + (1.0-x**2)*ddLi(i-1)               	! first derivative of lobatto polynomial order i

enddo

outp1 = Bi(P+1)							! P_th order of Lobatto polynomials
outp2 = dBi(P+1)						! derivative of P_th order Lobatto polynomial

END SUBROUTINE Lobatto_deriv_poly


! Legendre polynomials

SUBROUTINE Legendre_deriv_poly(P,x,outp1,outp2)

!=========================================================
! P	: Order of the polynomials
! x	: represent variable
! outp1	: P_th order of Legendre polynomial
! oupt2	: derivative of P_th order of Legendre polynomial
! 
!=========================================================

IMPLICIT NONE
INTEGER                           	      	:: i, P
DOUBLE PRECISION,INTENT(IN)                   	:: x
DOUBLE PRECISION,DIMENSION(0:P)             	:: Li, dLi
DOUBLE PRECISION, INTENT(OUT)                 	:: outp1, outp2

Li(0) 	= 1.0    		! lengendre polynomial order 0
dLi(0) 	= 0.0   		! first derivative of lengendre polynomial order 0
Li(1) 	= x      		! lengendre polynomial order 1
dLi(1) 	= 1.0   		! first derivative of lengendre polynomial order 1

do i = 2,P

Li(i) 	= ((2.0*i-1.0)/i)*x*Li(i-1) - ((i-1.0)/i)*Li(i-2)       ! lengendre polynomial order i

dLi(i)	= i*Li(i-1) + x*dLi(i-1)                               	! first derivative of lengendre polynomial order i


enddo

outp1 = Li(P)							! P_th order of Legendre polynomials
outp2 = dLi(P)							! derivative of P_th order Legendre polynomial

END SUBROUTINE Legendre_deriv_poly

! Lobatto points

SUBROUTINE Lobatto_points(P, xroots, xi)

!=========================================================
! P		: Order of the polynomials
! xi		: Lobatto points
! xroots	: Interpolation points
!=========================================================


IMPLICIT NONE

INTEGER   				:: i, k1, P
DOUBLE PRECISION 			:: outp1, outp2, er = 10.0e-10, Pi = 3.141
DOUBLE PRECISION 			:: xi0, xik, xikk
DOUBLE PRECISION, DIMENSION(0:P) 	:: xi
DOUBLE PRECISION, DIMENSION(0:P-1) 	:: xroots

! Lobatto points

do i = 0, P
   
   ! Chebchev points

   xi0 = cos(((2.0*i+1.0)/(2.0*P+2.0))*Pi)
   
   xik = xi0
   
   do k1 = 0,1000
      
      ! call Lobatto polynomials

      call Lobatto_deriv_poly(P, xik,outp1, outp2)
      
      ! Newton method

      xikk = xik - outp1/outp2

      if (abs(xikk-xik) < er) then
	
	exit
	
      endif
      xik = xikk
   enddo

   xi(i) = xikk
enddo

! Lobatto eta

do i = 0, P-1
   
   ! Chebchev points

   xi0 = cos(((2.0*i+1.0)/(2.0*(P-1)+2.0))*Pi)
   
   xik = xi0
   
   do k1 = 0,1000
      
      ! call Legendre and Lobatto polynomials

      call Lobatto_deriv_poly(P-1, xik,outp1, outp2)
      
      ! Newton method

      xikk = xik - outp1/outp2

      if (abs(xikk-xik) < er) then
	
	exit
	
      endif
      xik = xikk
   enddo

   xroots(i) = xikk
enddo



END SUBROUTINE Lobatto_points

! weight value

SUBROUTINE weight_LobattoPoints(P, xroots, xi, w)

!=========================================================
! P		: Order of the polynomials
! xi		: Lobatto points
! w		: weight values
! xroots	: Interpolation points
!=========================================================


IMPLICIT NONE
INTEGER				:: i, P
DOUBLE PRECISION			:: wi, outp1, outp2
DOUBLE PRECISION, DIMENSION(0:P) 	:: xi, w
DOUBLE PRECISION, DIMENSION(0:P-1) 	:: xroots

CALL Lobatto_points(P, xroots, xi)

do i  = 0,P
   
   ! call of Legendre and Lobatto polynomials and extract of Legendre output

   call Legendre_deriv_poly(P, xi(i),outp1, outp2)
   
   ! weight value using Legendre-Gauss-Lobatto

   wi = 2.0D0/(P*(P+1.0D0)*(outp1)**2)
   w(i) = wi
   
end do
END SUBROUTINE weight_LobattoPoints

END Module Weight_Leg_Lobat_mod
