Program Mass_matrix_DSS

USE Weight_Leg_Lobat_mod
USE Lagrange_mod
USE DSS_Mass_matrix_mod

INTEGER 				:: N = 1, Q = 2, u, v
INTEGER 				:: s,e, i, j, k1, Ne=2, Np=3
DOUBLE PRECISION 			:: p, q1, outp1, outp2
DOUBLE PRECISION, DIMENSION(0:2) 	:: w
DOUBLE PRECISION , DIMENSION(0:2,0:2) 	:: M
DOUBLE PRECISION , DIMENSION(0:1,0:1)	:: Me
INTEGER , DIMENSION(0:1) 		:: intm
DOUBLE PRECISION, DIMENSION(0:2) 	:: xi
DOUBLE PRECISION, DIMENSION(0:1) 	:: xroots

CALL weight_LobattoPoints(Q, xroots, xi, w)

CALL DSS_Mass_matrix(N, Ne, Q, Me, M)


! Matrix print

subroutine print_matrix(b,n,m)

IMPLICIT NONE
INTEGER					:: n,m,i
DOUBLE PRECISION			:: b(n,m) 

do i=1,n 
print '(20f6.2)',b(i,1:m) 
enddo

endsubroutine print_matrix

end program Mass_matrix_DSS
