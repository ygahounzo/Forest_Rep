MODULE Lagrange_mod

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE Lagrange_deriv(N, i, x, xroots, L, dL)
    
    !=======================================================
    ! N		: Order of Lagrange polynomials
    ! i  	: index Lagrange polynomial
    ! x		: sample point
    ! xroots	: Interpolation points
    ! L		: i_th Lagrange polynomial
    ! dL	: derivative of i_th Lagrange polynomial
    !=======================================================


    IMPLICIT NONE
    INTEGER                                     	:: i,j, k
    INTEGER, INTENT(IN)                         	:: N
    DOUBLE PRECISION, DIMENSION(0:N),INTENT(IN)   	:: xroots
    DOUBLE PRECISION, INTENT(IN)                	:: x
    DOUBLE PRECISION, INTENT(OUT)       		:: L, dL   
    DOUBLE PRECISION					:: prod
    
    
    !DO i=0,N
        L =1.0D0
	dL = 0.0
        DO j=0,N
	    prod = 1.0
            IF (i.ne.j) THEN

		L=L*(x-xroots(j))/(xroots(i)-xroots(j))

		do k = 0, N
		   IF ((k.ne.i) .and. (k.ne.j)) THEN 
			prod = prod*(x-xroots(k))/(xroots(i)-xroots(k))
		   END IF                
                END DO
		dL = dL + prod/(xroots(i)-xroots(j))
            END IF
        END DO!j
        
    !END DO!i

    RETURN

    END SUBROUTINE Lagrange_deriv

END MODULE Lagrange_mod

