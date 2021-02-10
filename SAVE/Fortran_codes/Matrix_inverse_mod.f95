MODULE Matrix_inverse_mod

IMPLICIT NONE

!====================================================================
!  Computing Inverse matrix Based on LU factorization method AX = B
!====================================================================

CONTAINS

! Inverse of a matrix

SUBROUTINE Matrix_inverse(n, M, Minv)

!====================================================================
! input		: n, M
! n		: dimension of the matrix
! M		: Matrix that we are computing the inverse
! 
! output	: Minv
! Minv		: Inverse of M
! 
!====================================================================

IMPLICIT NONE

INTEGER						:: n, i, j, k, N1
DOUBLE PRECISION, DIMENSION(0:n-1,0:n-1)	:: M, Minv
DOUBLE PRECISION, DIMENSION(0:n-1,0:n-1)	:: L, U
DOUBLE PRECISION, DIMENSION(0:n-1)		:: B, D, X
DOUBLE PRECISION				:: c 

N1 = n-1

! step 0: initialization for matrices L and U and B

L=0.0D0
U=0.0D0
b=0.0D0

! step 1: forward elimination

do k=0, N1-1
   do i=k+1,N1
      c=M(i,k)/M(k,k)
      L(i,k) = c
      do j=k+1,N1
         M(i,j) = M(i,j)-c*M(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 

! L matrix is a matrix of the elimination coefficient, the diagonal elements are 1.0D0

do i=0,N1

  L(i,i) = 1.0D0

end do

! U matrix is the upper triangular part of M

do j=0,N1
  do i=0,j

    U(i,j) = M(i,j)

  end do
end do

! Step 3: compute columns of the inverse matrix Minv

do k=0,N1

  B(k)=1.0
  D(0) = B(0)

! Step 3a: Solve LD=B using the forward substitution

  do i=0,N1
    D(i)=B(i)
    do j=0,i-1
      D(i) = D(i) - L(i,j)*D(j)
    end do
  end do

! Step 3b: Solve UX=D using the back substitution

  X(N1)=D(N1)/U(N1,N1)

  do i = N1-1,0,-1
    X(i) = D(i)
    do j=N1,i+1,-1
      X(i)=X(i)-U(i,j)*X(j)
    end do

    X(i) = X(i)/U(i,i)
  end do

! Step 3c: fill the solutions X(n) into column k of Minv

  do i=0,N1
    Minv(i,k) = X(i)
  end do

  B(k)=0.0D0

end do

END SUBROUTINE Matrix_inverse


END MODULE Matrix_inverse_mod
