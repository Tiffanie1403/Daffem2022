MODULE Interpolation

  USE Types
  USE Algebra
  USE Algebra2
  USE PublicVar
  USE Locate
  
  IMPLICIT NONE

CONTAINS

  !=====================================================!
  SUBROUTINE InterpolateSolution(X, interpSol, sol, ele)
  !=====================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: X !coordonnées du point a interpolé 
    REAL*8              , INTENT(OUT) :: interpSol ! valeur de la solution interpolé 
    REAL*8, DIMENSION(:), INTENT(IN)  :: sol ! vecteur de solution des noeuds 
    type(Element)       , INTENT(IN)  :: ele ! element sur lequel on se trouve 
    !----------------------------------------------
    REAL*8, DIMENSION(ele%nNodesPerElement) :: phi
    !-----------------------------------------------
    INTEGER :: i, nNodesPerElement, nDim
    !------------------------------------

    nNodesPerElement = ele%nNodesPerElement
    nDim = SIZE(X)    
    
    CALL computeBasisFunction(ele, X, phi)
   ! IF ( dbg ) PRINT*, "---- ", phi !doesn't work
    interpSol = 0.d0
    DO i = 1, nNodesPerElement
       interpSol = interpSol + phi(i)*sol(i)
    END DO
    
  END SUBROUTINE InterpolateSolution
  !=====================================================!

  !===========================================!
  SUBROUTINE ComputeBasisFunction(ele, X, phi)
  !===========================================!
  
  !interpolation linéaire avec les polynomes de lagrange 

    IMPLICIT NONE
    type(Element)       , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:), INTENT(IN)  :: X !coordonée ou on veut la solution interpolé 
    REAL*8, DIMENSION(:), INTENT(OUT) :: phi !
    !-----------------------------------------
    REAL*8, DIMENSION(ele%nNodesPerElement, ele%nNodesPerElement) :: A, Am1
    !---------------------------------------------------------------------------
    REAL*8, DIMENSION(ele%nNodesPerelement) :: Coeff, B
    !-------------------------------------------------------------------
    INTEGER :: nNodesPerElement, i, nDim, j
    !------------------------------------

    nNodesPerElement = ele%nNodesPerElement
    nDim = SIZE(X)
    
    DO i = 1, nNodesPerElement
       A(i,1:nDim) = ele%Coor(i,:)
       A(i,nNodesPerElement) = 1.d0
    END DO
    
    !CALL InverseMatrix_maxPivot(A, Am1)
    CALL Inverse3x3Matrix(A, Am1)
    CALL CheckInversion(A, Am1) 

    DO i = 1, nNodesPerElement
       B = 0.d0
       B(i) = 1.d0
       Coeff = MATMUL(Am1, B)
       phi(i) = 0.d0
       DO j = 1, nDim
          phi(i) = phi(i) + Coeff(j)*X(j)
       END DO
       phi(i) = phi(i) + coeff(nDim+1)
    END DO

!!$    PRINT*, phi
!!$    PRINT*, SUM(phi)
    
  END SUBROUTINE ComputeBasisFunction
  !===========================================!
  
END MODULE Interpolation
