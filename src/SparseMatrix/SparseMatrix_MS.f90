MODULE SparseMatrix_MS

  USE Types

  IMPLICIT NONE

CONTAINS

  !===============================================================!
  SUBROUTINE AssembleSparseMatrixMS(Mat, Aloc, elei, elej, Nv, Ne)
  !===============================================================!

  !*************************************************************************
  !   Adds the local matrix to the global matrix discretization for the Ms case
  !
  !  Parameters:
  !
  !    Input,Output, type(SparseMatrix) Mat, Structure for the global matrix
  !    Input, REAL*8, DIMENSION(:,:) Aloc, Local matrix
  !    Input, type(element) elej, element considerated
  !    Input, type(element) elei, element considerated
  !    Input, INTEGER  Nv, Number of nodes (vertices)
  !    Input, INTEGER  Ne, Number of elements
  !
  !*************************************************************************


    IMPLICIT NONE
    type(SparseMatrix)    , INTENT(INOUT) :: Mat
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: Aloc
    type(element)         , INTENT(IN)    :: elei, elej
    INTEGER               , INTENT(IN)    :: Nv, Ne
    !-------------------------------------------------
    INTEGER :: i, j, NU_i, NU_j, Pos_i, Pos_j, k, Pos_k
    INTEGER :: Nei, Nej, NvLoc, Nmax
    INTEGER :: l, m, iG, jG, iL, jL, ivar, jvar
    INTEGER :: nNodesPerElement, nVar
    !---------------------------------------------------
    LOGICAL :: found
    !------------------

    nMax = mat%nMax ; nNodesPerElement = elei%nNodesPerElement ; nVar = Mat%nVar

    DO i = 1, nNodesPerElement
       Pos_i = elei%Vertex(i)
       DO j = 1, nNodesPerElement
          Pos_j = elej%Vertex(j)
          IF ( Pos_i == Pos_j ) THEN
             DO ivar = 1, nVar
                DO jvar = 1, nVar
                   iG = Pos_i + (ivar - 1)*Nv
                   jG = 1 + (jvar - 1)*NMax
                   iL = i + (ivar - 1)*nNodesPerElement
                   jL = j + (jvar - 1)*nNodesPerElement
                   Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                END DO
             END DO
          ELSE
             NvLoc = Mat%Ind(Pos_i,1)
             DO k = 1, NvLoc
                Pos_k = Mat%Ind(Pos_i,1+k)
                IF ( Pos_j == Pos_k ) THEN
                   DO ivar = 1, nVar
                      DO jvar = 1, nVar
                         iG = Pos_i + (ivar-1)*Nv
                         jG = 1 + k + (jvar-1)*nMax
                         iL = i + (ivar-1)*nNodesPerElement
                         jL = j + (jvar-1)*nNodesPerElement
                         Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                      END DO
                   END DO
                   EXIT
                END IF
             END DO
          END IF
       END DO
    END DO

  END SUBROUTINE AssembleSparseMatrixMS
  !===================================================!

END MODULE SparseMatrix_MS
