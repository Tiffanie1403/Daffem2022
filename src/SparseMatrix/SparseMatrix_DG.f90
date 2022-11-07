MODULE SparseMatrix_DG

  USE Types

  IMPLICIT NONE

CONTAINS

  !===================================================!
  SUBROUTINE AssembleRHSDG(F, Floc, ele, Nv, Ne, nVar)
  !===================================================!

  !*************************************************************************
  !   Adds the local matrix to the global matrix discretization of the second member
  !   Discontinuous Galerkin case
  !
  !  Parameters:
  !
  !    Input,Output, REAL*8, DIMENSION(:) F, Structure for the global second member
  !    Input, REAL*8, DIMENSION(:,:) Floc, Local matrix for the second member
  !    Input, type(element) ele, element considerated
  !    Input, INTEGER  Nv, Number of nodes (vertices)
  !    Input, INTEGER  Ne, Number of elements
  !    Input, INTEGER  nVar, Number of variables
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: F
    REAL*8, DIMENSION(:), INTENT(IN)    :: Floc
    type(element)       , INTENT(IN)    :: ele
    INTEGER             , INTENT(IN)    :: Nv, Ne, nVar
    !----------------------------------------------------
    INTEGER :: nDim, nNodesPerElement
    INTEGER :: i, NU_i, ivar, iG, iL
    !-----------------------------------

    nNodesPerElement = ele%nNodesPerElement

    DO i = 1, nNodesPerElement
       NU_i = ele%surPos(i)
       DO ivar = 1, Nvar
          iL = i + (ivar-1)*nNodesPerElement
          iG = NU_i + (ivar-1)*nNodesPerElement*Ne
          F(iG) = F(iG) + Floc(iL)
       END DO
    END DO

  END SUBROUTINE AssembleRHSDG
  !========================================================!

  !===============================================================!
  SUBROUTINE AssembleSparseMatrixDG(Mat, Aloc, elei, elej, Nv, Ne)
  !===============================================================!

  !*************************************************************************
  !   Adds the local matrix to the global matrix discretization
  !   Discontinuous Galerkin case
  !
  !  Parameters:
  !
  !    Input,Output, type(SparseMatrix) Mat, Structure for the global matrix
  !    Input, REAL*8, DIMENSION(:,:) Aloc, Local matrix
  !    Input, type(element) elei,elej, elements considerated (share an edge)
  !    Input, INTEGER  Nv, Number of nodes (vertices)
  !    Input, INTEGER  Ne, Number of elements
  !
  !*************************************************************************

    IMPLICIT NONE
    type(SparseMatrix)    , INTENT(INOUT) :: Mat
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: Aloc
    type(element)         , INTENT(IN)    :: elei, elej
    INTEGER               , INTENT(IN)    :: Nv, Ne
    !-------------------------------------------------------
    INTEGER :: i, j, NU_i, NU_j, Pos_i, Pos_j, k, Pos_k
    INTEGER :: Nei, Nej, NvLoc, Nmax, nVar
    INTEGER :: l, m, iG, jG, iL, jL, nNodesPerElement
    !---------------------------------------------------
    LOGICAL :: found
    !-----------------

    Nei = elei%Num ; Nej = elej%Num
    NMax = Mat%NMax ; nVar = Mat%nVar
    nNodesPerElement = elei%nNodesPerElement

    DO i = 1, nNodesPerElement
       Pos_i = elei%surPos(i)
        DO j = 1, nNodesPerElement
          Pos_j = elej%surPos(j)
          IF ( Pos_i == Pos_j ) THEN
             DO l = 1, nVar
                DO m = 1, nVar
                   iG = Pos_i + (l-1)*nNodesPerElement*Ne
                   jG = 1 + (m-1)*Nmax
                   iL = i + (l-1)*nNodesPerElement
                   jL = j + (m-1)*nNodesPerElement
                   Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                END DO
             END DO
          ELSE
             NvLoc = Mat%Ind(Pos_i,1)
             DO k = 1, NvLoc
                Pos_k = Mat%Ind(Pos_i,1+k)
                IF ( Pos_j == Pos_k ) THEN
                   DO l = 1, nVar
                      DO m = 1, nVar
                         iG = Pos_i + (l-1)*nNodesPerElement*Ne
                         jG = 1 + k + (m-1)*Nmax
                         iL = i + (l-1)*nNodesPerElement
                         jL = j + (m-1)*nNodesPerElement
                         Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                      END DO
                   END DO
                   EXIT
                END IF
             END DO
          END IF
       END DO
    END DO

  END SUBROUTINE AssembleSparseMatrixDG
  !=================================================================!

  !==============================!
  SUBROUTINE FillIndDG(Mat, Mesh)
  !==============================!

  !*************************************************************************
  !   Fills Mat%Ind for the Storage of the discretized matrix
  !   Discontinuous version
  !
  !  Parameters:
  !
  !    Input,Output, type(MeshType) Mesh, Mesh structure
  !    Input,Output, type(SparseMatrix) Mat, Structure for the matrices of the resolution
  !
  !*************************************************************************



    IMPLICIT NONE
    type(SparseMatrix), INTENT(INOUT) :: Mat
    type(MeshType)    , INTENT(IN)    :: Mesh
    !------------------------------------------
    type(element) :: ele, Adj
    !---------------------------
    INTEGER :: NMax, NvLoc, nVar, nNodesPerElement
    INTEGER :: Nv, Ne, ie
    INTEGER :: AdjNum
    INTEGER :: i, j, NU_i, NU_j, k, NU_k
    INTEGER :: Pos_i, Pos_j, Pos_k
    !------------------------------------------------

    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nVar = Mat%nVar

    DO ie = 1, Ne

       ele = Mesh%ele(ie)
       nNodesPerElement = ele%nNodesPerElement

       IF ( ele%Solve ) THEN

          DO i = 1, nNodesPerElement

             Pos_i = ele%surPos(i)

             DO j = 1, nNodesPerElement

                Pos_j = ele%surPos(j)
                IF ( Pos_i /= Pos_j ) THEN
                   Mat%Ind(Pos_i,1) = Mat%Ind(Pos_i,1) + 1
                   NvLoc = Mat%Ind(Pos_i,1)
                   Mat%Ind(Pos_i,1+NvLoc) = Pos_j
                END IF
             END DO
             DO j = 1, nNodesPerElement
                IF ( ele%Adj(j) /= -1 ) THEN
                   AdjNum = ele%Adj(j)/nNodesPerElement
                   Adj = Mesh%Ele(AdjNum)
                   IF ( adj%solve ) THEN
                      DO k = 1, nNodesPerElement
                         Pos_k = Adj%surPos(k)
                         Mat%Ind(Pos_i,1) = Mat%Ind(Pos_i,1) + 1
                         NvLoc = Mat%Ind(Pos_i,1)
                         Mat%Ind(Pos_i,1+NvLoc) = Pos_k
                      END DO
                   END IF
                END IF
             END DO
          END DO
       END IF
    END DO

    Mat%NZ = 0
    DO i = 1, nNodesPerElement*mesh%NeActive
       Mat%NZ = Mat%NZ + Mat%Ind(i,1) + 1
    END DO
    Mat%NZ = nVar*nVar*Mat%NZ

  END SUBROUTINE FillIndDG
  !==============================!


  !===============================================================!
  SUBROUTINE AssembleSparseMatrixDG2(Mat, Aloc, elei, elej, Nv, Ne)
  !===============================================================!
  !*************************************************************************
  !   Assemble of the local matrices to the global one
  !   Subroutine for the interface terms
  !
  !  Parameters:
  !
  !    Input,Output, type(SparseMatrix) Mat, Structure for the matrices of the resolution
  !    Input, REAL*8, DIMENSION(:,:) Aloc, local matrices
  !    Input, type(element)  elei, elej, elements at the interface
  !    Input, INTEGER Nv, vertices number without the complementary nodes
  !    Input, INTEGER  Ne, elements number
  !
  !*************************************************************************


    IMPLICIT NONE
    type(SparseMatrix)    , INTENT(INOUT) :: Mat
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: Aloc
    type(element)         , INTENT(IN)    :: elei, elej
    INTEGER               , INTENT(IN)    :: Nv, Ne
    !-------------------------------------------------------
    INTEGER :: i, j, NU_i, NU_j, Pos_i, Pos_j, k, Pos_k
    INTEGER :: Nei, Nej, NvLoc, Nmax, nVar
    INTEGER :: l, m, iG, jG, iL, jL, nNodesPerElement
    !---------------------------------------------------
    LOGICAL :: found
    !-----------------

    Nei = elei%Num ; Nej = elej%Num
    NMax = Mat%NMax ; nVar = Mat%nVar
    nNodesPerElement = elei%nNodesPerElement

    DO i = 1, nNodesPerElement
       Pos_i = elei%SurVertex(i)

        DO j = 1, nNodesPerElement
          Pos_j = elej%SurVertex(j)

          IF ( Pos_i == Pos_j ) THEN
             DO l = 1, nVar
                DO m = 1, nVar

                   iG = Pos_i + (l - 1)*Nv
                   jG = 1 + (m - 1)*NMax
                   iL = i + (l - 1)*nNodesPerElement
                   jL = j + (m - 1)*nNodesPerElement

                   Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                END DO
             END DO
          ELSE
             NvLoc = Mat%Ind(Pos_i,1) ! nombre de voisin associé au noeud 
             DO k = 1, NvLoc
                Pos_k = Mat%Ind(Pos_i,1+k)
                IF ( Pos_j == Pos_k ) THEN
                   DO l = 1, nVar
                      DO m = 1, nVar

                         iG = Pos_i + (l - 1)*Nv
                         jG = 1 + k + (m - 1)*NMax
                         iL = i + (l - 1)*nNodesPerElement
                         jL = j + (m - 1)*nNodesPerElement

                         Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                      END DO
                   END DO
                   EXIT
                END IF
             END DO
          END IF
       END DO
    END DO

  END SUBROUTINE AssembleSparseMatrixDG2
  !=================================================================!


END MODULE SparseMatrix_DG
