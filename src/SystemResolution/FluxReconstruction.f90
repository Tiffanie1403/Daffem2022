MODULE FluxReconstruction

  USE Types
  USE PublicVar
  USE FluxReconstructionLHSRHS
  USE Algebra
  USE BCSelection
  USE Permeability

  IMPLICIT NONE

CONTAINS

  !================================================!
  SUBROUTINE FluxReconstruction_DG(Sol, Mesh, Data)
  !================================================!

  !*************************************************************************
  !   Builds flux in the case of a problem in primal form for discontinuous resolution
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************



    IMPLICIT NONE
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(DataStructure), INTENT(IN)    :: Data
    !------------------------------------------
    type(element) :: ele, Adj
    type(edge)    :: edge
    !-------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Aloc, Am1
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: pLoc, Floc, Bloc, pLocAdj
    !-----------------------------------------------------------------
    CHARACTER(len=20) :: BCType
    !-----------------------------
    INTEGER :: Nv, Ne, Ned, ied, i, j, iBC, id, ie
    INTEGER :: nNodesPerElement, NU_i, Pos_i, nDim, nVar
    !-----------------------------------------------------

    Nv = Mesh%Nv ; Ne = Mesh%Ne ; Ned = Mesh%Ned
    nDim = Data%nDim ; nVar = nDim + 1
    DO ie = 1, Ne
       ele = mesh%ele(ie)
       IF ( ele%solve ) THEN
          nNodesPerElement = ele%nNodesPerElement

          ALLOCATE ( Aloc(nDim*nNodesPerElement,nDim*nNodesPerElement) )
          ALLOCATE ( Am1(nDim*nNodesPerElement,nDim*nNodesPerElement) )
          ALLOCATE ( Floc(nDim*nNodesPerElement) )
          ALLOCATE ( Bloc(nDim*nNodesPerElement) )
          ALLOCATE ( pLoc(nNodesPerElement) )
          ALLOCATE ( pLocAdj(nDim*nNodesPerElement) )

          DO i = 1, nNodesPerElement
             Pos_i = ele%Pos(i)
             pLoc(i) = sol%pDG(Pos_i)
          END DO

          CALL ElementIntegral_GR(Aloc, Floc, pLoc, ele, Data)
          DO i = 1, nNodesPerElement
             ! Conformal BC
             IF ( ele%adj(i) == -1 ) THEN
                CALL RHS_edgInt(mesh%edg(ele%edg(i)), ele, ele, Floc, pLoc, pLoc, Data, i)
             ELSE
                adj = mesh%ele(ele%adj(i)/nNodesPerElement)
                ! Internal boundary edges
                IF ( adj%solve ) THEN
                   pLocAdj = sol%pDG(adj%Pos)
                   CALL RHS_edgInt(mesh%edg(ele%edg(i)), ele, Adj, Floc, pLoc, pLocAdj, Data, i)
                ! Embedded BC
                ELSE
                   CALL RHS_edgInt(mesh%edg(ele%edg(i)), ele, ele, Floc, pLoc, pLoc, Data, i)
                END IF
             END IF
          END DO

          CALL InverseMatrix(Aloc, Am1)
          CALL CheckInversion(Aloc,Am1)

          Bloc = MATMUL(Am1,Floc)
          DO i = 1, nNodesPerElement
             DO id = 1, nDim
                Sol%Bdg(ele%Pos(i),id) = Bloc(i + (id-1)*nNodesPerElement)
             END DO
          END DO
          DEALLOCATE ( Aloc, Am1, Floc, Bloc, pLoc, pLocAdj )
       END IF
    END DO

  END SUBROUTINE FluxReconstruction_DG
  !================================================!

  !========================================================!
  SUBROUTINE GreenGaussFluxReconstruction(Sol, Mesh, data)
  !========================================================!

  !*************************************************************************
  !   Builds flux in the case of a problem in primal form for continuous resolution
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************



    IMPLICIT NONE
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(DataStructure), INTENT(IN)    :: Data
    !-------------------------------------------
    type(element) :: ele
    !---------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: L
    REAL*8, DIMENSION(:),   ALLOCATABLE :: P1Grad, diag, GSF_i
    !----------------------------------------------------------
    INTEGER :: ie, i, NU_i
    INTEGER :: nDim, nNodesPerElement, Nv
    !-------------------------------------

    nDim = data%nDim ; Nv = mesh%Nv

    sol%Bcg = 0.d0

    ALLOCATE( P1Grad(nDim), diag(Nv), GSF_i(nDim) )
    ALLOCATE( L(nDim,nDim) )

    diag = 0.d0

    DO ie = 1, mesh%Ne
       ele = mesh%ele(ie)
       nNodesPerElement = ele%nNodesPerElement
       IF ( ele%solve ) THEN
          P1Grad = 0.d0
          CALL SetNodalMatrixPermeability(L, ele, data%icas)
          DO i = 1, nNodesPerElement
             NU_i = ele%vertex(i)
             diag(NU_i) = diag(NU_i) + ele%A
             CALL GSFEval(GSF_i, mesh%Vertex(NU_i)%X, i, ele)
             P1Grad = P1Grad + sol%pcg(NU_i)*GSF_i
          END DO
          DO i = 1, nNodesPerElement
             NU_i = ele%Vertex(i)
             Sol%Bcg(NU_i,:) = Sol%Bcg(NU_i,:) + MATMUL(L,P1Grad)*ele%A
          END DO
       END IF
    END DO

    DO i = 1, Nv
       IF ( diag(i) > 1e-10 ) THEN
          Sol%Bcg(i,:) = -Sol%Bcg(i,:)/diag(i)
       END IF
    END DO

    DEALLOCATE ( P1Grad, diag, GSF_i, L)

  END SUBROUTINE GreenGaussFluxReconstruction
  !========================================================!

END MODULE FluxReconstruction
