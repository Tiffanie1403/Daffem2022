MODULE Compute_DG

  USE Types
  USE PublicVar
  USE BCSelection
  USE SparseMatrix_DG
  USE ElementIntegral_Mixed
  USE EdgeIntegral_DG
  USE DirichletConformal_Mixed
  USE NeumannConformal_Mixed
  USE DirichletEmbedded_Mixed
  USE NeumannEmbedded_Mixed

  IMPLICIT NONE

CONTAINS

  !==================================================!
  SUBROUTINE Compute_LHS_RHS_DG(Data, Mesh, Sol, Mat)
  !==================================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    !--------------------------------------------
    type(element) :: NewElep, NewElem, Elep, Elem, ele, eleBC
    type(element) :: adj
    type(Edge)    :: Edg
    !----------------------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Aloc, Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Floc
    !------------------------------------------------------
    INTEGER :: Ne, Nv, nVar, ie, i, j
    INTEGER :: l, AdjNum
    !-----------------------------
    CHARACTER(len=20) :: BCType
    !-----------------------------
    INTEGER :: Nmax, ibc, ied, Ned, Nb, NeActive
    INTEGER :: nNodesPerElement
    !---------------------------------------
    LOGICAL :: eleSolve, eleID, eSm, eSp
    !-------------------------------------

    Ne = Mesh%Ne; Nv = Mesh%NvActive; Ned = Mesh%Ned ; Nb = Mesh%Nb ; nVar = Mat%nVar
    NeActive = mesh%NeActive
    DO ie = 1, Ne

       ele = Mesh%Ele(ie)
       lxx = ele%L(1,1) ; lxy = ele%L(1,2) ; lyx = ele%L(2,1) ; lyy = ele%L(2,2)
       nNodesPerElement = ele%nNodesPerElement

       ALLOCATE ( Aloc(nVar*nNodesPerElement,nVar*nNodesPerElement) )
       ALLOCATE ( Floc(nVar*nNodesPerElement) )

       CALL GetEmbeddedLogical(ele, Data%N_LS, eleSolve, eleID)

       IF ( ele%Solve ) THEN
          CALL Compute_ElementIntegral_Mixed(ele, Aloc, Floc, Data)
          CALL AssembleSparseMatrixDG(Mat, Aloc, ele, ele, Nv, NeActive)
          CALL AssembleRHSDG(Sol%RHS, Floc, ele, Nv, NeActive, nVar)

          DO i = 1, nNodesPerElement

             AdjNum = ele%Adj(i)/nNodesPerElement

             IF ( Ele%Adj(i) == -1 ) THEN

                CALL SelectBC(Ele, BCType, Data, i, iBC)

                SELECT CASE ( TRIM(ADJUSTL(BCType)) )
                CASE ( "Dirichlet" )
                   CALL DirichletConf_Mixed(Aloc, FLoc, ele, mesh%edg(ele%edg), Data, i, iBC)
                CASE ( "Neumann" )
                   CALL NeumannConf_Mixed(Aloc, FLoc, ele, mesh%edg(ele%edg), Data, i, iBC)
                END SELECT
                CALL AssembleSparseMatrixDG(Mat, Aloc, Ele, Ele, Nv, NeActive)
                CALL AssembleRHSDG(Sol%RHS, FLoc, Ele, Nv, NeActive, nVar)
             ELSE
                adj = mesh%ele(ele%adj(i)/nNodesPerElement)
                DO j = 1, Data%N_LS
                   IF ( adj%status(j) == 0 ) THEN
                      SELECT CASE ( TRIM(ADJUSTL(data%embeddedBCType(j))) )
                      CASE ( "Dirichlet","dirichlet" )
                         CALL DirichletEmb_Mixed(Aloc, Floc, ele, mesh%edg(ele%edg), &
                              sol%dist(j,ele%Vertex,:), data, i, j)
                      CASE ( "Neumann","neumann" )
                         CALL NeumannEmb_Mixed(Aloc, Floc, ele, mesh%edg(ele%edg), &
                              sol%dist(j,ele%Vertex,:), data, i, j)
                      END SELECT
                      CALL AssembleSparseMatrixDG(Mat, Aloc, Ele, Ele, Nv, NeActive)
                      CALL AssembleRHSDG(Sol%RHS, FLoc, Ele, Nv, NeActive, nVar)
                   END IF
                END DO
             END IF
          END DO
       END IF

       DEALLOCATE ( Aloc, Floc )

    END DO

    ! Internal edges integrals
    DO ied = 1, Ned
       Edg = Mesh%Edg(ied)
       IF ( Edg%Treated .eqv. .FALSE. ) THEN
          IF ( ( Edg%Tri(1) > 0 ) .AND. ( Edg%Tri(2) > 0 ) ) THEN
             EleM = Mesh%Ele(Edg%Tri(1))
             EleP = Mesh%Ele(Edg%Tri(2))
             CALL GetEmbeddedLogical(eleM, Data%N_LS, eSm, eleID)
             CALL GetEmbeddedLogical(eleP, Data%N_LS, eSp, eleID)
             IF ( eSm .AND. eSp ) THEN
                !PRINT*, ied
                nNodesPerElement = eleM%nNodesPerElement
                ALLOCATE ( Amm(nVar*nNodesPerElement,nVar*nNodesPerElement) )
                ALLOCATE ( Amp(nVar*nNodesPerElement,nVar*nNodesPerElement) )
                ALLOCATE ( Apm(nVar*nNodesPerElement,nVar*nNodesPerElement) )
                ALLOCATE ( App(nVar*nNodesPerElement,nVar*nNodesPerElement) )
                CALL ComputePMMatrices_DG(Edg, Mesh%Edg(eleM%Edg), Mesh%Edg(eleP%Edg), &
                     EleM, eleP, Amm, Amp, Apm, App, Data)

                CALL AssembleSparseMatrixDG(Mat, Amm, EleM, EleM, Nv, NeActive)
                CALL AssembleSparseMatrixDG(Mat, Amp, EleM, EleP, Nv, NeACtive)
                CALL AssembleSparseMatrixDG(Mat, Apm, EleP, EleM, Nv, NeActive)
                CALL AssembleSparseMatrixDG(Mat, App, EleP, EleP, Nv, NeActive)
                Mesh%Edg(ied)%Treated = .TRUE.
                DEALLOCATE ( Amm, Amp, Apm, App )

             END IF
          END IF
       END IF
    END DO


  END SUBROUTINE Compute_LHS_RHS_DG
  !==================================================!

END MODULE Compute_DG
