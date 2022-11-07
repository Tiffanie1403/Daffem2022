MODULE Compute_CG_Primal

  USE Types
  USE PublicVar
  USE ElementIntegral_Primal
  USE SparseMatrix_CG
  USE BCSelection
  USE DirichletConformal_Primal
  USE DirichletIcingConformal_Primal
  USE DirichletEmbedded_Primal
  USE NeumannConformal_Primal
  USE NeumannIcingConformal_Primal
  USE NeumannEmbedded_Primal
  USE NeumannEmbedded_Primal_AG
  USE NeumannJumpConformal_Primal
  USE NeumannJumpEmbedded_Primal
  USE Tools
  USE ImmersedElementIntegral
  USE ReadMesh_Mod

  IMPLICIT NONE

CONTAINS

  !=========================================================!
  SUBROUTINE Compute_LHS_RHS_CG_Primal(Data, Mesh, Sol, Mat)
  !=========================================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    !-------------------------------------------
    type(element) :: ele, eleM, eleP, eleBC, adj, ele1, ele2
    type(Edge)    :: Edg
    type(Vertex)  :: Node1, Node2, Node3
    !--------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Aloc
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Floc, Uloc
    !--------------------------------------------------
    CHARACTER(len=20) :: BCType
    !--------------------------
    INTEGER :: AdjNum, iFace, iFace2, iFace1
    INTEGER :: Nv, Ne, Ned, Nb, nVar, N, nDim
    INTEGER :: ie, ib, i, j, iBC, nNodesPerElement, k
    !------------------------------------------
    LOGICAL :: eleID, eleSolve
    !--------------------------

    Ne = Mesh%Ne  !nombre d'elements du maillage
    Nv = Mesh%NvActive
    Nb = Mesh%Nb
    Ned = Mesh%Ned
    nVar = Mat%NVar
    nDim = data%nDim !dimension du probleme
    Sol%RHS = 0.d0
    mat%Coeff = 0.d0

    !calcul et assemblage des (Gv,Gu)
    DO ie = 1, Ne
       ele = mesh%ele(ie) !on se place sur un element du maillage
       lxx = ele%L(1,1) ; lxy = ele%L(1,2) ; lyx = ele%L(2,1) ; lyy = ele%L(2,2) ! coefficient de la matrice de permeabilite associe a la zone de l'element
       nNodesPerElement = ele%nNodesPerElement
       ALLOCATE ( Aloc(nNodesPerElement,nNodesPerElement) ) !creation de la de matrice elementaire de l'element courant  que l'on nomme Aloc
       ALLOCATE ( Floc(nNodesPerElement), Uloc(nNodesPerElement) ) !deux vecteurs de la taille du nombre d'éléments

       CALL GetEmbeddedLogical(ele, data%N_LS, eleSolve, eleID)

        ! Recover local solution
        DO i = 1, nNodesPerElement
           Uloc(i) = Sol%pcg(ele%Vertex(i))
        END DO


        IF ( ele%solve ) THEN
           !PRINT*, ie, ele%L(1,1)
           CALL Compute_ElementIntegral_Primal(ele, Aloc, Uloc, Floc, Data)
           CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv, Ne)
           CALL AssembleRHSCG(Sol%RHS, Floc, ele, Nv, Ne, nVar)

        END IF

        DEALLOCATE ( Aloc, Floc, Uloc )
    END DO

    ! Conformal BC
    DO ib = 1, mesh%Nb
       edg = mesh%boundary(ib)
       ele = mesh%ele(edg%tri(1))
       IF ( ele%Solve ) THEN
          nNodesPerElement = ele%nNodesPerElement
          ALLOCATE ( Aloc(nNodesPerElement,nNodesPerElement) )
          ALLOCATE ( Floc(nNodesPerElement), Uloc(nNodesPerElement) )
          DO i = 1, nNodesPerElement
             Uloc(i) = Sol%pcg(ele%vertex(i))
          END DO
          CALL SelectBC2(edg, BCType, Data, iBC)
          CALL GetFace(ele, edg, iFace)
          SELECT CASE ( TRIM(ADJUSTL(bcType)) )
          CASE ( "Dirichlet" )
             IF ( .NOT. data%icing ) THEN
                CALL DirichletConf_Primal(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
             ELSE
                CALL DirichletIcingConf_Primal(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
             END IF
          CASE ( "Neumann" )
             IF ( .NOT. data%icing ) THEN
                CALL NeumannConf_Primal(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
             ELSE
                CALL NeumannIcingConf_Primal(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
             END IF
          CASE ( "NeumannJump" )
             ele1 = mesh%ele(edg%tri(1))
             ele2 = mesh%ele(edg%tri(2))
             CALL GetFace(ele1, edg, iFace1)
             CALL GetFace(ele2, edg, iFace2)
             CALL NeumannJumpConf_Primal(Aloc, Uloc, Floc, ele1, ele2, edg, iFace1, iFace2, iBC, data)
             ele = ele1
          CASE DEFAULT
             PRINT*, ib, TRIM(ADJUSTL(BCType))
             PRINT*, " ** ERROR in BC type... (Compute_CG_Primal.f90)"
             PRINT*, "    Available :"
             PRINT*, "      + Dirichlet/dirichlet"
             PRINT*, "      + Neumann/neumann"
             PRINT*, "      + NeumannJump"
             STOP
          END SELECT
          CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv, Ne)
          CALL AssembleRHSCG(Sol%RHS, Floc, ele, Nv, Ne, nVar)
          DEALLOCATE(Aloc, Floc, Uloc)
       END IF
    END DO

    ! Embedded BC
    IF ( data%embedded ) THEN
       DO ib = 1, mesh%nSurB
          edg = mesh%embBoundary(ib)
          nNodesPerElement = nDim + 1
          ALLOCATE ( Aloc(nNodesPerElement, nNodesPerElement) )
          ALLOCATE ( Floc(nNodesPerElement), Uloc(nNodesPerElement) )

          ! Revocer local solution on active element
          IF ( mesh%ele(mesh%EmbBoundary(ib)%Tri(1))%Solve ) THEN
             ele = mesh%ele(mesh%EmbBoundary(ib)%Tri(1))
          ELSE
             ele = mesh%ele(mesh%EmbBoundary(ib)%Tri(2))
          END IF
          DO j = 1, nNodesPerElement
             Uloc(j) = sol%pcg(ele%vertex(j))
          END DO

          CALL GetFace(ele, edg, iFace)

          SELECT CASE ( TRIM(ADJUSTL(mesh%EmbBoundary(ib)%BC_Type)) )
          CASE ( "Dirichlet" )
             CALL DirichletEmb_P(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), &
                  sol%dist(1,ele%vertex,:),data,iFace,mesh%embBoundary(ib)%embLS )
          CASE ( "Neumann" )

             CALL NeumannEmb_P(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), &
                  sol%dist(1,ele%vertex,:), data, iFace, mesh%embBoundary(ib)%embLS )
          CASE ( "NeumannJump" )

             ele1 = mesh%ele(edg%tri(1))
             ele2 = mesh%ele(edg%tri(2))
             CALL GetFace(ele1, edg, iFace1)
             CALL GetFace(ele2, edg, iFace2)
             CALL NeumannJumpEmb_P(Aloc, Uloc, Floc, ele1, ele2, edg, iFace1, iFace2, &
                  mesh%embBoundary(ib)%embLS, sol%dist(1,ele%vertex,:), data)
             ele = ele1

           CASE DEFAULT
              PRINT*, ie, TRIM(ADJUSTL(BCType))
              PRINT*, " ** ERROR in BC type... (Compute_CG_Primal.f90)"
              PRINT*, "    Available :"
              PRINT*, "      + Dirichlet/dirichlet"
              PRINT*, "      + Neumann/neumann"
              PRINT*, "      + NeumannJump"
              STOP
          END SELECT

          CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv, Ne)
          CALL AssembleRHSCG(Sol%RHS, Floc, ele, Nv, Ne, nVar)
          CALL AssembleRHSCG(sol%RHS, -MATMUL(Aloc,Uloc), ele, Nv, Ne, nVar)

          DEALLOCATE ( Aloc, Floc, Uloc )
       END DO
    END IF
    Sol%RHSTmp = Sol%RHS

  END SUBROUTINE Compute_LHS_RHS_CG_Primal
  !========================================================!

END MODULE Compute_CG_Primal
