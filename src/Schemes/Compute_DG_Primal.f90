MODULE Compute_DG_Primal

  USE Types
  USE PublicVar
  USE ElementIntegral_Primal
  USE EdgeIntegral_DG_Primal
  USE SparseMatrix_DG
  USE BCSelection
  USE DirichletConformal_Primal
  USE NeumannConformal_Primal
  USE DirichletEmbedded_Primal
  USE NeumannEmbedded_Primal
    
  IMPLICIT NONE

CONTAINS

  !=========================================================!
  SUBROUTINE Compute_LHS_RHS_DG_Primal(Data, Mesh, Sol, Mat)
  !=========================================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    !-------------------------------------------
    type(element) :: ele, eleM, eleP, eleBC, adj
    type(Edge)    :: Edg
    !-------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Aloc, Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Floc, Uloc, Ulocm, Ulocp, Fm, Fp
    !----------------------------------------------------------------------
    CHARACTER(len=20) :: BCType
    !---------------------------
    INTEGER :: AdjNum
    INTEGER :: Nv, Ne, Ned, Nb, nVar, nNodesPerElement
    INTEGER :: ie, i, j, iBC, ivar, it, N
    !-----------------------------------------------------
    LOGICAL :: PB
    !------------
    
    Ne = Mesh%Ne ; Nv = Mesh%Nv ; Nb = Mesh%Nb ; Ned = Mesh%Ned
    nVar = Mat%nVar

    Mat%Coeff = 0.d0 ; Sol%RHS = 0.d0

    it = 0
    
    DO ie = 1, Ne

       ele = Mesh%ele(ie)
       lxx = ele%L(1,1) ; lxy = ele%L(1,2) ; lyx = ele%L(2,1) ; lyy = ele%L(2,2)
       nNodesPerElement = ele%nNodesPerElement          
       ALLOCATE ( Aloc(nNodesPerElement,nNodesPerElement) )
       ALLOCATE ( Floc(nNodesPerElement), Uloc(nNodesPerElement) )
       
       IF ( ele%solve ) THEN

          ! Recover local solution
          DO i = 1, nNodesPerElement             
             Uloc(i) = Sol%pdg(ele%Pos(i))
          END DO
          
          CALL Compute_ElementIntegral_Primal(ele, Aloc, Uloc, Floc, Data)
          CALL AssembleSparseMatrixDG(Mat, Aloc, ele, ele, Nv, Ne)
          CALL AssembleRHSDG(Sol%RHS, Floc, ele, Nv, Ne, nVar)
                    
          DO i = 1, nNodesPerElement
             
             AdjNum = ele%Adj(i)/nNodesPerElement
             
             ! Conformal BC
             IF ( ele%Adj(i) == -1 ) THEN
                
                CALL SelectBC(ele, BCType, Data, i, iBC)
                
                SELECT CASE ( TRIM(ADJUSTL(BCType)) )
                CASE ( "Dirichlet", "dirichlet" )
                   CALL DirichletConf_Primal(Aloc, Uloc, Floc, ele, mesh%edg(ele%Edg), Data, i, iBC)
                CASE ( "Neumann" )
                   CALL NeumannConf_Primal(Aloc, Uloc, Floc, ele, mesh%edg(ele%Edg), Data, i, iBC)
                CASE DEFAULT
                   PRINT*, ie, TRIM(ADJUSTL(BCType))
                   PRINT*, " ** ERROR in BC type... (Compute_DG_Primal.f90)"                   
                   PRINT*, "    Available :"
                   PRINT*, "      + Dirichlet/dirichlet"
                   PRINT*, "      + Neumann/neumann"
                   STOP
                END SELECT
                
                CALL AssembleSparseMatrixDG(Mat, Aloc, ele, ele, Nv, Ne)
                CALL AssembleRHSDG(Sol%RHS, Floc, ele, Nv, Ne, nVar)
                                
                ! Embedded BC
             ELSE
                adj = mesh%ele(ele%adj(i)/nNodesPerElement)
                DO j = 1, data%N_LS
                   IF ( adj%status(j) == 0 ) THEN
                      SELECT CASE ( TRIM(ADJUSTL(data%embeddedBCType(j))) )
                      CASE ( "Dirichlet" )
                         CALL DirichletEmb_P(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), &
                              sol%dist(j,ele%vertex,:),  data, i, j)
                      CASE ( "Neumann" )
                         CALL NeumannEmb_P(Aloc, Uloc, Floc, ele, mesh%edg(ele%edg), &
                              sol%dist(j,ele%vertex,:), data, i, j)
                      END SELECT
                      CALL AssembleSparseMatrixDG(Mat, Aloc, ele, ele, Nv, Ne)
                      CALL AssembleRHSDG(Sol%RHS, Floc, ele, Nv, Ne, nVar)
                      CALL AssembleRHSDG(sol%RHS, -MATMUL(Aloc,Uloc), ele, Nv, Ne, nVar)
                   END IF
                END DO
             END IF
          END DO
       END IF
       DEALLOCATE ( Aloc, Floc, Uloc )
       
    END DO

    DO ie = 1, Mesh%Ned
       Edg = Mesh%Edg(ie)
       IF ( Edg%Treated .eqv. .FALSE. ) THEN
          IF ( ( Edg%Tri(1) > 0 ) .AND. ( Edg%Tri(2) > 0 ) ) THEN
             
             eleM = Mesh%Ele(Edg%Tri(1))
             eleP = Mesh%Ele(Edg%Tri(2))
             IF ( eleM%solve .AND. eleP%solve ) THEN
                nNodesPerElement = eleM%nNodesPerElement
                ALLOCATE ( Amm(nNodesPerElement,nNodesPerElement) )
                ALLOCATE ( Amp(nNodesPerElement,nNodesPerElement) )
                ALLOCATE ( Apm(nNodesPerElement,nNodesPerElement) )
                ALLOCATE ( App(nNodesPerElement,nNodesPerElement) )
                ALLOCATE ( Ulocm(nNodesPerElement), Ulocp(nNodesPerElement) )
                ALLOCATE ( Fm(nNodesPerElement), Fp(nNodesPerElement) )
                ! Recover local sol on each element
                DO i = 1, nNodesPerElement
                   Ulocm(i) = sol%pdg(eleM%Pos(i))
                   Ulocp(i) = sol%pdg(eleP%Pos(i))
                END DO

                CALL ComputePMMatrices_Primal(Edg, Mesh%Edg(eleM%edg), Mesh%Edg(eleP%edg), &
                     eleM, eleP, Amm, Amp, Apm, App, Ulocm, Ulocp, Fm, Fp, Data)
                
                CALL AssembleSparseMatrixDG(Mat, Amm, eleM, eleM, Nv, Ne)
                CALL AssembleSparseMatrixDG(Mat, Amp, eleM, eleP, Nv, Ne)
                CALL AssembleSparseMatrixDG(Mat, Apm, eleP, eleM, Nv, Ne)
                CALL AssembleSparseMatrixDG(Mat, App, eleP, eleP, Nv, Ne)

                CALL AssembleRHSDG(sol%rhs, -MATMUL(Amm,UlocM), eleM, Nv, Ne, nVar)
                CALL AssembleRHSDG(sol%rhs, -MATMUL(Amp,UlocP), eleM, Nv, Ne, nVar)
                CALL AssembleRHSDG(sol%rhs, -MATMUL(Apm,UlocM), eleP, Nv, Ne, nVar)
                CALL AssembleRHSDG(sol%rhs, -MATMUL(App,UlocP), eleP, Nv, Ne, nVar)

                !CALL AssembleRHSDG(Sol%RHS, Fm, eleM, Nv, Ne, nVar)
                !CALL AssembleRHSDG(Sol%RHS, Fp, eleP, Nv, Ne, nVar)
                
                
                Mesh%Edg(ie)%Treated = .TRUE.
                
                DEALLOCATE ( Amm, Amp, Apm, App, Ulocm, Ulocp, Fm, Fp )
             END IF
          END IF
       END IF
       
    END DO

    Sol%rhsTmp = sol%rhs
    
  END SUBROUTINE Compute_LHS_RHS_DG_Primal
  !============================================================!
  
END MODULE Compute_DG_Primal
