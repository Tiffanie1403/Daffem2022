MODULE Compute_MSEG
  
  USE Types
  USE PublicVar
  USE BCSelection
  USE SparseMatrix_CG
  USE SparseMatrix_MS
  USE ElementIntegral_EG
  USE EdgeIntegral_EG
  USE DirichletConformal_EG
  USE NeumannConformal_EG
  USE TransferOperator_MSEG
  USE Algebra
  
  IMPLICIT NONE

CONTAINS

  !====================================================!
  SUBROUTINE Compute_LHS_RHS_MSEG(Data, Mesh, Sol, Mat)
  !====================================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    !--------------------------------------------
    type(Element) :: Ele, Elem, Elep
    type(Edge)    :: Edg
    !---------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Mt, Mtm1, Mb, T, TT, MtM, MtMm1, MbM, TM
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: TTM, MtP, MtPm1, MbP, TP, TTP, MtBC, MtBCm1
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: MbBC, TBC, TTBC, Aloc, Amm, Amp, Apm, App
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: TTAlocT, TTAmmT, TTAmpT, TTApmT, TTAppT
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: FDirLoc, F, F2, Fp, Fm, F2p, F2m, TPF, TPFm, TPFp
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: TTFDirLoc, TTF, TTfp, TTFm, TTATPF, TTATPFp
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: TTATPFmm, TTApmTF, TTAppTF
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Phi, Phip, Phim, TTATF, TTAmmTF, TTAmpTF 
    !----------------------------------------------------------------------------------------
    INTEGER :: ie, i, j, NU_i, iBC, Ndim, Nv, Ne, nVar, nNodesPerElement, nDofPerElement
    !------------------------------------------------------------------------------------
    CHARACTER(len=20) :: BCType
    !---------------------------
    
    Ndim = Data%Ndim ; Nv = Mesh%Nv ; Ne = Mesh%Ne ; nVar = Mat%nVar
    
    DO ie = 1, Mesh%Ne

       ele = Mesh%Ele(ie)
       lxx = ele%L(1,1) ; lxy = ele%L(1,2) ; lyx = ele%L(2,1) ; lyy = ele%L(2,2)       
       nNodesPerElement = ele%nNodesPerElement ; nDofPerElement = nNodesPerElement + 1
       ALLOCATE ( Mt(nVar*nDofPerElement,nVar*nDofPerElement) )
       ALLOCATE ( Mtm1(nVar*nDofPerElement,nVar*nDofPerElement) )
       ALLOCATE ( Mb(nVar*nDofPerElement,nVar*nNodesPerElement) )

       ALLOCATE ( T(nVar*nDofPerElement,nVar*nNodesPerElement) )
       ALLOCATE ( TT(nVar*nNodesPerElement,nVar*nDofPerElement) )

       ALLOCATE ( Aloc(nVar*nDofPerElement,nVar*nDofPErElement) )
       ALLOCATE ( F(nVar*nDofPerElement) )
       ALLOCATE ( F2(nVar*nDofPerElement) )
       ALLOCATE ( TTF(nVar*nNodesPerElement) )
       ALLOCATE ( TPF(nVar*nNodesPerElement) )
       ALLOCATE ( TTATPF(nVar*nNodesPerElement) )
       
       ! Computation of Transfer matrix T and transpose TT
       CALL ComputeMt_MSEG(Ele, Mesh%Edg(Ele%Edg), Mt, Data)
       CALL ComputeMb_MSEG(Ele, Mb, Data)       
       CALL InverseMatrix(Mt, Mtm1)
       CALL CheckInversion(Mt, Mtm1)
       T = MATMUL(Mtm1, Mb)
       
       CALL Transpose(T, TT)
       
       CALL ComputeTransferSourceTerm_MSEG(F2, ele, Data)
       
       !** Element Integral **!
       ! Computation
       CALL Compute_ElementIntegral_EG(Ele, Aloc, F, Data)

       ! Identity Part
       TTF = MATMUL(TT,F)
       CALL AssembleRHSCG(Sol%RHS, TTF, ele, Nv, Ne, nVar)

       ! Transfered matrix
       TTAlocT = MATMUL(TT,MATMUL(Aloc,T))

       
       
       CALL AssembleSparseMatrixMS(Mat, TTAlocT, ele, ele, Nv, Ne)

       ! Source Term
       TPF = MATMUL(Mtm1, F2)
       TTATPF = MATMUL(TT,MATMUL(Aloc,TPF))       
       CALL AssembleRHSCG(Sol%RHS, -TTATPF, ele, Nv, Ne, nVar)
       
       DO i = 1, nNodesPerElement
          
          IF ( Ele%Adj(i) == -1 ) THEN

             CALL SelectBC(Ele, BCType, Data, i, iBC)
             
             SELECT CASE ( BCType )
             CASE ( "Dirichlet" )
                CALL DirichletConf_EG(Aloc, F, Ele, mesh%edg(ele%edg), Data, i, iBC)
             CASE ( "Neumann" )
                CALL NeumannConf_EG(Aloc, F, Ele, mesh%edg(ele%edg), Data, i, iBC)
             END SELECT
             
             TTAlocT = MATMUL(TT,MATMUL(Aloc,T))
             TTATF = MATMUL(TT,MATMUL(Aloc,TPF))
             
             CALL AssembleSparseMatrixMS(Mat, TTAlocT, Ele, Ele, Nv, Ne)
             CALL AssembleRHSCG(Sol%RHS, -TTATF, ele, Nv, Ne, nVar)
                          
             TTF = MATMUL(TT, F)
             CALL AssembleRHSCG(Sol%RHS, TTF, ele, Nv, Ne, nVar)
             
          END IF
          
       END DO

       DEALLOCATE ( Mt, Mtm1, Mb, T, TT, Aloc, F, F2, TTF, TPF, TTATPF )       
       
    END DO

    DO ie = 1, Mesh%Ned
       Edg = Mesh%Edg(ie)
       IF ( Edg%Treated .eqv. .FALSE. ) THEN
          !PRINT*, ie
          IF ( ( Edg%Tri(1) > 0 ) .AND. ( Edg%Tri(2) > 0 ) ) THEN
             eleM = Mesh%Ele(Edg%Tri(1))
             eleP = Mesh%Ele(Edg%Tri(2))
             nNodesPerElement = eleM%nNodesPerElement
             nDofPerElement = nNodesPerElement + 1             
             ALLOCATE ( MtM(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( MtMm1(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( MbM(nVar*nDofPerElement,nVar*nNodesPerElement) )

             ALLOCATE ( TM(nVar*nDofPerElement,nVar*nNodesPerElement) )
             ALLOCATE ( TTM(nVar*nNodesPerElement,nVar*nDofPerElement) )

             ALLOCATE ( MtP(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( MtPm1(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( MbP(nVar*nDofPerElement,nVar*nNodesPerElement) )

             ALLOCATE ( TP(nVar*nDofPerElement,nVar*nNodesPerElement) )
             ALLOCATE ( TTP(nVar*nNodesPerElement,nVar*nDofPerElement) )

             ALLOCATE ( F2M(nVar*nDofPerElement), F2P(nVar*nDofPerElement) )
             ALLOCATE ( TPFm(nVar*nDofPerElement), TPFp(nVar*nDofPerElement) )
             ALLOCATE ( Amm(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( Amp(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( Apm(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( App(nVar*nDofPerElement,nVar*nDofPerElement) )
             ALLOCATE ( TTAmmT(nVar*nNodesPerElement,nVar*nNodesPerElement) )
             ALLOCATE ( TTAmpT(nVar*nNodesPerElement,nVar*nNodesPerElement) )
             ALLOCATE ( TTApmT(nVar*nNodesPerElement,nVar*nNodesPerElement) )
             ALLOCATE ( TTAppT(nVar*nNodesPerElement,nVar*nNodesPerElement) )
             ALLOCATE ( TTAmmTF(nVar*nNodesPerElement), TTAmpTF(nVar*nNodesPerElement) )
             ALLOCATE ( TTApmTF(nVar*nNodesPerElement), TTAppTF(nVar*nNodesPerElement) )
             
             CALL ComputeMt_MSEG(EleM, Mesh%Edg(eleM%Edg), MtM, Data)
             CALL ComputeMb_MSEG(EleM, MbM, Data)             
             CALL InverseMatrix(MtM, MtMm1)             
             TM  = MATMUL(MtMm1, MbM)
             CALL Transpose(TM, TTM)

             CALL ComputeMt_MSEG(EleP, Mesh%Edg(eleP%Edg), MtP, Data)
             CALL ComputeMb_MSEG(EleP, MbP, Data)             
             CALL InverseMatrix(MtP, MtPm1)             
             TP  = MATMUL(MtPm1, MbP)

             CALL Transpose(TP, TTP)
             
             ! Source Terms
             CALL ComputeTransferSourceTerm_MSEG(F2M, eleM, Data)
             CALL ComputeTransferSourceTerm_MSEG(F2P, eleP, Data)
             
             TPFm = MATMUL(MtMm1, F2m)
             TPFp = MATMUL(MtPm1, F2p)             

             CALL ComputePMMatrices_EG(Edg, Mesh%Edg(eleM%Edg), Mesh%Edg(eleP%Edg), &
                  EleM, eleP, Amm, Amp, Apm, App, Data)
             
             TTAmmT = MATMUL(TTM,MATMUL(Amm,TM))
             TTAmpT = MATMUL(TTM,MATMUL(Amp,TP))
             TTApmT = MATMUL(TTP,MATMUL(Apm,TM))
             TTAppT = MATMUL(TTP,MATMUL(App,TP))
             
             TTAmmTF = MATMUL(TTM,MATMUL(Amm,TPFm))
             TTAmpTF = MATMUL(TTM,MATMUL(Amp,TPFp))
             TTApmTF = MATMUL(TTP,MATMUL(Apm,TPFm))
             TTAppTF = MATMUL(TTP,MATMUL(App,TPFp))
             
             CALL AssembleSparseMatrixMS(Mat, TTAmmT, EleM, EleM, Nv, Ne)
             CALL AssembleSparseMatrixMS(Mat, TTAmpT, EleM, EleP, Nv, Ne)
             CALL AssembleSparseMatrixMS(Mat, TTApmT, EleP, EleM, Nv, Ne)
             CALL AssembleSparseMatrixMS(Mat, TTAppT, EleP, EleP, Nv, Ne)
             
             CALL AssembleRHSCG(Sol%RHS, -TTAmmTF, EleM, Nv, Ne, nVar)
             CALL AssembleRHSCG(Sol%RHS, -TTAmpTF, EleM, Nv, Ne, nVar)
             CALL AssembleRHSCG(Sol%RHS, -TTApmTF, EleP, Nv, Ne, nVar)
             CALL AssembleRHSCG(Sol%RHS, -TTAppTF, EleP, Nv, Ne, nVar)
             
             Mesh%Edg(ie)%Treated = .TRUE.

             DEALLOCATE ( MtM, MtMm1, MbM, TM, TTM, MtP, MtPm1, MbP, TP, TTP, F2M, F2P)
             DEALLOCATE ( TPFm, TPFp, Amm, Amp, Apm, App, TTAmmT, TTAmpT, TTApmT)
             DEALLOCATE ( TTAppT, TTAmmTF, TTAmpTF, TTApmTF, TTAppTF )

          END IF
       END IF
              
    END DO
    
24  FORMAT(12(e10.3,3x))
25  FORMAT(9(e10.3,3x))
    
  END SUBROUTINE Compute_LHS_RHS_MSEG
  !====================================================!

  !================================!
  SUBROUTINE CG2EG(Data, Mesh, Sol)
  !================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    !--------------------------------------------
    type(element) :: ele
    !---------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Mt, Mtm1, Mb, T
    REAL*8, DIMENSION(:), ALLOCATABLE :: Ueg, Ucg, F2, TPF, Phi
    !-----------------------------------------------------------
    INTEGER :: ie, i, NU_i, id, iL
    INTEGER :: Ne, Nv, Ndim, nNodesPerElement
    INTEGER :: nDofPerElement, nVar
    !-------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = nDim + 1
    nDofPerElement = nNodesPerElement + 1
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nVar = nDim + 1
    
    DO ie = 1, Ne
       
       ele = Mesh%Ele(ie)
       nNodesPerElement = ele%nNodesPerElement ; nDofPerElement = nNodesPerElement + 1

       ALLOCATE ( Mt(nVar*nDofPerElement,nVar*nDofPerElement) )
       ALLOCATE ( Mb(nVar*nDofPerElement,nVar*nNodesPerElement) )
       ALLOCATE ( Mtm1(nVar*nDofPerElement,nVar*nDofPerElement) )
       ALLOCATE ( T(nVar*nNodesPerElement,nVar*nNodesPerElement) )
       ALLOCATE ( F2(nVar*nDofPerElement) )
       ALLOCATE ( TPF(nVar*nNodesPerElement) )
       ALLOCATE ( Ucg(nVar*nNodesPerElement) )
       ALLOCATE ( Ueg(nVar*nDofPerElement) )
       
       CALL ComputeMt_MSEG(ele, Mesh%Edg(ele%edg), Mt, Data)
       CALL ComputeMb_MSEG(ele, Mb, Data)
       CALL InverseMatrix(Mt, Mtm1)       
       T = MATMUL(Mtm1, Mb)
       
       CALL ComputeTransferSourceTerm_MSEG(F2, ele, Data)
       TPF = MATMUL(Mtm1,F2)
          
       DO i = 1, nNodesPerElement
          NU_i = ele%Vertex(i)
          DO id = 1, nDim
             iL = i + (id-1)*nNodesPerElement
             Ucg(iL) = Sol%Bcg(NU_i,id)
          END DO
          iL = i + nDim*nNodesPerElement
          Ucg(iL) = Sol%pcg(NU_i)
       END DO
       
       Ueg = MATMUL(T,Ucg)
       Ueg = Ueg + TPF
       DO id = 1, nDim
          iL = nDofPerElement + (id-1)*nDofPerElement
          Sol%Beg(ie,id) = Ueg(iL)
       END DO
       iL = nDofPerElement + nDim*nDofPerElement
       Sol%peg(ie) = Ueg(iL)
       
       DO i = 1, nNodesPerElement
          NU_i = ele%Vertex(i)
          DO id = 1, nDim
             iL = i + (id-1)*nDofPerElement
             Sol%Bcg(NU_i,id) = Ueg(iL)
          END DO
          iL = i + nDim*nDofPerElement
          Sol%pcg(NU_i) = Ueg(iL)
       END DO

       DEALLOCATE ( Mt, Mb, Mtm1, T, F2, TPF )
       DEALLOCATE ( Ucg, Ueg )
       
    END DO
    
  END SUBROUTINE CG2EG
  !================================!
  
END MODULE Compute_MSEG
