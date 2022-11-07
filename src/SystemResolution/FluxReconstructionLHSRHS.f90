MODULE FluxReconstructionLHSRHS

  USE Types
  USE PublicVar
  USE SourceTerm
  USE Quadrature
  USE Permeability
  USE ExactSolutiont

  IMPLICIT NONE

CONTAINS

  !=========================================================!
  SUBROUTINE ElementIntegral_GR(Aloc, Floc, pLoc, ele, Data)
  !=========================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(INOUT) :: Floc
    REAL*8, DIMENSION(:)  , INTENT(IN)    :: pLoc
    type(element)         , INTENT(IN)    :: ele
    type(DataStructure)   , INTENT(IN)    :: Data
    !----------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Lm1_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: X_q, Xq, GSF_i, GSF_j, phi_q
    !--------------------------------------------------------------------
    REAL*8 :: SF_i, SF_j, wq, Nrm
    !-----------------------------
    INTEGER :: iq, i, j, i1, j1, id, jd, nDim, nQuad
    INTEGER :: nNodesPerElement
    !-------------------------------------------------

    nDim = data%nDim ; nNodesPerElement = ele%nNodesPerElement

    ALLOCATE ( Lm1_q(nDim,nDim), phi_q(nDim+1) )
    ALLOCATE ( X_q(nDim), xq(nDim), GSF_i(nDim), GSF_j(nDim) )

    CALL GetNQuad(ele%eleType, nQuad)
    Aloc = 0.d0 ; Floc = 0.d0

    DO iq = 1, nQuad

       CALL CoorQuad(ele%eleType, Xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       X_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             X_q(id) = X_q(id) + SF_i*ele%Coor(i,id)
          END DO
       END DO

       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, Data%icas)

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             DO id = 1, nDim
                ! ( w , Lm1 B )
                DO jd = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement
                   j1 = j + (jd-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + Lm1_q(id,jd)*SF_i*SF_j*wq
                END DO
             END DO
             DO id = 1, nDim
                ! + ( div w , p )
                i1 = i + (id-1)*nNodesPerElement
                Floc(i1) = Floc(i1) + GSF_i(id)*SF_j*pLoc(j)*wq
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE ElementIntegral_GR
  !=========================================================!

  !====================================================================!
  SUBROUTINE RHS_edgInt(edg, ele, Adj, Floc, pLocE, pLocA, Data, iFace)
  !====================================================================!

    IMPLICIT NONE
    type(edge)          , INTENT(IN)    :: edg
    type(element)       , INTENT(IN)    :: ele, Adj
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Floc
    REAL*8, DIMENSION(:), INTENT(IN)    :: pLocE, pLocA
    type(DataStructure) , INTENT(IN)    :: Data
    INTEGER             , INTENT(IN)    :: iFace
    !-----------------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: xq_E, xq_A
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: N, Xq, xqE, xqA
    !----------------------------------------------------------------
    REAL*8 :: Nrm, L, h, p_qE, p_qA, SF_iE, SF_iA, wq
    !----------------------------------------------------
    INTEGER :: nQuad, nNodesPerElement, nDim
    INTEGER :: iq, id, jd, i, j, i1, j1, iFadj
    !------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement

    ALLOCATE ( N(nDim), xq_E(3,nDim), xq_A(3,nDim), xqE(nDim), xqA(nDim) )

    N = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm
    N = N/Nrm
    h = ele%A/L

    CALL IdentifyIFace(Adj, edg, iFadj, nDim)

    CALL QuadratureOnEdge(xq_E, ele, edg, iFace, nDim)
    CALL QuadratureOnEdge(xq_A, Adj, edg, iFadj, nDim)

    CALL GetNQuad(edg%eleType, nQuad)

    DO iq = 1, nQuad

       xqE = xq_E(iq,:) ; xqA = xq_A(iq,:)
       CALL WeightQuad(edg%eleType, wq, iq)
       wq = wq*L

       p_qE = 0.d0 ; p_qA = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_iE, xqE, i) ; CALL SFEval(Adj%eleType, SF_iA, xqA, i)
          p_qE = p_qE + pLocE(i)*SF_iE
          p_qA = p_qA + pLocA(i)*SF_iA
       END DO

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_iE, xqE, i)
          CALL SFEval(Adj%eleType, SF_iA, xqA, i)

          DO id = 1, nDim
             ! - < w.n , p >
             i1 = i + (id-1)*nNodesPerElement
             Floc(i1) = Floc(i1) - SF_iE*N(id)*0.5d0*(p_qE + p_qA)*wq
          END DO
       END DO

    END DO

  END SUBROUTINE RHS_edgInt
  !=====================================================================!

  !===========================================================!
  SUBROUTINE RHS_edgNeumann(edg, ele, Floc, pLoc, Data, iFace)
  !===========================================================!

    IMPLICIT NONE
    type(edge)          , INTENT(IN)    :: edg
    type(element)       , INTENT(IN)    :: ele
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Floc
    REAL*8, DIMENSION(:), INTENT(IN)    :: pLoc
    type(DataStructure) , INTENT(IN)    :: Data
    INTEGER             , INTENT(IN)    :: iFace
    !--------------------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: N, Xq
    !-------------------------------------------
    REAL*8 :: Nrm, L, h, p_q, SF_i, wq
    !-------------------------------------
    INTEGER :: nQuad, nNodesPerElement, nDim
    INTEGER :: iq, id, jd, i, j, i1, j1, iFadj
    !------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement

    ALLOCATE ( N(nDim),  xq(nDim) )

    CALL GetNQuad(ele%eleType, nQuad)

    N = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm

    CALL GetNQuad(edg%eleType, nQuad)

    DO iq = 1, nQuad

       CALL FaceCoorQuad(edg%eleType, xq, iq, iFace)
       CALL WeightQuad(edg%eleType, wq, iq)
       wq = wq*L

       p_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          p_q = p_q + pLoc(i)*SF_i
       END DO

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             Floc(i1) = Floc(i1) - SF_i*N(id)*p_q*wq
          END DO
       END DO

    END DO

  END SUBROUTINE RHS_edgNeumann
  !====================================================!

  !============================================================!
  SUBROUTINE RHS_edgDirichlet(edg, ele, Floc, pLoc, Data, iFace, iBC)
  !============================================================!

    IMPLICIT NONE
    type(edge)          , INTENT(IN)    :: edg
    type(element)       , INTENT(IN)    :: ele
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Floc, pLoc
    type(DataStructure) , INTENT(IN)    :: Data
    INTEGER             , INTENT(IN)    :: iFace, iBC
    !-------------------------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: N, Xq, x_q, Exa_q
    !------------------------------------------------------
    REAL*8 :: Nrm, L, h, p_q, pD_q, SF_i, wq
    !-----------------------------------------
    INTEGER :: nQuad, nNodesPerElement, nDim
    INTEGER :: iq, id, jd, i, j, i1, j1, iFadj
    !------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement

    ALLOCATE ( N(nDim),  xq(nDim), x_q(nDim), Exa_q(nDim+1) )

    CALL GetNQuad(ele%eleType, nQuad)

    N = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm

    CALL GetNQuad(edg%eleType, nQuad)

    DO iq = 1, nQuad

       CALL FaceCoorQuad(edg%eleType, xq, iq, iFace)
       CALL WeightQuad(edg%eleType, wq, iq)
       wq = wq*L

       x_q = 0.d0 ; p_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             x_q(id) = x_q(id) + ele%coor(i,id)*SF_i
          END DO
          p_q = p_q + pLoc(i)*SF_i
       END DO

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
          pD_q = Exa_q(nDim+1)
       ELSE
          pD_q = Data%BCVal(iBC)
       END IF

       PRINT*, pD_q, p_q

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)

          DO id = 1, nDim
             ! - < w.n , p >
             i1 = i + (id-1)*nNodesPerElement
             Floc(i1) = Floc(i1) - SF_i*N(id)*p_q*wq
          END DO
       END DO

    END DO

  END SUBROUTINE RHS_edgDirichlet
  !======================================================!

END MODULE FluxReconstructionLHSRHS
