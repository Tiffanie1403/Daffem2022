MODULE DirichletConformal_EG

  USE Types
  USE PublicVar
  USE ExactSolutiont
  USE Permeability
  USE Quadrature

  IMPLICIT NONE

CONTAINS

  !==================================================================!
  SUBROUTINE DirichletConf_EG(Aloc, F, Ele, EdgList, Data, iFace, iBC)
  !==================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    type(element)           , INTENT(IN)    :: Ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !-------------------------------------------------------
    type(edge) :: edg
    !------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Lm1_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, Xq, X_q
    !-----------------------------------------------------------------------
    REAL*8 :: L, Tau, pD_q, Nrm, a, h, bT
    REAL*8 :: SF_i, SF_j, mu_q, wq
    !---------------------------------------------------
    INTEGER :: i,id, i1, ii, j, k, jd, kd, ld, j1, k1, iq, nVar
    INTEGER :: nDim ,nNodesPerElement, nDofPerElement, nQuad
    !------------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = nDim + 1 ; nVar = nDim+1
    nDofPerElement = nNodesPerElement + 1

    ALLOCATE ( Exa_q(nDim+1), Lm1_q(nDim,nDim) )
    ALLOCATE ( N(nDim), GSF_i(nDim), GSF_j(nDim), Xq(nDim), X_q(nDim) )
    Aloc = 0.d0 ; F = 0.d0

    N = -Ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm)
    L = Nrm
    N = N/Nrm
    h = Ele%A/L
    edg = edgList(iFace)
    CALL GetNQuad(edg%eleType, nQuad)
    DO iq = 1, nQuad
       CALL FaceCoorQuad(edg%eleType, xq, iq, iFace)
       CALL WeightQuad(edg%eleType, wq, iq)
       wq = wq*L
       x_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             x_q(id) = x_q(id) + SF_i*ele%Coor(i,id)
          END DO
       END DO
       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, icas)

       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)
       a = Data%alphaD(iBC)/h*Nrm
       bT = Data%betaT(iBC)*h/Nrm

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, Exa_q, x_q,Data%T,ele, nDim)
          pD_q = Exa_q(nVar)
       ELSE
          pD_q = Data%BCVal(iBC)
       END IF

       ! - < w.nt , pD >
       DO i = 1, nDofPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             i1 = i + (id-1)*nDofPerElement
             F(i1) = F(i1) - N(id)*SF_i*pD_q*wq
          END DO
       END DO

       ! a < q , p  - pD >
       DO i = 1, nDofPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          i1 = i + nDim*nDofPerElement
          ! a < q , p >
          DO j = 1, nDofPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j)
             j1 = j + nDim*nDofPerElement
             Aloc(i1,j1) = Aloc(i1,j1) + a*SF_i*SF_j*wq
          END DO
          ! a < q , pD >
          F(i1) = F(i1) + a*SF_i*pD_q*wq
       END DO

    END DO

    DEALLOCATE ( Exa_q, Lm1_q, N, GSF_i, GSF_j, Xq )

  END SUBROUTINE DirichletConf_EG
  !=======================================================!

END MODULE DirichletConformal_EG
