MODULE NeumannConformal_EG

  USE Types
  USE PublicVar
  USE ExactSolutiont
  USE Permeability
  USE Quadrature

  IMPLICIT NONE

CONTAINS

  !=================================================================!
  SUBROUTINE NeumannConf_EG(Aloc, F, Ele, edgList, Data, iFace, iBC)
  !=================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    type(element)           , INTENT(IN)    :: Ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !------------------------------------------------------
    type(edge) :: edg
    !------------------
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Exa_q
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Lm1_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: N, GSF_i, GSF_j, Xq, X_q
    !---------------------------------------------------------------
    REAL*8 :: L, hN_q, Nrm, a, h, wq, SF_i, SF_j, mu_q
    !---------------------------------------------------
    INTEGER :: i, id, i1, ii, j, jd, j1, iq, nQuad, nVar
    INTEGER :: nDim ,nNodesPerElement, nDofPerElement
    !-----------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = nDim + 1 ; nVar = nDim + 1
    nDofPerElement = nNodesPerElement + 1

    ALLOCATE ( Exa_q(nVar), Lm1_q(nDim,nDim) )
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

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, Exa_q, x_q,Data%T, ele, nDim)
          hN_q = DOT_PRODUCT(Exa_q(1:nDim),N)
       ELSE
          hN_q = Data%BCVal(iBC)
       END IF

       DO i = 1, nDofPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO j = 1, nDofPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j)
             DO id = 1, nDim
                ! < w.n , p >
                i1 = i + (id-1)*nDofPerElement
                j1 = j + nDim*nDofPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + SF_i*SF_j*N(id)*wq
                ! - < q , B.n >
                i1 = i + nDim*nDofPerElement
                j1 = j + (id-1)*nDofPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - SF_i*SF_j*N(id)*wq
             END DO
          END DO
          ! < q , hN >
          i1 = i + nDim*nDofPerElement
          F(i1) = F(i1) - SF_i*hN_q*wq
       END DO

    END DO

    DEALLOCATE ( Exa_q, Lm1_q, N, GSF_i, GSF_j, Xq, X_q )

  END SUBROUTINE NeumannConf_EG
  !=======================================================!

END MODULE NeumannConformal_EG
