MODULE NeumannEmbedded_Mixed

  USE Types
  USE PublicVar
  USE Distance
  USE Quadrature
  USE Permeability
  USE ExactSolutiont
  USE ComputeNormalsAndTangents

  IMPLICIT NONE

CONTAINS

  !=========================================================================!
  SUBROUTINE NeumannEmb_Mixed(Aloc, F, ele, edgList, dist, data, iFace, iLS)
  !=========================================================================!

    IMPLICIT NONE
    REAL*8       , DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8       , DIMENSION(:)  , INTENT(OUT) :: F
    type(element)                , INTENT(IN)  :: ele
    type(edge)   , DIMENSION(:)  , INTENT(IN)  :: edgList
    REAL*8       , DIMENSION(:,:), INTENT(IN)  :: dist
    type(DataStructure)          , INTENT(IN)  :: data
    INTEGER                      , INTENT(IN)  :: iLS, iFace
    !-------------------------------------------------------
    type(edge) :: edg
    !-----------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, t_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: exa_q, nt, GSF_i, GSF_j, xq, x_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: d_q, tnt, n_q
    !--------------------------------------------------------------------------
    REAL*8 :: wq, L, Nrm, h, a, hN_q, SF_i, SF_j, LS_q, nnt, SF_k
    REAL*8 :: int, SF_i1, SF_i2, T
    !-------------------------------------------------------------
    INTEGER :: i, j, id, i1, j1, iq, jd, kd, ld
    INTEGER :: nDim, nNodesPerElement, nQuad, ii
    INTEGER :: nEdgesPerElement, k, k1, ii1, ii2, i2, iid
    !-------------------------------------------------

    nDim = data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; T=Data%T

    ALLOCATE ( exa_q(nDim+1), nt(nDim), GSF_i(nDim), GSF_j(nDim) )
    ALLOCATE ( xq(nDim), x_q(nDim), d_q(nDim), Lm1_q(nDim,nDim) )
    ALLOCATE ( t_q(nDim,nDim-1), tnt(nDim-1), n_q(nDim) )
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Xe(nNodesPerElement,nDim) )

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, icas)
    END DO
    Xe = ele%Coor

    Aloc = 0.d0 ; F = 0.d0

    nt = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + nt(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm ; h = ele%A/L
    nt = nt/Nrm

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

       CALL computeNodalDistance(x_q, d_q, LS_q, data, iLS)
       CALL getNormalAndTangents(d_q, n_q, t_q, nt, nnt, tnt, nDim)
       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, icas)

       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)
       a = Data%alphaD(iLS)/h*Nrm

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, exa_q, x_q + d_q,T, ele, nDim)
          hN_q = DOT_PRODUCT(exa_q(1:nDim),n_q)
       ELSE
          hN_q = Data%BCEmbVal(iLS)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          !- < q n.nt , hN >
          i1 = i + nDim*nNodesPerElement
          F(i1) = F(i1) - SF_i*hN_q*wq*nnt
          DO id = 1, nDim
             ! < w.nt , p >
             DO j = 1, nNodesPerElement
                CALL SFEval(ele%eleType, SF_j, xq, j)
                i1 = i + (id-1)*nNodesPerElement ; j1 = j + nDim*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + SF_i*SF_j*nt(id)*wq
             END DO
             ! < w.nt , p* >
             IF ( data%pressureEnrichment ) THEN
                DO ii = 1, nEdgesPerElement
                   j = permut(ii,1) ; k = permut(ii,2)
                   CALL SFEval(ele%eleType, SF_j, xq, j)
                   CALL SFEval(ele%eleType, SF_k, xq, k)
                   DO jd = 1, nDim
                      DO kd = 1, nDim
                         i1 = i + (id-1)*nNodesPerElement
                         j1 = j + (kd-1)*nNodesPerElement
                         k1 = k + (kd-1)*nNodesPerElement
                         int = 0.5d0*SF_i*SF_j*SF_k*nt(id)*wq
                         Aloc(i1,j1) = Aloc(i1,j1) + Lm1(j,jd,kd)*(Xe(j,jd) - Xe(k,jd))*int
                         Aloc(i1,k1) = Aloc(i1,k1) - Lm1(k,jd,kd)*(Xe(j,jd) - Xe(k,jd))*int
                      END DO
                   END DO
                END DO
             END  IF
             ! - < q n.nt , B.n + [(GB)d].n >
             ! - < q n.nt , B.n >
             DO j = 1, nNodesPerElement
                CALL SFEval(ele%eleType, SF_j, xq, j)
                CALL GSFEval(GSF_j, xq, j, ele)
                i1 = i + nDim*nNodesPerElement
                j1 = j + (id-1)*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - SF_i*SF_j*n_q(id)*wq*nnt
                ! - < q n.nt , [(GB)d].n >
                DO jd = 1, nDim
                   i1 = i + nDim*nNodesPerElement
                   j1 = j + (id-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) - SF_i*n_q(id)*GSF_j(jd)*d_q(jd)*wq*nnt
                END DO
             END DO
          END DO
       END DO

       IF ( data%PressureEnrichment .AND. data%testFunctionEnrichment ) THEN
          DO ii = 1, nEdgesPerElement
             ii1 = permut(ii,1) ; ii2 = permut(ii,2)
             CALL SFEval(ele%eleType, SF_i1, xq, ii1)
             CALL SFEval(ele%eleType, SF_i2, xq, ii2)
             DO id = 1, nDim
                DO iid = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   ! - < q* n.nt , hN >
                   int = -0.5d0*SF_i1*SF_i2*hN_q*nnt*wq
                   F(i1) = F(i1) + Lm1(ii1,iid,id)*(Xe(ii1,iid) - Xe(ii2,iid))*int
                   F(i2) = F(i2) - Lm1(ii2,iid,id)*(Xe(ii1,iid) - Xe(ii2,iid))*int
                   DO j = 1, nNodesPerElement
                      CALL SFEval(ele%eleType, SF_j, xq, j)
                      CALL GSFEval(GSF_j, xq, j, ele)
                      DO jd = 1, nDim
                         j1 = j + (jd-1)*nNodesPerElement
                         ! - < q* n.nt , B.n >
                         int = -0.5d0*SF_i1*SF_i2*nnt*SF_j*n_q(jd)*wq
                         Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,iid,id)*(Xe(ii1,iid) - Xe(ii2,iid))*int
                         Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,iid,id)*(Xe(ii1,iid) - Xe(ii2,iid))*int
                         DO kd = 1, nDim
                            ! - < q* n.nt , [(GB)d].n >
                            int = - 0.5d0*SF_i1*SF_i2*nnt*d_q(kd)*GSF_j(kd)*n_q(jd)*wq
                            Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,iid,id)*(Xe(ii1,iid) - Xe(ii2,iid))*int
                            Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,iid,id)*(Xe(ii1,iid) - Xe(ii2,iid))*int
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END IF

    END DO

    DEALLOCATE ( exa_q, nt, GSF_i, GSF_j, xq, x_q, d_q, Lm1_q )

  END SUBROUTINE NeumannEmb_Mixed
  !==============================================================================!


END MODULE NeumannEmbedded_Mixed
