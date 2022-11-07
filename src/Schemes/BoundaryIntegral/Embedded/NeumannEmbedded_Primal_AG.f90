MODULE NeumannEmbedded_Primal_AG

  USE Types
  USE Distance
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE ExactSolutiont
  USE ComputeNormalsAndTangents

  IMPLICIT NONE

CONTAINS

  !=======================================================================================!
  SUBROUTINE NeumannEmb_P_AG(Aloc, Uloc, F, ele, edgList, BRecons, dist, data, iFace, iLS)
  !=======================================================================================!

    IMPLICIT NONE
    REAL*8       , DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8       , DIMENSION(:)  , INTENT(OUT) :: F
    REAL*8       , DIMENSION(:)  , INTENT(IN)  :: Uloc
    type(element)                , INTENT(IN)  :: ele
    type(edge)   , DIMENSION(:)  , INTENT(IN)  :: edgList
    REAL*8       , DIMENSION(:,:), INTENT(IN)  :: dist
    REAL*8       , DIMENSION(:,:), INTENT(IN)  :: BRecons
    type(DataStructure)          , INTENT(IN)  :: data
    INTEGER                      , INTENT(IN)  :: iLS, iFace
    !-------------------------------------------------------
    type(edge) :: edg
    !-----------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: L_q, t_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: exa_q, nt, GSF_i, GSF_j, xq, x_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: d_q, tnt, n_q, G_q
    !--------------------------------------------------------------------------
    REAL*8 :: wq, L, Nrm, h, a, hN_q, SF_i, SF_j, LS_q, nnt, SF_k
    REAL*8 :: int, SF_i1, SF_i2
    !-------------------------------------------------------------
    INTEGER :: i, j, id, i1, j1, iq, jd, kd, ld
    INTEGER :: nDim, nNodesPerElement, nQuad, ii
    INTEGER :: nEdgesPerElement, k, k1, ii1, ii2, i2, iid
    INTEGER :: itangent
    !-------------------------------------------------

    nDim = data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement

    ALLOCATE ( exa_q(nDim+1), nt(nDim), GSF_i(nDim), GSF_j(nDim) )
    ALLOCATE ( xq(nDim), x_q(nDim), d_q(nDim), L_q(nDim,nDim) )
    ALLOCATE ( t_q(nDim,nDim-1), tnt(nDim-1), n_q(nDim) )
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Xe(nNodesPerElement,nDim) )
    ALLOCATE ( G_q(nDim) )
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
       CALL SetNodalMatrixPermeability(L_q, ele, icas)

       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, exa_q, x_q + d_q, Data%T, ele, nDim)
          hN_q = -DOT_PRODUCT(exa_q(1:nDim),n_q)
       ELSE
          hN_q = Data%BCVal(iLS)
       END IF

       ! AG SBM type
!!$       DO i = 1, nNodesPerElement
!!$          CALL SFEval(ele%eleType, SF_i, xq, i)
!!$          ! < q(n.nt) , hN >
!!$          F(i) = F(i) + SF_i*nnt*hN_q*wq
!!$          DO j = 1, nNodesPerElement
!!$             ! - < q(t.nt) , [L_q Gu].t >
!!$             CALL GSFEval(GSF_j, xq, j, ele)
!!$             DO id = 1, nDim - 1
!!$                DO jd = 1, nDim
!!$                   DO kd = 1, nDim
!!$                      int = SF_i*tnt(id)*L_q(kd,jd)*GSF_j(jd)*t_q(kd,id)*wq
!!$                      ! LHS
!!$                      Aloc(i,j) = Aloc(i,j) - int
!!$                   END DO
!!$                END DO
!!$             END DO
!!$             ! - < q(n.nt) , [G(Gu)d].n >
!!$             ! = < q(n.nt) , [GB]d.n >
!!$             DO id = 1, nDim
!!$                DO jd = 1, nDim
!!$                   ! Only RHS because lagged in Newton iterations!
!!$                   F(i) = F(i) + SF_i*nnt*GSF_j(jd)*d_q(jd)*n_q(id)*BRecons(j,id)*wq
!!$                END DO
!!$             END DO
!!$
!!$          END DO
!!$       END DO

       ! L test
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          ! < q , hN >
          F(i) = F(i) + SF_i*hN_q*wq
          DO j = 1, nNodesPerElement
             ! - < q , Gu.nt >
             CALL SFEval(ele%eleType, SF_j, xq, j)
             CALL GSFEval(GSF_j, xq, j, ele)
             DO id = 1, nDim
                !Aloc(i,j) = Aloc(i,j) - SF_i*GSF_j(id)*nt(id)*wq
                F(i) = F(i) + SF_i*GSF_j(id)*Uloc(j)*nt(id)*wq
             END DO
             ! - < q , B.n >
             DO id = 1, nDim
                F(i) = F(i) + SF_i*SF_j*BRecons(j,id)*n_q(id)*wq
             END DO
             ! < q , Gu.n >
!!$             DO id = 1, nDim
!!$                !Aloc(i,j) = Aloc(i,j) - SF_i*GSF_j(id)*n_q(id)*wq
!!$                F(i) = F(i) - SF_i*GSF_j(id)*Uloc(j)*n_q(id)*wq
!!$             END DO
             ! - < q , ([GB]d).n >
             DO id = 1, nDim
                DO jd = 1, nDim
                   F(i) = F(i) + SF_i*GSF_j(jd)*d_q(jd)*BRecons(j,id)*n_q(id)*wq
                END DO
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE NeumannEmb_P_AG
  !================================================================================!

END MODULE NeumannEmbedded_Primal_AG
