MODULE DirichletEmbedded_Mixed

  USE Types
  USE PublicVar
  USE Distance
  USE Quadrature
  USE Permeability
  USE ExactSolutiont
  USE ComputeNormalsAndTangents
  USE DirichletTangentialStabilization

  IMPLICIT NONE

CONTAINS

  !==============================================================================!
  SUBROUTINE DirichletEmb_Mixed(Aloc, Floc, ele, edgList, dist, data, iFace, iLS)
  !==============================================================================!

    IMPLICIT NONE
    REAL*8       , DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8       , DIMENSION(:)  , INTENT(OUT) :: Floc
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
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: exa_q, nt, GSF_i, GSF_j, xq, x_q, Bex
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: d_q, tnt, n_q, exFlux
    !--------------------------------------------------------------------------
    REAL*8 :: wq, L, Nrm, h, a, pD_q, SF_i, SF_j, LS_q, nnt, aT
    !--------------------------------------------------------
    INTEGER :: i, j, id, i1, j1, iq, jd, kd, ld
    INTEGER :: nDim, nNodesPerElement, nQuad
    !-----------------------------------------------

    nDim = data%nDim ; nNodesPerElement = ele%nNodesPerElement

    ALLOCATE ( exa_q(nDim+1), nt(nDim), GSF_i(nDim), GSF_j(nDim) )
    ALLOCATE ( xq(nDim), x_q(nDim), d_q(nDim), Lm1_q(nDim,nDim) )
    ALLOCATE ( t_q(nDim,nDim-1), tnt(nDim-1), n_q(nDim) )
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Xe(nNodesPerElement,nDim) )

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, icas)
    END DO
    Xe = ele%Coor

    Aloc = 0.d0 ; Floc = 0.d0

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
       a = Data%Emb_alphaD(iLS)/h*Nrm
       aT = data%Emb_betaT(iLS)*h/Nrm
       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, exa_q, x_q + d_q, Data%T, ele, nDim)
          pD_q = exa_q(nDim+1)
          Bex = exa_q(1:nDim)
       ELSE
          pD_q = Data%BCVal(iLS)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eletype, SF_i, xq, i)
          CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j)
             CALL GSFEval(GSF_j, xq, j, ele)
             ! < w.nt , Lm1 w.d >
             DO id = 1, nDim
                DO jd = 1, nDim
                   DO kd = 1, nDim
                      i1 = i + (id-1)*nNodesPerElement
                      j1 = j + (kd-1)*nNodesPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) + SF_i*SF_j*nt(id)*d_q(jd)*Lm1_q(jd,kd)*wq
                   END DO
                END DO
             END DO
          END DO
          ! - < w.nt , pD >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             Floc(i1) = Floc(i1) - SF_i*nt(id)*pD_q*wq
          END DO
       END DO

       ! < w.nt , 0.5 dT Lm1 G(B) d >
       IF ( data%pressureEnrichment ) THEN
          DO i = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_i, xq, i)
             CALL GSFEval(GSF_i, xq, i, ele)
             DO id = 1, nDim
                DO j = 1, nNodesPerElement
                   CALL GSFEval(GSF_j, xq, j, ele)
                   DO jd = 1, nDim
                      DO kd = 1, nDim
                         DO ld = 1, nDim
                            i1 = i + (id-1)*nNodesPerElement
                            j1 = j + (ld-1)*nNodesPerElement
                            Aloc(i1,j1) = Aloc(i1,j1) + &
                                 SF_i*nt(id)*0.5d0*Lm1_q(jd,ld)*d_q(jd)*GSF_j(kd)*d_q(kd)*wq
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END IF

       CALL Nitsche(Aloc, Floc, ele, d_q, Xq, wq, Lm1_q, nt, pD_q, a, nDim)
       IF ( Data%PressureEnrichment ) THEN
          CALL NitscheHessian(a, Aloc, Floc, xq, wq, Lm1_q, d_q, Lm1, Xe, n_q, ele, pD_q)
          CALL NitscheEnrichment(a, Aloc, Floc, xq, wq, Lm1_q, d_q, Lm1, Xe, n_q, ele, pD_q, data)
       END IF

       ! Tangential terms
       !CALL tangentialStabilization(ele, Aloc, Floc, xq, t_q, Lm1_q, d_q, Bex, aT, wq, nNodesPerElement, nDim)

    END DO

    DEALLOCATE ( exa_q, nt, GSF_i, GSF_j, xq, x_q, d_q, Lm1_q, t_q )

  END SUBROUTINE DirichletEmb_Mixed
  !==============================================================================!

  !========================================================================!
  SUBROUTINE Nitsche(Aloc, Floc, ele, d_q, Xq, wq, Lm1_q, N, pD_q, a, nDim)
  !========================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: Floc
    REAL*8, DIMENSION(:)    , INTENT(IN)    :: d_q, N, Xq
    REAL*8, DIMENSION(:,:)  , INTENT(IN)    :: Lm1_q
    REAL*8                  , INTENT(IN)    :: a, wq, pD_q
    type(element)           , INTENT(IN)    :: ele
    INTEGER                 , INTENT(IN)    :: nDim
    !-------------------------------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: GSF_i, GSF_j
    !-------------------------------------------------
    REAL*8 :: SF_i, SF_j, int
    !---------------------------------
    INTEGER :: i, j, id, jd, kd, i1, j1, ld, i2, ii1, ii2
    INTEGER :: nNodesPerElement, ii, jj, jj1, jj2
    !------------------------------------------------------------

    nNodesPerElement = ele%nNodesPerElement

    ALLOCATE ( GSF_i(nDim), GSF_j(nDim) )

    ! a < q - [Lm1 w].d , p + [Lm1 B].d - pD >
    DO i = 1, nNodesPerElement
       CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
       ! - a < q , pD >
       i1 = i + nDim*nNodesPerElement
       Floc(i1) = Floc(i1) + a*SF_i*pD_q*wq
       ! a <  [Lm1 w].d , pD >
       DO id = 1, nDim
          DO jd = 1, nDim
             i1 = i + (jd-1)*nNodesPerElement
             Floc(i1) = Floc(i1) - a*Lm1_q(id,jd)*d_q(id)*SF_i*pD_q*wq
          END DO
       END DO
       DO j = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
          ! a < q , p >
          i1 = i + nDim*nNodesPerElement ; j1 = j + nDim*nNodesPerElement
          Aloc(i1,j1) = Aloc(i1,j1) + a*SF_i*SF_j*wq
          ! - a < q , [Lm1 B].d >
          DO id = 1, nDim
             DO jd = 1, nDim
                i1 = i + nDim*nNodesPerElement ; j1 = j + (jd-1)*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - a*Lm1_q(id,jd)*d_q(id)*SF_i*SF_j*wq
             END DO
          END DO
          ! - a < [Lm1 w].d , p >
          DO id = 1, nDim
             DO jd = 1, nDim
                i1 = i + (jd-1)*nNodesPerElement ; j1 = j + nDim*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - a*Lm1_q(id,jd)*SF_i*SF_j*d_q(id)*wq
             END DO
          END DO
          ! a < [Lm1 w].d , [Lm1 B].d >
          DO id = 1, nDim
             DO jd = 1, nDim
                DO kd = 1, nDim
                   DO ld = 1, nDim
                      i1 = i + (jd-1)*nNodesPerElement ; j1 = j + (ld-1)*nNodesPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) &
                           + a*d_q(id)*Lm1_q(id,jd)*d_q(kd)*Lm1_q(kd,ld)*SF_i*SF_j*wq
                   END DO
                END DO
             END DO
          END DO
       END DO

    END DO

    DEALLOCATE ( GSF_i, GSF_j )

  END SUBROUTINE Nitsche
  !========================================================================!

  !===============================================================================!
  SUBROUTINE NitscheHessian(a, Aloc, F, xq, wq, Lm1_q, d_q, Lm1, Xe, N, ele, pD_q)
  !===============================================================================!

    IMPLICIT NONE
    type(element)           , INTENT(IN)    :: ele
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    REAL*8, DIMENSION(:,:,:), INTENT(IN)    :: Lm1
    REAL*8, DIMENSION(:,:)  , INTENT(IN)    :: Lm1_q, Xe
    REAL*8, DIMENSION(:)    , INTENT(IN)    :: xq, d_q, N
    REAL*8                  , INTENT(IN)    :: wq, pD_q, a
    !-------------------------------------------------------
    REAL*8 :: SF_i, SF_j, int1, int2
    !----------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: GSF_i, GSF_j
    !-------------------------------------------------
    INTEGER :: i, j, id, jd, kd, ld, i1, j1, md, nd, nDim
    INTEGER :: nNodesPerElement
    !-----------------------------------------------------

    nDim = SIZE(N) ; nNodesPerElement = ele%nNodesPerElement

    ALLOCATE ( GSF_i(nDim), GSF_j(nDim) )

    ! a < - 1/2 dT Lm1 G(w) d , p - [Lm1 B].d - dT Lm1 G(B) d - pD >
    DO i = 1, nNodesPerElement
       CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
       ! a < 1/2 dT Lm1 G(w) d , pD >
       DO id = 1, nDim
          DO jd = 1, nDim
             DO kd = 1, nDim
                i1 = i + (kd-1)*nNodesPerElement
                F(i1) = F(i1) - a*0.5d0*Lm1_q(id,kd)*d_q(id)*d_q(jd)*GSF_i(jd)*pD_q*wq
             END DO
          END DO
       END DO
       DO j = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
          ! a < - 1/2 dT Lm1 G(w) d , p >
          DO id = 1, nDim
             DO jd = 1, nDim
                DO kd = 1, nDim
                   i1 = i + (kd-1)*nNodesPerElement
                   j1 = j + nDim*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) &
                        - a*0.5d0*Lm1_q(id,kd)*d_q(id)*d_q(jd)*GSF_i(jd)*SF_j*wq
                END DO
             END DO
          END DO
          ! a < 1/2 dT Lm1 G(w) d , [Lm1 B].d >
          DO id = 1, nDim
             DO jd = 1, nDim
                DO kd = 1, nDim
                   DO ld = 1, nDim
                      DO md = 1, nDim
                         i1 = i + (kd-1)*nNodesPerElement
                         j1 = j + (md-1)*nNodesPerElement
                         Aloc(i1,j1) = Aloc(i1,j1) &
                              + a*0.5d0*Lm1_q(id,kd)*d_q(id)*d_q(jd)*d_q(ld)*Lm1_q(ld,md)*GSF_i(jd)*SF_j*wq
                      END DO
                   END DO
                END DO
             END DO
          END DO
          ! a < 1/2 dT Lm1 G(w) d , 1/2 dT Lm1 G(B) d >
          DO id = 1, nDim
             DO jd = 1, nDim
                DO kd = 1, nDim
                   DO ld = 1, nDim
                      DO md = 1, nDim
                         DO nd = 1, nDim
                            i1 = i + (kd-1)*nNodesPerElement
                            j1 = j + (nd-1)*nNodesPerElement
                            int1 = 0.5d0*d_q(id)*d_q(jd)*Lm1_q(id,kd)*GSF_i(jd)
                            int2 = 0.5d0*d_q(ld)*d_q(md)*Lm1_q(ld,nd)*GSF_j(md)
                            Aloc(i1,j1) = Aloc(i1,j1) + a*int1*int2*wq
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO

       END DO
    END DO

    ! a < q - [Lm1 w].d , - 1/2 dT Lm1 G(B) d >
    DO i = 1, nNodesPerElement
       CALL SFEval(ele%eleType, SF_i, xq, i)
       DO j = 1, nNodesPerElement
          CALL SFEVal(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
          ! - a < q , 1/2 dT Lm1 G(B) d >
          DO id = 1, nDim
             DO jd = 1, nDim
                DO kd = 1, nDim
                   i1 = i + nDim*nNodesPerElement
                   j1 = j + (kd-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) - &
                        a*SF_i*0.5d0*Lm1_q(id,kd)*d_q(id)*d_q(jd)*GSF_j(jd)*wq
                END DO
             END DO
          END DO
          ! a < Lm1 w.d , 1/2 dT Lm1 G(B) d >
          DO id = 1, nDim
             DO jd = 1, nDim
                DO kd = 1, nDim
                   DO ld = 1, nDim
                      DO md = 1, nDim
                         i1 = i + (jd-1)*nNodesPerElement
                         j1 = j + (md-1)*nNodesPerElement
                         Aloc(i1,j1) = Aloc(i1,j1) &
                              + a*d_q(id)*Lm1_q(id,jd)*SF_i*0.5d0*Lm1_q(kd,md)*d_q(kd)*d_q(ld)*GSF_j(ld)*wq
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE (GSF_i, GSF_j)

  END SUBROUTINE NitscheHessian
  !===============================================================================!

  !==================================================================================!
  SUBROUTINE NitscheEnrichment(a, Aloc, F, xq, wq, Lm1_q, d_q, Lm1, X, N, ele, pD_q, data)
  !==================================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    REAL*8, DIMENSION(:)    , INTENT(IN)    :: d_q, N, xq
    REAL*8, DIMENSION(:,:,:), INTENT(IN)    :: Lm1
    REAL*8, DIMENSION(:,:)  , INTENT(IN)    :: Lm1_q
    REAL*8, DIMENSION(:,:)  , INTENT(IN)    :: X
    REAL*8                  , INTENT(IN)    :: wq, a, pD_q
    type(element)           , INTENT(IN)    :: ele
    type(dataStructure)     , INTENT(IN)    :: data
    !-----------------------------------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: GSF_i, GSF_j
    !-------------------------------------------------
    REAL*8 :: SF_i, SF_j, SF_k, int, SF_i1, SF_j1, SF_i2, SF_j2
    !----------------------------------------------------------
    INTEGER :: i, j, k, id, jd, kd, i1, j1, k1, ii, ld, md, nd
    INTEGER :: nDim, nNodesPerElement, nEdgesPerElement
    INTEGER :: i2, j2, ii1, ii2, jj1, jj2, jj, iid, jjd
    !----------------------------------------------------------

    nDim = SIZE(N) ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement
    ALLOCATE ( GSF_i(nDim), GSF_j(nDim) )

    ! a < q - [Lm1 w].d - 1/2 dT Lm1 G(w) d , p* >
    DO i = 1, nNodesPerElement
       CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
       DO ii = 1, nEdgesPerElement
          j = Permut(ii,1) ; k = Permut(ii,2)
          CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)
          DO id = 1, nDim
             j1 = j + (id-1)*nNodesPerElement ; k1 = k + (id-1)*nNodesPerElement
             ! a < q , p* >
             i1 = i + nDim*nNodesPerElement
             int = a*SF_i*SF_j*SF_k*wq
             DO jd = 1, nDim
                Aloc(i1,j1) = Aloc(i1,j1) + 0.5d0*Lm1(j,jd,id)*(X(j,jd)-X(k,jd))*int
                Aloc(i1,k1) = Aloc(i1,k1) - 0.5d0*Lm1(k,jd,id)*(X(j,jd)-X(k,jd))*int
             END DO
             ! - a < [Lm1 w].d , p* >
             DO jd = 1, nDim
                DO kd = 1, nDim
                   DO ld = 1, nDim
                      i1 = i + (jd-1)*nNodesPerElement
                      int = a*d_q(kd)*Lm1_q(kd,jd)*SF_i*SF_j*SF_k*wq
                      Aloc(i1,j1) = Aloc(i1,j1) - 0.5d0*Lm1(j,ld,id)*(X(j,ld)-X(k,ld))*int
                      Aloc(i1,k1) = Aloc(i1,k1) + 0.5d0*Lm1(k,ld,id)*(X(j,ld)-X(k,ld))*int
                   END DO
                END DO
             END DO
             ! - a < 1/2 dT Lm1 G(w) d , p* >
             DO jd = 1, nDim
                DO kd = 1, nDim
                   DO ld = 1, nDim
                      DO md = 1, nDim
                         i1 = i + (ld-1)*nNodesPerElement
                         int = a*0.5d0*d_q(jd)*d_q(kd)*Lm1_q(jd,ld)*GSF_i(kd)*SF_j*SF_k*wq
                         Aloc(i1,j1) = Aloc(i1,j1) - 0.5d0*Lm1(j,md,id)*(X(j,md)-X(k,md))*int
                         Aloc(i1,k1) = Aloc(i1,k1) + 0.5d0*Lm1(k,md,id)*(X(j,md)-X(k,md))*int
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    ! a < q* , p + p* - [Lm1 B].d - 1/2 dT [Lm1(GB)] d - pD>
    IF ( data%testFunctionEnrichment ) THEN
       DO ii = 1, nEdgesPerElement
          ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
          CALL SFEval(ele%eleType, SF_i1, xq, ii1)
          CALL SFEval(ele%eleType, SF_i2, xq, ii2)
          DO id = 1, nDim
                i1 = ii1 + (id-1)*nNodesPerElement
                i2 = ii2 + (id-1)*nNodesPerElement
                ! a < q* , pD >
                int = 0.5d0*a*SF_i1*SF_i2*pD_q*wq
                DO iid = 1, nDim
                   F(i1) = F(i1) + Lm1(ii1,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                   F(i2) = F(i2) - Lm1(ii2,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                END DO
                !F(i1) = F(i1) + (X(ii1,id) - X(ii2,id))*int
                !F(i2) = F(i2) - (X(ii1,id) - X(ii2,id))*int
                DO j = 1, nNodesPerElement
                   CALL SFEval(ele%eleType, SF_j, xq, j)
                   CALL GSFEval(GSF_j, xq, j, ele)
                   ! a < q* , p >
                   j1 = j + nDim*nNodesPerElement
                   int = 0.5d0*a*SF_i1*SF_i2*SF_j*wq
                   DO iid = 1, nDim
                      Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                      Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                   END DO
                   DO jd = 1, nDim
                      ! - a < q* , [Lm1 B].d >
                      DO kd = 1, nDim
                         j1 = j + (jd-1)*nNodesPerElement
                         int = -a*0.5d0*SF_i1*SF_i2*Lm1_q(kd,jd)*d_q(kd)*SF_j*wq
                         DO iid = 1, nDim
                            Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                            Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                         END DO
                      END DO
                         ! - a < q* , 0.5 dT [Lm1 GB] d >
                      DO kd = 1, nDim
                         DO ld = 1, nDim
                            j1 = j + (ld-1)*nNodesPerElement
                            int = -a*0.25d0*SF_i1*SF_i2*d_q(kd)*d_q(jd)*Lm1_q(jd,ld)*GSF_j(kd)*wq
                            DO iid = 1, nDim
                               Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                               Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,iid,id)*(X(ii1,iid) - X(ii2,iid))*int
                            END DO
                         END DO
                      END DO
                   END DO
                END DO
                ! a < q* , p* >
                DO jj = 1, nEdgesPerElement
                   jj1 = Permut(jj,1) ; jj2 = Permut(jj,2)
                   CALL SFEval(ele%eleType, SF_j1, xq, jj1)
                   CALL SFEval(ele%eleType, SF_j2, xq, jj2)
                   DO jd = 1, nDim
                         j1 = jj1 + (jd-1)*nNodesPerElement
                         j2 = jj2 + (jd-1)*nNodesPerElement
                         int = 0.25d0*a*SF_i1*SF_i2*SF_j1*SF_j2*wq
                         DO iid = 1, nDim
                            DO jjd = 1, nDim
                         Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,iid,id)*(X(ii1,iid) - X(ii2,iid))*&
                                                     Lm1(jj1,jjd,jd)*(X(jj1,jjd) - X(jj2,jjd))*int
                         Aloc(i1,j2) = Aloc(i1,j2) - Lm1(ii1,iid,id)*(X(ii1,iid) - X(ii2,iid))*&
                                                     Lm1(jj2,jjd,jd)*(X(jj1,jjd) - X(jj2,jjd))*int
                         Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,iid,id)*(X(ii1,iid) - X(ii2,iid))*&
                                                     Lm1(jj1,jjd,jd)*(X(jj1,jjd) - X(jj2,jjd))*int
                         Aloc(i2,j2) = Aloc(i2,j2) + Lm1(ii2,iid,id)*(X(ii1,iid) - X(ii2,iid))*&
                                                     Lm1(jj2,jjd,jd)*(X(jj1,jjd) - X(jj2,jjd))*int
                      END DO
                   END DO
                END DO
             END DO
          END DO

       END DO
    END IF

    DEALLOCATE ( GSF_i, GSF_j)

  END SUBROUTINE NitscheEnrichment
  !===================================================================================!

END MODULE DirichletEmbedded_Mixed
