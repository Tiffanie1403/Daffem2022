MODULE DirichletEmbedded_Primal

  USE Types
  USE Distance
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE ExactSolutiont
  USE ComputeNormalsAndTangents

  IMPLICIT NONE

CONTAINS

  !====================================================================================!
  SUBROUTINE DirichletEmb_P2(Aloc, Uloc, F, ele, edg, dist, data, iFace, iLS)
  !====================================================================================!

    IMPLICIT NONE
    REAL*8       , DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8       , DIMENSION(:)  , INTENT(OUT) :: F
    REAL*8       , DIMENSION(:)  , INTENT(IN)  :: Uloc
    type(element)                , INTENT(IN)  :: ele
    type(edge)                   , INTENT(IN)  :: edg
    REAL*8       , DIMENSION(:,:), INTENT(IN)  :: dist
    type(DataStructure)          , INTENT(IN)  :: data
    INTEGER                      , INTENT(IN)  :: iLS, iFace
    !-------------------------------------------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: L_q, t_q, Xe, GG_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: exa_q, nt, GSF_i, GSF_j, xq, x_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: d_q, tnt, n_q, G_q
    !--------------------------------------------------------------------------
    REAL*8 :: wq, L, Nrm, h, a, hN_q, SF_i, SF_j, LS_q, nnt, SF_k
    REAL*8 :: int, SF_i1, SF_i2, pD_q
    !-------------------------------------------------------------
    INTEGER :: i, j, id, i1, j1, iq, jd, kd, ld
    INTEGER :: nDim, nNodesPerElement, nQuad, ii
    INTEGER :: nEdgesPerElement, k, k1, ii1, ii2, i2, iid
    INTEGER :: itangent
    !-------------------------------------------------


    PRINT*, " Je suis passe ici 1"
    nDim = data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement
!allocation dynamique des tableaux
    ALLOCATE ( exa_q(nDim+1), nt(nDim), GSF_i(nDim), GSF_j(nDim) )
    ALLOCATE ( xq(nDim), x_q(nDim), d_q(nDim), L_q(nDim,nDim) )
    ALLOCATE ( t_q(nDim,nDim-1), tnt(nDim-1), n_q(nDim) )
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Xe(nNodesPerElement,nDim) )
    ALLOCATE ( G_q(nDim), GG_q(nDim,nDim) )
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
       CALL getNormalAndTangents(d_q, n_q, t_q, nt, nnt, tnt, nDim) !t_q la tangente tableau
       CALL SetNodalMatrixPermeability(L_q, ele, icas)

       !d_q = 0.d0
       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)
       !a = Data%Emb_alphaD(iLS)/h*Nrm
       a = Data%Emb_alphaD(iLS)/h

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, exa_q, x_q + d_q, Data%T, ele, nDim)
          pD_q = exa_q(nDim+1)
       ELSE
          pD_q = Data%BCVal(iLS)
       END IF




       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)


             ! - < q , [L Gu].n >
             DO id = 1, nDim
                DO jd = 1, nDim
                   Aloc(i,j) = Aloc(i,j) - SF_i * L_q(jd,id) *GSF_j(id)* nt(jd)*wq
                END DO
             END DO

             !terme ajouté
             ! - < Gq.d , [L Gu].n >
             DO id = 1, nDim
                DO jd = 1, nDim
                   DO kd = 1, nDim
                      Aloc(i,j) = Aloc(i,j) - SF_i * L_q(jd,id) *GSF_j(id)* nt(jd) *GSF_j(kd)*d_q(kd)*wq
                   END DO
                END DO
             END DO

             DO id = 1, nDim
                DO jd = 1, nDim
                   ! - < [L Gq].n , u >
                   Aloc(i,j) = Aloc(i,j) - L_q(jd,id)*GSF_i(jd)*nt(id)*SF_j*wq
                   ! - < [L Gq].n , Gu.d >
                   DO kd = 1, nDim
                      int = L_q(id,jd)*GSF_i(jd)*nt(id)*GSF_j(kd)*d_q(kd)*wq
                      Aloc(i,j) = Aloc(i,j) - int
                   END DO
                END DO
             END DO
             ! a < q , u >
             Aloc(i,j) = Aloc(i,j) + a*SF_i*SF_j*wq

             DO id = 1, nDim
                ! a < q , Gu.d >
                Aloc(i,j) = Aloc(i,j) + a*SF_i*GSF_j(id)*d_q(id)*wq
                ! a < Gq.d , u >
                Aloc(i,j) = Aloc(i,j) + a*GSF_i(id)*d_q(id)*SF_j*wq
                DO jd = 1, nDim
                   ! a < Gq.d , Gu.d >
                   Aloc(i,j) = Aloc(i,j) + a*GSF_i(id)*d_q(id)*GSF_j(jd)*d_q(jd)*wq
                END DO
             END DO
          END DO
          ! a < q , uD >
          F(i) = F(i) + a*SF_i*pD_q*wq
          ! a < Gq.d , ud >
          DO id = 1, nDim
             F(i) = F(i) + a*GSF_i(id)*d_q(id)*pD_q*wq
          END DO
          ! + < [L Gq].n , uD >
          DO id = 1, nDim
             DO jd = 1, nDim
                F(i) = F(i) - L_q(jd,id)*GSF_i(id)*nt(jd)*pD_q*wq
             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE DirichletEmb_P2
  !============================================================================!


  !====================================================================================!
  SUBROUTINE DirichletEmb_P(Aloc, Uloc, F, ele, edgList, dist, data, iFace, iLS)
  !====================================================================================!

    IMPLICIT NONE
    REAL*8       , DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8       , DIMENSION(:)  , INTENT(OUT) :: F
    REAL*8       , DIMENSION(:)  , INTENT(IN)  :: Uloc
    type(element)                , INTENT(IN)  :: ele
    type(edge)   , DIMENSION(:)  , INTENT(IN)  :: edgList
    REAL*8       , DIMENSION(:,:), INTENT(IN)  :: dist
    type(DataStructure)          , INTENT(IN)  :: data
    INTEGER                      , INTENT(IN)  :: iLS, iFace
    !-------------------------------------------------------
    type(edge) :: edg
    !-----------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: L_q, t_q, Xe, GG_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: exa_q, nt, GSF_i, GSF_j, xq, x_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: d_q, tnt, n_q, G_q
    !--------------------------------------------------------------------------
    REAL*8 :: wq, L, Nrm, h, a, hN_q, SF_i, SF_j, LS_q, nnt, SF_k
    REAL*8 :: int, SF_i1, SF_i2, pD_q
    !-------------------------------------------------------------
    INTEGER :: i, j, id, i1, j1, iq, jd, kd, ld
    INTEGER :: nDim, nNodesPerElement, nQuad, ii
    INTEGER :: nEdgesPerElement, k, k1, ii1, ii2, i2, iid
    INTEGER :: itangent
    !-------------------------------------------------
    PRINT*, " Je suis passe ici 2"

    nDim = data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement

    ALLOCATE ( exa_q(nDim+1), nt(nDim), GSF_i(nDim), GSF_j(nDim) )
    ALLOCATE ( xq(nDim), x_q(nDim), d_q(nDim), L_q(nDim,nDim) )
    ALLOCATE ( t_q(nDim,nDim-1), tnt(nDim-1), n_q(nDim) )
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Xe(nNodesPerElement,nDim) )
    ALLOCATE ( G_q(nDim), GG_q(nDim,nDim) )
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
       !d_q = 0.d0
       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)
       !a = Data%Emb_alphaD(iLS)/h*Nrm
       a = Data%Emb_alphaD(iLS)/h

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, exa_q, x_q + d_q, Data%T, ele, nDim)
          pD_q = exa_q(nDim+1)
       ELSE
          pD_q = Data%BCVal(iLS)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)


             !terme ajoute
             ! - < q , [L Gu].n >
             DO id = 1, nDim
                DO jd = 1, nDim
                   Aloc(i,j) = Aloc(i,j) - SF_i*L_q(jd,id)*GSF_j(id)*nt(jd)*wq
                END DO
             END DO
             ! - < Gq.d , [L Gu].n >
             DO id = 1, nDim
                DO jd = 1, nDim
                   DO kd = 1, nDim
                      !DO ld = 1, nDim
                      Aloc(i,j) = Aloc(i,j) - SF_i * L_q(jd,id) *GSF_j(id)* nt(jd) *GSF_j(kd)*d_q(kd)
                      !END DO
                   END DO
                END DO
             END DO
             !fin terme ajoute



             DO id = 1, nDim
                DO jd = 1, nDim
                   ! - < [L Gq].n , u >
                   Aloc(i,j) = Aloc(i,j) - L_q(jd,id)*GSF_i(jd)*nt(id)*SF_j*wq
                   ! - < [L Gq].n , Gu.d >
                   DO kd = 1, nDim
                      int = L_q(id,jd)*GSF_i(jd)*nt(id)*GSF_j(kd)*d_q(kd)*wq
                      Aloc(i,j) = Aloc(i,j) - int
                   END DO
                END DO
             END DO
             ! a < q , u >
             Aloc(i,j) = Aloc(i,j) + a*SF_i*SF_j*wq

             DO id = 1, nDim
                ! a < q , Gu.d >
                Aloc(i,j) = Aloc(i,j) + a*SF_i*GSF_j(id)*d_q(id)*wq
                ! a < Gq.d , u >
                Aloc(i,j) = Aloc(i,j) + a*GSF_i(id)*d_q(id)*SF_j*wq
                DO jd = 1, nDim
                   ! a < Gq.d , Gu.d >
                   Aloc(i,j) = Aloc(i,j) + a*GSF_i(id)*d_q(id)*GSF_j(jd)*d_q(jd)*wq
                END DO
             END DO
          END DO
          ! a < q , uD >
          F(i) = F(i) + a*SF_i*pD_q*wq
          ! a < Gq.d , ud >
          DO id = 1, nDim
             F(i) = F(i) + a*GSF_i(id)*d_q(id)*pD_q*wq
          END DO
          ! + < [L Gq].n , uD >
          DO id = 1, nDim
             DO jd = 1, nDim
                F(i) = F(i) - L_q(jd,id)*GSF_i(id)*nt(jd)*pD_q*wq
             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE DirichletEmb_P
  !============================================================================!



END MODULE DirichletEmbedded_Primal
