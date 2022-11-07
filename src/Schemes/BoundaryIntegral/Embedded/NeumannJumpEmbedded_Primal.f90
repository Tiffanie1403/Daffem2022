MODULE NeumannJumpEmbedded_Primal

  USE Types
  USE Distance
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE ExactFunctionJumpt
  USE ComputeNormalsAndTangents

  IMPLICIT NONE

CONTAINS

  !===========================================================================!
  SUBROUTINE NeumannJumpEmb_P(Aloc, Uloc, F, ele1, ele2, edg, iFace1, iFace2, iLS, dist, data )
  !===========================================================================!

    IMPLICIT NONE
    REAL*8       , DIMENSION(:,:), INTENT(OUT)   :: Aloc
    REAL*8       , DIMENSION(:)  , INTENT(OUT)   :: F
    REAL*8       , DIMENSION(:)  , INTENT(IN)    :: Uloc
    type(element)                , INTENT(INOUT) :: ele1, ele2
    type(edge)                   , INTENT(IN)    :: edg
    REAL*8       , DIMENSION(:,:), INTENT(IN)    :: dist
    type(DataStructure)          , INTENT(IN)    :: data
    INTEGER                      , INTENT(IN)    :: iLS, iFace1, iFace2
    !--------------------------------------------------------------------
    type(Element) :: ele
    !---------------------
    REAL*8, DIMENSION(data%nDim,data%nDim)   :: L_q
    REAL*8, DIMENSION(data%nDim,data%nDim-1) :: t_q
    !----------------------------------------------------------
    REAL*8, DIMENSION(data%nDim+1) :: Exa1_q, exa2_q
    REAL*8, DIMENSION(data%nDim)   :: nt, GSF_i, GSF_j, xq, x_q, d_q, n_q
    REAL*8, DIMENSION(data%nDim-1) :: tnt
    !--------------------------------------------------------------------------------
    REAL*8 :: L, wq, hN_q, LS_q, nnt
    REAL*8 :: SF_i, SF_j, Nrm, h, uD_q, a, l1, l2
    !---------------------------------------------
    INTEGER :: i, j, iq, id, jd, nQuad, tmpFace, iFace
    INTEGER :: nDim, nNodesPerElement, nVar
    !---------------------------------------

    nDim = data%nDim ; nNodesPerElement = ele1%nNodesPerElement

    Aloc = 0.d0 ; F = 0.d0

    l1 = ele1%L(1,1) ; l2 = ele2%L(2,2)
    IF ( l1 < l2 ) THEN
       ele = ele2 ; iFace = iFace2
    ELSE
       ele = ele1 ; iFace = iFace1
    END IF

    nt = ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + nt(id)**2
    END DO
    Nrm = sqrt(Nrm)
    nt = nt/Nrm
    L = Nrm
    h = ele%A/L

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

       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)

       IF ( icas /= 0 ) THEN
          CALL ExactSolJumpt(ele1%L, exa1_q, x_q + d_q, Data%T, ele, nDim)
          CALL ExactSolJumpt(ele2%L, exa2_q, x_q + d_q, Data%T, ele, nDim)
          hN_q = -DOT_PRODUCT(ele2%L(1,1)*exa2_q(1:nDim) - ele2%L(1,1)*exa1_q(1:nDim),n_q)
          PRINT*, hN_q
       ELSE
          hN_q = Data%BCVal(iLS)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          ! - < q , hN >
          !PRINT*, l2 - l1
          F(i) = F(i) + SF_i*hN_q*wq
          !PRINT*, SF_i, ele1%vertex(i), edg%vertex
          !PRINT*, F(i)
       END DO
    END DO

    ele1 = ele



  END SUBROUTINE NeumannJumpEmb_P
  !======================================================================!

END MODULE NeumannJumpEmbedded_Primal
