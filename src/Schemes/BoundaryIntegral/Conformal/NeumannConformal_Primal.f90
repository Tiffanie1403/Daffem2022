MODULE NeumannConformal_Primal

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactSolutiont

  IMPLICIT NONE

CONTAINS

  !===========================================================================!
  SUBROUTINE NeumannConf_Primal(Aloc, Uloc, F, ele, edgList, Data, iFace, iBC)
  !===========================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(OUT) :: F
    REAL*8, DIMENSION(:)    , INTENT(IN)  :: Uloc
    type(element)           , INTENT(IN)  :: ele
    type(edge), DIMENSION(:), INTENT(IN)  :: edgList
    type(DataStructure)     , INTENT(IN)  :: Data
    INTEGER                 , INTENT(IN)  :: iBC, iFace
    !-------------------------------------------------
    type(Edge) :: edg
    !------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: L_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Exa_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: N, GSF_i, GSF_j, xq, x_q
    !----------------------------------------------------------------
    REAL*8 :: L, wq, hN_q
    REAL*8 :: SF_i, SF_j, Nrm, h, uD_q, a
    !---------------------------------------------
    INTEGER :: i, j, iq, id, jd, nQuad
    INTEGER :: nDim, nNodesPerElement, nVar
    !---------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement ; nVar = 1
    ALLOCATE ( L_q(nDim,nDim), Exa_q(nDim+1), N(nDim), GSF_i(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )

    Aloc = 0.d0 ; F = 0.d0

    N = - ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm)
    N = N/Nrm
    L = Nrm
    h = ele%A/L

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

       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)
       a = Data%alphaD(iBC)/h*Nrm

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
          hN_q = DOT_PRODUCT(Exa_q(1:nDim),N)
       ELSE
          hN_q = Data%BCVal(iBC)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          ! - < q , hN >
          F(i) = F(i) - SF_i*hN_q*wq
       END DO
    END DO

    DEALLOCATE ( L_q, Exa_q, N, GSF_i, GSF_j, xq )

  END SUBROUTINE NeumannConf_Primal
  !===========================================================================!

END MODULE NeumannConformal_Primal
