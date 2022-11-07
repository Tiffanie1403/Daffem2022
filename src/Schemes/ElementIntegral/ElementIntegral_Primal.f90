MODULE ElementIntegral_Primal

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE SourceTerm

  IMPLICIT NONE

CONTAINS

  !=========================================================================!
  SUBROUTINE Compute_ElementIntegral_Primal(ele, Aloc, Uloc, Floc, Data)
  !=========================================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)    :: ele
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(INOUT) :: Floc
    REAL*8, DIMENSION(:)  , INTENT(IN)    :: Uloc
    type(DataStructure)   , INTENT(IN)    :: Data
    !--------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: L_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: GSF_i, GSF_j, Xq, X_q, phi_q
    !--------------------------------------------------------------------
    REAL*8 :: SF_i, SF_j
    REAL*8 :: wq
    !--------------------
    INTEGER :: nNodesPerElement, nDim, nQuad
    INTEGER :: i, iq, j, id, jd, nVar
    !----------------------------------------
    CHARACTER(len=5) :: ET
    !------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    ET = ele%eleType ; nVar = 1

    ALLOCATE ( L_q(nDim,nDim), GSF_i(nDim), GSF_j(nDim) )
    ALLOCATE ( Xq(nDim), X_q(nDim), phi_q(nDim+1) )

    Aloc = 0.d0 ; Floc = 0.d0

    CALL GetNQuad(ele%eleType, nQuad)

    DO iq = 1, nQuad

       CALL CoorQuad(ele%eleType, Xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       X_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ET, SF_i, Xq, i)
          DO id = 1, nDim
             X_q(id) = X_q(id) + ele%Coor(i,id)*SF_i
          END DO
       END DO

       CALL SetNodalMatrixPermeability(L_q, ele, icas)

       DO i = 1, nNodesPerElement
          CALL SFEval(ET, SF_i, Xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ET, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             ! ( Gq , L Gu )
             DO id = 1, nDim
                DO jd = 1, nDim
                   ! LHS
                   Aloc(i,j) = Aloc(i,j) + GSF_i(id)*L_q(id,jd)*GSF_j(jd)*wq
                   ! RHS
                   Floc(i) = Floc(i) - GSF_i(id)*L_q(id,jd)*GSF_j(jd)*Uloc(j)*wq
                END DO
             END DO
          END DO

          ! ( q , phi )
          CALL Compute_SourceTerm(phi_q, X_q,Data%T, ele, Data%icas, nDim)          
          Floc(i) = Floc(i) + SF_i*phi_q(nDim+1)*wq
       END DO

    END DO

    DEALLOCATE ( L_q, GSF_i, GSF_j, Xq, X_q, phi_q )

  END SUBROUTINE Compute_ElementIntegral_Primal
  !======================================================================!

END MODULE ElementIntegral_Primal
