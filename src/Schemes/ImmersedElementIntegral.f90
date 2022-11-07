MODULE ImmersedElementIntegral

  USE Types
  USE PublicVar
  USE Quadrature
  USE ExactSolutiont
  USE Distance

  IMPLICIT NONE

CONTAINS

  SUBROUTINE IntegrateImmersedElement_CG_P(Aloc, Uloc, Floc, ele, LS, dist, Bloc, data)

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(IN)  :: Uloc
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: Floc
    type(element)         , INTENT(IN)  :: ele
    type(DataStructure)   , INTENT(IN)  :: data
    REAL*8, DIMENSION(:)  , INTENT(IN)  :: LS
    REAL*8, DIMENSION(:,:), INTENT(IN)  :: dist, Bloc
    !-----------------------------------------------
    REAL*8, DIMENSION(data%nDim)   :: xq, x_q, GSF_j, d_q
    REAL*8, DIMENSION(data%nDim+1) :: exa
    !-----------------------------------------------------
    REAL*8 :: wq, SF_i, SF_j, LS_q, LS_j, eta
    REAL*8 :: uD_j
    !------------------------------------------
    INTEGER :: nDim, nNodesPerElement, nQuad
    INTEGER :: iq, i, j, id, jd
    !--------------------------------------------

    Aloc = 0.d0 ; Floc = 0.d0
    nDim = data%nDim ; nNodesPerElement = ele%nNodesPerElement

    eta = 2.d0/ele%A

    CALL GetNQuad(ele%eleType, nQuad)

    DO iq = 1, nQuad

       CALL CoorQuad(ele%eleType, xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       x_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          x_q = x_q + SF_i*ele%Coor(i,:)
       END DO
       CALL ComputeNodalDistance(x_q, d_q, LS_q, data, 1)

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO j = 1, nNodesPerelement
             CALL SFEval(ele%eleType, SF_j, xq, j)
             IF ( LS(j) <= 0.d0 ) THEN
                CALL ExactSolt(ele%L, exa, ele%Coor(j,:),Data%T, ele, nDim)
             ELSE
                CALL ExactSolt(ele%L, exa, ele%Coor(j,:) + dist(j,:),Data%T, ele, nDim)
                !CALL ExactSol(exa, ele%Coor(j,:) , ele, nDim)
             END IF
             uD_j = exa(nDim+1)
             ! eta ( q , u - ud ) = eta sum (uj -uDj) ( qi , phi_j)
             Aloc(i,j) = Aloc(i,j) + eta*SF_i*SF_j*wq
             Floc(i) = Floc(i) + eta*SF_i*SF_j*uD_j*wq

             IF ( LS(j) > 0.d0 ) THEN
                ! eta ( qi , G phi_j . d ) u_j
                CALL GSFEval(GSF_j, xq, j, ele)
                DO id = 1, nDim
                   !Aloc(i,j) = Aloc(i,j) + eta*SF_i*GSF_j(id)*d_q(id)*wq
                   !Floc(i) = Floc(i) + eta*SF_i*Bloc(j,id)*d_q(id)*wq
                   !IF ( iq == 1 .AND. i == 1 .AND. id == 1 ) PRINT*, SF_i, GSF_j, d_q, wq
                END DO
             END IF

          END DO
       END DO

    END DO

  END SUBROUTINE IntegrateImmersedElement_CG_P



END MODULE ImmersedElementIntegral
