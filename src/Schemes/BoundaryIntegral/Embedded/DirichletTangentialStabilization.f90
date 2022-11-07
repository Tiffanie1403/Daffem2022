MODULE DirichletTangentialStabilization

  USE Types
  USE Quadrature
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE tangentialStabilization(ele, Aloc, Floc, xq, t_q, Lm1_q, d_q, Bex, aT, wq, nNodesPerElement, nDim)

    IMPLICIT NONE
    type(element)         , INTENT(IN)    :: ele
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: ALoc
    REAL*8, DIMENSION(:)  , INTENT(INOUT) :: Floc
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: t_q, Lm1_q
    REAL*8, DIMENSION(:)  , INTENT(IN)    :: d_q, Bex, xq    
    REAL*8                , INTENT(IN)    :: aT, wq
    INTEGER               , INTENT(IN)    :: nNodesPerElement, nDim 
    !-----------------------------------------------------
    REAL*8, DIMENSION(nDim) :: GSF_i, GSF_j, t
    !---------------------------------------
    REAL*8 :: SF_i, SF_j
    !--------------------
    INTEGER :: i, j, i1, j1, id, jd, kd, ld, md, nd
    INTEGER :: itangent, jtangent
    !-------------------------------------------------

    ! < ( w + [Gw]d ).t , aT (Lm1 B + [G(Lm1 B)]d).t) + dt pD >
    DO itangent = 1, nDim-1
       t = t_q(:,itangent)
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          ! aT < w.t , B.t >
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             DO id = 1, nDim
                DO jd = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement ; j1 = j + (jd-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + aT*SF_i*t(id)*t(jd)*SF_j*wq
                END DO
             END DO
          END DO
          ! aT < w.t , Bex.t >
          DO id = 1, nDim
             DO jd = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                Floc(i1) = Floc(i1) + aT*SF_i*t(id)*Bex(jd)*t(jd)*wq
             END DO
          END DO
          ! aT < w.t , ([(G B)d].t) >
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             DO id = 1, nDim
                DO jd = 1, nDim
                   DO kd = 1, nDim
                      i1 = i + (id-1)*nNodesPerElement
                      j1 = j + (jd-1)*nNodesPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) + aT*SF_i*t(id)*t(jd)*d_q(kd)*GSF_j(kd)*wq
                   END DO
                END DO
             END DO
          END DO

          ! aT < [(Gw)d].t , Bex.t >
          DO id = 1, nDim
             DO jd = 1, nDim
                DO kd = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement 
                   Floc(i1) = Floc(i1) + aT*t(id)*d_q(jd)*GSF_i(jd)*Bex(kd)*t(kd)*wq
                END DO
             END DO
          END DO

          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             DO id = 1, nDim
                DO jd = 1, nDim
                   DO kd = 1, nDim
                      ! aT < [(Gw)d].t , B.t >
                      i1 = i + (id-1)*nNodesPerelement ; j1 = j + (kd-1)*nNodesPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) + aT*t(id)*d_q(jd)*GSF_i(jd)*SF_j*t(kd)*wq
                      ! aT < [(Gw)d].t , [(GB)d].t >
                      DO ld = 1, nDim
                         i1 = i + (id-1)*nNodesPerElement ; j1 = j + (kd-1)*nNodesPerElement
                         Aloc(i1,j1) = Aloc(i1,j1) + &
                              aT*t(id)*d_q(jd)*GSF_i(jd)*t(kd)*d_q(ld)*GSF_j(ld)*wq
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO
    
  END SUBROUTINE tangentialStabilization

END MODULE DirichletTangentialStabilization
