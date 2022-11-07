MODULE ElementIntegral_EG

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE SourceTerm

  IMPLICIT NONE

CONTAINS

  !===========================================================!
  SUBROUTINE Compute_ElementIntegral_EG(ele, Aloc, Floc, Data)
  !===========================================================!

    IMPLICIT NONE
    type(DataStructure)   , INTENT(IN)  :: Data
    type(element)         , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: Floc
    !-------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Lm1_q, L_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: X_q,  phi_q, Xq, GSF_i, GSF_j
    !--------------------------------------------------------------------
    REAL*8 :: tauU, tauP, int, h
    REAL*8 :: alphaS, wq, SF_i, SF_j, mu_q
    !--------------------------------------------
    INTEGER :: i, j, id, k, m, i1, jd, kd, ld
    INTEGER :: k1, m1, j1, it, iq
    INTEGER :: nDim, nNodesPerElement
    INTEGER :: nDofPerElement, nQuad, nVar
    !------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement ; nVar = 1
    nDofPerElement = nNodesPerElement + 1

    ALLOCATE ( Lm1_q(nDim,nDim), L_q(nDim,nDim) )
    ALLOCATE ( X_q(nDim), phi_q(nDim+1), Xq(nDim), GSF_i(nDim), GSF_j(nDim) )

    alphaS = Data%alphaS
    h = ele%A**(1.d0/nDim)

    Aloc = 0.d0 ; Floc = 0.d0

    CALL GetNQuad(ele%eleType, nQuad)

    DO iq = 1, nQuad

       CALL CoorQuad(ele%eleType, Xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       X_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType,SF_i, xq, i)
          DO id = 1, nDim
             X_q(id) = X_q(id) + SF_i*ele%Coor(i,id)
          END DO
       END DO
       CALL SetNodalMatrixPermeability(L_q, ele, data%icas)
       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, Data%icas)
       mu_q = 0.d0
       DO id = 1, nDim
          mu_q = mu_q + ele%L(id,id)**2
       END DO
       mu_q = SQRT(mu_q)
       tauU = 0.5d0*mu_q
       tauP = alphaS*h*h/mu_q

       DO i = 1, nDofPerElement
          CALL SFEval(ele%eleType,SF_i, xq, i)
          CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nDofPerElement
             CALL SFEval(ele%eleType,SF_j, xq, j)
             CALL GSFEval(GSF_j, xq, j, ele)
             ! ( w , Lm1 B )
             DO id = 1, nDim
                DO jd = 1, nDim
                   i1 = i + (id-1)*nDofPerElement
                   j1 = j + (id-1)*nDofPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + Lm1_q(id,jd)*SF_i*SF_j*wq
                END DO
             END DO
             ! - ( div w , p )
             DO id = 1, nDim
                i1 = i + (id-1)*nDofPerElement
                j1 = j + nDim*nDofPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - GSF_i(id)*SF_j*wq
             END DO
             ! ( q , div B )
             DO id = 1, nDim
                i1 = i + nDim*nDofPerElement
                j1 = j + (id-1)*nDofPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + SF_i*GSF_j(id)*wq
             END DO

             ! Stabilization terms
             IF ( data%space_scheme == "MSEG-Stab" ) THEN
                DO id = 1, nDim
                   DO jd = 1, nDim
                      ! - 1/2 (Lm1 w , B)
                      i1 = i + (id-1)*nDofPerElement
                      j1 = j + (jd-1)*nDofPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) - 0.5d0*Lm1_q(jd,id)*SF_i*SF_j*wq
                      ! - 1/2 (Lm1 w , L Gp)
                      DO kd = 1, nDim
                         i1 = i + (id-1)*nDofPerElement
                         j1 = j + nDim*nDofPerElement
                         Aloc(i1,j1) = Aloc(i1,j1) - 0.5d0*Lm1_q(kd,id)*SF_i*L_q(kd,jd)*GSF_j(jd)*wq
                      END DO
                      ! 0.5d0 ( Gq , L Gp)
                      i1 = i + nDim*nDofPerElement
                      j1 = j + nDim*nDofPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) + 0.5d0*GSF_i(id)*L_q(id,jd)*GSF_j(jd)*wq
                   END DO
                   ! 0.5 ( Gq , B )
                   i1 = i + nDim*nDofPerElement
                   j1 = j + (id-1)*nDofPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + 0.5d0*GSF_i(id)*SF_j*wq

                END DO
             END IF
          END DO

          ! ( q , phi )
          CALL Compute_SourceTerm(phi_q, X_q,Data%T, ele, Data%icas, nDim)
          i1 = i + nDim*nDofPerElement
          Floc(i1) = Floc(i1) + SF_i*phi_q(nDim+1)*wq

       END DO

    END DO

    DEALLOCATE ( Lm1_q, L_q, X_q, phi_q, Xq, GSF_i, GSF_j )

  END SUBROUTINE Compute_ElementIntegral_EG
  !===============================================================!

END MODULE ElementIntegral_EG
