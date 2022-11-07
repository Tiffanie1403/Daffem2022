MODULE ElementIntegral_Mixed

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE SourceTerm

  IMPLICIT NONE

CONTAINS

  !==============================================================!
  SUBROUTINE Compute_ElementIntegral_Mixed(ele, Aloc, Floc, Data)
  !==============================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: Floc
    type(DataStructure)   , INTENT(IN)  :: Data
    !----------------------------------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, Xe, L_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_i, GSF_j, GSF_k
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_i1, GSF_i2, GSF_j1, GSF_j2
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Xq, X_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: phi_q
    !------------------------------------------------------------
    REAL*8 :: SF_i, SF_j, h, SF_k
    REAL*8 :: SF_i1, SF_i2, SF_j1, SF_j2
    REAL*8 :: wq, mu_q
    REAL*8 :: tauU, tauP, int
    !---------------------------------
    INTEGER :: nNodesPerElement, nDim, nQuad
    INTEGER :: i, iq, j, id, jd, i1, i2, j1, j2, kd, md, nd
    INTEGER :: nVar, ii, k, k1, ld, nEdgesPerElement
    INTEGER :: ii1, ii2, jj, jj1, jj2
    !-----------------------------------------------

    nDim = Data%nDim ; nVar = nDim + 1 ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Lm1_q(nDim,nDim), L_q(nDim,nDim) )
    ALLOCATE ( GSF_i(nDim), GSF_j(nDim), GSF_k(nDim), Xq(nDim), X_q(nDim) )
    ALLOCATE ( phi_q(nVar), Xe(nNodesPerElement,nDim) )
    ALLOCATE ( GSF_i1(nDim), GSF_i2(nDim), GSF_j1(nDim), GSF_j2(nDim) )

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, Data%icas)
    END DO

    Xe = ele%Coor

    Aloc = 0.d0 ; Floc = 0.d0
    h = ele%A**(1./nDim)

    CALL GetNQuad(ele%eleType, nQuad)

    DO iq = 1, nQuad

       CALL CoorQuad(ele%eleType, Xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       X_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, Xq, i)
          DO id = 1, nDim
             X_q(id) = X_q(id) + ele%Coor(i,id)*SF_i
          END DO
       END DO

       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, Data%icas)
       CALL SetNodalMatrixPermeability(L_q, ele, data%icas)
       mu_q = 0.d0
       DO id = 1, nDim
          mu_q = ele%L(id,id)**2
       END DO
       mu_q = sqrt(mu_q)
       tauU = 0.5d0
       tauP = Data%alphaS*h*h/mu_q*0.5d0

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             ! ( w , Lm1 B )
             DO id = 1, nDim
                DO jd = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement
                   j1 = j + (jd-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + Lm1_q(id,jd)*SF_i*SF_j*wq
                END DO
             END DO
             ! - ( div w , p )
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - GSF_i(id)*SF_j*wq
             END DO
             ! ( q , div B )
             DO id = 1, nDim
                i1 = i + nDim*nNodesPerElement
                j1 = j + (id-1)*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + SF_i*GSF_j(id)*wq
             END DO
          END DO

          ! ( q , phi )
          CALL Compute_SourceTerm(phi_q, X_q,Data%T, ele, Data%icas, nDim)
          !PRINT*, phi_q(nDim+1)
          i1 = i + nDim*nNodesPerElement
          Floc(i1) = Floc(i1) + SF_i*phi_q(nDim+1)*wq

       END DO

       ! ( q* , phi )
       IF ( data%testFunctionEnrichment ) THEN
          DO ii = 1, nEdgesPerElement
             ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
             CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)
             DO id = 1, nDim !decalement pour l'indice dands la matrice
                DO jd = 1, nDim !vecteur ej,k
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   int = 0.5d0*SF_i1*SF_i2*phi_q(nDim+1)*wq
                   Floc(i1) = Floc(i1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                   Floc(i2) = Floc(i2) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                END DO
             END DO
          END DO
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEVal(ele%eleType, SF_j, xq, j)
             CALL GSFEval(GSF_j, xq, j, ele)
             DO id = 1, nDim
                DO jd = 1, nDim
                   ! tauU ( -Lm1 w , B )
                   i1 = i + (id-1)*nNodesPerElement
                   j1 = j + (jd-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) - tauU*Lm1_q(jd,id)*SF_i*SF_j*wq
                   ! tau U ( -Lm1 w , L Gp )
                   DO kd = 1, nDim
                      i1 = i + (id-1)*nNodesPerElement
                      j1 = j + nDim*nNodesPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) - tauU*Lm1_q(kd,id)*SF_i*L_q(kd,jd)*GSF_j(jd)*wq
                   END DO
                   ! tauU ( Gq , L Gp )
                   i1 = i + nDim*nNodesPerElement
                   j1 = j + nDim*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + tauU*GSF_i(id)*L_q(id,jd)*GSF_j(jd)*wq
                END DO
                ! tauU ( Gq , B )
                i1 = i + nDim*nNodesPerElement
                j1 = j + (id-1)*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + tauU*GSF_i(id)*SF_j*wq
             END DO

             ! tauP ( div w , div u )
             DO id = 1, nDim
                DO jd = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement
                   j1 = j + (jd-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + tauP*GSF_i(id)*GSF_j(jd)*wq
                END DO
             END DO

          END DO
          ! tauP ( div w , phi )
          CALL Compute_SourceTerm(phi_q, X_q,Data%T, ele, Data%icas, nDim)
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             Floc(i1) = Floc(i1) + tauP*GSF_i(id)*phi_q(nDim+1)*wq
          END DO
       END DO

       ! Pressure enrichment components
       IF ( Data%PressureEnrichment ) THEN
          ! Test function enrichment
          IF ( Data%TestFunctionEnrichment ) THEN
             DO ii = 1, nEdgesPerElement
                ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
                CALL SFEval(ele%eleType, SF_i1, xq, ii1)
                CALL SFEval(ele%eleType, SF_i2, xq, ii2)
                CALL GSFEval(GSF_i1, xq, ii1, ele)
                CALL GSFEval(GSF_i2, xq, ii2, ele)
                DO j = 1, nNodesPerElement
                   CALL SFEval(ele%eleType, SF_j, xq, j)
                   CALL GSFEval(GSF_j, xq, j, ele)
                   ! ( q* , div B )
                   DO id = 1, nDim
                      DO jd = 1, nDim
                         DO kd = 1, nDim
                            i1 = ii1 + (id-1)*nNodesPerElement
                            i2 = ii2 + (id-1)*nNodesPerElement
                            j1 = j + (kd-1)*nNodesPerElement
                            int = 0.5d0*SF_i1*SF_i2*GSF_j(kd)*wq
                            Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                            Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                         END DO
                      END DO
                   END DO

                   ! tauU ( Gq* , B )
                   DO id = 1, nDim
                      DO jd = 1, nDim
                         DO kd = 1, nDim
                            i1 = ii1 + (id-1)*nNodesPerElement
                            i2 = ii2 + (id-1)*nNodesPerElement
                            j1 = j + (kd-1)*nNodesPerElement
                            int = tauU*0.5d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*SF_j*wq
                            Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                            Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                         END DO
                      END DO
                   END DO

                   ! tauU ( Gq* , L Gp )
                   DO id = 1, nDim
                      DO jd = 1, nDim
                         DO kd = 1, nDim
                            DO ld = 1, nDim
                               i1 = ii1 + (id-1)*nNodesPerElement
                               i2 = ii2 + (id-1)*nNodesPerElement
                               j1 = j + nDim*nNodesPerElement
                               int = tauU*0.5d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*L_q(kd,ld)*GSF_j(ld)*wq
                               Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                               Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                            END DO
                         END DO
                      END DO
                   END DO
                END DO
                ! tauU( Gq* , L Gp* )
                DO jj = 1, nEdgesPerElement
                   jj1 = Permut(jj,1) ; jj2 = Permut(jj,2)
                   CALL SFEval(ele%eleType, SF_j1, xq, jj1)
                   CALL SFEval(ele%eleType, SF_j2, xq, jj2)
                   CALL GSFEval(GSF_j1, xq, jj1, ele)
                   CALL GSFEval(GSF_j2, xq, jj2, ele)
                   DO id = 1, nDim
                      DO jd = 1, nDim
                         DO kd = 1, nDim
                            DO ld = 1, nDim
                               DO md = 1, nDim
                                  DO nd = 1, nDim
                                     i1 = ii1 + (id-1)*nNodesPerElement
                                     i2 = ii2 + (id-1)*nNodesPerElement
                                     j1 = jj1 + (kd-1)*nNodesPerElement
                                     j2 = jj2 + (kd-1)*nNodesPerElement
                                     int = tauU*0.25d0*(SF_i1*GSF_i2(md) + GSF_i1(md)*SF_i2)*&
                                          L_q(md,nd)*(SF_j1*GSF_j2(nd) + GSF_j1(nd)*SF_j2)*wq
                                     Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                          Lm1(jj1,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                                     Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                          Lm1(jj1,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                                     Aloc(i1,j2) = Aloc(i1,j2) - Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                          Lm1(jj2,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                                     Aloc(i2,j2) = Aloc(i2,j2) + Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                          Lm1(jj2,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                                  END DO
                               END DO
                            END DO
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END IF

          DO i = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_i, xq, i)
             CALL GSFEval(GSF_i, xq, i, ele)
             DO ii = 1, nEdgesPerElement
                j = Permut(ii,1) ; k = Permut(ii,2)
                CALL SFEval(ele%eleType, SF_j, xq, j)
                CALL SFEval(ele%eleType, SF_k, xq, k)
                CALL GSFEval(GSF_j, xq, j, ele)
                CALL GSFEval(GSF_k, xq, k, ele)

                DO id = 1, nDim
                   DO jd = 1, nDim
                      DO kd = 1, nDim
                         ! - ( div w , T* )
                         i1 = i + (id-1)*nNodesPerElement
                         j1 = j + (kd-1)*nNodesPerElement
                         k1 = k + (kd-1)*nNodesPerElement
                         int = 0.5d0*GSF_i(id)*SF_j*SF_k*wq
                         Aloc(i1,j1) = Aloc(i1,j1) - Lm1(j,jd,kd)*(Xe(j,jd) - Xe(k,jd))*int
                         Aloc(i1,k1) = Aloc(i1,k1) + Lm1(k,jd,kd)*(Xe(j,jd) - Xe(k,jd))*int

                         ! tauU ( - Lm1 w , L Gp* )
                         DO ld = 1, nDim
                            DO md = 1, nDim
                               i1 = i + (id-1)*nNodesPerElement
                               j1 = j + (ld-1)*nNodesPerElement
                               k1 = k + (ld-1)*nNodesPerElement
                               int = tauU*0.5d0*Lm1_q(jd,id)*SF_i*L_q(jd,md)*(SF_j*GSF_k(md) + SF_k*GSF_j(md))*wq
                               Aloc(i1,j1) = Aloc(i1,j1) - Lm1(j,kd,ld)*(Xe(j,kd) - Xe(k,kd))*int
                               Aloc(i1,k1) = Aloc(i1,k1) + Lm1(k,kd,ld)*(Xe(j,kd) - Xe(k,kd))*int
                            END DO
                         END DO

                         ! tauU ( Gq , LGp* )
                         DO ld = 1, nDim
                            i1 = i + nDim*nNodesPerElement
                            j1 = j + (kd-1)*nNodesPerElement
                            k1 = k + (kd-1)*nNodesPerElement
                            int = 0.5d0*tauU*GSF_i(id)*L_q(id,ld)*(SF_k*GSF_j(ld) + SF_j*GSF_k(ld))*wq
                            Aloc(i1,j1) = Aloc(i1,j1) + Lm1(j,jd,kd)*(Xe(j,jd) - Xe(k,jd))*int
                            Aloc(i1,k1) = Aloc(i1,k1) - Lm1(k,jd,kd)*(Xe(j,jd) - Xe(k,jd))*int
                         END DO
                      END DO
                   END DO
                END DO
             END DO

          END DO

       END IF

    END DO

    DEALLOCATE ( Lm1_q, GSF_i, GSF_j, GSF_k, Xe, Xq, X_q, phi_q, L_q )
    DEALLOCATE ( GSF_i1, GSF_i2, GSF_j1, GSF_j2 )

  END SUBROUTINE Compute_ElementIntegral_Mixed
  !================================================================!

END MODULE ElementIntegral_Mixed
