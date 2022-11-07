MODULE ElementIntegral_Mixed_Icing_Nicolson

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE SourceTerm

  IMPLICIT NONE

CONTAINS

  !==============================================================!
  SUBROUTINE Compute_ElementIntegralIcing_Mixed_Nico(ele, Aloc, Floc, Floc_prec, Data)
  !==============================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: Floc
    REAL*8, DIMENSION(:)  , INTENT(IN)  :: Floc_prec
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
    REAL*8 :: wq, mu_q, lambda, int1, int2
    REAL*8 :: tauU, tauP, int, Dt, T, a2
    !---------------------------------
    INTEGER :: nNodesPerElement, nDim, nQuad
    INTEGER :: i, iq, j, id, jd, i1, i2, j1, j2, kd, md, nd
    INTEGER :: nVar, ii, k, k1, ld, nEdgesPerElement
    INTEGER :: ii1, ii2, jj, jj1, jj2
    !-----------------------------------------------

    nDim = Data%nDim ; nVar = nDim + 1 ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; Dt = Data%Dt ; T=Data%T
    int1=0.d0 ; int2=0.d0
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Lm1_q(nDim,nDim), L_q(nDim,nDim) )
    ALLOCATE ( GSF_i(nDim), GSF_j(nDim), GSF_k(nDim), Xq(nDim), X_q(nDim) )
    ALLOCATE ( phi_q(nVar), Xe(nNodesPerElement,nDim) )
    ALLOCATE ( GSF_i1(nDim), GSF_i2(nDim), GSF_j1(nDim), GSF_j2(nDim) )


    a2 = 2.d0

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, Data%icas)
    END DO

    Xe = ele%Coor ; lambda=ele%lambda
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
          mu_q = ele%lambda**2
       END DO
       mu_q = sqrt(mu_q)
       tauU = 0.5d0
       tauP = Data%alphaS*h*h/mu_q

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             ! ( w , Lm1 B^n+1 )
             DO id = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement
                   j1 = j + (id-1)*nNodesPerElement
                   Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*SF_i*SF_j*wq
             END DO

             ! -( w , Lm1 B^n )
             DO id = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement
                   j1 = j + (id-1)*nNodesPerElement
                   Floc(i1) = Floc(i1) - (1.d0/lambda)*SF_i*SF_j*Floc_prec(j1)*wq
             END DO

             ! - ( div w , T^n+1 )
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - GSF_i(id)*SF_j*wq
             END DO

             !  ( div w , T^n )
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Floc(i1) = Floc(i1) + GSF_i(id)*SF_j*Floc_prec(j1)*wq
             END DO

             ! tauP a2/dt ( div w , T^n+1 )
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) +(a2/Dt)*tauP*GSF_i(id)*SF_j*wq
             END DO


             ! ( q , div B^n+1 )
             DO id = 1, nDim
                i1 = i + nDim*nNodesPerElement
                j1 = j + (id-1)*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + SF_i*GSF_j(id)*wq
             END DO

             ! -( q , div B^n )
    !         DO id = 1, nDim
    !            i1 = i + nDim*nNodesPerElement
    !            j1 = j + (id-1)*nNodesPerElement
    !            Floc(i1) = Floc(i1) - SF_i*GSF_j(id)*Floc_prec(j1)*wq
    !         END DO

             !ajout du terme avec la dependance en temps
             ! a2/Dt ( q , T^n+1 )
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Aloc(i1,j1) = Aloc(i1,j1) + (a2/Dt)*SF_i*SF_j*wq

             ! a2/dT ( q, T^n )
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Floc(i1) = Floc(i1) + (a2/Dt)*SF_i*SF_j*Floc_prec(j1)*wq

          END DO

          ! ( q , phi^n+1 )
          CALL Compute_SourceTerm(phi_q, X_q, T, ele, Data%icas, nDim)
          i1 = i + nDim*nNodesPerElement
          Floc(i1) = Floc(i1) + SF_i*phi_q(nDim+1)*wq

          ! ( q , phi^n )
          CALL Compute_SourceTerm(phi_q, X_q, T-Dt, ele, Data%icas, nDim)
          i1 = i + nDim*nNodesPerElement
          Floc(i1) = Floc(i1) + SF_i*phi_q(nDim+1)*wq

       END DO

       IF ( data%TestFunctionEnrichment ) THEN

          DO ii = 1, nEdgesPerElement
             ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
             CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)

             ! ( q* , phi^n+1 )
             CALL Compute_SourceTerm(phi_q, X_q, T, ele, Data%icas, nDim)
             DO id = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   int = 0.5d0*SF_i1*SF_i2*phi_q(nDim+1)*wq
                   Floc(i1) = Floc(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                   Floc(i2) = Floc(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
             END DO

             ! ( q* , phi^n )
             CALL Compute_SourceTerm(phi_q, X_q, T-Dt, ele, Data%icas, nDim)
             DO id = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   int = 0.5d0*SF_i1*SF_i2*phi_q(nDim+1)*wq
                   Floc(i1) = Floc(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                   Floc(i2) = Floc(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
             END DO

          END DO
       END IF

       !ajout des termes de stabilisation
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEVal(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

             ! tauU ( -Lm1 w , B^n+1 )
             DO id = 1, nDim
                 i1 = i + (id-1)*nNodesPerElement
                 j1 = j + (id-1)*nNodesPerElement
                 Aloc(i1,j1) = Aloc(i1,j1) - tauU*(1.d0/lambda)*SF_i*SF_j*wq
             END DO

             ! tauU ( Lm1 w , B^n )
             DO id = 1, nDim
                 i1 = i + (id-1)*nNodesPerElement
                 j1 = j + (id-1)*nNodesPerElement
                 Floc(i1) = Floc(i1) + tauU*(1.d0/lambda)*SF_i*SF_j*Floc_prec(j1)*wq
             END DO

            ! tau U ( -Lm1 w , L GT^n+1 )
            DO id = 1, nDim
                  i1 = i + (id-1)*nNodesPerElement
                  j1 = j + nDim*nNodesPerElement
                  Aloc(i1,j1) = Aloc(i1,j1) - tauU*SF_i*GSF_j(id)*wq
            END DO

            ! tau U ( Lm1 w , L GT^n )
            DO id = 1, nDim
              DO jd=1,nDim
                  i1 = i + (id-1)*nNodesPerElement
                  j1 = j + (jd-1)*nNodesPerElement
                  Floc(i1) = Floc(i1) - tauU*SF_i*GSF_j(id)*(1.d0/lambda)*Floc_prec(j1)*wq
                END DO
            END DO

           ! tauU ( Gq , L GT^n+1 )
           DO id = 1, nDim
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Aloc(i1,j1) = Aloc(i1,j1) + tauU*GSF_i(id)*lambda*GSF_j(id)*wq
           END DO

           ! -tauU ( Gq , L GT^n )
           DO id = 1, nDim
             DO jd=1,Ndim
               i1 = i + nDim*nNodesPerElement
               j1 = j + (jd-1)*nNodesPerElement
               Floc(i1) = Floc(i1) + tauU*GSF_i(id)*(1.d0/lambda)*Floc_prec(j1)*lambda*GSF_j(id)*wq
            END DO
           END DO

           ! tauU ( Gq , B^n+1 )
           DO id = 1, nDim
            i1 = i + nDim*nNodesPerElement
            j1 = j + (id-1)*nNodesPerElement
            Aloc(i1,j1) = Aloc(i1,j1) + tauU*GSF_i(id)*SF_j*wq
           END DO

           ! -tauU ( Gq , B^n )
           DO id = 1, nDim
            i1 = i + nDim*nNodesPerElement
            j1 = j + (id-1)*nNodesPerElement
            Floc(i1) = Floc(i1) - tauU*GSF_i(id)*Floc_prec(j1)*SF_j*wq
           END DO

           ! tauP ( div w , div B^n+1 )
           DO id = 1, nDim
            DO jd = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             j1 = j + (jd-1)*nNodesPerElement
             Aloc(i1,j1) = Aloc(i1,j1) + tauP*GSF_i(id)*GSF_j(jd)*wq
            END DO
           END DO

           ! -tauP ( div w , div B^n )
      !     DO id = 1, nDim
      !      DO jd = 1, nDim
      !       i1 = i + (id-1)*nNodesPerElement
      !       j1 = j + (jd-1)*nNodesPerElement
      !       Floc(i1) = Floc(i1) - tauP*GSF_i(id)*GSF_j(jd)*Floc_prec(j1)*wq
      !      END DO
      !     END DO


           ! tauP  a2/dT ( div w , T^n )
           DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Floc(i1) = Floc(i1) + tauP*(a2/Dt)*GSF_i(id)*SF_j*Floc_prec(j1)*wq
           END DO

        END DO

        CALL Compute_SourceTerm(phi_q, X_q, T, ele, Data%icas, nDim)
        ! tauP ( div w , phi ) !associé au terme de stabilisation div-div
        DO id = 1, nDim
         i1 = i + (id-1)*nNodesPerElement
         Floc(i1) = Floc(i1) + tauP*GSF_i(id)*phi_q(nDim+1)*wq
        END DO

        CALL Compute_SourceTerm(phi_q, X_q, T-Dt, ele, Data%icas, nDim)
        ! tauP ( div w , phi ) !associé au terme de stabilisation div-div
        DO id = 1, nDim
         i1 = i + (id-1)*nNodesPerElement
         Floc(i1) = Floc(i1) + tauP*GSF_i(id)*phi_q(nDim+1)*wq
        END DO

       END DO

       ! Pressure enrichment components
       IF ( Data%PressureEnrichment ) THEN
          DO i = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
             DO ii = 1, nEdgesPerElement
                j = Permut(ii,1) ; k = Permut(ii,2)
                CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
                CALL SFEval(ele%eleType, SF_k, xq, k) ; CALL GSFEval(GSF_k, xq, k, ele)


                !  - ( div w , T^n+1* )

                DO id = 1, nDim
                   DO jd = 1, nDim
                     i1 = i + (id-1)*nNodesPerElement
                     j1 = j + (jd-1)*nNodesPerElement
                     k1 = k + (jd-1)*nNodesPerElement
                     int = -0.5d0*GSF_i(id)*SF_j*SF_k*wq
                     Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                     Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                   END DO
                END DO


                !  ( div w , T^n* )

                DO id = 1, nDim
                   DO jd = 1, nDim
                     i1 = i + (id-1)*nNodesPerElement
                     j1 = j + (jd-1)*nNodesPerElement
                     k1 = k + (jd-1)*nNodesPerElement
                     int = 0.5d0*GSF_i(id)*SF_j*SF_k*wq
                     Floc(i1) = Floc(i1) + Floc_prec(j1)*(1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                     Floc(i1) = Floc(i1) - Floc_prec(k1)*(1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                   END DO
                END DO



                !tauP a2/Dt ( div w , T^n+1* )
                DO id = 1, nDim
                   DO jd = 1, nDim
                     i1 = i + (id-1)*nNodesPerElement
                     j1 = j + (jd-1)*nNodesPerElement
                     k1 = k + (jd-1)*nNodesPerElement
                     int = (a2/Dt)*tauP*0.5d0*GSF_i(id)*SF_j*SF_k*wq
                     Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                     Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                   END DO
                END DO

               ! tauU ( - Lm1 w , L G(T^n+1)* )
               DO id = 1, nDim
                DO jd = 1, nDim
                 i1 = i + (id-1)*nNodesPerElement
                 j1 = j + (jd-1)*nNodesPerElement
                 k1 = k + (jd-1)*nNodesPerElement
                 int = -tauU*0.5d0*SF_i*(SF_j*GSF_k(id) + SF_k*GSF_j(id))*wq
                 Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                 Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                END DO
               END DO

               ! -tauU ( - Lm1 w , L GT^n* )
               DO id = 1, nDim
                DO jd = 1, nDim
                 i1 = i + (id-1)*nNodesPerElement
                 j1 = j + (jd-1)*nNodesPerElement
                 k1 = k + (jd-1)*nNodesPerElement
                 int = tauU*0.5d0*SF_i*(SF_j*GSF_k(id) + SF_k*GSF_j(id))*wq
                 Floc(i1) = Floc(i1) + Floc_prec(j1)*(1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                 Floc(i1) = Floc(i1) - Floc_prec(k1)*(1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                END DO
               END DO

               ! tauU ( Gq , L G(T^n+1)* )
               DO id = 1, nDim
                DO jd = 1, nDim
                  i1 = i + nDim*nNodesPerElement
                  j1 = j + (jd-1)*nNodesPerElement
                  k1 = k + (jd-1)*nNodesPerElement
                  int = 0.5d0*tauU*GSF_i(id)*(SF_k*GSF_j(id) + SF_j*GSF_k(id))*wq
                  Aloc(i1,j1) = Aloc(i1,j1) + (Xe(j,jd) - Xe(k,jd))*int
                  Aloc(i1,k1) = Aloc(i1,k1) - (Xe(j,jd) - Xe(k,jd))*int
                END DO
               END DO

               ! -tauU ( Gq , L G(T^n)* )
               DO id = 1, nDim
                DO jd = 1, nDim
                  i1 = i + nDim*nNodesPerElement
                  j1 = j + (jd-1)*nNodesPerElement
                  k1 = k + (jd-1)*nNodesPerElement
                  int = 0.5d0*tauU*GSF_i(id)*(SF_k*GSF_j(id) + SF_j*GSF_k(id))*wq
                  Floc(i1) = Floc(i1) + Floc_prec(j1)*(Xe(j,jd) - Xe(k,jd))*int
                  Floc(i1) = Floc(i1) - Floc_prec(k1)*(Xe(j,jd) - Xe(k,jd))*int
                END DO
               END DO


               ! a2/Dt ( q , T^(n+1)* )
               DO id = 1, nDim
                  i1 = i + nDim*nNodesPerElement
                  j1 = j + (id-1)*nNodesPerElement
                  k1 = k + (id-1)*nNodesPerElement
                  int = 0.5d0*(a2/Dt)*SF_i*SF_j*SF_k*wq
                  Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,id) - Xe(k,id))*int
                  Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,id) - Xe(k,id))*int
               END DO

               ! a2/Dt ( q , T^n* )
               DO id = 1, nDim
                  i1 = i + nDim*nNodesPerElement
                  j1 = j + (id-1)*nNodesPerElement
                  k1 = k + (id-1)*nNodesPerElement
                  int = (1.d0/lambda)*0.5d0*(a2/Dt)*SF_i*SF_j*SF_k*wq
                  Floc(i1) = Floc(i1) + Floc_prec(j1)*(Xe(j,id) - Xe(k,id))*int
                  Floc(i1) = Floc(i1) - Floc_prec(k1)*(Xe(j,id) - Xe(k,id))*int
               END DO

               ! - tauP  a2/dT ( div w , T^n* )
               DO id = 1, nDim
                 DO jd= 1,nDim
                  i1 = i + (id-1)*nNodesPerElement
                  j1 = j + (jd-1)*nNodesPerElement
                  k1 = k + (jd-1)*nNodesPerElement
                  int = (a2/lambda)*0.5d0*(1.d0/Dt)*tauP*SF_j*SF_k*wq
                  Floc(i1) = Floc(i1) + GSF_i(id)*Floc_prec(j1)*(Xe(j,jd) - Xe(k,jd))*int
                  Floc(i1) = Floc(i1) - GSF_i(id)*Floc_prec(k1)*(Xe(j,jd) - Xe(k,jd))*int
              END DO
             END DO

           END DO
        END DO
       END IF

       ! Test function enrichment
       IF ( Data%TestFunctionEnrichment ) THEN
          DO ii = 1, nEdgesPerElement
            ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
            CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL GSFEval(GSF_i1, xq, ii1, ele)
            CALL SFEval(ele%eleType, SF_i2, xq, ii2) ; CALL GSFEval(GSF_i2, xq, ii2, ele)
               DO j = 1, nNodesPerElement
                  CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

                  ! ( q* , div B^n+1 )
                  DO id = 1, nDim
                     DO jd = 1, nDim
                       i1 = ii1 + (id-1)*nNodesPerElement
                       i2 = ii2 + (id-1)*nNodesPerElement
                       j1 = j   + (jd-1)*nNodesPerElement
                       int = 0.5d0*SF_i1*SF_i2*GSF_j(jd)*wq
                       Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                       Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                      END DO
                  END DO

                  ! -( q* , div B^n1 )
                  DO id = 1, nDim
                     DO jd = 1, nDim
                       i1 = ii1 + (id-1)*nNodesPerElement
                       i2 = ii2 + (id-1)*nNodesPerElement
                       j1 = j   + (jd-1)*nNodesPerElement
                       int = -0.5d0*SF_i1*SF_i2*GSF_j(jd)*wq
                       Floc(i1) = Floc(i1) + Floc_prec(j1)*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                       Floc(i2) = Floc(i2) - Floc_prec(j1)*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                      END DO
                  END DO


                  ! a2/Dt ( q* , T^(n+1) )
                  DO id = 1, nDim
                     i1 = ii1 + (id-1)*nNodesPerElement
                     i2 = ii2 + (id-1)*nNodesPerElement
                     j1 = j   + nDim*nNodesPerElement

                     int = 0.5d0*(a2/Dt)*SF_i1*SF_i2*SF_j*wq
                     Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                     Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                  END DO

                  ! a2/DT ( q* , T^n)
                  DO id = 1, nDim
                        i1 = ii1 + (id-1)*nNodesPerElement
                        i2 = ii2 + (id-1)*nNodesPerElement
                        j1 = j   + nDim*nNodesPerElement

                        int = 0.5d0*(a2/Dt)*SF_i1*SF_i2*SF_j*Floc_prec(j1)*wq
                        Floc(i1) = Floc(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                        Floc(i2) = Floc(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                  END DO


                 !terme de stabilisation enrichie
                 ! tauU ( Gq* , B^n+1 )
                 DO id = 1, nDim
                   DO kd = 1, nDim
                    i1 = ii1 + (id-1)*nNodesPerElement
                    i2 = ii2 + (id-1)*nNodesPerElement
                    j1 = j   + (kd-1)*nNodesPerElement
                    int = tauU*0.5d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*SF_j*wq
                    Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                    Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                  END DO
                 END DO

                 !terme de stabilisation enrichie
                 ! -tauU ( Gq* , B^n )
                 DO id = 1, nDim
                   DO kd = 1, nDim
                    i1 = ii1 + (id-1)*nNodesPerElement
                    i2 = ii2 + (id-1)*nNodesPerElement
                    j1 = j   + (kd-1)*nNodesPerElement
                    int = -tauU*0.5d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*SF_j*wq
                    Floc(i1) = Floc(i1) + Floc_prec(j1)*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                    Floc(i2) = Floc(i2) - Floc_prec(j1)*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                  END DO
                 END DO


                 ! tauU ( Gq* , L GT^n+1)
                 DO id = 1, nDim
                    DO jd = 1, nDim
                       i1 = ii1 + (id-1)*nNodesPerElement
                       i2 = ii2 + (id-1)*nNodesPerElement
                       j1 = j   + nDim*nNodesPerElement
                       int = tauU*0.5d0*(SF_i1*GSF_i2(jd) + GSF_i1(jd)*SF_i2)*GSF_j(jd)*wq
                       Aloc(i1,j1) = Aloc(i1,j1) + (Xe(ii1,id) - Xe(ii2,id))*int
                       Aloc(i2,j1) = Aloc(i2,j1) - (Xe(ii1,id) - Xe(ii2,id))*int
                    END DO
                 END DO

                 ! tauU ( Gq* , L GT^n)
                 DO id = 1, nDim
                    DO jd = 1, nDim
                       i1 = ii1 + (id-1)*nNodesPerElement
                       i2 = ii2 + (id-1)*nNodesPerElement
                       j1 = j   + nDim*nNodesPerElement
                       int = -tauU*0.5d0*(SF_i1*GSF_i2(jd) + GSF_i1(jd)*SF_i2)*GSF_j(jd)*wq
                       Floc(i1) = Floc(i1) + Floc_prec(j1)*(Xe(ii1,id) - Xe(ii2,id))*int
                       Floc(i2) = Floc(i2) - Floc_prec(j1)*(Xe(ii1,id) - Xe(ii2,id))*int
                    END DO
                 END DO


                END DO
              END DO
            END IF

            IF ( Data%TestFunctionEnrichment .AND. Data%PressureEnrichment ) THEN
              DO ii = 1, nEdgesPerElement
                ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
                CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL GSFEval(GSF_i1, xq, ii1, ele)
                CALL SFEval(ele%eleType, SF_i2, xq, ii2) ; CALL GSFEval(GSF_i2, xq, ii2, ele)
                DO jj = 1, nEdgesPerElement
                  jj1 = Permut(jj,1) ; jj2 = Permut(jj,2)
                  CALL SFEval(ele%eleType, SF_j1, xq, jj1) ; CALL GSFEval(GSF_j1, xq, jj1, ele)
                  CALL SFEval(ele%eleType, SF_j2, xq, jj2) ; CALL GSFEval(GSF_j2, xq, jj2, ele)

                   ! tauU( Gq* , L G(T^n+1)* )
                   DO id = 1, nDim
                      DO jd = 1, nDim
                         DO kd = 1, nDim
                           i1 = ii1 + (id-1)*nNodesPerElement
                           i2 = ii2 + (id-1)*nNodesPerElement
                           j1 = jj1 + (jd-1)*nNodesPerElement
                           j2 = jj2 + (jd-1)*nNodesPerElement

                           int = tauU*0.25d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*&
                                (SF_j1*GSF_j2(kd) + GSF_j1(kd)*SF_j2)*wq

                           Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
                          (Xe(jj1,jd) - Xe(jj2,jd))*int
                           Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
                          (Xe(jj1,jd) - Xe(jj2,jd))*int
                          Aloc(i1,j2) = Aloc(i1,j2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
                          (Xe(jj1,jd) - Xe(jj2,jd))*int
                           Aloc(i2,j2) = Aloc(i2,j2) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
                          (Xe(jj1,jd) - Xe(jj2,jd))*int
                       END DO
                    END DO
                 END DO

                 ! tauU( Gq* , L G(T^n)* )
                 DO id = 1, nDim
                    DO jd = 1, nDim
                       DO kd = 1, nDim
                         i1 = ii1 + (id-1)*nNodesPerElement
                         i2 = ii2 + (id-1)*nNodesPerElement
                         j1 = jj1 + (jd-1)*nNodesPerElement
                         j2 = jj2 + (jd-1)*nNodesPerElement

                         int = -(1.d0/lambda)*tauU*0.25d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*&
                              (SF_j1*GSF_j2(kd) + GSF_j1(kd)*SF_j2)*wq

                         Floc(i1) = Floc(i1) + Floc_prec(j1)*(Xe(ii1,id) - Xe(ii2,id))*&
                        (Xe(jj1,jd) - Xe(jj2,jd))*int
                         Floc(i2) = Floc(i2) -Floc_prec(j1)* (Xe(ii1,id) - Xe(ii2,id))*&
                        (Xe(jj1,jd) - Xe(jj2,jd))*int
                        Floc(i1) = Floc(i1) - Floc_prec(j2)*(Xe(ii1,id) - Xe(ii2,id))*&
                        (Xe(jj1,jd) - Xe(jj2,jd))*int
                         Floc(i2) = Floc(i2) + Floc_prec(j2)*(Xe(ii1,id) - Xe(ii2,id))*&
                        (Xe(jj1,jd) - Xe(jj2,jd))*int
                     END DO
                  END DO
               END DO

                 ! a2/Dt ( q* , T^(n+1)* )
                 DO id = 1, nDim
                    DO jd = 1, nDim
                       i1 = ii1 + (id-1)*nNodesPerElement
                       i2 = ii2 + (id-1)*nNodesPerElement
                       j1 = jj1 + (jd-1)*nNodesPerElement
                       j2 = jj2 + (jd-1)*nNodesPerElement

                       int = (1.d0/(lambda**2))*0.25d0*(a2/Dt)*SF_i1*SF_i2*SF_j1*SF_j2*wq
                       Aloc(i1,j1) = Aloc(i1,j1) + (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                       Aloc(i1,j2) = Aloc(i1,j2) - (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                       Aloc(i2,j1) = Aloc(i2,j1) - (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                       Aloc(i2,j2) = Aloc(i2,j2) + (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                    END DO
                 END DO

                 ! a2/Dt ( q* , T^n* )
                 DO id = 1, nDim
                    DO jd = 1, nDim
                       i1 = ii1 + (id-1)*nNodesPerElement
                       i2 = ii2 + (id-1)*nNodesPerElement
                       j1 = jj1 + (jd-1)*nNodesPerElement
                       j2 = jj2 + (jd-1)*nNodesPerElement

                       int = (1.d0/(lambda**2))*0.25d0*(a2/Dt)*SF_i1*SF_i2*SF_j1*SF_j2*wq
                       Floc(i1) = Floc(i1) + Floc_prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                       Floc(i1) = Floc(i1) - Floc_prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                       Floc(i2) = Floc(i2) - Floc_prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                       Floc(i2) = Floc(i2) + Floc_prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
                    END DO
                 END DO

             END DO
          END DO
        END IF
    END DO

    Floc= 0.5d0*Floc ; Aloc = 0.5d0*Aloc
    DEALLOCATE ( Lm1_q, GSF_i, GSF_j, GSF_k, Xe, Xq, X_q, phi_q, L_q )
    DEALLOCATE ( GSF_i1, GSF_i2, GSF_j1, GSF_j2 )

  END SUBROUTINE Compute_ElementIntegralIcing_Mixed_Nico
  !================================================================!

END MODULE ElementIntegral_Mixed_Icing_Nicolson
