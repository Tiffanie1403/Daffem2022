MODULE NeumannConformal_Mixed_Nicolson

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactSolutiont
  USE ExactFunctionJumpt


  IMPLICIT NONE

CONTAINS

  !====================================================================!
  SUBROUTINE NeumannConf_Mixed_Nicolson(Aloc, F, F_prec, ele, edgList, Data, iFace, iBC)
  !====================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    REAL*8, DIMENSION(:)    , INTENT(IN)    :: F_prec
    type(element)           , INTENT(IN)    :: ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !-------------------------------------------------------
    type(Edge) :: Edg
    !------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, L_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, xq, x_q
    !-------------------------------------------------------------------------
    REAL*8 :: L, wq, x_iq, y_iq, z_iq, SF_i1, SF_i2, SF_k
    REAL*8 :: SF_i, SF_j, Nrm, h, hN_q, a, int, lambda, hN_q_prec
    !-----------------------------------------------------
    INTEGER :: i, j, iq, id, jd, nDim, i1, j1, ii1, ii2, k
    INTEGER :: nNodesPerElement, nQuad, nVar, ii, i2, k1, kd
    INTEGER :: nEdgesPerElement
    !------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement
    nVar = nDim + 1

    ALLOCATE ( Lm1(nNodesPerElement, nDim, nDim) )
    ALLOCATE ( Lm1_q(nDim,nDim), L_q(nDim,nDim), Xe(nNodesPerElement, nDim) )
    ALLOCATE ( Exa_q(nVar), N(nDim), GSF_i(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, icas)
    END DO

    Xe = ele%Coor

    Aloc = 0.d0 ; F = 0.d0
    N = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm)
    N = N/Nrm ; L = Nrm
    h = ele%A/L
    a = alphaIP2*h/Nrm

    edg = edgList(iFace)
    CALL GetNQuad(edg%eleType, nQuad)

    DO iq = 1, nQuad

       CALL FaceCoorQuad(edg%eleType, xq, iq, iFace)
       CALL WeightQuad(edg%EleType, wq, iq)
       wq = wq*L

       x_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             x_q(id) = x_q(id) + SF_i*ele%Coor(i,id)
          END DO
       END DO

      lambda=ele%lambda

       IF ( icas /= 0 ) THEN
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
          hN_q = DOT_PRODUCT(Exa_q(1:nDim),N)
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T-Data%Dt, ele, nDim) !valeur de la condition au pas de temps precedent
          hN_q_prec = DOT_PRODUCT(Exa_q(1:nDim),N)
       ELSE
          hN_q = Data%BCVal(iBC) ;   hN_q_prec = Data%BCVal(iBC)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

             ! < w.nt , T^n+1 >
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + SF_i*SF_j*N(id)*wq
             END DO


            ! -< w.nt , T^n >
            DO id = 1, nDim
               i1 = i + (id-1)*nNodesPerElement
               j1 = j + nDim*nNodesPerElement
               F(i1) = F(i1) - SF_i*SF_j*F_prec(j1)*N(id)*wq
            END DO

            !terme de Nitsche
            ! - < q , (B^n+1).n >
            DO id = 1, nDim
              i1 = i + nDim*nNodesPerElement
              j1 = j + (id-1)*nNodesPerElement
              Aloc(i1,j1) = Aloc(i1,j1) - SF_i*SF_j*N(id)*wq
            END DO

            !terme de Nitsche
            !  < q , (B^n).n >
            DO id = 1, nDim
              i1 = i + nDim*nNodesPerElement
              j1 = j + (id-1)*nNodesPerElement
              F(i1) = F(i1) + SF_i*SF_j*F_prec(j1)*N(id)*wq
            END DO

          END DO

          ! -< q , hN >
          i1 = i + nDim*nNodesPerElement
          F(i1) = F(i1) - SF_i*hN_q*wq

          ! -< q , hN_prec >
          i1 = i + nDim*nNodesPerElement
          F(i1) = F(i1) - SF_i*hN_q_prec*wq


          IF ( data%pressureEnrichment ) THEN
              DO ii = 1, nEdgesPerElement
                 j = permut(ii,1) ; k = permut(ii,2)
                 CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)

                 ! < w.nt , (T^n+1)* >
                 DO id = 1, nDim
                   DO jd = 1, nDim
                         i1 = i + (id-1)*nNodesPerElement
                         j1 = j + (jd-1)*nNodesPerElement
                         k1 = k + (jd-1)*nNodesPerElement
                         int = 0.5d0*SF_i*SF_j*SF_k*N(id)*wq
                         Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                         Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                   END DO
                END DO

                ! -< w.nt , (T^n)* >
                DO id = 1, nDim
                  DO jd = 1, nDim
                        i1 = i + (id-1)*nNodesPerElement
                        j1 = j + (jd-1)*nNodesPerElement
                        k1 = k + (jd-1)*nNodesPerElement
                        int = -(1.d0/lambda)*0.5d0*SF_i*SF_j*SF_k*N(id)*wq
                        F(i1) = F(i1) + F_prec(j1)*(Xe(j,jd) - Xe(k,jd))*int
                        F(i1) = F(i1) - F_prec(k1)*(Xe(j,jd) - Xe(k,jd))*int
                  END DO
               END DO

             END DO
          END IF

       END DO
       IF ( data%TestFunctionEnrichment ) THEN
           DO ii = 1, nEdgesPerElement
              ii1 = permut(ii,1) ; ii2 = permut(ii,2)
              CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)

              ! -< q* , hN >
              DO id=1,nDim
                i1 = ii1 + (id-1)*nNodesPerElement
                i2 = ii2 + (id-1)*nNodesPerElement
                int = -0.5d0*SF_i1*SF_i2*hN_q*wq
                F(i1) = F(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                F(i2) = F(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
              END DO

              ! -< q* , hN_prec >
              DO id=1,nDim
                i1 = ii1 + (id-1)*nNodesPerElement
                i2 = ii2 + (id-1)*nNodesPerElement
                int = -0.5d0*SF_i1*SF_i2*hN_q_prec*wq
                F(i1) = F(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                F(i2) = F(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
              END DO

              DO j = 1, nNodesPerElement
                 CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

                 ! - < q* , (B^n+1).n >
                 DO id = 1, nDim
                   DO jd=1,nDim
                     i1 = ii1 + (id-1)*nNodesPerElement
                     i2 = ii2 + (id-1)*nNodesPerElement
                     j1 = j   + (jd-1)*nNodesPerElement
                     int = -0.5d0*SF_i1*SF_i2*SF_j*N(jd)*wq
                     Aloc(i1,j1) = Aloc(i1,j1) + (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                     Aloc(i2,j1) = Aloc(i2,j1) - (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                  END DO
                 END DO

                 !  < q* , (B^n).n >
                 DO id = 1, nDim
                   DO jd=1,nDim
                     i1 = ii1 + (id-1)*nNodesPerElement
                     i2 = ii2 + (id-1)*nNodesPerElement
                     j1 = j   + (jd-1)*nNodesPerElement
                     int = (1/lambda)*0.5d0*SF_i1*SF_i2*SF_j*N(jd)*wq
                     F(i1) = F(i1) + F_prec(j1)*(Xe(ii1,id) - Xe(ii2,id))*int
                     F(i2) = F(i2) - F_prec(j1)*(Xe(ii1,id) - Xe(ii2,id))*int
                  END DO
                 END DO

              END DO

          END DO
       END IF

    END DO

    Aloc = 0.5d0*Aloc ; F = 0.5d0*F
    DEALLOCATE ( Lm1_q, L_q, Exa_q, N, GSF_i, GSF_j, xq, x_q )

  END SUBROUTINE NeumannConf_Mixed_Nicolson
  !=============================================================!

END MODULE NeumannConformal_Mixed_Nicolson
