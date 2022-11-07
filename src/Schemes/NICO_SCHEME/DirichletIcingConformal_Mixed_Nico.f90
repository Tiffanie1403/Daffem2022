MODULE DirichletIcingConformal_Mixed_Nicolson

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactFunctionJumpt


  IMPLICIT NONE

CONTAINS

  !======================================================================!
  SUBROUTINE DirichletIcingConf_Mixed_Nico(Aloc, F,F_prec, ele, edgList, Data, iFace, iBC)
  !======================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    REAL*8, DIMENSION(:)    , INTENT(IN)    :: F_prec
    type(element)           , INTENT(IN)    :: ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !-----------------------------------------------------
    type(Edge) :: Edg
    !------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, L_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, xq, x_q, GSF_i1, GSF_i2
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_j1, GSF_j2, GSF_k
    !-----------------------------------------------------------------------
    REAL*8 :: L, wq, x_iq, y_iq, z_iq, int, lambda,hn_q, hn_q_prec, pD_q_prec
    REAL*8 :: SF_i, SF_j, Nrm, h, pD_q, a, SF_k, SF_i1, SF_j1, SF_i2, SF_j2
    !----------------------------------------------------------------------
    INTEGER :: i, j, iq, id, jd, nDim, i1, j1, ii, ii1, ii2, i2, k, k1, j2
    INTEGER :: jj, nNodesPerElement, nQuad, nVar, nEdgesPerElement
    !----------------------------------------------------------

    nDim = Data%nDim; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; nVar = nDim + 1

    ALLOCATE ( Lm1_q(nDim,nDim), L_q(nDim,nDim), Lm1(nNodesPerElement,nDim,nDim) )
    ALLOCATE ( Exa_q(nVar), N(nDim), GSF_i(nDim), GSF_i1(nDim), GSF_k(nDim),  &
    GSF_i2(nDim),GSF_j1(nDim), GSF_j2(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )
    ALLOCATE ( Xe(nNodesPerElement,nDim) )
    !PRINT*,ele%vertex,ele%ref,ele%lambda
    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, icas)
    END DO
    Xe = ele%Coor !coordonnes des points definissant l'element de reference
    Aloc = 0.d0 ; F = 0.d0
    N = - ele%N(iFace,:) ! recuperation de la normale
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm ; h = ele%A/L
    N = N/Nrm ! normalisation

    lambda=ele%lambda

    edg = edgList(iFace)
    CALL GetNQuad(edg%eleType, nQuad)

    DO iq = 1, nQuad

       CALL FaceCoorQuad(edg%eleType, xq, iq, iFace)
       !xq resort les points de quadrature definit sur le triangle unité
       CALL WeightQuad(edg%EleType, wq, iq)
       wq = wq*L
       x_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             x_q(id) = x_q(id) + SF_i*ele%Coor(i,id) !passage a l'element courant car calcul fait sur l'element de reference
          END DO
       END DO

       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%lambda**2
       END DO
       Nrm = sqrt(Nrm)
       a = alphaIP1*Nrm/h

       IF ( icas /= 0 ) THEN
          Exa_q = 0d0
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
          pD_q = Exa_q(nDim+1) ;  Exa_q = 0d0
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T-Data%Dt, ele, nDim)
          pD_q_prec = Exa_q(nDim+1)
       ELSE
          pD_q = Data%BCVal(iBC) !prend la condition definit dans le fichier de donnees
          pD_q_prec = Data%BCVal(iBC) !prend la condition definit dans le fichier de donnees
       END IF

       IF ( icas /= 0 ) THEN
          Exa_q = 0d0
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
          hN_q = DOT_PRODUCT(Exa_q(1:nDim),N) ; Exa_q = 0d0 
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T-Data%Dt, ele, nDim)
          hN_q_prec = DOT_PRODUCT(Exa_q(1:nDim),N)
       ELSE
          hN_q = Data%BCVal(iBC) ; hN_q_prec = Data%BCVal(iBC)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)

          ! - < w.nt , pd_q >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             F(i1) = F(i1) - N(id)*SF_i*pD_q*wq
          END DO

          ! - < w.nt , pd_q_prec >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             F(i1) = F(i1) - N(id)*SF_i*pD_q_prec*wq
          END DO

          !stabilization term for the dirichlet condition
          DO j = 1, nNodesPerElement
            CALL SFEval(ele%eleType, SF_j, xq, j)


             ! a < q , T^n+1 >
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Aloc(i1,j1) = Aloc(i1,j1) + a*SF_i*SF_j*wq

             ! -a < q , T^n >
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             F(i1) = F(i1) - a*SF_i*SF_j*F_prec(j1)*wq



          END DO

          ! a < q , pD >
          i1 = i + nDim*nNodesPerElement
          F(i1) = F(i1) + a*SF_i*pD_q*wq

          ! a < q , pD_prec >
          i1 = i + nDim*nNodesPerElement
          F(i1) = F(i1) + a*SF_i*pD_q_prec*wq

          IF ( Data%PressureEnrichment ) THEN
             DO ii = 1, nEdgesPerElement
                j = Permut(ii,1) ; k = Permut(ii,2)
                CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)
                ! a < q , (T^n+1)* >
                DO id = 1, nDim
                      i1 = i + nDim*nNodesPerElement
                      j1 = j + (id-1)*nNodesPerElement
                      k1 = k + (id-1)*nNodesPerElement
                      int = 0.5d0*a*SF_i*SF_j*SF_k*wq
                      Aloc(i1,j1) = Aloc(i1,j1) + (1/lambda)*(Xe(j,id) - Xe(k,id))*int
                      Aloc(i1,k1) = Aloc(i1,k1) - (1/lambda)*(Xe(j,id) - Xe(k,id))*int
                END DO

                ! -a < q , (T^n)* >
                DO id = 1, nDim
                      i1 = i + nDim*nNodesPerElement
                      j1 = j + (id-1)*nNodesPerElement
                      k1 = k + (id-1)*nNodesPerElement
                      int = -(1/lambda)*0.5d0*a*SF_i*SF_j*SF_k*wq
                      F(i1) = F(i1) + F_prec(j1)*(Xe(j,id) - Xe(k,id))*int
                      F(i1) = F(i1) - F_prec(k1)*(Xe(j,id) - Xe(k,id))*int
                END DO


             END DO
          END IF
       END DO
       IF(Data%TestFunctionEnrichment) THEN
         DO ii = 1, nEdgesPerElement
            ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
            CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)

            ! a < q* , pD >
            DO id=1,nDim
              i1 = ii1 + (id-1)*nNodesPerElement
              i2 = ii2 + (id-1)*nNodesPerElement
              int = a*0.5d0*SF_i1*SF_i2*pD_q*wq
              F(i1) = F(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
              F(i2) = F(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
            END DO

            ! a < q* , pD_prec >
            DO id=1,nDim
              i1 = ii1 + (id-1)*nNodesPerElement
              i2 = ii2 + (id-1)*nNodesPerElement
              int = a*0.5d0*SF_i1*SF_i2*pD_q_prec*wq
              F(i1) = F(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
              F(i2) = F(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
            END DO

            DO j = 1, nNodesPerElement
              CALL SFEval(ele%eleType, SF_j, xq, j)
              ! a < q* , T^n+1 >
              DO id = 1, nDim
                    i1 = ii1 + (id-1)*nNodesPerElement
                    i2 = ii2 + (id-1)*nNodesPerElement
                    j1 = j   + nDim*nNodesPerElement
                    int = 0.5d0*a*SF_i1*SF_j*SF_i2*wq
                    Aloc(i1,j1) = Aloc(i1,j1) + (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                    Aloc(i2,j1) = Aloc(i2,j1) - (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
              END DO

              ! -a < q* , T^n >
              DO id = 1, nDim
                    i1 = ii1 + (id-1)*nNodesPerElement
                    i2 = ii2 + (id-1)*nNodesPerElement
                    j1 = j   + nDim*nNodesPerElement
                    int = -(1/lambda)*0.5d0*a*SF_i1*SF_j*SF_i2*wq
                    F(i1) = F(i1) + F_prec(j1)*(Xe(ii1,id) - Xe(ii2,id))*int
                    F(i2) = F(i2) - F_prec(j1)*(Xe(ii1,id) - Xe(ii2,id))*int
              END DO

           END DO

         END DO
       END IF
       IF ( Data%PressureEnrichment .AND. Data%TestFunctionEnrichment ) THEN
         DO ii = 1, nEdgesPerElement
            ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
            CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)
            DO jj = 1, nEdgesPerElement
               j = Permut(jj,1) ; k = Permut(jj,2)
               CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)

               ! a < q* , (T^n+1)* >
               DO id = 1, nDim
                 DO jd=1,nDim
                     i1 = ii1 + (id-1)*nNodesPerElement
                     i2 = ii2 + (id-1)*nNodesPerElement
                     j1 = j   + (jd-1)*nNodesPerElement
                     k1 = k   + (jd-1)*nNodesPerElement
                     int = 0.25d0*a*SF_i1*SF_i2*SF_j*SF_k*wq
                     Aloc(i1,j1) = Aloc(i1,j1) + (1/lambda**2)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                     Aloc(i1,k1) = Aloc(i1,k1) - (1/lambda**2)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                     Aloc(i2,j1) = Aloc(i2,j1) - (1/lambda**2)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                     Aloc(i2,k1) = Aloc(i2,k1) + (1/lambda**2)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                   END DO
               END DO


               ! -a < q* , (T^n)* >
               DO id = 1, nDim
                 DO jd=1,nDim
                     i1 = ii1 + (id-1)*nNodesPerElement
                     i2 = ii2 + (id-1)*nNodesPerElement
                     j1 = j   + (jd-1)*nNodesPerElement
                     k1 = k   + (jd-1)*nNodesPerElement
                     int = -(1/lambda**2)*0.25d0*a*SF_i1*SF_i2*SF_j*SF_k*wq
                     F(i1) = F(i1) + F_prec(j1)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                     F(i1) = F(i1) - F_prec(k1)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                     F(i2) = F(i2) - F_prec(j1)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                     F(i2) = F(i2) + F_prec(k1)*(Xe(j,jd) - Xe(k,jd))&
                     *(Xe(ii1,id) - Xe(ii2,id))*int
                   END DO
               END DO

             END DO
           END DO

       END IF
    END DO

    F = 0.5d0*F ; Aloc = 0.5d0*Aloc
    DEALLOCATE ( Lm1_q, L_q, Exa_q, N, GSF_i, GSF_j, xq, x_q )

  END SUBROUTINE DirichletIcingConf_Mixed_Nico
  !=============================================================!


END MODULE DirichletIcingConformal_Mixed_Nicolson
