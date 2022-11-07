MODULE DirichletConformal_Mixed

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactSolutiont
  USE ExactFunctionJumpt


  IMPLICIT NONE

CONTAINS

  !======================================================================!
  SUBROUTINE DirichletConf_Mixed(Aloc, F, ele, edgList, Data, iFace, iBC)
  !======================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    type(element)           , INTENT(IN)    :: ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !-----------------------------------------------------
    type(Edge) :: Edg
    !------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, L_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, xq, x_q
    !-----------------------------------------------------------------------
    REAL*8 :: L, wq, x_iq, y_iq, z_iq, int
    REAL*8 :: SF_i, SF_j, Nrm, h, pD_q, a, SF_k, SF_i1, SF_j1, SF_i2, SF_j2
    !----------------------------------------------------------------------
    INTEGER :: i, j, iq, id, jd, nDim, i1, j1, ii, ii1, ii2, i2, k, k1, j2
    INTEGER :: jj, jj1, jj2, kd, ld
    INTEGER :: nNodesPerElement, nQuad, nVar, nEdgesPerElement
    !----------------------------------------------------------

    nDim = Data%nDim; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; nVar = nDim + 1

    ALLOCATE ( Lm1_q(nDim,nDim), L_q(nDim,nDim), Lm1(nNodesPerElement,nDim,nDim) )
    ALLOCATE ( Exa_q(nVar), N(nDim), GSF_i(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )
    ALLOCATE ( Xe(nNodesPerElement,nDim) )

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, icas)
    END DO
    Xe = ele%Coor
    Aloc = 0.d0 ; F = 0.d0
    N = - ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm ; h = ele%A/L
    N = N/Nrm

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

       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, icas)
       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)
       a = Data%alphaD(iBC)/h!*Nrm
      ! IF ( icas /= 0 ) THEN
      !    CALL ExactSol(ele%L, Exa_q, x_q, ele, nDim)
      !    pD_q = Exa_q(nDim+1)
       !ELSE
        !  pD_q = Data%BCVal(iBC)
       !END IF

       IF ( icas /= 0 ) THEN
          CALL ExactSolt(ele%L, Exa_q, x_q,Data%T, ele, nDim)
          pD_q = Exa_q(nDim+1)
       ELSE
          pD_q = Data%BCVal(iBC) !prend la condition definit dans le fichier de donnees
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          ! - < w.nt , pD >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             F(i1) = F(i1) - N(id)*SF_i*pD_q*wq
          END DO
          ! a < q , p - pD >
          ! a < q , p >
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Aloc(i1,j1) = Aloc(i1,j1) + a*SF_i*SF_j*wq
          END DO

          IF ( Data%PressureEnrichment ) THEN
             ! a < q , p* >
             DO ii = 1, nEdgesPerElement
                j = Permut(ii,1) ; k = Permut(ii,2)
                CALL SFEval(ele%eleType, SF_j, xq, j)
                CALL SFEval(ele%eleType, SF_k, xq, k)
                int = a*SF_i*SF_j*SF_k*wq
                DO id = 1, nDim
                   DO jd = 1, nDim
                      i1 = i + nDim*nNodesPerElement
                      j1 = j + (jd-1)*nNodesPerElement
                      k1 = k + (jd-1)*nNodesPerElement
                      Aloc(i1,j1) = Aloc(i1,j1) + 0.5d0*Lm1(j,id,jd)*(Xe(j,id) - Xe(k,id))*int
                      Aloc(i1,k1) = Aloc(i1,k1) - 0.5d0*Lm1(k,id,jd)*(Xe(j,id) - Xe(k,id))*int
                   END DO
                END DO
             END DO
          END IF

          ! a < q , pD >
          i1 = i + nDim*nNodesPerElement
          F(i1) = F(i1) + a*SF_i*pD_q*wq
       END DO

       IF ( data%pressureEnrichment .AND. data%testFunctionEnrichment ) THEN

          DO ii = 1, nEdgesPerElement
             ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
             CALL SFEval(ele%eleType, SF_i1, xq, ii1)
             CALL SFEval(ele%eleType, SF_i2, xq, ii2)
             ! a < q* , pD >
             DO id = 1, nDim
                DO jd = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   int = a*0.5d0*SF_i1*SF_i2*pD_q*wq
                   F(i1) = F(i1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                   F(i2) = F(i2) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                END DO
             END DO

             DO j = 1, nNodesPerElement
                CALL SFEval(ele%eleType, SF_j, xq, j)
                ! a < q* , p >
                DO id = 1, nDim
                   DO jd = 1, nDim
                      i1 = ii1 + (id-1)*nNodesPerElement
                      i2 = ii2 + (id-1)*nNodesPerElement
                      j1 = j + nDim*nNodesPerElement
                      int = a*0.5d0*SF_i1*SF_i2*SF_j*wq
                      Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                      Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*int
                   END DO
                END DO
             END DO

             ! a < q* , p* >
             DO jj = 1, nEdgesPerElement
                jj1 = permut(jj,1) ; jj2 = permut(jj,2)
                CALL SFEval(ele%eleType, SF_j1, xq, jj1)
                CALL SFEval(ele%eleType, SF_j2, xq, jj2)
                DO id = 1, nDim
                   DO jd = 1, nDim
                      DO kd = 1, nDim
                         DO ld = 1, nDim
                            int = a*0.25d0*SF_i1*SF_i2*SF_j1*SF_j2*wq
                            i1 = ii1 + (id-1)*nNodesPerElement
                            i2 = ii2 + (id-1)*nNodesPerElement
                            j1 = jj1 + (kd-1)*nNodesPerElement
                            j2 = jj2 + (kd-1)*nNodesPerElement
                            Aloc(i1,j1) = Aloc(i1,j1) + Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                                        Lm1(jj1,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                            Aloc(i1,j2) = Aloc(i1,j2) - Lm1(ii1,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                                        Lm1(jj2,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                            Aloc(i2,j1) = Aloc(i2,j1) - Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                                        Lm1(jj1,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                            Aloc(i2,j2) = Aloc(i2,j2) + Lm1(ii2,jd,id)*(Xe(ii1,jd) - Xe(ii2,jd))*&
                                                        Lm1(jj2,ld,kd)*(Xe(jj1,ld) - Xe(jj2,ld))*int
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO

       END IF

    END DO

    DEALLOCATE ( Lm1_q, L_q, Exa_q, N, GSF_i, GSF_j, xq, x_q )

  END SUBROUTINE DirichletConf_Mixed
  !=============================================================!

END MODULE DirichletConformal_Mixed
