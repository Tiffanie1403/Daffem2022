MODULE NeumannIcingConformal_Mixed

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactFunctionJumpt
  USE ExactSolutiont

  IMPLICIT NONE

CONTAINS

  !====================================================================!
  SUBROUTINE NeumannIcingConf_Mixed(Aloc, F, ele, edgList, Data, iFace, iBC)
    !====================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    type(element)           , INTENT(IN)    :: ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !-------------------------------------------------------
    type(Edge) :: Edg
    !------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, xq
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_i, GSF_j
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Xe

    REAL*8, DIMENSION(data%nDim) :: x_q
    !-------------------------------------------------------------------------
    REAL*8 :: L, wq, hN_q, SF_k
    REAL*8 :: SF_i, SF_j, Nrm, int, lambda, SF_i1, SF_i2
    !-----------------------------------------------------
    INTEGER :: i, j, iq, id, jd, nDim, i1, j1, k1, k
    INTEGER :: nNodesPerElement, nQuad, nVar, i2
    INTEGER :: nEdgesPerElement, ii, ii1, ii2
    !------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement
    nVar = nDim + 1

    ALLOCATE ( Exa_q(nVar), N(nDim), xq(nDim))
    ALLOCATE ( GSF_i(nDim), GSF_j(nDim))
    ALLOCATE ( Xe(nNodesPerElement,nDim) )


    Aloc = 0.d0 ; F = 0.d0
    N = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm)
    N = N/Nrm ; L = Nrm
    Xe = ele%Coor !coordonnes des points definissant l'element de reference
    lambda=ele%lambda
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

       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + ele%lambda**2
       END DO
       Nrm = sqrt(Nrm)

       IF ( icas /= 0 ) THEN
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
          hN_q = DOT_PRODUCT(Exa_q(1:nDim),N)
       ELSE
          hN_q = Data%BCVal(iBC)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j); CALL GSFEval(GSF_j, xq, j, ele)

             ! < w.nt , T >
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) + SF_i*SF_j*N(id)*wq
             END DO

             !-<q, B.n>
             DO id = 1, nDim
                i1 = i + nDim*nNodesPerElement
                j1 = j + (id-1)*nNodesPerElement
                Aloc(i1,j1) = Aloc(i1,j1) - SF_j*N(id)*SF_i*wq ! ici on a une divergence donc une seule boucle
             END DO

          END DO

          ! -< q , hN >
          i1 = i + nDim*nNodesPerElement
          F(i1) = F(i1) - SF_i*hN_q*wq
       END DO

       IF ( data%pressureEnrichment ) THEN
         DO i = 1, nNodesPerElement
           CALL SFEval(ele%eleType, SF_i, xq, i); CALL GSFEval(GSF_i, xq, i, ele)
          DO ii = 1, nEdgesPerElement
             j = permut(ii,1) ; k = permut(ii,2)
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)

             DO id = 1, nDim
               DO jd = 1,nDim
                ! < w.nt , T* > morceaux non pris en compte par le dirichlet
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + (jd-1)*nNodesPerElement
                k1 = k + (jd-1)*nNodesPerElement
                int = SF_k*SF_j*SF_i*wq*N(id)
                Aloc(i1,j1) = Aloc(i1,j1) + 0.5d0*(1/lambda)*(Xe(j,jd) - Xe(k,jd))*int
                Aloc(i1,k1) = Aloc(i1,k1) - 0.5d0*(1/lambda)*(Xe(j,jd) - Xe(k,jd))*int
              END DO
             END DO


           END DO
         END DO
         IF(Data%TestFunctionEnrichment) THEN
           DO ii = 1, nEdgesPerElement
              ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
              CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)
              DO j = 1, nNodesPerElement
                 CALL SFEval(ele%eleType, SF_j, xq, j)

                 !-<q*, B.n>
                 DO id = 1, nDim
                   DO jd= 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   j1 = j   + (jd-1)*nNodesPerElement

                   int =  -SF_i1* SF_i2*SF_j*wq*N(jd)
                   Aloc(i1,j1) = Aloc(i1,j1) + 0.5d0*(1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                   Aloc(i2,j1) = Aloc(i2,j1) - 0.5d0*(1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                 END DO
               END DO
             END DO

             ! -< q* , hN >
             DO id = 1, nDim
             i1 = ii1 + (id-1)*nNodesPerElement
             i2 = ii2 + (id-1)*nNodesPerElement
             int = -SF_i1* SF_i2*hN_q*wq
             F(i1) = F(i1) + 0.5d0*int*1/lambda*(Xe(ii1,id) - Xe(ii2,id))
             F(i2) = F(i2) - 0.5d0*int*1/lambda*(Xe(ii1,id) - Xe(ii2,id))
            END DO

           END DO
         END IF
       END IF

    END DO

    DEALLOCATE ( Exa_q, N, xq, Xe )

  END SUBROUTINE NeumannIcingConf_Mixed
  !=============================================================!

END MODULE NeumannIcingConformal_Mixed
