MODULE DirichletIcingConformal_Primal

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactFunctionJumpt

  IMPLICIT NONE

CONTAINS

  !=============================================================================!
  SUBROUTINE DirichletIcingConf_Primal(Aloc, Uloc, F, ele, edgList, Data, iFace, iBC)
  !=============================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:)  , INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)    , INTENT(OUT) :: F
    REAL*8, DIMENSION(:)    , INTENT(IN)  :: Uloc
    type(element)           , INTENT(IN)  :: ele
    type(edge), DIMENSION(:), INTENT(IN)  :: edgList
    type(DataStructure)     , INTENT(IN)  :: Data
    INTEGER                 , INTENT(IN)  :: iBC, iFace
    !-------------------------------------------------
    type(Edge) :: edg
    !------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: L_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Exa_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: N, GSF_i, GSF_j, xq, x_q
    !-----------------------------------------------------------------
    REAL*8 :: L, wq, pD_q
    REAL*8 :: SF_i, SF_j, Nrm, h, a
    !---------------------------------------------
    INTEGER :: i, j, iq, id, jd, nQuad
    INTEGER :: nDim, nNodesPerElement, nVar
    !---------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement ; nVar = 1

    ALLOCATE ( L_q(nDim,nDim), Exa_q(nDim+1), N(nDim), GSF_i(nDim) )
    ALLOCATE ( GSF_j(nDim), xq(nDim), x_q(nDim) )

    Aloc = 0.d0 ; F = 0.d0

    N = - ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm)
    N = N/Nrm
    L = Nrm
    h = ele%A/L

    edg = edgList(iFace)
    CALL GetNQuad(edg%eleType, nQuad)

    DO iq = 1, nQuad

       CALL FaceCoorQuad(edg%eleType, xq, iq, iFace)
       CALL WeightQuad(edg%eleType, wq, iq)
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
          Nrm = Nrm + ele%L(id,id)**2
       END DO
       Nrm = sqrt(Nrm)
       a = 2.d0/h*Nrm!Data%alphaD(iBC)/h*Nrm
       CALL SetNodalMatrixPermeability(L_q, ele, Data%icas)

       IF ( icas /= 0 ) THEN
          CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
          pD_q = Exa_q(nDim+1)
       ELSE
          pD_q = Data%BCVal(iBC)
       END IF

       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
          DO j = 1, nNodesPerElement
             CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

             ! - < q , [L Gu].n >
             DO id = 1, nDim
                DO jd = 1, nDim
                   ! LHS
                   Aloc(i,j) = Aloc(i,j) - SF_i*L_q(id,jd)*GSF_j(jd)*N(id)*wq
                   ! RHS
                   F(i) = F(i) + SF_i*L_q(id,jd)*GSF_j(jd)*N(id)*Uloc(j)*wq
                END DO
             END DO
             ! + < [L Gq].n , u >
             DO id = 1, nDim
                DO jd = 1, nDim
                   ! LHS
                   !Aloc(i,j) = Aloc(i,j) + L_q(id,jd)*GSF_i(jd)*N(id)*SF_j*wq
                   ! RHS
                   !F(i) = F(i) - L_q(id,jd)*GSF_i(jd)*N(id)*SF_j*Uloc(j)*wq
                END DO
             END DO
             ! a < q , u >
             ! LHS
             Aloc(i,j) = Aloc(i,j) + a*SF_i*SF_j*wq
             ! RHS
             F(i) = F(i) - a*SF_i*SF_j*Uloc(j)*wq
          END DO
          ! a < q , uD >
          F(i) = F(i) + a*SF_i*pD_q*wq

          ! + < [L Gq].n , uD >
          DO id = 1, nDim
             DO jd = 1, nDim
                !F(i) = F(i) + L_q(id,jd)*GSF_i(jd)*N(id)*pD_q*wq
             END DO
          END DO
       END DO

    END DO


    DEALLOCATE ( L_q, Exa_q, N, GSF_i, GSF_j, xq )

  END SUBROUTINE DirichletIcingConf_Primal
  !================================================================!

END MODULE DirichletIcingConformal_Primal
