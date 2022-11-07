MODULE TransferOperator_MSEG

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE SourceTerm

  IMPLICIT NONE

CONTAINS

  !============================================!
  SUBROUTINE ComputeMt_MSEG(ele, Edg, Mt, Data)
  !============================================!

    IMPLICIT NONE
    type(element)           , INTENT(IN)  :: ele
    type(Edge), DIMENSION(:), INTENT(IN)  :: Edg
    REAL*8, DIMENSION(:,:)  , INTENT(OUT) :: Mt
    type(DataStructure)     , INTENT(IN)  :: Data
    !----------------------------------------------
    type(Edge) :: EdgInt
    !----------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Id_Mat, Lm1_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: N, Xq, X_q
    !--------------------------------------------------
    REAL*8 :: aP, aB, hOrt, NL, Nrm, wq, mu_q, L, SF_i
    !---------------------------------------------------
    INTEGER :: i, j, i1, j1, iq, ii, iCheck, iVar, nVar, nQuad
    INTEGER :: id, jd, nDim, nNodesPerElement, nDofPerElement
    !------------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nDofPerElement = nNodesPerElement + 1 ; nVar = nDim + 1

    ALLOCATE ( Id_Mat(nDim,nDim), Lm1_q(nDim,nDim) )
    ALLOCATE ( N(nDim), Xq(nDim), x_q(nDim) )

    Mt = 0.d0
    DO ivar = 1, nVar
       DO i = 1, nNodesPerElement
          i1 = i + (ivar-1)*nDofPerElement
          Mt(i1,i1) = 1.d0
       END DO
    END DO

    CALL GetNQuad(ele%eleType, nQuad)

    DO iq = 1, nQuad

       ! Element Integral
       CALL CoorQuad(ele%eleType, Xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       X_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             X_q(id) = X_q(id) + SF_i*ele%Coor(i,id)
          END DO
       END DO
       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, Data%icas)

!!$       ! ( wk , Lm1 Bk )
!!$       DO id = 1, nDim
!!$          DO jd = 1, nDim
!!$             i1 = nDofPerElement + (id-1)*nDofPerElement
!!$             j1 = nDofPerElement + (jd-1)*nDofPerElement
!!$             Mt(i1,j1) = Mt(i1,j1) + Lm1_q(id,jd)*wq
!!$          END DO
!!$       END DO
       ! 1/2 ( wk , Lm1 Bk )
       DO id = 1, nDim
          DO jd = 1, nDim
             i1 = nDofPerElement + (id-1)*nDofPerElement
             j1 = nDofPerElement + (jd-1)*nDofPerElement
             Mt(i1,j1) = Mt(i1,j1) + 0.5d0*Lm1_q(id,jd)*wq
          END DO
       END DO
    END DO

    ! Integrals boundary
    DO ii = 1, nNodesPerElement

       EdgInt = Edg(ii)
       N = -ele%N(ii,:)
       Nrm = 0.d0
       DO id = 1, nDim
          Nrm = Nrm + N(id)**2
       END DO
       Nrm = sqrt(Nrm)
       N = N/Nrm ; L = Nrm
       hOrt = ele%A/L

       CALL GetNQuad(edgInt%eleType, nQuad)

       DO iq = 1, nQuad

          CALL FaceCoorQuad(edgInt%eleType, xq, iq, ii)
          CALL WeightQuad(edgInt%eleType, wq, iq)
          wq = wq*L

          X_q = 0.d0
          DO i = 1, nNodesPerElement
             CALL SFEval(edgInt%eleType, SF_i, xq, i)
             DO id = 1, nDim
                X_q(id) = X_q(id) + ele%Coor(i,id)*SF_i
             END DO
          END DO
          mu_q = 0.d0
          DO id = 1, nDim
             mu_q = mu_q + ele%L(id,id)**2
          END DO
          mu_q = sqrt(mu_q)

          aP = alphaP*hOrt/mu_q
          aB = alphaB/hort*mu_q

          DO id = 1, nDim
             ! - < wK.n , pK >
             i1 = nDofPerElement + (id-1)*nDofPerElement
             j1 = nDofPerElement + nDim*nDofPerElement
             Mt(i1,j1) = Mt(i1,j1) - N(id)*wq

             ! aP < wK.n , Bk.n >
             DO jd = 1, nDim
                i1 = nDofPerElement + (id-1)*nDofPerElement
                j1 = nDofPerElement + (jd-1)*nDofPerElement
                Mt(i1,j1) = Mt(i1,j1) + aP*N(id)*N(jd)*wq
             END DO


             ! - < qK , BK.n >
             i1 = nDofPerElement + nDim*nDofPerElement
             j1 = nDofPerElement + (id-1)*nDofPerElement
             Mt(i1,j1) = Mt(i1,j1) - N(id)*wq
          END DO

          ! aB < qK , pK >
          i1 = nDofPerElement + nDim*nDofPerElement
          j1 = nDofPerElement + nDim*nDofPerElement
          Mt(i1,j1) = Mt(i1,j1) + aB*wq

       END DO
    END DO

    DEALLOCATE ( Id_Mat, Lm1_q, N, Xq, x_q )

  END SUBROUTINE ComputeMt_MSEG
  !===============================================!

  !=======================================!
  SUBROUTINE ComputeMb_MSEG(Ele, Mb, Data)
  !=======================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)    :: ele
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: Mb
    type(DataStructure)   , INTENT(IN)    :: Data
    !---------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Lm1_q
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: N, Xq, GSF_j
    !-----------------------------------------------------
    REAL*8 :: L, Nrm, wq, SF_j
    REAL*8 :: aP, aB, h, mu_q
    !---------------------------
    INTEGER :: nDim, nNodesPerElement, nDofPerElement
    INTEGER :: i, j, id, jd, iq, nVar, ivar, i1, j1, nQuad
    !-------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nDofPerElement = nNodesPerElement + 1 ; nVar = nDim + 1

    ALLOCATE ( Lm1_q(nDim,nDim), N(nDim), Xq(nDim), GSF_j(nDim) )

    Mb = 0.d0

    DO ivar = 1, nVar
       DO i = 1, nNodesPerElement
          i1 = i + (ivar-1)*nDofPerElement
          j1 = i + (ivar-1)*nNodesPerElement
          Mb(i1,j1) = 1.d0
       END DO
    END DO

    ! Element Integral
    CALL GetNQuad(ele%eleType, nQuad)

    DO iq = 1, nQuad

       CALL CoorQuad(ele%eleType, Xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, Data%icas)
       mu_q = 0.d0
       DO id = 1, nDim
          mu_q = mu_q + ele%L(id,id)**2
       END DO
       mu_q = sqrt(mu_q)

       DO j = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
          DO id = 1, nDim
!!$             ! - ( wK , Lm1 B )
!!$             DO jd = 1, nDim
!!$                i1 = nDofPerElement + (id-1)*nDofPerElement
!!$                j1 = j + (id-1)*nNodesPerElement
!!$                Mb(i1,j1) = Mb(i1,j1) - SF_j*wq
!!$             END DO
!!$             ! - ( wk , Gp )
!!$             i1 = nDofPerElement + (id-1)*nDofPerElement
!!$             j1 = j + nDim*nNodesPerElement
!!$             Mb(i1,j1) = Mb(i1,j1) - GSF_j(id)*wq
             ! - 1/2 ( wK , Lm1 B )
             DO jd = 1, nDim
                i1 = nDofPerElement + (id-1)*nDofPerElement
                j1 = j + (id-1)*nNodesPerElement
                Mb(i1,j1) = Mb(i1,j1) - 0.5d0*SF_j*wq
             END DO
             ! - 1/2 ( wk , Gp )
             i1 = nDofPerElement + (id-1)*nDofPerElement
             j1 = j + nDim*nNodesPerElement
             Mb(i1,j1) = Mb(i1,j1) - 0.5d0*GSF_j(id)*wq

             ! - ( qK , div B )
             i1 = nDofPerElement + nDim*nDofPerElement
             j1 = j + (id-1)*nNodesPerElement
             Mb(i1,j1) = Mb(i1,j1) - GSF_j(id)*wq
          END DO
       END DO

    END DO

  END SUBROUTINE ComputeMb_MSEG
  !=======================================!

  !======================================================!
  SUBROUTINE ComputeTransferSourceTerm_MSEG(F, ele, Data)
  !======================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: F
    type(element)       , INTENT(IN)    :: ele
    type(DataStructure) , INTENT(IN)    :: Data
    !--------------------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: phi_q, Xq, X_q
    !----------------------------------------------------
    REAL*8 :: wq, SF_i
    !-------------------
    INTEGER :: nNodesPerElement, nDim, nDofPerElement
    INTEGER :: nVar, iq, i, i1, id, ivar, nQuad
    !-------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nDofPerElement = nNodesPerElement + 1 ; nVar = nDim + 1

    ALLOCATE ( phi_q(nDim+1), Xq(nDim), x_q(nDim) )

    F = 0.d0

    CALL GetNQuad(ele%eleType, nQuad)

    DO iq = 1, nQuad

       CALL CoorQuad(ele%eleType, Xq, iq)
       CALL WeightQuad(ele%eleType, wq, iq)
       wq = wq*ele%A

       X_q = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(ele%eleType, SF_i, xq, i)
          DO id = 1, nDim
             X_q(id) = X_q(id) + SF_i*ele%Coor(i,id)
          END DO
       END DO

       CALL Compute_SourceTerm(phi_q, X_q,Data%T, ele, Data%icas, nDim)
       ! ( qK , phi )
       DO ivar = 1, nVar
          i1 = nDofPerElement + (ivar-1)*nDofPerElement
          F(i1) = F(i1) + phi_q(ivar)*wq
       END DO
    END DO

    DEALLOCATE ( phi_q, Xq, x_q )

  END SUBROUTINE ComputeTransferSourceTerm_MSEG
  !==============================================!

END MODULE TransferOperator_MSEG
