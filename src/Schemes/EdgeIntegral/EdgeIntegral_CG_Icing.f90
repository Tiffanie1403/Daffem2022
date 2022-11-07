MODULE EdgeIntegral_CG_Icing

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE ExactSolutiont
  USE ExactFunctionJumpt
  USE ReadMesh_mod
  USE Distance
  USE ComputeNormalsAndTangents
  USE IcingTest_mod
  USE Velocity

  IMPLICIT NONE

CONTAINS

  !=====================================================================================!
  SUBROUTINE ComputePMMatrices_CG(mesh, iBC, edg, edgM, edgP, eleM, eleP, Amm, Amp, Apm, App, F1, F2, Data)
    !=====================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)    :: Edg
    type(Edge), DIMENSION(:), INTENT(IN)    :: edgM, edgP
    type(element)           , INTENT(IN)    :: eleM, eleP
    REAL*8, DIMENSION(:,:)  , INTENT(OUT)   :: Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F1,F2
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC
    type(MeshType)          , INTENT(IN)    :: Mesh
    !----------------------------------------------------------
    type(edge) :: edgIntM, edgIntP
    !----------------------------------------------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Np, Nm, Xq_M, Xq_P, Exa_p, Exa_m
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: XeM, XeP
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: XxqP, XxqM
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1_M, Lm1_P
    REAL*8, DIMENSION(data%nDim+1)        :: exa1_q, exa2_q
    REAL*8, DIMENSION(data%nDim)          :: x_qP, x_qM
    LOGICAL 						      :: testJumpT
    !-------------------------------------------------------------
    REAL*8  :: Nrm, wq, L, HORT,lambdaM,lambdaP
    REAL*8  :: SF_iM, SF_iP, SF_jM, SF_jP, a
    REAL*8  :: tauP, tauU, h, jn2, jt_bis, jn2_bis
    REAL*8  :: SF_i1M, SF_i1P, SF_i2P, SF_i2M, SF_j1M
    REAL*8  :: SF_j2P, SF_j1P,SF_j2M
    REAL*8  :: Sigma_bP, Sigma_bM, jt
    REAL*8  :: SF_iMM, SF_iMP, SF_iPM, SF_iPP, SF_km,SF_kp
    REAL*8  :: SF_jMM, SF_jMP, SF_jPM, SF_jPP, int, SF_i1, SF_i2
    !------------------------------------------------------------
    INTEGER :: iM, iP, i, jM, jP,iq, i1, j1, j, id, jd
    INTEGER :: ii,ii1,ii2,k,jj,jj1,jj2, i2, k1,kD,ld,j2
    INTEGER :: nDim, nNodesPerElement, nQuad, nEdgesPerElement, nVar
    !--------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nVar = nDim + 1 ; testJumpT = .FALSE. 
    nEdgesPerElement = eleM%nEdgesPerElement
    ALLOCATE ( Np(nDim), Nm(nDim), Xq_M(nDim), Xq_P(nDim))
    ALLOCATE ( XxqP(nDim), XxqM(nDim) )
    ALLOCATE ( Exa_p(nVar), Exa_M(nVar)  )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( XeM(nNodesPerElement,nDim), XeP(nNodesPerElement,nDim) )
    ALLOCATE ( Lm1_M(nNodesPerElement, nDim, nDim), Lm1_P(nNodesPerElement, nDim, nDim) )

    XeM = eleM%Coor ; XeP = eleP%Coor
    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1_M(i,:,:), eleM, Data%icas)
       CALL SetNodalInverseMatrixPermeability(Lm1_P(i,:,:), eleP, Data%icas)
    END DO
    CALL IdentifyIFace(eleM, edg, iM, nDim)
    CALL IdentifyIFace(eleP, edg, iP, nDim)
    
	
    lambdaM=eleM%lambda ; lambdaP=eleP%lambda
    edgIntM = edgM(iM)  ; edgIntP = edgP(iP)

    Np = -eleP%N(iP,:)
    Nrm = 0.d0

    DO id = 1, nDim
       Nrm = Nrm + Np(id)**2
    END DO
    Nrm  = sqrt(Nrm) ; L = Nrm
    Np   = Np/Nrm ! on normalise pour avoir des normales unitaires
    Nm   = -Np    !normale oppossé en m
    hOrt = (Elem%A + Elep%A)/(2.d0*L)
    h    = hOrt
    tauU = 0.d0 !alphaIP1*Nrm
    tauP = alphaIP2*h*h/Nrm
    CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim) ! iM numerotation of the edge for eleM
    CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)

    Amm = 0.d0 ; Amp = 0.d0 ; Apm = 0.d0 ; App = 0.d0
    F1 = 0.d0 ; F2 = 0.d0
    CALL GetNQuad(edg%eleType, nQuad)
    DO iq = 1, nQuad

       xq_M = xqM(iq,:)
       xq_P = xqP(iq,:)

       x_qP = 0.d0
       x_qM = 0.d0
       !position of the quadrature node from the reference element to the current element 
       DO i = 1, nNodesPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i)
          CALL SFEval(eleP%eleType, SF_iP, xq_P, i)
          DO id = 1, nDim
             x_qP(id) = x_qP(id) + SF_iP*eleP%Coor(i,id)
             x_qM(id) = x_qM(id) + SF_iM*eleM%Coor(i,id)
          END DO
       END DO
	
       IF ( icas /=  0 ) THEN
          CALL ExactSolJumpt(eleP%L, exa1_q, x_qP, Data%T,eleP, nDim)
          CALL ExactSolJumpt(eleM%L, exa2_q, x_qM,Data%T, eleM, nDim)
          jn2 = DOT_PRODUCT(exa1_q(1:nDim) - exa2_q(1:nDim), Np)
       ELSE
          jn2 = data%BCVal(iBC)
       END IF

       IF ( icas /= 0 ) THEN
          CALL ExactSolJumpt(eleP%L, exa1_q, x_qP,Data%T, eleP, nDim)
          CALL ExactSolJumpt(eleM%L, exa2_q, x_qM,Data%T, eleM, nDim)
          jt = exa1_q(nDim+1) - exa2_q(nDim+1)
       ELSE
          jt =  Data%BCVal(iBC)
       END IF
       CALL WeightQuad(edg%eleType, wq, iq)
       wq = wq*L

       ! terme a l'interface
       DO i = 1, nNodesPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL SFEval(eleP%eleType, SF_iP, xq_P, i)
          DO j = 1, nNodesPerElement
             CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL SFEval(eleP%eleType, SF_jP, xq_P, j)

             ! < [ w ] , { T } >
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Amm(i1,j1) = Amm(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jM*wq
                Amp(i1,j1) = Amp(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jP*wq
                Apm(i1,j1) = Apm(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jM*wq
                App(i1,j1) = App(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jP*wq
             END DO

             ! - < { q } , [ B.n ] >
             DO id = 1, nDim
                i1 = i + nDim*nNodesPerElement
                j1 = j + (id-1)*nNodesPerElement
                Amm(i1,j1) = Amm(i1,j1) - 0.5d0*SF_iM*SF_jM*Nm(id)*wq
                Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*SF_jP*Np(id)*wq
                Apm(i1,j1) = Apm(i1,j1) - 0.5d0*SF_iP*SF_jM*Nm(id)*wq
                App(i1,j1) = App(i1,j1) - 0.5d0*SF_iP*SF_jP*Np(id)*wq
             END DO

             ! tauU/h < [ q ] , [ T ] >
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Amm(i1,j1) = Amm(i1,j1) + tauU/h*SF_iM*SF_jM*wq
             Amp(i1,j1) = Amp(i1,j1) - tauU/h*SF_iM*SF_jP*wq
             Apm(i1,j1) = Apm(i1,j1) - tauU/h*SF_iP*SF_jM*wq
             App(i1,j1) = App(i1,j1) + tauU/h*SF_iP*SF_jP*wq

             !terme de stabilisation supplémentaire
             ! a2 < [w] , [B] >
             DO id = 1, nDim
                DO jd = 1, nDim
                   i1 = i + (id-1)*nNodesPerElement
                   j1 = j + (jd-1)*nNodesPerElement
                   Amm(i1,j1) = Amm(i1,j1) + tauP/h*SF_iM*SF_jM*Nm(id)*Nm(jd)*wq
                   Amp(i1,j1) = Amp(i1,j1) + tauP/h*SF_iM*SF_jP*Nm(id)*Np(jd)*wq
                   Apm(i1,j1) = Apm(i1,j1) + tauP/h*SF_iP*SF_jM*Np(id)*Nm(jd)*wq
                   App(i1,j1) = App(i1,j1) + tauP/h*SF_iP*SF_jP*Np(id)*Np(jd)*wq
                END DO
             END DO

          END DO
          ! - < {q} , j2 >
          i1 = i + nDim*nNodesPerElement
          F1(i1) = F1(i1) - 0.5d0*SF_iM*jn2*wq
          F2(i1) = F2(i1) - 0.5d0*SF_iP*jn2*wq

          ! a < [q] , jt >
          i1 = i + nDim*nNodesPerElement
          F1(i1) = F1(i1) - tauU/h*SF_iM*jt*wq
          F2(i1) = F2(i1) + tauU/h*SF_iP*jt*wq

          ! - < {w} , jt >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             F1(i1) = F1(i1) + 0.5d0*Nm(id)*SF_iM*jt*wq
             F2(i1) = F2(i1) - 0.5d0*Np(id)*SF_iP*jt*wq
          END DO

          ! a2 < [w] , j2 >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             F1(i1) = F1(i1) + tauP/h*SF_iM*Nm(id)*jn2*wq
             F2(i1) = F2(i1) + tauP/h*SF_iP*Np(id)*jn2*wq
          END DO
       END DO

       IF ( data%pressureEnrichment ) THEN
          DO i = 1, nNodesPerElement
             CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL SFEval(eleP%eleType, SF_iP, xq_P, i)
             DO ii = 1, nEdgesPerElement
                j = permut(ii,1) ; k = permut(ii,2)
                CALL SFEval(eleP%eleType, SF_jp, xq_P, j) ; CALL SFEval(eleP%eleType, SF_kp, xq_P, k)
                CALL SFEval(eleM%eleType, SF_jm, xq_M, j) ; CALL SFEval(eleM%eleType, SF_km, xq_M, k)

                ! < [ w ] , { T* } >
                DO id = 1, nDim
                   DO jd = 1, nDim
                      i1 = i + (id-1)*nNodesPerElement
                      j1 = j + (jd-1)*nNodesPerElement
                      k1 = k + (jd-1)*nNodesPerElement
                      int = 0.25d0*SF_iM*SF_jM*SF_kM*Nm(id)*wq
                      Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                      Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                      int = 0.25d0*SF_iM*SF_jP*SF_kP*Nm(id)*wq
                      Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                      Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int

                      int = 0.25d0*SF_iP*SF_jM*SF_kM*Np(id)*wq
                      Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                      Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                      int = 0.25d0*SF_iP*SF_jP*SF_kP*Np(id)*wq
                      App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                      App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                   END DO
                END DO

                ! tauU/h < [ q ] , [ T* ] >
                DO id = 1, nDim
                   i1 = i + nDim*nNodesPerElement
                   j1 = j + (id-1)*nNodesPerElement
                   k1 = k + (id-1)*nNodesPerElement
                   int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*wq
                   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
                   Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

                   int = -0.5d0*tauU/h*SF_iM*SF_jP*SF_kP*wq
                   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
                   Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int

                   int = -0.5d0*tauU/h*SF_iP*SF_jM*SF_kM*wq
                   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
                   Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

                   int = 0.5d0*tauU/h*SF_iP*SF_jP*SF_kP*wq
                   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
                   App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
                END DO

             END DO
          END DO
       END IF

       IF ( Data%TestFunctionEnrichment ) THEN
          DO ii = 1, nEdgesPerElement
             ii1 = permut(ii,1) ; ii2 = permut(ii,2)
             CALL SFEval(eleP%eleType, SF_i1p, xq_P, ii1) ; CALL SFEval(eleP%eleType, SF_i2p, xq_P, ii2)
             CALL SFEval(eleM%eleType, SF_i1m, xq_M, ii1) ; CALL SFEval(eleM%eleType, SF_i2m, xq_M, ii2)
             DO j = 1, nNodesPerElement
                CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL SFEval(eleM%eleType, SF_jM, xq_M, j)

                ! - < { q* } , [ B ] >
                DO id = 1, nDim
                   DO kd = 1, nDim
                      i1 = ii1 + (id-1)*nNodesPerElement
                      i2 = ii2 + (id-1)*nNodesPerElement
                      j1 = j + (kd-1)*nNodesPerElement
                      int = -0.25d0*SF_i1m*SF_i2m*SF_jm*Nm(kd)*wq
                      Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      int = -0.25d0*SF_i1m*SF_i2m*SF_jp*Np(kd)*wq
                      Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      int = -0.25d0*SF_i1p*SF_i2p*SF_jm*Nm(kd)*wq
                      Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      int = -0.25d0*SF_i1p*SF_i2p*SF_jp*Np(kd)*wq
                      App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                   END DO
                END DO

                ! tauU/h < [ q* ] , [ T ] >
                DO id = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   j1 = j + nDim*nNodesPerElement

                   int = 0.5d0*tauU/h*SF_i1m*SF_i2m*SF_jm*wq
                   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                   Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                   int = -0.5d0*tauU/h*SF_i1m*SF_i2m*SF_jp*wq
                   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                   Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                   int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*wq
                   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                   Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

                   int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*wq
                   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                   App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                END DO
             END DO

             ! - < {q*} , j2 >
             DO id= 1, nDim
                i1 = ii1 + (id-1)*nNodesPerElement
                i2 = ii2 + (id-1)*nNodesPerElement
                int=-0.25d0*SF_i1M*SF_i2m*jn2*wq
                F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                int=-0.25d0*SF_i1p*SF_i2p*jn2*wq
                F2(i1) = F2(i1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                F2(i2) = F2(i2) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
             END DO

             ! a < [q*] , jt >
             DO id= 1, nDim
                i1 = ii1 + (id-1)*nNodesPerElement
                i2 = ii2 + (id-1)*nNodesPerElement

                int=-0.5d0*tauU/h*SF_i1M*SF_i2M*jt*wq
                F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                int=0.5d0*tauU/h*SF_i1P*SF_i2P*jt*wq
                F2(i1) = F2(i1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                F2(i2) = F2(i2) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
             END DO

          END DO
       END IF

       IF ( Data%TestFunctionEnrichment .AND. Data%PressureEnrichment ) THEN
          DO ii = 1, nEdgesPerElement
             ii1 = permut(ii,1) ; ii2 = permut(ii,2)
             CALL SFEval(eleP%eleType, SF_i1p, xq_P, ii1) ; CALL SFEval(eleP%eleType, SF_i2p, xq_P, ii2)
             CALL SFEval(eleM%eleType, SF_i1m, xq_M, ii1) ; CALL SFEval(eleM%eleType, SF_i2m, xq_M, ii2)
             DO jj = 1, nEdgesPerElement
                jj1 = permut(jj,1) ; jj2 = permut(jj,2)
                CALL SFEval(eleM%eleType, SF_j1M, xq_M, jj1) ; CALL SFEval(eleM%eleType, SF_j2M, xq_M, jj2)
                CALL SFEval(eleP%eleType, SF_j1P, xq_P, jj1) ; CALL SFEval(eleP%eleType, SF_j2P, xq_P, jj2)

                ! tauU/h < [ q* ] , [ T* ] >
                DO id = 1, nDim
                   DO jd = 1, nDim
                      i1 = ii1 + (id-1)*nNodesPerElement
                      i2 = ii2 + (id-1)*nNodesPerElement
                      j1 = jj1 + (jd-1)*nNodesPerElement
                      j2 = jj2 + (jd-1)*nNodesPerElement

                      int = 0.25d0*tauU/h*SF_i1m*SF_i2m*SF_j1m*SF_j2m*wq
                      Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM**2)*(XeM(ii1,id)-XeM(ii2,id)) &
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int
                      Amm(i1,j2) = Amm(i1,j2) - (1.d0/lambdaM**2)*(XeM(ii1,id)-XeM(ii2,id)) &
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int
                      Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM**2)*(XeM(ii1,id)-XeM(ii2,id)) &
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int
                      Amm(i2,j2) = Amm(i2,j2) + (1.d0/lambdaM**2)*(XeM(ii1,id)-XeM(ii2,id)) &
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int

                      int = -0.25d0*tauU/h*SF_i1m*SF_i2m*SF_j1p*SF_j2p*wq
                      Amp(i1,j1) = Amp(i1,j1) + (1.d0/(lambdaM*lambdaP))*(XeM(ii1,id)-XeM(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int
                      Amp(i1,j2) = Amp(i1,j2) - (1.d0/(lambdaM*lambdaP))*(XeM(ii1,id)-XeM(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int
                      Amp(i2,j1) = Amp(i2,j1) - (1.d0/(lambdaM*lambdaP))*(XeM(ii1,id)-XeM(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int
                      Amp(i2,j2) = Amp(i2,j2) + (1.d0/(lambdaM*lambdaP))*(XeM(ii1,id)-XeM(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int

                      int = -0.25d0*tauU/h*SF_i1p*SF_i2p*SF_j1m*SF_j2m*wq
                      Apm(i1,j1) = Apm(i1,j1) + (1.d0/(lambdaM*lambdaP))*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int
                      Apm(i1,j2) = Apm(i1,j2) - (1.d0/(lambdaM*lambdaP))*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int
                      Apm(i2,j1) = Apm(i2,j1) - (1.d0/(lambdaM*lambdaP))*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int
                      Apm(i2,j2) = Apm(i2,j2) + (1.d0/(lambdaM*lambdaP))*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeM(jj1,jd)-XeM(jj2,jd))*int

                      int = 0.25d0*tauU/h*SF_i1p*SF_i2p*SF_j1p*SF_j2p*wq
                      App(i1,j1) = App(i1,j1) + (1.d0/lambdaP**2)*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int
                      App(i1,j2) = App(i1,j2) - (1.d0/lambdaP**2)*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int
                      App(i2,j1) = App(i2,j1) - (1.d0/lambdaP**2)*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int
                      App(i2,j2) = App(i2,j2) + (1.d0/lambdaP**2)*(XeP(ii1,id)-XeP(ii2,id))&
                           *(XeP(jj1,jd)-XeP(jj2,jd))*int
                   END DO
                END DO
             END DO
          END DO
       END IF

    END DO !fin de la boucle sur les noeuds de quadrature


    DEALLOCATE ( Np, Nm, Xq_M, Xq_P, XqM, XqP, X_M, X_P, Exa_p, Exa_m)
    DEALLOCATE (XxqP,XxqM,XeM,XeP,Lm1_M,Lm1_P)

  END SUBROUTINE ComputePMMatrices_CG
  !======================================================================!


  !=====================================================================================!
  SUBROUTINE ComputePMMatrices_CG_Embedded(PhysicalInterface, mesh, edg, edgM, edgP, eleM, eleP, &
  Amm, Amp, Apm, App, F1, F2, Data, iLS, Sol_prec,Sol)
    !=====================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)    :: Edg
    type(Edge), DIMENSION(:), INTENT(IN)    :: edgM, edgP
    type(element)           , INTENT(IN)    :: eleM, eleP
    REAL*8, DIMENSION(:,:)  , INTENT(OUT)   :: Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)    , INTENT(OUT)   :: F1,F2
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: ILS
    type(MeshType)          , INTENT(IN)    :: Mesh
    type(SolStructure) 		, INTENT(IN)    :: Sol_prec,Sol
    type(VectorInterface)   , INTENT(IN)    :: PhysicalInterface 
    !----------------------------------------------------------
    type(edge) :: edgIntM, edgIntP
    type(element) :: ELE_inter
    !----------------------------------------------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE    :: Np, Nm, Xq_M, Xq_P, Exa_p, Exa_m
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE    :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE    :: XeM, XeP
    REAL*8, DIMENSION(:)    , ALLOCATABLE    :: XxqP, XxqM
    REAL*8, DIMENSION(data%nDim)             :: d_qP, n_qP, d_qM, n_qM
    REAL*8, DIMENSION(data%nDim-1)           :: tntM, tntP
    REAL*8, DIMENSION(data%nDim,data%nDim-1) :: t_qP, t_qM
    REAL*8, DIMENSION(data%nDim+1)           :: exa1_q, exa2_q, exa1_q_prec, exa2_q_prec
    REAL*8, DIMENSION(data%nDim)             :: x_qP, x_qM
    REAL*8, DIMENSION(:)    , ALLOCATABLE    :: GSF_iM, GSF_jM,GSF_iP, GSF_jP
    REAL*8, DIMENSION(:)    , ALLOCATABLE    :: GSF_kM, GSF_kp, Exa_p_prec, Exa_M_prec
	LOGICAL								     :: testjumpT
    !-------------------------------------------------------------
    REAL*8  :: Nrm, wq, jt, pD_M, pD_P, nnt,jt_bis, jn2_bis
    REAL*8  :: SF_iM, SF_iP, SF_jM, SF_jP, a, ttM, ttP,res, mu_qM, mu_qP
    REAL*8  :: SF_i1M, SF_i1P, SF_i2P, SF_i2M, SF_j1M, SF_j2P, SF_j1P,SF_j2M,alphaM,alphaP
    REAL*8  :: tauU, h, LS_qP, LS_qM, lambdaM, lambdaP,res1,pD_M_prec, pD_P_prec
    REAL*8  :: Sigma_bP, Sigma_bM, jn2, L, int, SF_kM, SF_kP,res2, jump,Jumpfx,JumpFy,Tm, meanT 
    !------------------------------------------------------------
    INTEGER :: iM, iP, i, jM, jP,iq, i1, j1, j, id, jd, ii
    INTEGER :: ii1,ii2, jj,jj1,jj2,k,i2,k1, kd ,ld,j2,Ns_1, Ns_2
    INTEGER :: nDim, nNodesPerElement, nQuad, nEdgesPerElement, nVar
    !--------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nVar = nDim + 1  ; nEdgesPerElement = eleM%nEdgesPerElement
    ALLOCATE ( Np(nDim), Nm(nDim), Xq_M(nDim), Xq_P(nDim))
    ALLOCATE (  XxqP(nDim), XxqM(nDim) ) ; ALLOCATE ( Exa_p(nVar), Exa_M(nVar)  )
    ALLOCATE ( Exa_p_prec(nVar), Exa_M_prec(nVar)  )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( XeM(nNodesPerElement,nDim), XeP(nNodesPerElement,nDim) )
    ALLOCATE ( GSF_iM(nDim), GSF_jM(nDim), GSF_iP(nDim), GSF_jP(nDim) )
    ALLOCATE ( GSF_kP(nDim), GSF_kM(nDim) )

	Ns_1 = edg%Vertex(1) ; Ns_2 = edg%Vertex(2)
	XeM = eleM%Coor ; XeP = eleP%Coor
	CALL IdentifyIFace(eleM, edg, iM, nDim) ; CALL IdentifyIFace(eleP, edg, iP, nDim)

	res=0.d0 ; res1=0.d0 ; res2=0.d0
	edgIntM = edgM(iM)  ; edgIntP = edgP(iP)
	lambdaM = eleM%lambda ; lambdaP = eleP%lambda
	
	Np  = -eleP%N(iP,:)
	Nrm = 0.d0
	DO id = 1, nDim
		Nrm = Nrm + Np(id)**2
	END DO
	Nrm = sqrt(Nrm) ; L  = Nrm
	Np  = Np/Nrm    ; Nm = -Np
	
	!value of the characteristic length 
	h = (eleM%A + eleP%A)/(2.d0*L)
	
	!try this with a 0 value 
	tauU = alphaIP1*Nrm !stab from the Nitsche's method coefficient which is problematic in 1D 
	
	CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim)
	CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)
	!xqM,xqP coordinates on the reference triangle 
	
	!Initialization 
	Amm = 0.d0 ; Amp = 0.d0 ; Apm = 0.d0 ; App = 0.d0
	F1  = 0.d0 ; F2 = 0.d0  ; nnt = 0.d0
	
	!number of quadrature points 
	CALL GetNQuad(edg%eleType, nQuad)

	DO iq = 1, nQuad !boucle sur le nombre de points de quadrature

		!Coordonnées des neouds de quadrature sur l'éléments de références 
		xq_M = xqM(iq,:)  ;  xq_P = xqP(iq,:)

		x_qP = 0.d0  ;  x_qM = 0.d0
		
		DO i = 1, nNodesPerElement
			CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL SFEval(eleP%eleType, SF_iP, xq_P, i)
			DO id = 1, nDim
				x_qP(id) = x_qP(id) + SF_iP*eleP%Coor(i,id)
				x_qM(id) = x_qM(id) + SF_iM*eleM%Coor(i,id)
			END DO
		END DO
     
		CALL ComputeNodalDistancePhysicalinter(PhysicalInterface, x_qP, d_qP, LS_qP, Data, iLS)
		CALL getNormalAndTangents(d_qP, n_qP, t_qP, nP, nnt, tntP, nDim)

		IF (nnt>1 .OR. nnt<0) THEN
			PRINT*,"Produit des deux normales", eleP%num,eleP%ref,eleM%Num, eleM%ref,"value", nnt
			PRINT*, "Problème sur le calcul des normales" 
			!Stop
		END IF

		n_qM = -n_qP
		d_qM = -d_qP
		
		!condition sur le flux, condition de stefan 
		!IF(icas==001) THEN 
		!	jn2 = Data%Jump 
		!	IF(DATA%TSCHEME=="NICOLSON") THEN
		!		!Reconstruction on the jump for the edg with the new surrogate at time n+1
		!		jn2 = (jn2+Data%PrevJumpSurrogate)*0.5d0			!real reconstruction without the analytical solution 
		!	END IF
		!ELSE 
		IF (icas /=  0) THEN
			CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP ,Data%T, eleP, nDim)
			CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM ,Data%T, eleM, nDim)
			jn2 = DOT_PRODUCT(exa1_q(1:nDim) - exa2_q(1:nDim), n_qP)

			IF(DATA%TSCHEME=="NICOLSON") THEN
				CALL ExactSolJumpt(eleP%L, exa1_q_prec, x_qP + d_qP ,Data%T-Data%Dt, eleP, nDim)
				CALL ExactSolJumpt(eleM%L, exa2_q_prec, x_qM - d_qM ,Data%T-Data%Dt, eleM, nDim)
				
				jn2 = 0.5d0*(jn2+DOT_PRODUCT(exa1_q_prec(1:nDim) - exa2_q_prec(1:nDim), n_qP))
				!CALL JumpPrevNewInter(Data, Mesh, Sol, Sol_prec,JumpFx,JumpFy) 	
				!jn2 = DOT_PRODUCT(exa1_q(1:nDim) - exa2_q(1:nDim), n_qP) 
				!jn2 = (jn2 + JumpFx*n_qP(1)+JumpFy*n_qP(2))*0.5
			END IF
		ELSE
			jn2 = Data%BCEmbVal(iLS)
		END IF		
		 
		!probleme ici sur le calcul du saut au temps précédent 
		!on a bien le saut de temperature nul pour le cas 001
		!IF(icas==001) THEN 
		!	jt = 0.d0 !value at time n+1 continuity throught the two phases 
		!	IF(Data%TSCHEME=="NICOLSON") THEN 
		!		!computation for time n on the interface n+1 
		!		jump = 0.d0 
		!		CALL ValJumpT(Data, Mesh, Sol_prec, edgIntP, jump,MeanT, d_qP, iq,lambdaM,lambdaP) 
		!		jt = (jt + jump)*0.5d0
		!	
		!	END IF 
        !ELSE  
        
        IF (icas /= 0) THEN
			CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T, eleP, nDim)
			CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T, eleM, nDim)
			pD_P = exa1_q(nDim+1) ; pD_M = exa2_q(nDim+1)
			jt =pD_P - pD_M
	
			IF(DATA%TSCHEME=="NICOLSON") THEN

				CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T-Data%Dt, eleP, nDim)
				CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T-Data%Dt, eleM, nDim)
				pD_P_prec = exa1_q(nDim+1) ; pD_M_prec = exa2_q(nDim+1)
				jt_bis = pD_P_prec - pD_M_prec
				jt = 0.5d0*(jt+ jt_bis )
				
				!computation for time n on the interface n+1 
				!jump = 0.d0
				!CALL ValJumpT(Data, Mesh, Sol_prec, edgIntP, jump,meanT, d_qP, iq,lambdaM,lambdaP) 
				!jt = (jt + jump)*0.5d0
			END IF
		ELSE
			jt = Data%BCEmbVal(iLS)
		END IF

		CALL WeightQuad(edg%eleType, wq, iq)
		wq = wq*L
		
		! terme a l'interface
		DO i = 1, nNodesPerElement
			CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL GSFEval(GSF_iM, xq_M, i,eleM)
			CALL SFEval(eleP%eleType, SF_iP, xq_P, i) ; CALL GSFEval(GSF_iP, xq_P, i,eleP)

			DO j = 1, nNodesPerElement
				CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL GSFEval(GSF_jM, xq_M, j,eleM)
				CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL GSFEval(GSF_jP, xq_P, j,eleP)

				! < [ w ].~n , { T } > 
				DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Amm(i1,j1) = Amm(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jM*wq
					Amp(i1,j1) = Amp(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jP*wq
					Apm(i1,j1) = Apm(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jM*wq
					App(i1,j1) = App(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jP*wq
				END DO
				
				!< {w}.~n , [Lm B.d ] >
				DO id = 1, nDim
					DO jd = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement
						Amm(i1,j1) = Amm(i1,j1) - 0.5d0*nM(id)*SF_iM*(1.d0/lambdaM)*SF_jM*d_qM(jd)*wq
						Amp(i1,j1) = Amp(i1,j1) - 0.5d0*nM(id)*SF_iM*(1.d0/lambdaP)*SF_jP*d_qP(jd)*wq
						Apm(i1,j1) = Apm(i1,j1) + 0.5d0*nP(id)*SF_iP*(1.d0/lambdaM)*SF_jM*d_qM(jd)*wq
						App(i1,j1) = App(i1,j1) + 0.5d0*nP(id)*SF_iP*(1.d0/lambdaP)*SF_jP*d_qP(jd)*wq
					END DO
				END DO
				
				! - < { q } , [ B ].n(n.~n) >
				DO id = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (id-1)*nNodesPerElement
					Amm(i1,j1) = Amm(i1,j1) - 0.5d0*SF_iM*SF_jM*n_qM(id)*nnt*wq
					Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*SF_jP*n_qP(id)*nnt*wq
					Apm(i1,j1) = Apm(i1,j1) - 0.5d0*SF_iP*SF_jM*n_qM(id)*nnt*wq
					App(i1,j1) = App(i1,j1) - 0.5d0*SF_iP*SF_jP*n_qP(id)*nnt*wq
				END DO

				! - < { q } , [ DBd ].n(n.~n) >
				DO id = 1, nDim
					DO jd=1, nDim
						i1 = i + nDim*nNodesPerElement
						j1 = j + (id-1)*nNodesPerElement
						Amm(i1,j1) = Amm(i1,j1) + 0.5d0*SF_iM*N_qm(id)*GSF_jM(jd)*d_qM(jd)*nnt*wq
						Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*N_qp(id)*GSF_jP(jd)*d_qP(jd)*nnt*wq
						Apm(i1,j1) = Apm(i1,j1) + 0.5d0*SF_iP*N_qm(id)*GSF_jM(jd)*d_qM(jd)*nnt*wq
						App(i1,j1) = App(i1,j1) - 0.5d0*SF_iP*N_qp(id)*GSF_jP(jd)*d_qP(jd)*nnt*wq
					END DO
				END DO


				!ajout d'un terme de stabilisation sur la temperature
				!a<[q + Dq . d],[T + DT.d]>
				!on decompose en 4 parties

				!a<[q],[T]>
				i1 = i + nDim*nNodesPerElement
				j1 = j + nDim*nNodesPerElement
				Amm(i1,j1) = Amm(i1,j1) + tauU/h*SF_iM*SF_jM*wq
				Amp(i1,j1) = Amp(i1,j1) - tauU/h*SF_iM*SF_jP*wq
				Apm(i1,j1) = Apm(i1,j1) - tauU/h*SF_iP*SF_jM*wq
				App(i1,j1) = App(i1,j1) + tauU/h*SF_iP*SF_jP*wq
				
				!- a* < [q] ,[ LM B.d] >
				DO jd = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (jd-1)*nNodesPerElement

					int= tauU/h*(1.d0/lambdaM)*SF_iM*wq
					Amm(i1,j1) = Amm(i1,j1) + SF_jM*d_qM(jd)*int
					int= tauU/h*(1.d0/lambdaP)*SF_iM*wq
					Amp(i1,j1) = Amp(i1,j1) + SF_jP*d_qP(jd)*int

					int= tauU/h*(1.d0/lambdaM)*SF_iP*wq
					Apm(i1,j1) = Apm(i1,j1) - SF_jM*d_qM(jd)*int

					int= tauU/h*(1.d0/lambdaP)*SF_iP*wq
					App(i1,j1) = App(i1,j1) - SF_jP*d_qP(jd)*int
				END DO

				!- a* < [ LM w.d] ,[T] >
				DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + nDim*nNodesPerElement

					int= tauU/h*(1.d0/lambdaM)*SF_jM*wq
					Amm(i1,j1) = Amm(i1,j1) + SF_iM*d_qM(id)*int

					int= tauU/h*(1.d0/lambdaM)*SF_jP*wq
					Amp(i1,j1) = Amp(i1,j1) - SF_iM*d_qM(id)*int

					int= tauU/h*(1.d0/lambdaP)*SF_jM*wq
					Apm(i1,j1) = Apm(i1,j1) + SF_iP*d_qP(id)*int

					int= tauU/h*(1.d0/lambdaP)*SF_jP*wq
					App(i1,j1) = App(i1,j1) - SF_iP*d_qP(id)*int
				END DO

				! tauU/h* < [LM w.d], [LM  B.D ] >
				DO id = 1,nDim
					DO jd = 1,nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement

						int =tauU/h*(1.d0/(lambdaM**2))*wq
						Amm(i1,j1) = Amm(i1,j1) + SF_iM*d_qM(id)*SF_jM*d_qM(jd)*int

						int =tauU/h*(1.d0/(lambdaM*lambdaP))*wq
						Amp(i1,j1) = Amp(i1,j1) + SF_iM*d_qM(id)*SF_jP*d_qP(jd)*int

						int =tauU/h*(1.d0/(lambdaP*lambdaM))*wq
						Apm(i1,j1) = Apm(i1,j1) + SF_iP*d_qP(id)*SF_jM*d_qM(jd)*int

						int =tauU/h*(1.d0/(lambdaP**2))*wq
						App(i1,j1) = App(i1,j1) + SF_iP*d_qP(id)*SF_jP*d_qP(jd)*int
					END DO
				END DO

				!terme supplémentaire issu du ddl
				!on ajoute tout de meme au non enrichie pour comparaison
				! 0.5* < { w }.~n , [ d^t ( LM * DB) d ] > 
				DO id = 1, nDim
					DO jd=1, nDim
						DO kd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) + 0.25d0*SF_iM*Nm(id)*(1.d0/lambdaM)* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							Amp(i1,j1) = Amp(i1,j1) - 0.25d0*SF_iM*Nm(id)*(1.d0/lambdaP)* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
							Apm(i1,j1) = Apm(i1,j1) - 0.25d0*SF_iP*Np(id)*(1.d0/lambdaM)* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							App(i1,j1) = App(i1,j1) + 0.25d0*SF_iP*Np(id)*(1.d0/lambdaP)* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
						END DO
					END DO
				END DO


				!- 0.5*a* < [q] , d^t[ LM * DB ]d >
				DO jd=1, nDim
					DO kd=1, nDim
						i1 = i + nDim*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement

						Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_iM*(1.d0/lambdaM)* &
						GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
						Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_iM*(1.d0/lambdaP)* &
						GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
						Apm(i1,j1) = Apm(i1,j1) + tauU/h*0.5d0*SF_iP*(1.d0/lambdaM)* &
						GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
						App(i1,j1) = App(i1,j1) - tauU/h*0.5d0*SF_iP*(1.d0/lambdaP)* &
						GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
					END DO
				END DO

				! 0.5*a* < [Lm w.d] , [d^t (LM * DB) d ] >
				DO id=1, nDim
					DO jd=1, nDim
						DO kd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM**2))* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM*lambdaP))* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
							Apm(i1,j1) = Apm(i1,j1) - tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP*lambdaM))* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							App(i1,j1) = App(i1,j1) + tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP**2))* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
						END DO
					END DO
				END DO

				!- 0.5*a* < [d^t (LM * Dw) d ],[T] >
				DO id = 1, nDim
					DO kd = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement

						Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_jM*(1.d0/lambdaM)* &
						GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

						Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_jP*(1.d0/lambdaM)* &
						GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

						Apm(i1,j1) = Apm(i1,j1) + tauU/h*0.5d0*SF_jM*(1.d0/lambdaP)* &
						GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq

						App(i1,j1) = App(i1,j1) - tauU/h*0.5d0*SF_jP*(1.d0/lambdaP)* &
						GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
					END DO
				END DO

				! 0.5*a* <[d^t (LM * Dw ) d],[Lm B.d] >
				DO id=1, nDim
					DO kd=1,nDim
						DO jd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM**2))* &
							GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
							Amp(i1,j1) = Amp(i1,j1) - tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP*lambdaM))* &
							GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
							Apm(i1,j1) = Apm(i1,j1) + tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM*lambdaP))* &
							GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
							App(i1,j1) = App(i1,j1) + tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP**2))* &
							GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
						END DO
					END DO
				END DO

				! 0.25*a* <[d^t (LM * Dw) d], [d^t (Lm DB) d] >
				DO id=1, nDim
					DO kd=1, nDim
						DO jd=1, nDim
							DO ld=1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement

								Amm(i1,j1) = Amm(i1,j1) + tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
								(1.d0/(lambdaM**2))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
								Amp(i1,j1) = Amp(i1,j1) - tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
								(1.d0/(lambdaM*lambdaP))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
								Apm(i1,j1) = Apm(i1,j1) - tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
								(1.d0/(lambdaP*lambdaM))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
								App(i1,j1) = App(i1,j1) + tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
							   (1.d0/(lambdaP**2))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
							END DO
						END DO
					END DO
				END DO
            
			END DO

			! - < {q} , j2(n.~n)>
			i1 = i + nDim*nNodesPerElement
			F1(i1) = F1(i1) - 0.5d0*SF_iM*jn2*nnt*wq
			F2(i1) = F2(i1) - 0.5d0*SF_iP*jn2*nnt*wq

			! - < {w}.n , jt >
			DO id = 1, nDim
				i1 = i + (id-1)*nNodesPerElement
				F1(i1) = F1(i1) + 0.5d0*nM(id)*SF_iM*jt*wq
				F2(i1) = F2(i1) - 0.5d0*nP(id)*SF_iP*jt*wq
			END DO

				
			!ajout du terme sur le second membre issu de la stabilisation
			!a<[q+Dq.d],jT>

			!a<[q],jT>
			i1 = i + nDim*nNodesPerElement
			F1(i1) = F1(i1) - tauU/h*SF_iM*jt*wq
			F2(i1) = F2(i1) + tauU/h*SF_iP*jt*wq

			!terme pour le second membre
            !- a* < [ LM w.d] ,jT >
            DO id=1, nDim
               i1 = i + (id-1)*nNodesPerElement

               int= tauU/h*(1.d0/lambdaM)*jT*wq
               F1(i1) = F1(i1) - SF_iM*d_qM(id)*int

               int= tauU/h*(1.d0/lambdaP)*jT*wq
               F2(i1) = F2(i1) - SF_iP*d_qP(id)*int
            END DO

            !- 0.5*a* < d^t[ LM * Dw d],jT >
            DO id=1, nDim
				DO kd=1,nDim
					i1 = i + (id-1)*nNodesPerElement

					res= GSF_iM(kd)*d_qM(kd)*d_qM(id)
					F1(i1) = F1(i1) + tauU/h*0.5d0*(1.d0/lambdaM)*res*jT*wq
					res= GSF_iP(kd)*d_qP(kd)*d_qP(id)
					F2(i1) = F2(i1) - tauU/h*0.5d0*(1.d0/lambdaP)*res*jT*wq
				END DO
            END DO
            
		END DO

		IF ( Data%PressureEnrichment ) THEN
			DO i = 1, nNodesPerElement
				CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL GSFEval(GSF_iM, xq_M, i,eleM)
				CALL SFEval(eleP%eleType, SF_iP, xq_P, i) ; CALL GSFEval(GSF_iP, xq_P, i,eleP)

				DO ii = 1, nEdgesPerElement
					j = permut(ii,1) ; k = permut(ii,2)
					CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL SFEval(eleP%eleType, SF_kP, xq_P, k)
					CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL SFEval(eleM%eleType, SF_kM, xq_M, k)
					
					CALL GSFEval(GSF_jM, xq_M, j,eleM) ; CALL GSFEval(GSF_kM, xq_M, k,eleM) !!new term 
					CALL GSFEval(GSF_jP, xq_P, j,eleP) ; CALL GSFEval(GSF_kP, xq_P, k,eleP) !!new term
					

					! < [ w ].~n , { T* } >
					DO id = 1, nDim
					   DO jd = 1, nDim
						  i1 = i + (id-1)*nNodesPerElement
						  j1 = j + (jd-1)*nNodesPerElement
						  k1 = k + (jd-1)*nNodesPerElement
						  int = 0.25d0*SF_iM*Nm(id)*SF_jM*SF_kM*wq
						  Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
						  Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

						  int = 0.25d0*SF_iM*Nm(id)*SF_jP*SF_kP*wq
						  Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
						  Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int

						  int = 0.25d0*SF_iP*Np(id)*SF_jM*SF_kM*wq
						  Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
						  Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

						  int = 0.25d0*SF_iP*Np(id)*SF_jP*SF_kP*wq
						  App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
						  App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
					   END DO
					END DO
				

					! tauU/h < [ q ] , [ T* ] >
					DO id = 1, nDim
					   i1 = i + nDim*nNodesPerElement
					   j1 = j + (id-1)*nNodesPerElement
					   k1 = k + (id-1)*nNodesPerElement
					   int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*wq
					   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
					   Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

					   int = -0.5d0*tauU/h*SF_iM*SF_jP*SF_kP*wq
					   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					   Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int

					   int = -0.5d0*tauU/h*SF_iP*SF_jM*SF_kM*wq
					   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
					   Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

					   int = 0.5d0*tauU/h*SF_iP*SF_jP*SF_kP*wq
					   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					   App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					END DO

					! -tauU/h* < LM [w.d],[ T* ] >
					DO id=1,nDim
						DO jd=1,nDim
						
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement
							k1 = k + (jd-1)*nNodesPerElement
							int =0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jM*SF_kM*wq
							Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							int =-0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jP*SF_kP*wq
							Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int


							int =0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jM*SF_kM*wq
							Apm(i1,j1) = Apm(i1,j1) +  (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							Apm(i1,k1) = Apm(i1,k1) -  (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							int =-0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jP*SF_kP*wq
							App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							  
						END DO
					END DO

					! -tauU/h*0.5 < LM [d^t Dw.d],[ T* ] >
					DO id = 1, nDim
						DO kd = 1, nDim
							DO jd = 1, nDim
							
							   i1 = i + (id-1)*nNodesPerElement
							   j1 = j + (jd-1)*nNodesPerElement
							   k1 = k + (jd-1)*nNodesPerElement

							   res= (1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
							   int= -tauU/h*0.25d0*res*SF_jM*SF_kM*wq
							   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							   Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							   res= (1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
							   int= tauU/h*0.25d0*res*SF_jP*SF_kP*wq
							   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int

							   res= (1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
							   int= tauU/h*0.25d0*res*SF_jM*SF_kM*wq
							   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							   Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							   res= (1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
							   int= -tauU/h*0.25d0*res*SF_jP*SF_kP*wq
							   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   
							END DO
						END DO
					END DO

				END DO
			END DO
		END IF

		IF ( Data%TestFunctionEnrichment ) THEN
			DO ii = 1, nEdgesPerElement
				ii1 = permut(ii,1) ; ii2 = permut(ii,2)
				CALL SFEval(eleP%eleType, SF_i1P, xq_P, ii1) ; CALL SFEval(eleP%eleType, SF_i2P, xq_P, ii2)
				CALL SFEval(eleM%eleType, SF_i1M, xq_M, ii1) ; CALL SFEval(eleM%eleType, SF_i2M, xq_M, ii2)
				DO j = 1, nNodesPerElement
					CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL GSFEval(GSF_jM, xq_M, j,eleM)
					CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL GSFEval(GSF_jP, xq_P, j,eleP)

					! - < { q* } , [ B ].n(n.~n) >
					DO id = 1, nDim
					    DO kd = 1, nDim
						  i1 = ii1 + (id-1)*nNodesPerElement
						  i2 = ii2 + (id-1)*nNodesPerElement
						  j1 = j   + (kd-1)*nNodesPerElement

						  int = -0.25d0*SF_i1M*SF_i2M*SF_jM*N_qM(kd)*nnt*wq
						  Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						  Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						  int = -0.25d0*SF_i1M*SF_i2M*SF_jP*N_qP(kd)*nnt*wq
						  Amp(i1,j1)  = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						  Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						  int = -0.25d0*SF_i1P*SF_i2P*SF_jM*N_qM(kd)*nnt*wq
						  Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						  Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						  int = -0.25d0*SF_i1P*SF_i2P*SF_jP*N_qP(kd)*nnt*wq
						  App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						  App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					    END DO
					END DO

					! - < { q* } , [ DB.d ].n(n.~n) >
					DO id = 1, nDim
					   DO kd = 1, nDim
						 DO jd= 1,nDim
							i1 = ii1 + (id-1)*nNodesPerElement
							i2 = ii2 + (id-1)*nNodesPerElement
							j1 = j   + (kd-1)*nNodesPerElement
							int = 0.25d0*SF_i1M*SF_i2M*GSF_jM(jd)*d_qM(jd)*N_qM(kd)*nnt*wq
							Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
							Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
							int = -0.25d0*SF_i1M*SF_i2M*GSF_jP(jd)*N_qP(jd)*d_qP(kd)*nnt*wq
							Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
							Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
							int = 0.25d0*SF_i1P*SF_i2P*GSF_jM(jd)*N_qM(jd)*d_qM(kd)*nnt*wq
							Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
							Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
							int = -0.25d0*SF_i1P*SF_i2P*GSF_jP(jd)*N_qP(jd)*d_qP(kd)*nnt*wq
							App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
							App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						 END DO
					   END DO
					END DO

					! tauU/h < [ q* ] , [ T ] >
					DO id = 1, nDim
					   i1 = ii1 + (id-1)*nNodesPerElement
					   i2 = ii2 + (id-1)*nNodesPerElement
					   j1 = j   + nDim*nNodesPerElement

					   int = 0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*wq
					   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					   Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					   int = -0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*wq
					   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					   Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					   int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*wq
					   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					   Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

					   int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*wq
					   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					   App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					END DO

					! -tauU/h < [ q* ] , [ LM B.d ] >
					DO id = 1, nDim
					  DO jd= 1, nDim
						 i1 = ii1 + (id-1)*nNodesPerElement
						 i2 = ii2 + (id-1)*nNodesPerElement
						 j1 = j   + (jd-1)*nNodesPerElement

						 int = 0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
						 Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						 Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

						 int = 0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
						 Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						 Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

						 int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
						 Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						 Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

						 int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
						 App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						 App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					  END DO
					END DO

					! -tauU/h*0.5 < [ q* ],LM [d^t DB.d] >
					DO id =1, nDim
						DO jd =1, nDim
							DO kd = 1,nDim
							
								i1 = ii1 + (id-1)*nNodesPerElement
								i2 = ii2 + (id-1)*nNodesPerElement
								j1 = j   + (jd-1)*nNodesPerElement

								res=(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
								int=-tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
								Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
								Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

								res=(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
								int= tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
								Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
								Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

								res=(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
								int= tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
								Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

								res=(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
								int= -tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
								App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								
							END DO
						END DO
					END DO


				END DO

				! - < {q*} , j2(n.~n) >
				DO id= 1, nDim
					i1 = ii1 + (id-1)*nNodesPerElement
					i2 = ii2 + (id-1)*nNodesPerElement

					int=-0.25d0*SF_i1M*SF_i2M*jn2*nnt*wq
					F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					int=-0.25d0*SF_i1P*SF_i2P*jn2*nnt*wq
					F2(i1) = F2(i1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					F2(i2) = F2(i2) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
				END DO

				! a < [q*] , jt >
				DO id= 1, nDim
					i1 = ii1 + (id-1)*nNodesPerElement
					i2 = ii2 + (id-1)*nNodesPerElement

					int=-0.5d0*tauU/h*SF_i1M*SF_i2M*jt*wq
					F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					int=0.5d0*tauU/h*SF_i1P*SF_i2P*jt*wq
					F2(i1) = F2(i1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					F2(i2) = F2(i2) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
				END DO

			END DO
		END IF

		IF ( Data%TestFunctionEnrichment .AND. Data%PressureEnrichment ) THEN
			DO ii = 1, nEdgesPerElement
			
				ii1 = permut(ii,1) ; ii2 = permut(ii,2)
				CALL SFEval(eleP%eleType, SF_i1P, xq_P, ii1) ; CALL SFEval(eleP%eleType, SF_i2P, xq_P, ii2)
				CALL SFEval(eleM%eleType, SF_i1M, xq_M, ii1) ; CALL SFEval(eleM%eleType, SF_i2M, xq_M, ii2)
				DO jj = 1, nEdgesPerElement
				
					jj1 = permut(jj,1) ; jj2 = permut(jj,2)
					CALL SFEval(eleM%eleType, SF_j1M, xq_M, jj1) ; CALL SFEval(eleM%eleType, SF_j2M, xq_M, jj2)
					CALL SFEval(eleP%eleType, SF_j1P, xq_P, jj1) ; CALL SFEval(eleP%eleType, SF_j2P, xq_P, jj2)

					! tauU/h < [ q* ] , [ T* ] >
					DO id = 1, nDim
						DO jd = 1, nDim
							i1 = ii1 + (id-1)*nNodesPerElement
							i2 = ii2 + (id-1)*nNodesPerElement
							j1 = jj1 + (jd-1)*nNodesPerElement
							j2 = jj2 + (jd-1)*nNodesPerElement

							int = 0.25d0*(1.d0/(lambdaM**2))*tauU/h*SF_i1M*SF_i2M*SF_j1M*SF_j2M*wq
							Amm(i1,j1) = Amm(i1,j1)+(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i1,j2) = Amm(i1,j2)-(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i2,j1) = Amm(i2,j1)-(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i2,j2) = Amm(i2,j2)+(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

							int = -0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1M*SF_i2M*SF_j1P*SF_j2P*wq
							Amp(i1,j1) = Amp(i1,j1) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i1,j2) = Amp(i1,j2) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i2,j1) = Amp(i2,j1) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i2,j2) = Amp(i2,j2) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int

							int = -0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1P*SF_i2P*SF_j1M*SF_j2M*wq
							Apm(i1,j1) = Apm(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i1,j2) = Apm(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i2,j1) = Apm(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i2,j2) = Apm(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

							int = 0.25d0*(1.d0/(lambdaP**2))*tauU/h*SF_i1P*SF_i2P*SF_j1P*SF_j2P*wq
							App(i1,j1) = App(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i1,j2) = App(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i2,j1) = App(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i2,j2) = App(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
						END DO
					END DO

				END DO
			END DO
		END IF

    END DO
    
	
    DEALLOCATE ( Np, Nm, Xq_M, Xq_P, XqM, XqP, X_M, X_P, Exa_p, Exa_m)
    DEALLOCATE ( XeP, XeM )
    DEALLOCATE (GSF_iM, GSF_jM, GSF_iP, GSF_jP, GSF_kM, GSF_kP)

 END SUBROUTINE ComputePMMatrices_CG_Embedded
 !======================================================================!

 !=====================================================================================!
  SUBROUTINE ComputePMMatrices_CG_Embedded2(PhysicalInterface, Mesh, edg, edgM, edgP, eleM, eleP,&
   Amm, Amp, Apm, App, F1, F2, Data, iLS, Sol_prec, Sol)
 !=====================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)    :: Edg
    type(Edge), DIMENSION(:), INTENT(IN)    :: edgM, edgP
    type(element)           , INTENT(IN)    :: eleM, eleP
    REAL*8, DIMENSION(:,:)  , INTENT(OUT)   :: Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)    , INTENT(OUT)   :: F1,F2
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: ILS
    type(MeshType)          , INTENT(IN)    :: Mesh
    type(SolStructure) 		, INTENT(IN)    :: Sol_prec,Sol
    type(VectorInterface)   , INTENT(IN)    :: PhysicalInterface 
    !----------------------------------------------------------
    type(edge)    :: edgIntM, edgIntP
    type(element) :: ELE_inter
    !----------------------------------------------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE 	 :: Np, Nm, Xq_M, Xq_P, Exa_p, Exa_m
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE 	 :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE 	 :: XeM, XeP
    REAL*8, DIMENSION(:)    , ALLOCATABLE 	 :: XxqP, XxqM
    REAL*8, DIMENSION(data%nDim)  		  	 :: d_qP, n_qP, d_qM, n_qM
    REAL*8, DIMENSION(data%nDim-1) 		  	 :: tntM, tntP
    REAL*8, DIMENSION(data%nDim,data%nDim-1) :: t_qP, t_qM
    REAL*8, DIMENSION(data%nDim+1) 			 :: exa1_q, exa2_q, exa1_q_prec, exa2_q_prec
    REAL*8, DIMENSION(data%nDim)   			 :: x_qP, x_qM
    REAL*8, DIMENSION(:)    , ALLOCATABLE 	 :: GSF_iM, GSF_jM,GSF_iP, GSF_jP
    REAL*8, DIMENSION(:)    , ALLOCATABLE 	 :: GSF_kM, GSF_kp, Exa_p_prec, Exa_M_prec
	LOGICAL								  	 :: testjumpT
    !-------------------------------------------------------------
    REAL*8  :: Nrm, wq, jt, pD_M, pD_P, nnt,jt_bis, jn2_bis
    REAL*8  :: SF_iM, SF_iP, SF_jM, SF_jP, a, ttM, ttP,res, mu_qM, mu_qP
    REAL*8  :: SF_i1M, SF_i1P, SF_i2P, SF_i2M, SF_j1M, SF_j2P, SF_j1P,SF_j2M,alphaM,alphaP
    REAL*8  :: tauU, h, LS_qP, LS_qM, lambdaM, lambdaP,res1,pD_M_prec, pD_P_prec
    REAL*8  :: Sigma_bP, Sigma_bM, jn2, L, int, SF_kM, SF_kP,res2, jump,Jumpfx,JumpFy,Tm, MeanT
    !------------------------------------------------------------
    INTEGER :: iM, iP, i, jM, jP,iq, i1, j1, j, id, jd, ii
    INTEGER :: ii1,ii2, jj,jj1,jj2,k,i2,k1, kd ,ld,j2,Ns_1, Ns_2
    INTEGER :: nDim, nNodesPerElement, nQuad, nEdgesPerElement, nVar
    !--------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nVar = nDim + 1  ; nEdgesPerElement = eleM%nEdgesPerElement
    
    ALLOCATE (Np(nDim),Nm(nDim), Xq_M(nDim), Xq_P(nDim))
    ALLOCATE (XxqP(nDim), XxqM(nDim)) ; ALLOCATE (Exa_p(nVar), Exa_M(nVar))
    ALLOCATE (Exa_p_prec(nVar), Exa_M_prec(nVar))
    ALLOCATE (XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim))
    ALLOCATE (X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim))
    ALLOCATE (XeM(nNodesPerElement,nDim), XeP(nNodesPerElement,nDim))
    ALLOCATE (GSF_iM(nDim), GSF_jM(nDim), GSF_iP(nDim), GSF_jP(nDim))
    ALLOCATE (GSF_kP(nDim), GSF_kM(nDim))

	Ns_1 = edg%Vertex(1) ; Ns_2 = edg%Vertex(2)
	XeM = eleM%Coor      ; XeP = eleP%Coor
	CALL IdentifyIFace(eleM, edg, iM, nDim) ; CALL IdentifyIFace(eleP, edg, iP, nDim)

	res=0.d0 ; res1=0.d0  ; res2    = 0.d0
	edgIntM = edgM(iM)    ; edgIntP = edgP(iP)
	lambdaM = eleM%lambda ; lambdaP = eleP%lambda
	
	Np  = -eleP%N(iP,:)
	Nrm = 0.d0
	DO id = 1, nDim
		Nrm = Nrm + Np(id)**2
	END DO
	Nrm = sqrt(Nrm) ;  L  = Nrm
	Np  = Np/Nrm    ;  Nm = -Np
	
	!value of the characteristic length 
	h = (eleM%A + eleP%A)/(2.d0*L)
	
	!try this with a 0 value 
	tauU = alphaIP1!*Nrm !stab from the Nitsche's method coefficient which is problematic in 1D 
	
	CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim)
	CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)
	!xqM,xqP coordinates on the reference triangle 
	
	!Initialization 
	Amm = 0.d0 ; Amp = 0.d0  ; Apm = 0.d0 ; App = 0.d0
	F1  = 0.d0 ; F2  = 0.d0  ; nnt = 0.d0
	
	!number of quadrature points 
	CALL GetNQuad(edg%eleType, nQuad)

	DO iq = 1, nQuad !boucle sur le nombre de points de quadrature

		!Coordonnées des neouds de quadrature sur l'éléments de références 
		xq_M = xqM(iq,:)  ;  xq_P = xqP(iq,:)

		x_qP = 0.d0  ;  x_qM = 0.d0
		
		DO i = 1, nNodesPerElement
			CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL SFEval(eleP%eleType, SF_iP, xq_P, i)
			DO id = 1, nDim
				x_qP(id) = x_qP(id) + SF_iP*eleP%Coor(i,id)
				x_qM(id) = x_qM(id) + SF_iM*eleM%Coor(i,id)
			END DO
		END DO
     
		CALL ComputeNodalDistancePhysicalinter(PhysicalInterface, x_qP, d_qP, LS_qP, Data, iLS)
		CALL getNormalAndTangents(d_qP, n_qP, t_qP, nP, nnt, tntP, nDim)

		IF (nnt>1 .OR. nnt<0) THEN
			PRINT*,"Produit des deux normales", eleP%num,eleP%ref,eleM%Num, eleM%ref,"value", nnt
			PRINT*, "Problème sur le calcul des normales" 
			!STOP
		END IF

		n_qM = -n_qP
		d_qM = -d_qP

        IF ( icas /= 0 ) THEN
        	!variable for the jump temperature 
			CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T, eleP, nDim)
			CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T, eleM, nDim)
			pD_P = exa1_q(nDim+1) ; pD_M = exa2_q(nDim+1)
			jt = pD_P - pD_M
			Tm = 0.5*(pD_P + pD_M)

			IF(DATA%TSCHEME=="NICOLSON") THEN

				CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T-Data%Dt, eleP, nDim)
				CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T-Data%Dt, eleM, nDim)
				
				pD_P_prec = exa1_q(nDim+1) ; pD_M_prec = exa2_q(nDim+1)
				!jt_bis = pD_P_prec - pD_M_prec
				
				CALL ValJumpT(Data, Mesh, Sol_prec, edgIntP, jt_bis, MeanT, d_qP, iq, lambdaM, lambdaP) 
				jt_bis = pD_P_prec - pD_M_prec
				jt = (jt + jt_bis)*0.5d0
				
				!PRINT*,"val Mean", MeanT, (pD_P_prec+pD_M_prec)*0.5d0, iq
				Tm = 0.5d0*(Tm+(pD_P_prec+pD_M_prec)*0.5d0)
				!Tm = 0.5d0*(Tm+MeanT)
				PRINT*,"Difference Jump",MeanT, (pD_P_prec+pD_M_prec)*0.5d0, iq
			END IF
		ELSE
			jt = Data%BCEmbVal(iLS)
			Tm = Icing_Tm !a changer pour ce cas précis 
		END IF
				
		CALL WeightQuad(edg%eleType, wq, iq)
		wq = wq*L
		
		! terme a l'interface
		DO i = 1, nNodesPerElement
			CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL GSFEval(GSF_iM, xq_M, i,eleM)
			CALL SFEval(eleP%eleType, SF_iP, xq_P, i) ; CALL GSFEval(GSF_iP, xq_P, i,eleP)

			DO j = 1, nNodesPerElement
				CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL GSFEval(GSF_jM, xq_M, j,eleM)
				CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL GSFEval(GSF_jP, xq_P, j,eleP)

				!  < [w].n, { Lm B.d } > !new term !! new term form with B 
				! < w+.n+ + w-.n- , (Lm+B+.d+ - Lm-B-d-)*0.5 >
				DO id = 1, nDim
					DO jd = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement

						Amm(i1,j1) = Amm(i1,j1) - 1.d0/lambdaM*SF_iM*nM(id)*SF_jM*d_qM(jd)*wq*0.5d0
						Amp(i1,j1) = Amp(i1,j1) + 1.d0/lambdaP*SF_iM*nM(id)*SF_jP*d_qP(jd)*wq*0.5d0
						Apm(i1,j1) = Apm(i1,j1) - 1.d0/lambdaM*SF_iP*nP(id)*SF_jM*d_qM(jd)*wq*0.5d0
						App(i1,j1) = App(i1,j1) + 1.d0/lambdaP*SF_iP*nP(id)*SF_jP*d_qP(jd)*wq*0.5d0
												
					END DO
				END DO 
				
				!  0.5d0*<[w].n, {d^t (Lm DB) d} > !new term !! be careful not the divergence operator but the tensor 
				!  0.5d0*<w+n+ + w_n_, 0.5*(d^t+ LmDB+ d+ + d^t- LmDB- d-)>
				DO id = 1, nDim
					DO jd = 1, nDim
						DO kd = 1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) + 0.25d0*nM(id)*SF_iM*GSF_jM(jd)*d_qM(jd)*d_qM(kd)*1.d0/lambdaM*wq			
							Amp(i1,j1) = Amp(i1,j1) + 0.25d0*nM(id)*SF_iM*GSF_jP(jd)*d_qP(jd)*d_qP(kd)*1.d0/lambdaP*wq
							Apm(i1,j1) = Apm(i1,j1) + 0.25d0*nP(id)*SF_iP*GSF_jM(jd)*d_qM(jd)*d_qM(kd)*1.d0/lambdaM*wq
							App(i1,j1) = App(i1,j1) + 0.25d0*nP(id)*SF_iP*GSF_jP(jd)*d_qP(jd)*d_qP(kd)*1.d0/lambdaP*wq
									
						END DO
					END DO
				END DO

				!< {w}.~n , [Lm B.d ] > 
				! < (w+n+ - w-n-)*0.5, LmP B+d+ + LmM B-d- > 
				DO id = 1, nDim
					DO jd = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement !vecteur distance correcte 
						Amm(i1,j1) = Amm(i1,j1) - 0.5d0*SF_iM*nM(id)*1.d0/lambdaM*SF_jM*d_qM(jd)*wq
						Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*nM(id)*1.d0/lambdaP*SF_jP*d_qP(jd)*wq
						Apm(i1,j1) = Apm(i1,j1) + 0.5d0*SF_iP*nP(id)*1.d0/lambdaM*SF_jM*d_qM(jd)*wq
						App(i1,j1) = App(i1,j1) + 0.5d0*SF_iP*nP(id)*1.d0/lambdaP*SF_jP*d_qP(jd)*wq
					END DO
				END DO
				
				
				!terme supplémentaire issu du ddl
				! 0.5* < { w }.~n , [ d^t ( LM * DB) d ] > 
				! 0.5* < (w+n+ - w-n-)*0.5 , d^t+ LmDB+ d+ - d^t- lmDB- d->
				DO id = 1, nDim
					DO jd=1, nDim
						DO kd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) + 0.25d0*SF_iM*Nm(id)*(1.d0/lambdaM)*GSF_jM(jd)*d_qM(jd)*d_qM(kd)*wq
							Amp(i1,j1) = Amp(i1,j1) - 0.25d0*SF_iM*Nm(id)*(1.d0/lambdaP)*GSF_jP(jd)*d_qP(jd)*d_qP(kd)*wq
							Apm(i1,j1) = Apm(i1,j1) - 0.25d0*SF_iP*Np(id)*(1.d0/lambdaM)*GSF_jM(jd)*d_qM(jd)*d_qM(kd)*wq
							App(i1,j1) = App(i1,j1) + 0.25d0*SF_iP*Np(id)*(1.d0/lambdaP)*GSF_jP(jd)*d_qP(jd)*d_qP(kd)*wq
						END DO
					END DO
				END DO


				!  < { q } , [ B ].n > !test witht his term added 
				!DO id = 1, nDim
				!	i1 = i + nDim*nNodesPerElement
				!	j1 = j + (id-1)*nNodesPerElement
				!	Amm(i1,j1) = Amm(i1,j1) + 0.5d0*SF_iM*SF_jM*nM(id)*wq
				!	Amp(i1,j1) = Amp(i1,j1) + 0.5d0*SF_iM*SF_jP*nP(id)*wq
				!	Apm(i1,j1) = Apm(i1,j1) + 0.5d0*SF_iP*SF_jM*nM(id)*wq
				!	App(i1,j1) = App(i1,j1) + 0.5d0*SF_iP*SF_jP*nP(id)*wq
				!END DO
				
				!  < [ q ] , { B }.n > !test with this term added 
				!DO id = 1, nDim
				!	i1 = i + nDim*nNodesPerElement
				!	j1 = j + (id-1)*nNodesPerElement
				!	Amm(i1,j1) = Amm(i1,j1) + 0.5d0*SF_iM*SF_jM*nM(id)*wq
				!	Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*SF_jP*nP(id)*wq
				!	Apm(i1,j1) = Apm(i1,j1) - 0.5d0*SF_iP*SF_jM*nM(id)*wq
				!	App(i1,j1) = App(i1,j1) + 0.5d0*SF_iP*SF_jP*nP(id)*wq
				!END DO

				!ajout d'un terme de stabilisation sur la temperature
				!a<[q + Dq . d],[T + DT.d]>
				!on decompose en 4 parties

				!a<[q],[T]>
				i1 = i + nDim*nNodesPerElement
				j1 = j + nDim*nNodesPerElement
				Amm(i1,j1) = Amm(i1,j1) + tauU/h*SF_iM*SF_jM*wq
				Amp(i1,j1) = Amp(i1,j1) - tauU/h*SF_iM*SF_jP*wq
				Apm(i1,j1) = Apm(i1,j1) - tauU/h*SF_iP*SF_jM*wq
				App(i1,j1) = App(i1,j1) + tauU/h*SF_iP*SF_jP*wq
				
				!a< {q} , {T} >
				i1 = i + nDim*nNodesPerElement
				j1 = j + nDim*nNodesPerElement
				Amm(i1,j1) = Amm(i1,j1) + 0.25d0*tauU/h*SF_iM*SF_jM*wq
				Amp(i1,j1) = Amp(i1,j1) + 0.25d0*tauU/h*SF_iM*SF_jP*wq
				Apm(i1,j1) = Apm(i1,j1) + 0.25d0*tauU/h*SF_iP*SF_jM*wq
				App(i1,j1) = App(i1,j1) + 0.25d0*tauU/h*SF_iP*SF_jP*wq
				
				!- a* < [q] ,[ LM B.d] >
				DO jd = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (jd-1)*nNodesPerElement

					int= tauU/h*(1.d0/lambdaM)*SF_iM*wq
					Amm(i1,j1) = Amm(i1,j1) + SF_jM*d_qM(jd)*int
					int= tauU/h*(1.d0/lambdaP)*SF_iM*wq
					Amp(i1,j1) = Amp(i1,j1) + SF_jP*d_qP(jd)*int

					int= tauU/h*(1.d0/lambdaM)*SF_iP*wq
					Apm(i1,j1) = Apm(i1,j1) - SF_jM*d_qM(jd)*int

					int= tauU/h*(1.d0/lambdaP)*SF_iP*wq
					App(i1,j1) = App(i1,j1) - SF_jP*d_qP(jd)*int
				END DO
				
				!- a* < {q} , {LM B.d} >
				DO jd = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (jd-1)*nNodesPerElement

					int= -0.25d0*tauU/h*(1.d0/lambdaM)*SF_iM*wq
					Amm(i1,j1) = Amm(i1,j1) - SF_jM*d_qM(jd)*int
					int= -0.25d0*tauU/h*(1.d0/lambdaP)*SF_iM*wq
					Amp(i1,j1) = Amp(i1,j1) + SF_jP*d_qP(jd)*int

					int= -0.25d0*tauU/h*(1.d0/lambdaM)*SF_iP*wq
					Apm(i1,j1) = Apm(i1,j1) - SF_jM*d_qM(jd)*int

					int= -0.25d0*tauU/h*(1.d0/lambdaP)*SF_iP*wq
					App(i1,j1) = App(i1,j1) + SF_jP*d_qP(jd)*int
				END DO


				!- a* < [ LM w.d] ,[T] >
				DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + nDim*nNodesPerElement

					int= tauU/h*(1.d0/lambdaM)*SF_jM*wq
					Amm(i1,j1) = Amm(i1,j1) + SF_iM*d_qM(id)*int

					int= tauU/h*(1.d0/lambdaM)*SF_jP*wq
					Amp(i1,j1) = Amp(i1,j1) - SF_iM*d_qM(id)*int

					int= tauU/h*(1.d0/lambdaP)*SF_jM*wq
					Apm(i1,j1) = Apm(i1,j1) + SF_iP*d_qP(id)*int

					int= tauU/h*(1.d0/lambdaP)*SF_jP*wq
					App(i1,j1) = App(i1,j1) - SF_iP*d_qP(id)*int
				END DO
				
				!- a* < { LM w.d} , {T} >
				DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + nDim*nNodesPerElement

					int= -0.25d0*tauU/h*(1.d0/lambdaM)*SF_jM*wq
					Amm(i1,j1) = Amm(i1,j1) - SF_iM*d_qM(id)*int

					int= -0.25d0*tauU/h*(1.d0/lambdaM)*SF_jP*wq
					Amp(i1,j1) = Amp(i1,j1) - SF_iM*d_qM(id)*int

					int= -0.25d0*tauU/h*(1.d0/lambdaP)*SF_jM*wq
					Apm(i1,j1) = Apm(i1,j1) + SF_iP*d_qP(id)*int

					int= -0.25d0*tauU/h*(1.d0/lambdaP)*SF_jP*wq
					App(i1,j1) = App(i1,j1) + SF_iP*d_qP(id)*int
				END DO
				

				! tauU/h* < [LM w.d], [LM  B.D ] >
				DO id = 1,nDim
					DO jd = 1,nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement

						int =tauU/h*(1.d0/(lambdaM**2))*wq
						Amm(i1,j1) = Amm(i1,j1) + SF_iM*d_qM(id)*SF_jM*d_qM(jd)*int

						int =tauU/h*(1.d0/(lambdaM*lambdaP))*wq
						Amp(i1,j1) = Amp(i1,j1) + SF_iM*d_qM(id)*SF_jP*d_qP(jd)*int

						int =tauU/h*(1.d0/(lambdaP*lambdaM))*wq
						Apm(i1,j1) = Apm(i1,j1) + SF_iP*d_qP(id)*SF_jM*d_qM(jd)*int

						int =tauU/h*(1.d0/(lambdaP**2))*wq
						App(i1,j1) = App(i1,j1) + SF_iP*d_qP(id)*SF_jP*d_qP(jd)*int
					END DO
				END DO

				! tauU/h* < {LM w.d}, {LM  B.D} >
				DO id = 1,nDim
					DO jd = 1,nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement

						int = 0.25d0*tauU/h*(1.d0/(lambdaM**2))*wq
						Amm(i1,j1) = Amm(i1,j1) + SF_iM*d_qM(id)*SF_jM*d_qM(jd)*int

						int = 0.25d0*tauU/h*(1.d0/(lambdaM*lambdaP))*wq
						Amp(i1,j1) = Amp(i1,j1) - SF_iM*d_qM(id)*SF_jP*d_qP(jd)*int

						int = 0.25d0*tauU/h*(1.d0/(lambdaP*lambdaM))*wq
						Apm(i1,j1) = Apm(i1,j1) - SF_iP*d_qP(id)*SF_jM*d_qM(jd)*int

						int = 0.25d0*tauU/h*(1.d0/(lambdaP**2))*wq
						App(i1,j1) = App(i1,j1) + SF_iP*d_qP(id)*SF_jP*d_qP(jd)*int
					END DO
				END DO

				!- 0.5*a* < [q] , d^t[ LM * DB ]d >
				DO jd=1, nDim
					DO kd=1, nDim
						i1 = i + nDim*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement

						Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_iM*(1.d0/lambdaM)* &
						GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
						Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_iM*(1.d0/lambdaP)* &
						GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
						Apm(i1,j1) = Apm(i1,j1) + tauU/h*0.5d0*SF_iP*(1.d0/lambdaM)* &
						GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
						App(i1,j1) = App(i1,j1) - tauU/h*0.5d0*SF_iP*(1.d0/lambdaP)* &
						GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
					END DO
				END DO

				!- 0.5*a* < {q} , d^t{ LM * DB }d >
				DO jd=1, nDim
					DO kd=1, nDim
						i1 = i + nDim*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement

						Amm(i1,j1) = Amm(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_iM*(1.d0/lambdaM)* &
						GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
						Amp(i1,j1) = Amp(i1,j1) -  0.25d0*tauU/h*0.5d0*SF_iM*(1.d0/lambdaP)* &
						GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
						Apm(i1,j1) = Apm(i1,j1) -  0.25d0*tauU/h*0.5d0*SF_iP*(1.d0/lambdaM)* &
						GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
						App(i1,j1) = App(i1,j1) -  0.25d0*tauU/h*0.5d0*SF_iP*(1.d0/lambdaP)* &
						GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
					END DO
				END DO
				
				! 0.5*a* < [Lm w.d] , [d^t (LM * DB) d ] >
				DO id=1, nDim
					DO jd=1, nDim
						DO kd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM**2))* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM*lambdaP))* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
							Apm(i1,j1) = Apm(i1,j1) - tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP*lambdaM))* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							App(i1,j1) = App(i1,j1) + tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP**2))* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
						END DO
					END DO
				END DO
				
				! 0.5*a* < {Lm w.d} , {d^t (LM * DB) d} >
				DO id=1, nDim
					DO jd=1, nDim
						DO kd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM**2))* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							Amp(i1,j1) = Amp(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM*lambdaP))* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
							Apm(i1,j1) = Apm(i1,j1) + 0.25d0*tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP*lambdaM))* &
							GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
							App(i1,j1) = App(i1,j1) + 0.25d0*tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP**2))* &
							GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
						END DO
					END DO
				END DO

				!- 0.5*a* < [d^t (LM * Dw) d ],[T] >
				DO id = 1, nDim
					DO kd = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement

						Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_jM*(1.d0/lambdaM)* &
						GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

						Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_jP*(1.d0/lambdaM)* &
						GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

						Apm(i1,j1) = Apm(i1,j1) + tauU/h*0.5d0*SF_jM*(1.d0/lambdaP)* &
						GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq

						App(i1,j1) = App(i1,j1) - tauU/h*0.5d0*SF_jP*(1.d0/lambdaP)* &
						GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
					END DO
				END DO
				
				!- 0.5*a* < {d^t (LM * Dw) d } , {T} >
				DO id = 1, nDim
					DO kd = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement

						Amm(i1,j1) = Amm(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_jM*(1.d0/lambdaM)* &
						GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

						Amp(i1,j1) = Amp(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_jP*(1.d0/lambdaM)* &
						GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

						Apm(i1,j1) = Apm(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_jM*(1.d0/lambdaP)* &
						GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq

						App(i1,j1) = App(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_jP*(1.d0/lambdaP)* &
						GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
					END DO
				END DO

				! 0.5*a* <[d^t (LM * Dw ) d],[Lm B.d] >
				DO id=1, nDim
					DO kd=1,nDim
						DO jd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM**2))* &
							GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
							Amp(i1,j1) = Amp(i1,j1) - tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP*lambdaM))* &
							GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
							Apm(i1,j1) = Apm(i1,j1) + tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM*lambdaP))* &
							GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
							App(i1,j1) = App(i1,j1) + tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP**2))* &
							GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
						END DO
					END DO
				END DO
				
				! 0.5*a* < {d^t (LM * Dw ) d} , {Lm B.d} >
				DO id=1, nDim
					DO kd=1,nDim
						DO jd=1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement

							Amm(i1,j1) = Amm(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM**2))* &
							GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
							Amp(i1,j1) = Amp(i1,j1) + 0.25d0*tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP*lambdaM))* &
							GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
							Apm(i1,j1) = Apm(i1,j1) - 0.25d0*tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM*lambdaP))* &
							GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
							App(i1,j1) = App(i1,j1) + 0.25d0*tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP**2))* &
							GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
						END DO
					END DO
				END DO

				! 0.25*a* <[d^t (LM * Dw) d], [d^t (Lm DB) d] >
				DO id=1, nDim
					DO kd=1, nDim
						DO jd=1, nDim
							DO ld=1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement

								Amm(i1,j1) = Amm(i1,j1) + tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
								(1.d0/(lambdaM**2))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
								Amp(i1,j1) = Amp(i1,j1) - tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
								(1.d0/(lambdaM*lambdaP))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
								Apm(i1,j1) = Apm(i1,j1) - tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
								(1.d0/(lambdaP*lambdaM))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
								App(i1,j1) = App(i1,j1) + tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
							   (1.d0/(lambdaP**2))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
							END DO
						END DO
					END DO
				END DO
				
				! 0.25*a* < {d^t (LM * Dw) d} , {d^t (Lm DB) d} >
				DO id=1, nDim
					DO kd=1, nDim
						DO jd=1, nDim
							DO ld=1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement

								Amm(i1,j1) = Amm(i1,j1) + 0.25d0*tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
								(1.d0/(lambdaM**2))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
								Amp(i1,j1) = Amp(i1,j1) + 0.25d0*tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
								(1.d0/(lambdaM*lambdaP))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
								Apm(i1,j1) = Apm(i1,j1) + 0.25d0*tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
								(1.d0/(lambdaP*lambdaM))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
								App(i1,j1) = App(i1,j1) + 0.25d0*tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
							   (1.d0/(lambdaP**2))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
							END DO
						END DO
					END DO
				END DO
            
			END DO

			! - < {w}.n , jt >
			! - < (w+.n+ - w-.n-)*0.5 , jt >
			DO id = 1, nDim
				i1 = i + (id-1)*nNodesPerElement
				F1(i1) = F1(i1) + 0.5d0*nM(id)*SF_iM*jt*wq
				F2(i1) = F2(i1) - 0.5d0*nP(id)*SF_iP*jt*wq
			END DO

			!- < [ w ].~n , Tm > !new term !!!!
			! - < w+n+ + w-n-, Tm >
			DO id = 1, nDim
				i1 = i + (id-1)*nNodesPerElement
				F1(i1) = F1(i1) - nM(id)*SF_iM*Tm*wq
				F2(i1) = F2(i1) - nP(id)*SF_iP*Tm*wq
			END DO
				
			!ajout du terme sur le second membre issu de la stabilisation
			!a<[q+Dq.d],jT>

			!a<[q],jT>
			i1 = i + nDim*nNodesPerElement
			F1(i1) = F1(i1) - tauU/h*SF_iM*jt*wq
			F2(i1) = F2(i1) + tauU/h*SF_iP*jt*wq
			
			!a<{q},Tm>
			i1 = i + nDim*nNodesPerElement
			F1(i1) = F1(i1) + 0.5d0*tauU/h*SF_iM*Tm*wq
			F2(i1) = F2(i1) + 0.5d0*tauU/h*SF_iP*Tm*wq

			!terme pour le second membre
            !- a* < [ LM w.d] ,jT >
            DO id=1, nDim
               i1 = i + (id-1)*nNodesPerElement

               int= tauU/h*(1.d0/lambdaM)*jT*wq
               F1(i1) = F1(i1) - SF_iM*d_qM(id)*int

               int= tauU/h*(1.d0/lambdaP)*jT*wq
               F2(i1) = F2(i1) - SF_iP*d_qP(id)*int
            END DO
            
            !terme pour le second membre
            !- a* < { LM w.d} ,jT >
            DO id=1, nDim
               i1 = i + (id-1)*nNodesPerElement

               int= -0.5d0*tauU/h*(1.d0/lambdaM)*Tm*wq
               F1(i1) = F1(i1) - SF_iM*d_qM(id)*int

               int= -0.5d0*tauU/h*(1.d0/lambdaP)*Tm*wq
               F2(i1) = F2(i1) + SF_iP*d_qP(id)*int
            END DO

            !- 0.5*a* < d^t[ LM * Dw d],jT >
            DO id=1, nDim
				DO kd=1,nDim
					i1 = i + (id-1)*nNodesPerElement

					res= GSF_iM(kd)*d_qM(kd)*d_qM(id)
					F1(i1) = F1(i1) + tauU/h*0.5d0*(1.d0/lambdaM)*res*jT*wq
					res= GSF_iP(kd)*d_qP(kd)*d_qP(id)
					F2(i1) = F2(i1) - tauU/h*0.5d0*(1.d0/lambdaP)*res*jT*wq
				END DO
            END DO
            
            !- 0.5*a* < d^t{LM * Dw d}, Tm >
            DO id=1, nDim
				DO kd=1,nDim
					i1 = i + (id-1)*nNodesPerElement

					res= GSF_iM(kd)*d_qM(kd)*d_qM(id)
					F1(i1) = F1(i1) - 0.5d0*tauU/h*0.5d0*(1.d0/lambdaM)*res*Tm*wq
					res= GSF_iP(kd)*d_qP(kd)*d_qP(id)
					F2(i1) = F2(i1) - 0.5d0*tauU/h*0.5d0*(1.d0/lambdaP)*res*Tm*wq
				END DO
            END DO
            
		END DO

		IF ( Data%PressureEnrichment ) THEN
			DO i = 1, nNodesPerElement
				CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL GSFEval(GSF_iM, xq_M, i,eleM)
				CALL SFEval(eleP%eleType, SF_iP, xq_P, i) ; CALL GSFEval(GSF_iP, xq_P, i,eleP)
				
				DO ii = 1, nEdgesPerElement
					j = permut(ii,1) ; k = permut(ii,2)
					CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL SFEval(eleP%eleType, SF_kP, xq_P, k)
					CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL SFEval(eleM%eleType, SF_kM, xq_M, k)
					
					CALL GSFEval(GSF_jM, xq_M, j,eleM) ; CALL GSFEval(GSF_kM, xq_M, k,eleM) !!new term 
					CALL GSFEval(GSF_jP, xq_P, j,eleP) ; CALL GSFEval(GSF_kP, xq_P, k,eleP) !!new term		
			
					! tauU/h < [ q ] , [ T* ] >
					DO id = 1, nDim
					   i1 = i + nDim*nNodesPerElement
					   j1 = j + (id-1)*nNodesPerElement
					   k1 = k + (id-1)*nNodesPerElement
					   int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*wq
					   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
					   Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

					   int = -0.5d0*tauU/h*SF_iM*SF_jP*SF_kP*wq
					   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					   Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int

					   int = -0.5d0*tauU/h*SF_iP*SF_jM*SF_kM*wq
					   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
					   Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

					   int = 0.5d0*tauU/h*SF_iP*SF_jP*SF_kP*wq
					   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					   App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					END DO

					! tauU/h < {q} , {T*} >
					DO id = 1, nDim
					   i1 = i + nDim*nNodesPerElement
					   j1 = j + (id-1)*nNodesPerElement
					   k1 = k + (id-1)*nNodesPerElement
					   int = 0.25d0*0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*wq
					   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
					   Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

					   int = 0.25d0*0.5d0*tauU/h*SF_iM*SF_jP*SF_kP*wq
					   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					   Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int

					   int = 0.25d0*0.5d0*tauU/h*SF_iP*SF_jM*SF_kM*wq
					   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
					   Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

					   int = 0.25d0*0.5d0*tauU/h*SF_iP*SF_jP*SF_kP*wq
					   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					   App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
					END DO
					
					! -tauU/h* < LM [w.d],[ T* ] >
					DO id=1,nDim
						DO jd=1,nDim
						
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement
							k1 = k + (jd-1)*nNodesPerElement
							int =0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jM*SF_kM*wq
							Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							int =-0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jP*SF_kP*wq
							Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int


							int =0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jM*SF_kM*wq
							Apm(i1,j1) = Apm(i1,j1) +  (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							Apm(i1,k1) = Apm(i1,k1) -  (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							int =-0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jP*SF_kP*wq
							App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							  
						END DO
					END DO
					
					! -tauU/h* < LM {w.d} ,{ T*} >
					DO id=1,nDim
						DO jd=1,nDim
						
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement
							k1 = k + (jd-1)*nNodesPerElement
							int = 0.25d0*0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jM*SF_kM*wq
							Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							int = 0.25d0*0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jP*SF_kP*wq
							Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int


							int = -0.25d0*0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jM*SF_kM*wq
							Apm(i1,j1) = Apm(i1,j1) +  (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							Apm(i1,k1) = Apm(i1,k1) -  (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							int = -0.25d0*0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jP*SF_kP*wq
							App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							  
						END DO
					END DO

					! -tauU/h*0.5 < LM [d^t Dw.d],[ T* ] >
					DO id = 1, nDim
						DO kd = 1, nDim
							DO jd = 1, nDim
							
							   i1 = i + (id-1)*nNodesPerElement
							   j1 = j + (jd-1)*nNodesPerElement
							   k1 = k + (jd-1)*nNodesPerElement

							   res= (1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
							   int= -tauU/h*0.25d0*res*SF_jM*SF_kM*wq
							   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							   Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							   res= (1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
							   int= tauU/h*0.25d0*res*SF_jP*SF_kP*wq
							   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int

							   res= (1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
							   int= tauU/h*0.25d0*res*SF_jM*SF_kM*wq
							   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							   Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							   res= (1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
							   int= -tauU/h*0.25d0*res*SF_jP*SF_kP*wq
							   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   
							END DO
						END DO
					END DO
					
					! -tauU/h*0.5 < LM {d^t Dw.d} , {T*} >
					DO id = 1, nDim
						DO kd = 1, nDim
							DO jd = 1, nDim
							
							   i1 = i + (id-1)*nNodesPerElement
							   j1 = j + (jd-1)*nNodesPerElement
							   k1 = k + (jd-1)*nNodesPerElement

							   res= 0.25d0*(1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
							   int= -tauU/h*0.25d0*res*SF_jM*SF_kM*wq
							   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							   Amm(i1,k1) = Amm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							   res= 0.25d0*(1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
							   int= -tauU/h*0.25d0*res*SF_jP*SF_kP*wq
							   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   Amp(i1,k1) = Amp(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int

							   res= 0.25d0*(1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
							   int= -tauU/h*0.25d0*res*SF_jM*SF_kM*wq
							   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
							   Apm(i1,k1) = Apm(i1,k1) - (1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

							   res= 0.25d0*(1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
							   int= -tauU/h*0.25d0*res*SF_jP*SF_kP*wq
							   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   App(i1,k1) = App(i1,k1) - (1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
							   
							END DO
						END DO
					END DO

				END DO
			END DO
		END IF

		IF ( Data%TestFunctionEnrichment ) THEN
			DO ii = 1, nEdgesPerElement
				ii1 = permut(ii,1) ; ii2 = permut(ii,2)
				CALL SFEval(eleP%eleType, SF_i1P, xq_P, ii1) ; CALL SFEval(eleP%eleType, SF_i2P, xq_P, ii2)
				CALL SFEval(eleM%eleType, SF_i1M, xq_M, ii1) ; CALL SFEval(eleM%eleType, SF_i2M, xq_M, ii2)
				DO j = 1, nNodesPerElement
					CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL GSFEval(GSF_jM, xq_M, j,eleM)
					CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL GSFEval(GSF_jP, xq_P, j,eleP)

					! tauU/h < [ q* ] , [ T ] >
					DO id = 1, nDim
					   i1 = ii1 + (id-1)*nNodesPerElement
					   i2 = ii2 + (id-1)*nNodesPerElement
					   j1 = j   + nDim*nNodesPerElement

					   int = 0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*wq
					   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					   Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					   int = -0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*wq
					   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					   Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					   int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*wq
					   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					   Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

					   int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*wq
					   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					   App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					END DO
					
					! tauU/h < {q*} , {T} >
					DO id = 1, nDim
					   i1 = ii1 + (id-1)*nNodesPerElement
					   i2 = ii2 + (id-1)*nNodesPerElement
					   j1 = j   + nDim*nNodesPerElement

					   int = 0.25d0*0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*wq
					   Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					   Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					   int = 0.25d0*0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*wq
					   Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					   Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					   int = 0.25d0*0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*wq
					   Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					   Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

					   int = 0.25d0*0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*wq
					   App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					   App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					END DO
					

					! -tauU/h < [ q* ] , [ LM B.d ] >
					DO id = 1, nDim
					  DO jd= 1, nDim
						 i1 = ii1 + (id-1)*nNodesPerElement
						 i2 = ii2 + (id-1)*nNodesPerElement
						 j1 = j   + (jd-1)*nNodesPerElement

						 int = 0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
						 Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						 Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

						 int = 0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
						 Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						 Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

						 int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
						 Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						 Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

						 int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
						 App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						 App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					  END DO
					END DO


					! -tauU/h < {q*} , {LM B.d} >
					DO id = 1, nDim
					  DO jd= 1, nDim
						 i1 = ii1 + (id-1)*nNodesPerElement
						 i2 = ii2 + (id-1)*nNodesPerElement
						 j1 = j   + (jd-1)*nNodesPerElement

						 int = 0.25d0*0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
						 Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						 Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

						 int = -0.25d0*0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
						 Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
						 Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

						 int = 0.25d0*0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
						 Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						 Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

						 int = -0.25d0*0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
						 App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
						 App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					  END DO
					END DO
					
					! -tauU/h*0.5 < [ q* ],LM [d^t DB.d] >
					DO id =1, nDim
						DO jd =1, nDim
							DO kd = 1,nDim
							
								i1 = ii1 + (id-1)*nNodesPerElement
								i2 = ii2 + (id-1)*nNodesPerElement
								j1 = j   + (jd-1)*nNodesPerElement

								res=(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
								int=-tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
								Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
								Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

								res=(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
								int= tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
								Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
								Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

								res=(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
								int= tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
								Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

								res=(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
								int= -tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
								App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								
							END DO
						END DO
					END DO
					
					! -tauU/h*0.5 < {q*} , LM {d^t DB.d} >
					DO id =1, nDim
						DO jd =1, nDim
							DO kd = 1,nDim
							
								i1 = ii1 + (id-1)*nNodesPerElement
								i2 = ii2 + (id-1)*nNodesPerElement
								j1 = j   + (jd-1)*nNodesPerElement

								res = 0.25d0*(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
								int=-tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
								Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
								Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

								res = 0.25d0*(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
								int= -tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
								Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
								Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

								res = 0.25d0*(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
								int= -tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
								Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

								res = 0.25d0*(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
								int= -tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
								App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
								
							END DO
						END DO
					END DO


				END DO

				! a < [q*] , jt >
				DO id= 1, nDim
					i1 = ii1 + (id-1)*nNodesPerElement
					i2 = ii2 + (id-1)*nNodesPerElement

					int=-0.5d0*tauU/h*SF_i1M*SF_i2M*jt*wq
					F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					int=0.5d0*tauU/h*SF_i1P*SF_i2P*jt*wq
					F2(i1) = F2(i1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					F2(i2) = F2(i2) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
				END DO
				
				! a < {q*} , jt >
				DO id= 1, nDim
					i1 = ii1 + (id-1)*nNodesPerElement
					i2 = ii2 + (id-1)*nNodesPerElement

					int=0.5*0.5d0*tauU/h*SF_i1M*SF_i2M*Tm*wq
					F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
					F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

					int=0.5*0.5d0*tauU/h*SF_i1P*SF_i2P*Tm*wq
					F2(i1) = F2(i1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
					F2(i2) = F2(i2) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
				END DO

			END DO
		END IF

		IF ( Data%TestFunctionEnrichment .AND. Data%PressureEnrichment ) THEN
			DO ii = 1, nEdgesPerElement
			
				ii1 = permut(ii,1) ; ii2 = permut(ii,2)
				CALL SFEval(eleP%eleType, SF_i1P, xq_P, ii1) ; CALL SFEval(eleP%eleType, SF_i2P, xq_P, ii2)
				CALL SFEval(eleM%eleType, SF_i1M, xq_M, ii1) ; CALL SFEval(eleM%eleType, SF_i2M, xq_M, ii2)
				DO jj = 1, nEdgesPerElement
				
					jj1 = permut(jj,1) ; jj2 = permut(jj,2)
					CALL SFEval(eleM%eleType, SF_j1M, xq_M, jj1) ; CALL SFEval(eleM%eleType, SF_j2M, xq_M, jj2)
					CALL SFEval(eleP%eleType, SF_j1P, xq_P, jj1) ; CALL SFEval(eleP%eleType, SF_j2P, xq_P, jj2)

					! tauU/h < [ q* ] , [ T* ] >
					DO id = 1, nDim
						DO jd = 1, nDim
							i1 = ii1 + (id-1)*nNodesPerElement
							i2 = ii2 + (id-1)*nNodesPerElement
							j1 = jj1 + (jd-1)*nNodesPerElement
							j2 = jj2 + (jd-1)*nNodesPerElement

							int = 0.25d0*(1.d0/(lambdaM**2))*tauU/h*SF_i1M*SF_i2M*SF_j1M*SF_j2M*wq
							Amm(i1,j1) = Amm(i1,j1)+(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i1,j2) = Amm(i1,j2)-(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i2,j1) = Amm(i2,j1)-(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i2,j2) = Amm(i2,j2)+(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

							int = -0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1M*SF_i2M*SF_j1P*SF_j2P*wq
							Amp(i1,j1) = Amp(i1,j1) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i1,j2) = Amp(i1,j2) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i2,j1) = Amp(i2,j1) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i2,j2) = Amp(i2,j2) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int

							int = -0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1P*SF_i2P*SF_j1M*SF_j2M*wq
							Apm(i1,j1) = Apm(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i1,j2) = Apm(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i2,j1) = Apm(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i2,j2) = Apm(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

							int = 0.25d0*(1.d0/(lambdaP**2))*tauU/h*SF_i1P*SF_i2P*SF_j1P*SF_j2P*wq
							App(i1,j1) = App(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i1,j2) = App(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i2,j1) = App(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i2,j2) = App(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
						END DO
					END DO
					
					
					! tauU/h < {q} , {T*} >
					DO id = 1, nDim
						DO jd = 1, nDim
							i1 = ii1 + (id-1)*nNodesPerElement
							i2 = ii2 + (id-1)*nNodesPerElement
							j1 = jj1 + (jd-1)*nNodesPerElement
							j2 = jj2 + (jd-1)*nNodesPerElement

							int = 0.25d0*0.25d0*(1.d0/(lambdaM**2))*tauU/h*SF_i1M*SF_i2M*SF_j1M*SF_j2M*wq
							Amm(i1,j1) = Amm(i1,j1)+(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i1,j2) = Amm(i1,j2)-(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i2,j1) = Amm(i2,j1)-(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Amm(i2,j2) = Amm(i2,j2)+(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

							int = 0.25d0*0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1M*SF_i2M*SF_j1P*SF_j2P*wq
							Amp(i1,j1) = Amp(i1,j1) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i1,j2) = Amp(i1,j2) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i2,j1) = Amp(i2,j1) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							Amp(i2,j2) = Amp(i2,j2) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int

							int = 0.25d0*0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1P*SF_i2P*SF_j1M*SF_j2M*wq
							Apm(i1,j1) = Apm(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i1,j2) = Apm(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i2,j1) = Apm(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
							Apm(i2,j2) = Apm(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

							int = 0.25d0*0.25d0*(1.d0/(lambdaP**2))*tauU/h*SF_i1P*SF_i2P*SF_j1P*SF_j2P*wq
							App(i1,j1) = App(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i1,j2) = App(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i2,j1) = App(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
							App(i2,j2) = App(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
						END DO
					END DO


				END DO
			END DO
		END IF

    END DO
    
	
    DEALLOCATE ( Np, Nm, Xq_M, Xq_P, XqM, XqP, X_M, X_P, Exa_p, Exa_m)
    DEALLOCATE ( XeP, XeM )
    DEALLOCATE (GSF_iM, GSF_jM, GSF_iP, GSF_jP, GSF_kM, GSF_kP)

  END SUBROUTINE ComputePMMatrices_CG_Embedded2
  !======================================================================!
  
   !===================================
	SUBROUTINE MeanTemperatureRes(Data, Mesh, Sol,MeanT)
  	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: MeanT

    !-------------------------------------------------- 
    INTEGER 			  				:: Ne, Nv, nDim, nNodesPerElement, i
    INTEGER 			 				:: Nu_i, Nu_j, Ns 
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               				:: nele ! numerotation of the element countaining the node on the true interface
    REAL*8, ALLOCATABLE, DIMENSION(:) 	:: n_interface
    !-------------------------------------------------- 


    ALLOCATE(X_cur(2),n_interface(2))
     n_interface = 0.d0 ; MeanT = 0.d0 
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 

    !analytical velocity 
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
		X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! Sol_prec est definis via Sol en post processing 
		!d = Icing_L1 - mesh%Vertex(Nu_i)%X(1) 
		
		MeanT = MeanT + (Sol%Pcg(Nu_i)+Sol%Pcg(Nu_j))*0.5d0 
		MeanT = MeanT - (Icing_lambda1*Sol%Bcg(Nu_i,1)*Sol%dist(1,Nu_i,1)+Icing_lambda2*Sol%Bcg(Nu_j,1)*Sol%dist(1,Nu_i,1))*0.5d0
		MeanT = meanT - (Icing_lambda1*Sol%Bcg(Nu_i,2)*Sol%dist(1,Nu_i,2)+Icing_lambda2*Sol%Bcg(Nu_j,2)*Sol%dist(1,Nu_i,2))*0.5d0
		
		!PRINT*,"look",MeanT
		
    END DO
    MeanT = MeanT/Ns 
    !PRINT*,"trueMEan",MeanT
	!STOP 
  END SUBROUTINE MeanTemperatureRes
  !===================================  

  END MODULE EdgeIntegral_CG_Icing
