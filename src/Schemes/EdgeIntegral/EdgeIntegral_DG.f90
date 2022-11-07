MODULE EdgeIntegral_DG

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability

  IMPLICIT NONE

CONTAINS

  !=====================================================================================!
  SUBROUTINE ComputePMMatrices_DG(edg, edgM, edgP, eleM, eleP, Amm, Amp, Apm, App, Data)
  !=====================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)  :: Edg
    type(Edge), DIMENSION(:), INTENT(IN)  :: edgM, edgP
    type(element)           , INTENT(IN)  :: eleM, eleP
    REAL*8, DIMENSION(:,:)  , INTENT(OUT) :: Amm, Amp, Apm, App
    type(DataStructure)     , INTENT(IN)  :: Data
    !----------------------------------------------------------
    type(edge) :: edgIntM, edgIntP
    !--------------------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1_M, Lm1_P
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Np, Nm, Xq_M, Xq_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: LM_q, LP_q, XeM, XeP
    !-------------------------------------------------------------
    REAL*8 :: Nrm, wq, L, mu_q, hOrt, NrmM, NrmP
    REAL*8 :: SF_iM, SF_iP, SF_jM, SF_jP, a, SF_kM, SF_kP
    REAL*8 :: tauP, tauU, h, int, SF_i1m, SF_i2m, SF_i1p, SF_i2p
    REAL*8 :: SF_j1p, SF_j2p, SF_j1m, SF_j2m
    !------------------------------------------------------------
    INTEGER :: iM, iP, i, jM, jP, ii, i2, ii1, ii2
    INTEGER :: iq, i1, j1, j, id, k, k1, jd, jj, jj1, jj2
    INTEGER :: i1m, i1p, j1m, j1p, kd, j2, ld
    INTEGER :: nDim, nNodesPerElement, nQuad
    INTEGER :: nEdgesPerElement
    !--------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nEdgesPerElement = eleM%nEdgesPerElement
    ALLOCATE ( Np(nDim), Nm(nDim), Xq_M(nDim), Xq_P(nDim) )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( LM_q(nDim,nDim), LP_q(nDim,nDim) )
    ALLOCATE ( Lm1_M(nNodesPerElement, nDim, nDim), Lm1_P(nNodesPerElement, nDim, nDim) )
    ALLOCATE ( XeM(nNodesPerElement,nDim), XeP(nNodesPerElement,nDim) )

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1_M(i,:,:), eleM, Data%icas)
       CALL SetNodalInverseMatrixPermeability(Lm1_P(i,:,:), eleP, Data%icas)
    END DO
    XeM = eleM%Coor ; XeP = eleP%Coor
    
    CALL IdentifyIFace(eleM, edg, iM, nDim)
    CALL IdentifyIFace(eleP, edg, iP, nDim)
    
    edgIntM = edgM(iM) ; edgIntP = edgP(iP)
    
    Np = -eleP%N(iP,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + Np(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm
    Np = Np/Nrm
    Nm = -Np

    CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim)
    CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)
    
    Amm = 0.d0 ; Amp = 0.d0 ; Apm = 0.d0 ; App = 0.d0

    CALL GetNQuad(edg%eleType, nQuad)
    
    DO iq = 1, nQuad

       xq_M = xqM(iq,:)
       xq_P = xqP(iq,:)
       
       CALL WeightQuad(edg%eleType, wq, iq)
       wq = wq*L
       
       CALL SetNodalMatrixPermeability(LM_q, eleM, icas)
       CALL SetNodalMatrixPermeability(LP_q, eleP, icas)
       NrmM = 0.d0 ; NrmP = 0.d0
       DO id = 1, nDim
          NrmM = NrmM + eleM%L(id,id)**2
          NrmP = NrmP + eleP%L(id,id)**2
       END DO
       Nrm = MAX(sqrt(NrmM), sqrt(NrmP))
       hOrt = (Elem%A + Elep%A)/(2.d0*L)
       h = hOrt
       tauU = alphaIP1*Nrm
       tauP = alphaIP2*h*h/Nrm
       
       DO i = 1, nNodesPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i)
          CALL SFEval(eleM%eleType, SF_iP, xq_P, i)

          DO j = 1, nNodesPerElement
             CALL SFEval(eleM%eleType, SF_jM, xq_M, j)
             CALL SFEval(eleM%eleType, SF_jP, xq_P, j)

             ! < [ w ] , { p } >
             DO id = 1, nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + nDim*nNodesPerElement
                Amm(i1,j1) = Amm(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jM*wq
                Amp(i1,j1) = Amp(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jP*wq
                Apm(i1,j1) = Apm(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jM*wq
                App(i1,j1) = App(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jP*wq
             END DO

             ! - < { q } , [ B ] >
             DO id = 1, nDim
                i1 = i + nDim*nNodesPerElement
                j1 = j + (id-1)*nNodesPerElement
                Amm(i1,j1) = Amm(i1,j1) - 0.5d0*SF_iM*SF_jM*Nm(id)*wq
                Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*SF_jP*Np(id)*wq
                Apm(i1,j1) = Apm(i1,j1) - 0.5d0*SF_iP*SF_jM*Nm(id)*wq
                App(i1,j1) = App(i1,j1) - 0.5d0*SF_iP*SF_jP*Np(id)*wq                
             END DO
          END DO
       END DO

       DO i = 1, nNodesPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i)
          CALL SFEval(eleP%eleType, SF_iP, xq_P, i)

          DO j = 1, nNodesPerElement
             CALL SFEval(eleM%eleType, SF_jM, xq_M, j)
             CALL SFEval(eleP%eleType, SF_jP, xq_P, j)

             ! tauP/h < [ w ] , [ B ] >
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
             ! tauU/h < [ q ] , [ p ] >
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Amm(i1,j1) = Amm(i1,j1) + tauU/h*SF_iM*SF_jM*wq
             Amp(i1,j1) = Amp(i1,j1) - tauU/h*SF_iM*SF_jP*wq
             Apm(i1,j1) = Apm(i1,j1) - tauU/h*SF_iP*SF_jM*wq
             App(i1,j1) = App(i1,j1) + tauU/h*SF_iP*SF_jP*wq

          END DO
       
          IF ( data%pressureEnrichment ) THEN
             DO ii = 1, nEdgesPerElement
                j = permut(ii,1) ; k = permut(ii,2)
                CALL SFEval(eleP%eleType, SF_jp, xq_P, j)
                CALL SFEval(eleP%eleType, SF_kp, xq_P, k)
                CALL SFEval(eleM%eleType, SF_jm, xq_M, j)
                CALL SFEval(eleM%eleType, SF_km, xq_M, k)                
                DO id = 1, nDim
                   DO jd = 1, nDim
                      DO kd = 1, nDim
                         ! < [ w ] , { p } >
                         i1 = i + (id-1)*nNodesPerElement
                         j1 = j + (kd-1)*nNodesPerElement
                         k1 = k + (kd-1)*nNodesPerElement
                         int = 0.25d0*SF_iM*SF_jM*SF_kM*Nm(id)*wq
                         Amm(i1,j1) = Amm(i1,j1) + Lm1_M(j,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         Amm(i1,k1) = Amm(i1,k1) - Lm1_M(k,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         int = 0.25d0*SF_iM*SF_jP*SF_kP*Nm(id)*wq
                         Amp(i1,j1) = Amp(i1,j1) + Lm1_P(j,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                         Amp(i1,k1) = Amp(i1,k1) - Lm1_P(k,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                         int = 0.25d0*SF_iP*SF_jM*SF_kM*Np(id)*wq
                         Apm(i1,j1) = Apm(i1,j1) + Lm1_M(j,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         Apm(i1,k1) = Apm(i1,k1) - Lm1_M(k,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         int = 0.25d0*SF_iP*SF_jP*SF_kP*Np(id)*wq
                         App(i1,j1) = App(i1,j1) + Lm1_P(j,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                         App(i1,k1) = App(i1,k1) - Lm1_P(k,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                         
                         
                         ! tauU/h < [ q ] , [ p ] >
                         i1 = i + nDim*nNodesPerElement
                         int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*Nm(id)*Nm(id)*wq
                         Amm(i1,j1) = Amm(i1,j1) + Lm1_M(j,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         Amm(i1,j1) = Amm(i1,j1) - Lm1_M(k,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*Nm(id)*Np(id)*wq
                         Amp(i1,j1) = Amp(i1,j1) + Lm1_P(j,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                         Amp(i1,j1) = Amp(i1,j1) - Lm1_P(k,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                         int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*Np(id)*Nm(id)*wq
                         Apm(i1,j1) = Apm(i1,j1) + Lm1_M(j,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         Apm(i1,j1) = Apm(i1,j1) - Lm1_M(k,jd,kd)*(XeM(j,jd) - XeM(k,jd))*int
                         int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*Np(id)*Np(id)*wq
                         App(i1,j1) = App(i1,j1) + Lm1_P(j,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                         App(i1,j1) = App(i1,j1) - Lm1_P(k,jd,kd)*(XeP(j,jd) - XeP(k,jd))*int
                      END DO
                   END DO
                END DO
             END DO
          END IF
       END DO

       IF ( data%pressureEnrichment .AND. data%testFunctionEnrichment ) THEN
          DO ii = 1, nEdgesPerElement
             ii1 = permut(ii,1) ; ii2 = permut(ii,2)
             CALL SFEval(eleP%eleType, SF_i1p, xq_P, ii1)
             CALL SFEval(eleP%eleType, SF_i2p, xq_P, ii2)          
             CALL SFEval(eleM%eleType, SF_i1m, xq_M, ii1)
             CALL SFEval(eleM%eleType, SF_i2m, xq_M, ii2)
             DO j = 1, nNodesPerElement
                CALL SFEval(eleP%eleType, SF_jP, xq_P, j)
                CALL SFEval(eleM%eleType, SF_jM, xq_M, j)
                
                DO id = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   DO jd = 1, nDim
                      ! - < { q* } , [ B ] >
                      j1 = j + (jd-1)*nNodesPerElement
                      int = -0.25d0*SF_i1m*SF_i2m*SF_jm*Nm(jd)*wq
                      Amm(i1,j1) = Amm(i1,j1) + (XeM(ii1,id) - XeM(ii2,id))*int
                      Amm(i2,j1) = Amm(i2,j1) - (XeM(ii1,id) - XeM(ii2,id))*int
                      int = -0.25d0*SF_i1m*SF_i2m*SF_jp*Np(jd)*wq
                      Amp(i1,j1) = Amp(i1,j1) + (XeM(ii1,id) - XeM(ii2,id))*int
                      Amp(i2,j1) = Amp(i2,j1) - (XeM(ii1,id) - XeM(ii2,id))*int
                      int = -0.25d0*SF_i1p*SF_i2p*SF_jm*Nm(jd)*wq
                      Apm(i1,j1) = Apm(i1,j1) + (XeP(ii1,id) - XeP(ii2,id))*int
                      Apm(i2,j1) = Apm(i2,j1) - (XeP(ii1,id) - XeP(ii2,id))*int
                      int = 0.25d0*SF_i1p*SF_i2p*SF_jp*Nm(jd)*wq
                      App(i1,j1) = App(i1,j1) + (XeP(ii1,id) - XeP(ii2,id))*int
                      App(i2,j1) = App(i2,j1) - (XeP(ii1,id) - XeP(ii2,id))*int
                      DO kd = 1, nDim
                         ! tauU/h < [ q* ] , [ p ] >
                         j1 = j + nDim*nNodesPerElement
                         int = 0.5d0*tauU/h*SF_i1m*SF_i2m*SF_jm*Nm(jd)*Nm(kd)*wq
                         Amm(i1,j1) = Amm(i1,j1) + (XeM(ii1,id) - XeM(ii2,id))*int
                         Amm(i2,j1) = Amm(i2,j1) - (XeM(ii1,id) - XeM(ii2,id))*int
                         int = 0.5d0*tauU/h*SF_i1m*SF_i2m*SF_jp*Nm(jd)*Np(kd)*wq
                         Amp(i1,j1) = Amp(i1,j1) + (XeM(ii1,id) - XeM(ii2,id))*int
                         Amp(i2,j1) = Amp(i2,j1) - (XeM(ii1,id) - XeM(ii2,id))*int
                         int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*Np(jd)*Nm(kd)*wq
                         Apm(i1,j1) = Apm(i1,j1) + (XeP(ii1,id) - XeP(ii2,id))*int
                         Apm(i2,j1) = Apm(i2,j1) - (XeP(ii1,id) - XeP(ii2,id))*int
                         int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*Np(jd)*Np(kd)*wq
                         App(i1,j1) = App(i1,j1) + (XeP(ii1,id) - XeP(ii2,id))*int
                         App(i2,j1) = App(i2,j1) - (XeP(ii1,id) - XeP(ii2,id))*int                      
                      END DO
                   END DO
                END DO
             END DO
             
             ! tauU/h < [ q* ] , [ p* ] >
             DO jj = 1, nEdgesPerElement
                DO id = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   jj1 = permut(jj,1) ; jj2 = permut(jj,2)
                   CALL SFEval(eleM%eleType, SF_j1M, xq_M, jj1)
                   CALL SFEval(eleM%eleType, SF_j2M, xq_M, jj2)
                   CALL SFEval(eleP%eleType, SF_j1P, xq_P, jj1)
                   CALL SFEval(eleP%eleType, SF_j2P, xq_P, jj2)
                   DO jd = 1, nDim
                      j1 = jj1 + (jd - 1)*nNodesPerElement
                      j2 = jj2 + (jd - 1)*nNodesPerElement
                      DO kd = 1, nDim
                         DO ld = 1, nDim
                            int = 0.25d0*tauU/h*SF_i1m*SF_i2m*SF_j1m*SF_j2m*Nm(kd)*Nm(ld)*wq
                            Amm(i1,j1) = Amm(i1,j1) + (XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            Amm(i1,j2) = Amm(i1,j2) - (XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            Amm(i2,j1) = Amm(i2,j1) - (XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            Amm(i2,j2) = Amm(i2,j2) + (XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            int = 0.25d0*tauU/h*SF_i1m*SF_i2m*SF_j1p*SF_j2p*Nm(kd)*Np(ld)*wq
                            Amp(i1,j1) = Amp(i1,j1) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                            Amp(i1,j2) = Amp(i1,j2) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                            Amp(i2,j1) = Amp(i2,j1) - (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                            Amp(i2,j2) = Amp(i2,j2) + (XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                            int = 0.25d0*tauU/h*SF_i1p*SF_i2p*SF_j1m*SF_j2m*Np(kd)*Nm(ld)*wq
                            Apm(i1,j1) = Apm(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            Apm(i1,j2) = Apm(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            Apm(i2,j1) = Apm(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            Apm(i2,j2) = Apm(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                            int = 0.25d0*tauU/h*SF_i1p*SF_i2p*SF_j1p*SF_j2p*Np(kd)*Np(ld)*wq
                            App(i1,j1) = App(i1,j1) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                            App(i1,j2) = App(i1,j2) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                            App(i2,j1) = App(i2,j1) - (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                            App(i2,j2) = App(i2,j2) + (XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END IF
       
       END DO
    
    DEALLOCATE ( Np, Nm, Xq_M, Xq_P, XqM, XqP, X_M, X_P, LM_q, LP_q )
    DEALLOCATE ( XeP, XeM, Lm1_P, Lm1_M )
    
  END SUBROUTINE ComputePMMatrices_DG
  !======================================================================!

END MODULE EdgeIntegral_DG
