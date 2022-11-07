MODULE EdgeIntegral_EG

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability

  IMPLICIT NONE

CONTAINS

  !=====================================================================================!
  SUBROUTINE ComputePMMatrices_EG(edg, edgM, edgP, eleM, eleP, Amm, Amp, Apm, App, Data)
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
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Np, Nm, Xq_M, Xq_P
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: LM_q, LP_q
    !-----------------------------------------------
    REAL*8 :: Nrm, wq, L, mu_q, hOrt, NrmM, NrmP
    REAL*8 :: SF_iM, SF_iP, SF_jM, SF_jP, a, SF_kM, SF_kP
    REAL*8 :: tauP, tauU, h
    !--------------------------------------------------
    INTEGER :: iM, iP, i, jM, jP, ii
    INTEGER :: iq, i1, j1, j, id, k, k1, jd
    INTEGER :: i1m, i1p, j1m, j1p, kd, nDofPerElement
    INTEGER :: nDim, nNodesPerElement, nQuad
    !-------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nDofPerElement = nNodesPerElement + 1
    
    ALLOCATE ( Np(nDim), Nm(nDim), Xq_M(nDim), Xq_P(nDim) )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( LM_q(nDim,nDim), LP_q(nDim,nDim) )
    
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
       
       DO i = 1, nDofPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i)
          CALL SFEval(eleM%eleType, SF_iP, xq_P, i)

          DO j = 1, nDofPerElement
             CALL SFEval(eleM%eleType, SF_jM, xq_M, j)
             CALL SFEval(eleM%eleType, SF_jP, xq_P, j)

             ! < [ w ] , { p } >
             DO id = 1, nDim
                i1 = i + (id-1)*nDofPerElement
                j1 = j + nDim*nDofPerElement
                Amm(i1,j1) = Amm(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jM*wq
                Amp(i1,j1) = Amp(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jP*wq
                Apm(i1,j1) = Apm(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jM*wq
                App(i1,j1) = App(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jP*wq
             END DO

             ! - < { q } , [ B ] >
             DO id = 1, nDim
                i1 = i + nDim*nDofPerElement
                j1 = j + (id-1)*nDofPerElement
                Amm(i1,j1) = Amm(i1,j1) - 0.5d0*SF_iM*SF_jM*Nm(id)*wq
                Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*SF_jP*Np(id)*wq
                Apm(i1,j1) = Apm(i1,j1) - 0.5d0*SF_iP*SF_jM*Nm(id)*wq
                App(i1,j1) = App(i1,j1) - 0.5d0*SF_iP*SF_jP*Np(id)*wq
             END DO             
          END DO
       END DO

       DO i = 1, nDofPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i)
          CALL SFEval(eleM%eleType, SF_iP, xq_P, i)

          DO j = 1, nDofPerElement
             CALL SFEval(eleM%eleType, SF_jM, xq_M, j)
             CALL SFEval(eleM%eleType, SF_jP, xq_P, j)

             ! tauP/h < [ w ] , [ B ] >
             DO id = 1, nDim
                DO jd = 1, nDim
                   i1 = i + (id-1)*nDofPerElement
                   j1 = j + (jd-1)*nDofPerElement
                   Amm(i1,j1) = Amm(i1,j1) + tauP/h*SF_iM*SF_jM*Nm(id)*Nm(jd)*wq
                   Amp(i1,j1) = Amp(i1,j1) + tauP/h*SF_iM*SF_jP*Nm(id)*Np(jd)*wq
                   Apm(i1,j1) = Apm(i1,j1) + tauP/h*SF_iP*SF_jM*Np(id)*Nm(jd)*wq
                   App(i1,j1) = App(i1,j1) + tauP/h*SF_iP*SF_jP*Np(id)*Np(jd)*wq
                END DO
             END DO
             ! tauU/h < [ q ] , [ p ] >
             i1 = i + nDim*nDofPerElement
             j1 = j + nDim*nDofPerElement
             Amm(i1,j1) = Amm(i1,j1) + tauU/h*SF_iM*SF_jM*wq
             Amp(i1,j1) = Amp(i1,j1) - tauU/h*SF_iM*SF_jP*wq
             Apm(i1,j1) = Apm(i1,j1) - tauU/h*SF_iP*SF_jM*wq
             App(i1,j1) = App(i1,j1) + tauU/h*SF_iP*SF_jP*wq

          END DO
       END DO
    END DO


  END SUBROUTINE ComputePMMatrices_EG
  !======================================================================!
  
END MODULE EdgeIntegral_EG
