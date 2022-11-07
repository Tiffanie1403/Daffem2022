MODULE EdgeIntegral_DG_Primal

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability

  IMPLICIT NONE

CONTAINS

  !=========================================================================================!
  SUBROUTINE ComputePMMatrices_Primal(Edg, edgM, edgP, eleM, eleP, Amm, Amp, Apm, App, Ulocm, Ulocp, Fm, Fp, Data)
  !=========================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)  :: Edg
    type(Edge), DIMENSION(:), INTENT(IN) :: edgM, edgP
    type(Element)           , INTENT(IN)  :: eleM, eleP
    REAL*8, DIMENSION(:,:)  , INTENT(INOUT) :: Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)    , INTENT(IN)    :: Ulocm, Ulocp
    REAL*8, DIMENSION(:)    , INTENT(OUT)   :: Fm, Fp
    type(DataStructure)     , INTENT(IN)    :: Data
    !---------------------------------------------------------
    type(Edge) :: edgIntM, edgIntP
    !--------------------------------
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Np, Nm, Xq_M, Xq_P
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: xM, yM, zM, xP, yP, zP
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: LM_q, LP_q
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: GSF_jM, GSF_jP, GSF_iM, GSF_iP
    !----------------------------------------------------------------------
    REAL*8 :: Nrm, wq, L, mu_q_p, mu_q_m, mu_q, hOrt
    REAL*8 :: SF_iM, SF_iP, SF_jM, SF_jP, a, SF_kM, SF_kP
    REAL*8 :: x_iq_M, y_iq_M, x_iq_P, y_iq_P
    REAL*8 :: x1, y1, x2, y2, int
    REAL*8 :: tauP, tauU, h, NrmP, NrmM
    !--------------------------------------------------
    INTEGER :: nDim, nNodesPerElement
    INTEGER :: iM, iP, i, jM, jP, ii
    INTEGER :: iq, i1, j1, j, id, k, k1, jd
    INTEGER :: i1m, i1p, j1m, j1p, nQuad
    !------------------------------------------
    
    nDim = data%nDim ; nNodesPerElement = eleM%nNodesPerElement

    ALLOCATE ( Np(nDim), Nm(nDim), Xq_M(nDim), Xq_P(nDim) )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( LM_q(nDim,nDim), LP_q(nDim,nDim) )
    ALLOCATE ( GSF_iM(nDim), GSF_iP(nDim), GSF_jM(nDim), GSF_jP(nDim) )

    CALL IdentifyIFace(eleM, edg, iM, nDim)
    CALL IdentifyIFace(eleP, edg, iP, nDim)
    
    edgIntM = edgM(iM) ; edgIntP = edgP(iP)
    
    Np = -eleP%N(iP,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + Np(id)**2
    END DO
    Nrm = sqrt(Nrm)
    Np = Np/Nrm
    Nm = -Np
    L = Nrm
    
    Amm = 0.d0 ; Amp = 0.d0 ; Apm = 0.d0 ; App = 0.d0
    Fp = 0.d0 ; Fm = 0.d0
    CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim)
    CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)

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
       a = alphaIP1/h*Nrm
       
       DO i = 1, nNodesPerElement
          CALL SFEval(EleM%EleType,SF_iM, xq_M, i)
          CALL SFEval(EleM%EleType,SF_iP, xq_P, i)
          CALL GSFEval(GSF_iM, xq_M, i, eleM) ; CALL GSFEval(GSF_iP, xq_P, i, eleP)
          DO j = 1, nNodesPerElement
             CALL SFEval(EleM%EleType,SF_jM, xq_M, j)
             CALL SFEval(EleM%EleType,SF_jP, xq_P, j)
             CALL GSFEval(GSF_jM, xq_M, j, eleM) ; CALL GSFEval(GSF_jP, xq_P, j, eleP)
             
             ! - < [ q ] , { L Gp } >
             DO id = 1, nDim
                DO jd = 1, nDim
                   Amm(i,j) = Amm(i,j) - 0.5d0*Nm(id)*SF_iM*LM_q(id,jd)*GSF_jM(jd)*wq
                   Fm(i)    = Fm(i)    + 0.5d0*Nm(id)*SF_iM*LM_q(id,jd)*GSF_jM(jd)*wq*Ulocm(j)

                   Amp(i,j) = Amp(i,j) - 0.5d0*Nm(id)*SF_iM*LP_q(id,jd)*GSF_jP(jd)*wq
                   Fm(i)    = Fm(i)    + 0.5d0*Nm(id)*SF_iM*LP_q(id,jd)*GSF_jP(jd)*wq*Ulocp(j)

                   Apm(i,j) = Apm(i,j) - 0.5d0*Np(id)*SF_iP*LM_q(id,jd)*GSF_jM(jd)*wq
                   Fp(i)    = Fp(i)    + 0.5d0*Np(id)*SF_iP*LM_q(id,jd)*GSF_jM(jd)*wq*Ulocm(j)

                   App(i,j) = App(i,j) - 0.5d0*Np(id)*SF_iP*LP_q(id,jd)*GSF_jP(jd)*wq
                   Fp(i)    = Fp(j)    + 0.5d0*Np(id)*SF_iP*LP_q(id,jd)*GSF_jP(jd)*wq*Ulocp(j)
                END DO
             END DO
             
             ! - < { L Gq } , [ p ] >
             DO id = 1, nDim
                DO jd = 1, nDim
                   Amm(i,j) = Amm(i,j) - 0.5d0*Nm(id)*SF_jM*LM_q(id,jd)*GSF_iM(jd)*wq
                   Fm(i)    =    Fm(i) + 0.5d0*Nm(id)*SF_jM*LM_q(id,jd)*GSF_iM(jd)*wq*Ulocm(j)

                   Amp(i,j) = Amp(i,j) - 0.5d0*Np(id)*SF_jP*LM_q(id,jd)*GSF_iM(jd)*wq
                   Fm(i)    = Fm(i)    + 0.5d0*Np(id)*SF_jP*LM_q(id,jd)*GSF_iM(jd)*wq*Ulocp(j)
                   
                   Apm(i,j) = Apm(i,j) - 0.5d0*Nm(id)*SF_jM*LP_q(id,jd)*GSF_iP(jd)*wq
                   Fp(i)    =    Fp(i) + 0.5d0*Nm(id)*SF_jM*LP_q(id,jd)*GSF_iP(jd)*wq*Ulocm(j)

                   App(i,j) = App(i,j) - 0.5d0*Np(id)*SF_jP*LP_q(id,jd)*GSF_iP(jd)*wq
                   Fp(i)    =    Fp(i) + 0.5d0*Np(id)*SF_jP*LP_q(id,jd)*GSF_iP(jd)*wq*Ulocp(j)
                END DO
             END DO
             
             ! a < [ q ] , [ u ] >
             Amm(i,j) = Amm(i,j) + a*SF_iM*SF_jM*DOT_PRODUCT(Nm,Nm)*wq
             Fm(i)    =    Fm(i) - a*SF_iM*SF_jM*DOT_PRODUCT(Nm,Nm)*wq*Ulocm(j)

             Amp(i,j) = Amp(i,j) + a*SF_iM*SF_jP*DOT_PRODUCT(Nm,Np)*wq
             Fm(i)    = Fm(i)    - a*SF_iM*SF_jP*DOT_PRODUCT(Nm,Np)*wq*Ulocp(j)
             
             Apm(i,j) = Apm(i,j) + a*SF_iP*SF_jM*DOT_PRODUCT(Np,Nm)*wq
             Fp(i)    = Fp(i)    - a*SF_iP*SF_jM*DOT_PRODUCT(Np,Nm)*wq*Ulocm(j)
             
             App(i,j) = App(i,j) + a*SF_iP*SF_jP*DOT_PRODUCT(Np,Np)*wq
             Fp(i)    =    Fp(i) - a*SF_iP*SF_jP*DOT_PRODUCT(Np,Np)*wq*Ulocp(j)
          END DO
       END DO
       
    END DO

    DEALLOCATE ( Np, Nm, Xq_M, Xq_P, XqM, XqP, X_M, X_P, LM_q, LP_q )
    DEALLOCATE ( GSF_iM, GSF_jM, GSF_iP, GSF_jP )
    
  END SUBROUTINE ComputePMMatrices_Primal
  !=======================================================================!

END MODULE EdgeIntegral_DG_Primal
