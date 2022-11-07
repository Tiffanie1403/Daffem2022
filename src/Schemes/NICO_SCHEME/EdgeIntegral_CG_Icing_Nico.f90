MODULE EdgeIntegral_CG_Icing_Nicolson

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE ExactSolutiont
  USE ExactFunctionJumpt
  USE ReadMesh_mod
  USE Distance
  USE ComputeNormalsAndTangents

  IMPLICIT NONE

CONTAINS

  !=====================================================================================!
  SUBROUTINE ComputePMMatrices_CG_Embedded_Nico(mesh, edg, edgM, edgP, eleM, eleP, Amm, Amp, Apm, App, F1, &
     F2, Data, iLS,FM_prec,FP_prec)
    !=====================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)    :: Edg
    type(Edge), DIMENSION(:), INTENT(IN)    :: edgM, edgP
    type(element)           , INTENT(IN)    :: eleM, eleP
    REAL*8, DIMENSION(:,:)  , INTENT(OUT)   :: Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F1,F2
    REAL*8, DIMENSION(:)    , INTENT(IN)    :: FM_prec,FP_prec
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: ILS
    type(MeshType)          , INTENT(IN)    :: Mesh
    !----------------------------------------------------------
    type(edge) :: edgIntM, edgIntP
    type(element) :: ELE_inter
    !----------------------------------------------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Np, Nm, Xq_M, Xq_P, Exa_p, Exa_m
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: XeM, XeP
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: XxqP, XxqM
    REAL*8, DIMENSION(data%nDim)   :: d_qP, n_qP, d_qM, n_qM
    REAL*8, DIMENSION(data%nDim-1) :: tntM, tntP
    REAL*8, DIMENSION(data%nDim,data%nDim-1) :: t_qP, t_qM
    REAL*8, DIMENSION(data%nDim+1) :: exa1_q, exa2_q
    REAL*8, DIMENSION(data%nDim)   :: x_qP, x_qM
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_iM, GSF_jM,GSF_iP, GSF_jP
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_kM, GSF_kp

    !-------------------------------------------------------------
    REAL*8  :: Nrm, wq, HORT, jt, pD_M, pD_P, nnt, pd_q_bis, jt_bis, hn_bis
    REAL*8  :: SF_iM, SF_iP, SF_jM, SF_jP, a, ttM, ttP,res,jn2_bis
    REAL*8  :: SF_i1M, SF_i1P, SF_i2P, SF_i2M, SF_j1M, SF_j2P, SF_j1P,SF_j2M
    REAL*8  :: tauP, tauU, h, LS_qP, LS_qM, lambdaM, lambdaP,res1
    REAL*8  :: Sigma_bP, Sigma_bM, jn2, L, int, SF_kM, SF_kP,res2
    !------------------------------------------------------------
    INTEGER :: iM, iP, i, jM, jP,iq, i1, j1, j, id, jd, ii
    INTEGER :: ii1,ii2, jj,jj1,jj2,k,i2,k1, kd ,ld,j2
    INTEGER :: nDim, nNodesPerElement, nQuad, nEdgesPerElement, nVar
    !--------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nVar = nDim + 1 ; nEdgesPerElement = eleM%nEdgesPerElement
    ALLOCATE ( Np(nDim), Nm(nDim), Xq_M(nDim), Xq_P(nDim))
    ALLOCATE (  XxqP(nDim), XxqM(nDim) ) ; ALLOCATE ( Exa_p(nVar), Exa_M(nVar)  )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( XeM(nNodesPerElement,nDim), XeP(nNodesPerElement,nDim) )
    ALLOCATE ( GSF_iM(nDim), GSF_jM(nDim), GSF_iP(nDim), GSF_jP(nDim) )
    ALLOCATE ( GSF_kP(nDim), GSF_kM(nDim) )

    XeM = eleM%Coor ; XeP = eleP%Coor

    CALL IdentifyIFace(eleM, edg, iM, nDim) ; CALL IdentifyIFace(eleP, edg, iP, nDim)

    res=0.d0 ; res1=0.d0 ; res2=0.d0
    edgIntM = edgM(iM) ; edgIntP = edgP(iP)
    lambdaM=eleM%lambda ; lambdaP= eleP%lambda
    Np = -eleP%N(iP,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + Np(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L=Nrm
    Np = Np/Nrm
    Nm = -Np
    hOrt = (eleM%A + eleP%A)/(2.d0*L)
    h = hOrt
    tauU = 0 !alphaIP1*Nrm
    tauP = 0 !alphaIP2*h*h/Nrm
    CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim)
    CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)

    Amm = 0.d0 ; Amp = 0.d0 ; Apm = 0.d0 ; App = 0.d0
    F1 = 0.d0 ; F2 = 0.d0 ; nnt = 0.d0
    CALL GetNQuad(edg%eleType, nQuad)
    DO iq = 1, nQuad !boucle sur le nombre de points de quadrature

       xq_M = xqM(iq,:)
       xq_P = xqP(iq,:)

       x_qP = 0.d0
       x_qM = 0.d0
       DO i = 1, nNodesPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i)
          CALL SFEval(eleP%eleType, SF_iP, xq_P, i)
          DO id = 1, nDim
             x_qP(id) = x_qP(id) + SF_iP*eleP%Coor(i,id)
             x_qM(id) = x_qM(id) + SF_iM*eleM%Coor(i,id)
          END DO
       END DO
       CALL computeNodalDistance(x_qP, d_qP, LS_qP, Data, iLS)
       CALL getNormalAndTangents(d_qP, n_qP, t_qP, nP, nnt, tntP, nDim)

       IF (nnt>1 .OR. nnt<0) THEN
         STOP "Problème sur le calcul des normales"
       END IF

       n_qM = -n_qP
       d_qM = -d_qP

       IF ( icas /=  0 ) THEN

          CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP ,Data%T, eleP, nDim)
          CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM ,Data%T, eleM, nDim)
          jn2 = DOT_PRODUCT(exa1_q(1:nDim) - exa2_q(1:nDim), n_qP)

          CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP ,Data%T-Data%Dt, eleP, nDim)
          CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM ,Data%T-Data%Dt, eleM, nDim)
          jn2_bis = DOT_PRODUCT(exa1_q(1:nDim) - exa2_q(1:nDim), n_qP)
       ELSE
          jn2 = Data%BCEmbVal(iLS)
       END IF

       IF ( icas /= 0 ) THEN
          CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T-Data%Dt, eleP, nDim)
          CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T-Data%Dt, eleM, nDim)
          pD_P = exa1_q(nDim+1)
          pD_M = exa2_q(nDim+1)
          jt_bis = pD_P - pD_M

          CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T, eleP, nDim)
          CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T, eleM, nDim)
          pD_P = exa1_q(nDim+1)
          pD_M = exa2_q(nDim+1)
          jt = pD_P - pD_M

       ELSE
          jt = Data%BCEmbVal(iLS)
       END IF
       !PRINT*,"jt & j2n",jT,jn2

       CALL WeightQuad(edg%eleType, wq, iq)
       wq = wq*L
       ! terme a l'interface
       DO i = 1, nNodesPerElement
          CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL GSFEval(GSF_iM, xq_M, i,eleM)
          CALL SFEval(eleP%eleType, SF_iP, xq_P, i) ; CALL GSFEval(GSF_iP, xq_P, i,eleP)

          DO j = 1, nNodesPerElement
             CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL GSFEval(GSF_jM, xq_M, j,eleM)
             CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL GSFEval(GSF_jP, xq_P, j,eleP)

            ! < [ w ].~n , { T^n+1 } >
            DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Amm(i1,j1) = Amm(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jM*wq
             Amp(i1,j1) = Amp(i1,j1) + 0.5d0*Nm(id)*SF_iM*SF_jP*wq
             Apm(i1,j1) = Apm(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jM*wq
             App(i1,j1) = App(i1,j1) + 0.5d0*Np(id)*SF_iP*SF_jP*wq
            END DO

            ! -< [ w ].~n , { T^n } >
            DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             F1(i1) = F1(i1) - 0.5d0*Nm(id)*SF_iM*SF_jM*FM_prec(j1)*wq
             F1(i1) = F1(i1) - 0.5d0*Nm(id)*SF_iM*SF_jP*FP_prec(j1)*wq
             F2(i1) = F2(i1) - 0.5d0*Np(id)*SF_iP*SF_jM*FM_prec(j1)*wq
             F2(i1) = F2(i1) - 0.5d0*Np(id)*SF_iP*SF_jP*FP_prec(j1)*wq
            END DO

        !< {w}.~n , [Lm (B^n+1).d ] >
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

        !-< {w}.~n , [Lm B^n.d ] >
        DO id = 1, nDim
          DO jd = 1, nDim
            i1 = i + (id-1)*nNodesPerElement
            j1 = j + (jd-1)*nNodesPerElement
            F1(i1) = F1(i1) + 0.5d0*nM(id)*SF_iM*(1.d0/lambdaM)*SF_jM*FM_prec(j1)*d_qM(jd)*wq
            F1(i1) = F1(i1) + 0.5d0*nM(id)*SF_iM*(1.d0/lambdaP)*SF_jP*FP_prec(j1)*d_qP(jd)*wq
            F2(i1) = F2(i1) - 0.5d0*nP(id)*SF_iP*(1.d0/lambdaM)*SF_jM*FM_prec(j1)*d_qM(jd)*wq
            F2(i1) = F2(i1) - 0.5d0*nP(id)*SF_iP*(1.d0/lambdaP)*SF_jP*FP_prec(j1)*d_qP(jd)*wq
          END DO
        END DO


        ! - < { q } , [ B^n+1 ].n(n.~n) >
        DO id = 1, nDim
          i1 = i + nDim*nNodesPerElement
          j1 = j + (id-1)*nNodesPerElement
          Amm(i1,j1) = Amm(i1,j1) - 0.5d0*SF_iM*SF_jM*n_qM(id)*nnt*wq
          Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*SF_jP*n_qP(id)*nnt*wq
          Apm(i1,j1) = Apm(i1,j1) - 0.5d0*SF_iP*SF_jM*n_qM(id)*nnt*wq
          App(i1,j1) = App(i1,j1) - 0.5d0*SF_iP*SF_jP*n_qP(id)*nnt*wq
        END DO

        ! < { q } , [ B^n ].n(n.~n) >
        DO id = 1, nDim
          i1 = i + nDim*nNodesPerElement
          j1 = j + (id-1)*nNodesPerElement
          F1(i1) = F1(i1) + 0.5d0*SF_iM*SF_jM*FM_prec(j1)*n_qM(id)*nnt*wq
          F1(i1) = F1(i1) + 0.5d0*SF_iM*SF_jP*FP_prec(j1)*n_qP(id)*nnt*wq
          F2(i1) = F2(i1) + 0.5d0*SF_iP*SF_jM*FM_prec(j1)*n_qM(id)*nnt*wq
          F2(i1) = F2(i1) + 0.5d0*SF_iP*SF_jP*FP_prec(j1)*n_qP(id)*nnt*wq
        END DO


          ! - < { q } , [ D(B^n+1)d ].n(n.~n) >
  !        DO id = 1, nDim
  !          DO jd=1, nDim
  !            i1 = i + nDim*nNodesPerElement
  !            j1 = j + (id-1)*nNodesPerElement
  !            Amm(i1,j1) = Amm(i1,j1) + 0.5d0*SF_iM*N_qm(id)*GSF_jM(jd)*d_qM(jd)*nnt*wq
  !            Amp(i1,j1) = Amp(i1,j1) - 0.5d0*SF_iM*N_qp(id)*GSF_jP(jd)*d_qP(jd)*nnt*wq
  !            Apm(i1,j1) = Apm(i1,j1) + 0.5d0*SF_iP*N_qm(id)*GSF_jM(jd)*d_qM(jd)*nnt*wq
  !            App(i1,j1) = App(i1,j1) - 0.5d0*SF_iP*N_qp(id)*GSF_jP(jd)*d_qP(jd)*nnt*wq
  !          END DO
  !        END DO

          !  < { q } , [ D(B^n)d ].n(n.~n) >
  !        DO id = 1, nDim
  !          DO jd=1, nDim
  !            i1 = i + nDim*nNodesPerElement
  !            j1 = j + (id-1)*nNodesPerElement
  !            F1(i1) = F1(i1) - 0.5d0*SF_iM*N_qm(id)*FM_prec(j1)*GSF_jM(jd)*d_qM(jd)*nnt*wq
  !            F1(i1) = F1(i1) + 0.5d0*SF_iM*N_qp(id)*FP_prec(j1)*GSF_jP(jd)*d_qP(jd)*nnt*wq
  !            F2(i1) = F2(i1) - 0.5d0*SF_iP*N_qm(id)*FM_prec(j1)*GSF_jM(jd)*d_qM(jd)*nnt*wq
  !            F2(i1) = F2(i1) + 0.5d0*SF_iP*N_qp(id)*FP_prec(j1)*GSF_jP(jd)*d_qP(jd)*nnt*wq
  !          END DO
  !        END DO


          !ajout d'un terme de stabilisation sur la temperature
          !a<[q + Dq . d],[T + DT.d]>
          !on decompose en 4 parties

            !a<[q],[T^n+1]>
             i1 = i + nDim*nNodesPerElement
             j1 = j + nDim*nNodesPerElement
             Amm(i1,j1) = Amm(i1,j1) + tauU/h*SF_iM*SF_jM*wq
             Amp(i1,j1) = Amp(i1,j1) - tauU/h*SF_iM*SF_jP*wq
             Apm(i1,j1) = Apm(i1,j1) - tauU/h*SF_iP*SF_jM*wq
             App(i1,j1) = App(i1,j1) + tauU/h*SF_iP*SF_jP*wq

             !-a<[q],[T^n]>
              i1 = i + nDim*nNodesPerElement
              j1 = j + nDim*nNodesPerElement
              F1(i1) = F1(i1) - tauU/h*SF_iM*FM_prec(j1)*SF_jM*wq
              F1(i1) = F1(i1) + tauU/h*SF_iM*FP_prec(j1)*SF_jP*wq
              F2(i1) = F2(i1) + tauU/h*SF_iP*FM_prec(j1)*SF_jM*wq
              F2(i1) = F2(i1) - tauU/h*SF_iP*FP_prec(j1)*SF_jP*wq


        !- a* < [q] ,[ LM (B^n+1).d] >
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

          !a* < [q] ,[ LM (B^n).d] >
            DO jd = 1, nDim
             i1 = i + nDim*nNodesPerElement
             j1 = j + (jd-1)*nNodesPerElement

              int= -tauU/h*(1.d0/lambdaM)*SF_iM*wq
              F1(i1) = F1(i1) + SF_jM*FM_prec(j1)*d_qM(jd)*int

              int= -tauU/h*(1.d0/lambdaP)*SF_iM*wq
              F1(i1) = F1(i1) + SF_jP*FP_prec(j1)*d_qP(jd)*int

              int= -tauU/h*(1.d0/lambdaM)*SF_iP*wq
              F2(i1) = F2(i1) - SF_jM*FM_prec(j1)*d_qM(jd)*int

              int= -tauU/h*(1.d0/lambdaP)*SF_iP*wq
              F2(i1) = F2(i1) - SF_jP*FP_prec(j1)*d_qP(jd)*int
            END DO


          !- a* < [ LM w.d] ,[T^n+1] >
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

          ! a* < [ LM w.d] ,[T^n] >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             j1 = j + nDim*nNodesPerElement

             int= -tauU/h*(1.d0/lambdaM)*SF_jM*wq
             F1(i1) = F1(i1) + FM_prec(j1)*SF_iM*d_qM(id)*int

             int= -tauU/h*(1.d0/lambdaM)*SF_jP*wq
             F1(i1) = F1(i1) - FP_prec(j1)*SF_iM*d_qM(id)*int

             int= -tauU/h*(1.d0/lambdaP)*SF_jM*wq
             F2(i1) = F2(i1) + FM_prec(j1)*SF_iP*d_qP(id)*int

             int= -tauU/h*(1.d0/lambdaP)*SF_jP*wq
             F2(i1) = F2(i1) - FP_prec(j1)*SF_iP*d_qP(id)*int
          END DO

          ! tauU/h* < [LM w.d], [LM  (B^n+1).D ] >
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

          ! -tauU/h* < [LM w.d], [LM  (B^n).D ] >
          DO id = 1,nDim
             DO jd = 1,nDim
                i1 = i + (id-1)*nNodesPerElement
                j1 = j + (jd-1)*nNodesPerElement

                int =-tauU/h*(1.d0/(lambdaM**2))*wq
                F1(i1) = F1(i1) + FM_prec(j1)*SF_iM*d_qM(id)*SF_jM*d_qM(jd)*int

                int =-tauU/h*(1.d0/(lambdaM*lambdaP))*wq
                F1(i1) = F1(i1) + FP_prec(j1)*SF_iM*d_qM(id)*SF_jP*d_qP(jd)*int

                int =-tauU/h*(1.d0/(lambdaP*lambdaM))*wq
                F2(i1) = F2(i1) + FM_prec(j1)*SF_iP*d_qP(id)*SF_jM*d_qM(jd)*int

                int =-tauU/h*(1.d0/(lambdaP**2))*wq
                F2(i1) = F2(i1) + FP_prec(j1)*SF_iP*d_qP(id)*SF_jP*d_qP(jd)*int
             END DO
          END DO


          !terme supplémentaire issu du ddl
          !on ajoute tout de meme au non enrichie pour comparaison
          ! 0.5* < { w }.~n , [ d^t ( LM * D(B^n+1) ) d ] >
    ! !     DO id = 1, nDim
        !    DO jd=1, nDim
        !      DO kd=1, nDim
        !        i1 = i + (id-1)*nNodesPerElement
!                j1 = j + (jd-1)*nNodesPerElement
!
!                Amm(i1,j1) = Amm(i1,j1) + 0.25d0*SF_iM*Nm(id)*(1.d0/lambdaM)* &
!                GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                Amp(i1,j1) = Amp(i1,j1) - 0.25d0*SF_iM*Nm(id)*(1.d0/lambdaP)* &
!                GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                Apm(i1,j1) = Apm(i1,j1) - 0.25d0*SF_iP*Np(id)*(1.d0/lambdaM)* &
!                GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                App(i1,j1) = App(i1,j1) + 0.25d0*SF_iP*Np(id)*(1.d0/lambdaP)* &
!                GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!              END DO
!            END DO
!          END DO

          ! -0.5* < { w }.~n , [ d^t ( LM * D(B^n)) d ] >
    !      DO id = 1, nDim
    !        DO jd=1, nDim
    !          DO kd=1, nDim
    !            i1 = i + (id-1)*nNodesPerElement
    !            j1 = j + (jd-1)*nNodesPerElement
!
!                F1(i1) = F1(i1) - FM_prec(j1)*0.25d0*SF_iM*Nm(id)*(1.d0/lambdaM)* &
!                GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                F1(i1) = F1(i1) + FP_prec(j1)*0.25d0*SF_iM*Nm(id)*(1.d0/lambdaP)* &
!                GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                F2(i1) = F2(i1) + FM_prec(j1)*0.25d0*SF_iP*Np(id)*(1.d0/lambdaM)* &
!                GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                F2(i1) = F2(i1) - FP_prec(j1)*0.25d0*SF_iP*Np(id)*(1.d0/lambdaP)* &
!                GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!              END DO
!            END DO
!          END DO


          !- 0.5*a* < [q] , d^t[ LM * D(B^n+1) ]d >
  !          DO jd=1, nDim
  !            DO kd=1, nDim
  !              i1 = i + nDim*nNodesPerElement
  !              j1 = j + (jd-1)*nNodesPerElement
!
!                Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_iM*(1.d0/lambdaM)* &
!                GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_iM*(1.d0/lambdaP)* &
!                GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                Apm(i1,j1) = Apm(i1,j1) + tauU/h*0.5d0*SF_iP*(1.d0/lambdaM)* &
!                GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                App(i1,j1) = App(i1,j1) - tauU/h*0.5d0*SF_iP*(1.d0/lambdaP)* &
!                GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!              END DO
!            END DO

            ! 0.5*a* < [q] , d^t[ LM * D(B^n) ]d >
        !      DO jd=1, nDim
        !        DO kd=1, nDim
        !          i1 = i + nDim*nNodesPerElement
        !          j1 = j + (jd-1)*nNodesPerElement
!
!                  F1(i1) = F1(i1) + FM_prec(j1)*tauU/h*0.5d0*SF_iM*(1.d0/lambdaM)* &
!                  GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                  F1(i1) = F1(i1) - FP_prec(j1)*tauU/h*0.5d0*SF_iM*(1.d0/lambdaP)* &
!                  GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                  F2(i1) = F2(i1) - FM_prec(j1)*tauU/h*0.5d0*SF_iP*(1.d0/lambdaM)* &
!                  GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                  F2(i1) = F2(i1) + FP_prec(j1)*tauU/h*0.5d0*SF_iP*(1.d0/lambdaP)* &
!                  GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                END DO
!              END DO

            ! 0.5*a* < [Lm w.d] , [d^t (LM * D(B^n+1)) d ] >
  !          DO id=1, nDim
  !            DO jd=1, nDim
  !              DO kd=1, nDim
  !                i1 = i + (id-1)*nNodesPerElement
  !                j1 = j + (jd-1)*nNodesPerElement
!
!                  Amm(i1,j1) = Amm(i1,j1) - tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM**2))* &
!                  GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                  Amp(i1,j1) = Amp(i1,j1) + tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM*lambdaP))* &
!                  GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                  Apm(i1,j1) = Apm(i1,j1) - tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP*lambdaM))* &
!                  GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                  App(i1,j1) = App(i1,j1) + tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP**2))* &
!                  GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                END DO
!              END DO
!            END DO

            ! -0.5*a* < [Lm w.d] , [d^t (LM * D(B^n)) d ] >
        !    DO id=1, nDim
        !      DO jd=1, nDim
        !        DO kd=1, nDim
        !          i1 = i + (id-1)*nNodesPerElement
        !          j1 = j + (jd-1)*nNodesPerElement
!
!                  F1(i1) = F1(i1) + FM_prec(j1)*tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM**2))* &
!                  GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                  F1(i1) = F1(i1) - FP_prec(j1)*tauU/h*0.5d0*SF_iM*d_qM(id)*(1.d0/(lambdaM*lambdaP))* &
!                  GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                  F2(i1) = F2(i1) + FM_prec(j1)*tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP*lambdaM))* &
!                  GSF_jM(kd)*d_qM(kd)*d_qM(jd)*wq
!                  F2(i1) = F2(i1) - FP_prec(j1)*tauU/h*0.5d0*SF_iP*d_qP(id)*(1.d0/(lambdaP**2))* &
!                  GSF_jP(kd)*d_qP(kd)*d_qP(jd)*wq
!                END DO
!              END DO
!            END DO

            !- 0.5*a* < [d^t (LM * Dw) d ],[T^n+1] >
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

           ! 0.5*a* < [d^t (LM * Dw) d ],[T^n] >
             DO id = 1, nDim
               DO kd = 1, nDim
                 i1 = i + (id-1)*nNodesPerElement
                 j1 = j + nDim*nNodesPerElement

                  F1(i1) = F1(i1) + FM_prec(j1)*tauU/h*0.5d0*SF_jM*(1.d0/lambdaM)* &
                 GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

                 F1(i1) = F1(i1) - FP_prec(j1)*tauU/h*0.5d0*SF_jP*(1.d0/lambdaM)* &
                 GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq

                 F2(i1) = F2(i1) - FM_prec(j1)*tauU/h*0.5d0*SF_jM*(1.d0/lambdaP)* &
                 GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq

                 F2(i1) = F2(i1) + FP_prec(j1)*tauU/h*0.5d0*SF_jP*(1.d0/lambdaP)* &
                 GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
            END DO
          END DO

          ! 0.5*a* <[d^t (LM * Dw ) d],[Lm (B^n+1).d] >
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

           ! -0.5*a* <[d^t (LM * Dw ) d],[Lm (B^n).d] >
           DO id=1, nDim
            DO kd=1,nDim
              DO jd=1, nDim
                 i1 = i + (id-1)*nNodesPerElement
                 j1 = j + (jd-1)*nNodesPerElement

                 F1(i1) = F1(i1) + FM_prec(j1)*tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM**2))* &
                 GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
                 F1(i1) = F1(i1) + FP_prec(j1)*tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP*lambdaM))* &
                 GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
                 F2(i1) = F2(i1) - FM_prec(j1)*tauU/h*0.5d0*SF_jM*d_qM(jd)*(1.d0/(lambdaM*lambdaP))* &
                 GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
                 F2(i1) = F2(i1) - FP_prec(j1)*tauU/h*0.5d0*SF_jP*d_qP(jd)*(1.d0/(lambdaP**2))* &
                 GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
              END DO
             END DO
            END DO

          ! 0.25*a* <[d^t (LM * Dw) d], [d^t (Lm D(B^n+1)) d] >
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

        ! END DO

         ! 0.25*a* <[d^t (LM * Dw) d], [d^t (Lm D(B^n) d] >
  !       DO id=1, nDim
  !         DO kd=1, nDim
  !           DO jd=1, nDim
  !              DO ld=1, nDim
  !                 i1 = i + (id-1)*nNodesPerElement
  !                 j1 = j + (jd-1)*nNodesPerElement
!
!!                   F1(i1) = F1(i1) -  FM_prec(j1)*tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
  !                 (1.d0/(lambdaM**2))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
  !                 F1(i1) = F1(i1) +  FP_prec(j1)*tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
  !                 (1.d0/(lambdaM*lambdaP))*GSF_iM(kd)*d_qM(kd)*d_qM(id)*wq
  !                 F2(i1) = F2(i1) +  FM_prec(j1)*tauU/h*0.25d0*GSF_jM(ld)*d_qM(ld)*d_qM(jd)* &
  !                 (1.d0/(lambdaP*lambdaM))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
  !                 F2(i1) = F2(i1) -  FP_prec(j1)*tauU/h*0.25d0*GSF_jP(ld)*d_qP(ld)*d_qP(jd)* &
  !                 (1.d0/(lambdaP**2))*GSF_iP(kd)*d_qP(kd)*d_qP(id)*wq
  !               END DO
  !             END DO
  !           END DO
!           END DO
!
        END DO

          ! - < {q} , j2(n.~n)>
    !      i1 = i + nDim*nNodesPerElement
    !      F1(i1) = F1(i1) - 0.5d0*SF_iM*jn2*nnt*wq
    !      F2(i1) = F2(i1) - 0.5d0*SF_iP*jn2*nnt*wq
!
!          ! - < {q} , j2_bis(n.~n)>
!          i1 = i + nDim*nNodesPerElement
!          F1(i1) = F1(i1) - 0.5d0*SF_iM*jn2_bis*nnt*wq
!          F2(i1) = F2(i1) - 0.5d0*SF_iP*jn2_bis*nnt*wq

          ! - < {w}.n , jt >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             F1(i1) = F1(i1) + 0.5d0*nM(id)*SF_iM*jt*wq
             F2(i1) = F2(i1) - 0.5d0*nP(id)*SF_iP*jt*wq
          END DO

          ! - < {w}.n , jt_bis >
          DO id = 1, nDim
             i1 = i + (id-1)*nNodesPerElement
             F1(i1) = F1(i1) + 0.5d0*nM(id)*SF_iM*jt_bis*wq
             F2(i1) = F2(i1) - 0.5d0*nP(id)*SF_iP*jt_bis*wq
          END DO

          !ajout du terme sur le second membre issu de la stabilisation
          !a<[q+Dq.d],jT>
          i1 = i + nDim*nNodesPerElement
          F1(i1) = F1(i1) - tauU/h*SF_iM*jt*wq
          F2(i1) = F2(i1) + tauU/h*SF_iP*jt*wq

          !a<[q+Dq.d],jT_bis>
          i1 = i + nDim*nNodesPerElement
          F1(i1) = F1(i1) - tauU/h*SF_iM*jt_bis*wq
          F2(i1) = F2(i1) + tauU/h*SF_iP*jt_bis*wq

          !terme pour le second membre
             !- a* < [ LM w.d] ,jT >
             DO id=1, nDim
                i1 = i + (id-1)*nNodesPerElement

                int= tauU/h*(1.d0/lambdaM)*jT*wq
                F1(i1) = F1(i1) - SF_iM*d_qM(id)*int

                int= tauU/h*(1.d0/lambdaP)*jT*wq
                F2(i1) = F2(i1) - SF_iP*d_qP(id)*int
             END DO

             !- a* < [ LM w.d] ,jT_bis >
             DO id=1, nDim
                i1 = i + (id-1)*nNodesPerElement

                int= tauU/h*(1.d0/lambdaM)*jT_bis*wq
                F1(i1) = F1(i1) - SF_iM*d_qM(id)*int

                int= tauU/h*(1.d0/lambdaP)*jT_bis*wq
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

             !- 0.5*a* < d^t[ LM * Dw d],jT_bis >
             DO id=1, nDim
               DO kd=1,nDim
                i1 = i + (id-1)*nNodesPerElement

                res= GSF_iM(kd)*d_qM(kd)*d_qM(id)
                F1(i1) = F1(i1) + tauU/h*0.5d0*(1.d0/lambdaM)*res*jT_bis*wq
                res= GSF_iP(kd)*d_qP(kd)*d_qP(id)
                F2(i1) = F2(i1) - tauU/h*0.5d0*(1.d0/lambdaP)*res*jT_bis*wq
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

                ! < [ w ].~n , { (T^n+1)* } >
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

                ! -< [ w ].~n , { T^n* } >
                DO id = 1, nDim
                   DO jd = 1, nDim
                      i1 = i + (id-1)*nNodesPerElement
                      j1 = j + (jd-1)*nNodesPerElement
                      k1 = k + (jd-1)*nNodesPerElement
                      int = -0.25d0*SF_iM*Nm(id)*SF_jM*SF_kM*wq
                      F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                      F1(i1) = F1(i1) - FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                      int = -0.25d0*SF_iM*Nm(id)*SF_jP*SF_kP*wq
                      F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                      F1(i1) = F1(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int

                      int = -0.25d0*SF_iP*Np(id)*SF_jM*SF_kM*wq
                      F2(i1) = F2(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                      F2(i1) = F2(i1) - FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                      int = -0.25d0*SF_iP*Np(id)*SF_jP*SF_kP*wq
                      F2(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                      F2(i1) = F2(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                   END DO
                END DO


                ! -tauU/h < [ q ] , [ T^n* ] >
                DO id = 1, nDim
                   i1 = i + nDim*nNodesPerElement
                   j1 = j + (id-1)*nNodesPerElement
                   k1 = k + (id-1)*nNodesPerElement
                   int = 0.5d0*tauU/h*SF_iM*SF_jM*SF_kM*wq
                   F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
                   F1(i1) = F1(i1) - FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

                   int = -0.5d0*tauU/h*SF_iM*SF_jP*SF_kP*wq
                   F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
                   F1(i1) = F1(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int

                   int = -0.5d0*tauU/h*SF_iP*SF_jM*SF_kM*wq
                   F2(i1) = F2(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int
                   F2(i1) = F2(i1) - FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,id) - XeM(k,id))*int

                   int = 0.5d0*tauU/h*SF_iP*SF_jP*SF_kP*wq
                   F2(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
                   F2(i1) = F2(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,id) - XeP(k,id))*int
                END DO

                ! tauU/h < [ q ] , [ T^n+1* ] >
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

                ! -tauU/h* < LM [w.d],[ T^n+1* ] >
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

                ! tauU/h* < LM [w.d],[ T^n* ] >
                DO id=1,nDim
                   DO jd=1,nDim
                      i1 = i + (id-1)*nNodesPerElement
                      j1 = j + (jd-1)*nNodesPerElement
                      k1 = k + (jd-1)*nNodesPerElement
                      int =-0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jM*SF_kM*wq
                      F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                      F1(i1) = F1(i1) - FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                      int =0.5d0*tauU/h*(1.d0/lambdaM)*SF_iM*d_qM(id)*SF_jP*SF_kP*wq
                      F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                      F1(i1) = F1(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int


                      int =-0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jM*SF_kM*wq
                      F2(i1) = F2(i1) +  FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                      F2(i1) = F2(i1) -  FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                      int =0.5d0*tauU/h*(1.d0/lambdaP)*SF_iP*d_qP(id)*SF_jP*SF_kP*wq
                      F2(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                      F2(i1) = F2(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                   END DO
                END DO

                ! -tauU/h*0.5 < LM [d^t Dw.d],[ T^n+1* ] >
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

                ! tauU/h*0.5 < LM [d^t Dw.d],[ T^n* ] >
                DO id = 1, nDim
                  DO kd = 1, nDim
                    DO jd = 1, nDim
                       i1 = i + (id-1)*nNodesPerElement
                       j1 = j + (jd-1)*nNodesPerElement
                       k1 = k + (jd-1)*nNodesPerElement

                       res= (1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
                       int= tauU/h*0.25d0*res*SF_jM*SF_kM*wq
                       F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                       F1(i1) = F1(i1) - FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                       res= (1.d0/lambdaM)*GSF_iM(kd)*d_qM(kd)*d_qM(id)
                       int= -tauU/h*0.25d0*res*SF_jP*SF_kP*wq
                       F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                       F1(i1) = F1(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int

                       res= (1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
                       int= -tauU/h*0.25d0*res*SF_jM*SF_kM*wq
                       F2(i1) = F2(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int
                       F2(i1) = F2(i1) - FM_prec(k1)*(1.d0/lambdaM)*(XeM(j,jd) - XeM(k,jd))*int

                       res= (1.d0/lambdaP)*GSF_iP(kd)*d_qP(kd)*d_qP(id)
                       int= tauU/h*0.25d0*res*SF_jP*SF_kP*wq
                       F2(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
                       F2(i1) = F2(i1) - FP_prec(k1)*(1.d0/lambdaP)*(XeP(j,jd) - XeP(k,jd))*int
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

                ! - < { q* } , [ B^n+1 ].n(n.~n) >
                DO id = 1, nDim
                   DO kd = 1, nDim
                      i1 = ii1 + (id-1)*nNodesPerElement
                      i2 = ii2 + (id-1)*nNodesPerElement
                      j1 = j   + (kd-1)*nNodesPerElement

                      int = -0.25d0*SF_i1M*SF_i2M*SF_jM*N_qM(kd)*nnt*wq
                      Amm(i1,j1) = Amm(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      Amm(i2,j1) = Amm(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      int = -0.25d0*SF_i1M*SF_i2M*SF_jP*N_qP(kd)*nnt*wq
                      Amp(i1,j1) = Amp(i1,j1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      Amp(i2,j1) = Amp(i2,j1) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      int = -0.25d0*SF_i1P*SF_i2P*SF_jM*N_qM(kd)*nnt*wq
                      Apm(i1,j1) = Apm(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      Apm(i2,j1) = Apm(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      int = -0.25d0*SF_i1P*SF_i2P*SF_jP*N_qP(kd)*nnt*wq
                      App(i1,j1) = App(i1,j1) + (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      App(i2,j1) = App(i2,j1) - (1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                   END DO
                END DO

                !  < { q* } , [ B^n ].n(n.~n) >
                DO id = 1, nDim
                   DO kd = 1, nDim
                      i1 = ii1 + (id-1)*nNodesPerElement
                      i2 = ii2 + (id-1)*nNodesPerElement
                      j1 = j   + (kd-1)*nNodesPerElement

                      int = 0.25d0*SF_i1M*SF_i2M*SF_jM*N_qM(kd)*nnt*wq
                      F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      F1(i2) = F1(i2) - FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      int = 0.25d0*SF_i1M*SF_i2M*SF_jP*N_qP(kd)*nnt*wq
                      F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      F1(i2) = F1(i2) - FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                      int = 0.25d0*SF_i1P*SF_i2P*SF_jM*N_qM(kd)*nnt*wq
                      F2(i1) = F2(i1) + FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      F2(i2) = F2(i2) - FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      int = 0.25d0*SF_i1P*SF_i2P*SF_jP*N_qP(kd)*nnt*wq
                      F2(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                      F2(i2) = F2(i2) - FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                   END DO
                END DO

                ! - < { q* } , [ D(B^n+1).d ].n(n.~n) >
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

                !  < { q* } , [ DB.d ].n(n.~n) >
                DO id = 1, nDim
                   DO kd = 1, nDim
                     DO jd= 1,nDim
                        i1 = ii1 + (id-1)*nNodesPerElement
                        i2 = ii2 + (id-1)*nNodesPerElement
                        j1 = j   + (kd-1)*nNodesPerElement
                        int = -0.25d0*SF_i1M*SF_i2M*GSF_jM(jd)*d_qM(jd)*N_qM(kd)*nnt*wq
                        F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                        F1(i2) = F1(i2) - FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                        int = 0.25d0*SF_i1M*SF_i2M*GSF_jP(jd)*N_qP(jd)*d_qP(kd)*nnt*wq
                        F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                        F1(i2) = F1(i2) - FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                        int = -0.25d0*SF_i1P*SF_i2P*GSF_jM(jd)*N_qM(jd)*d_qM(kd)*nnt*wq
                        F2(i1) = F2(i1) + FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                        F2(i2) = F2(i2) - FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                        int = 0.25d0*SF_i1P*SF_i2P*GSF_jP(jd)*N_qP(jd)*d_qP(kd)*nnt*wq
                        F1(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                        F1(i2) = F2(i2) - FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                     END DO
                   END DO
                END DO

                ! tauU/h < [ q* ] , [ T^n+1 ] >
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

                ! -tauU/h < [ q* ] , [ T^n ] >
                DO id = 1, nDim
                   i1 = ii1 + (id-1)*nNodesPerElement
                   i2 = ii2 + (id-1)*nNodesPerElement
                   j1 = j   + nDim*nNodesPerElement

                   int = -0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*wq
                   F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                   F1(i2) = F1(i2) - FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                   int = 0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*wq
                   F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                   F1(i2) = F1(i2) - FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                   int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*wq
                   F2(i1) = F2(i1) + FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                   F2(i2) = F2(i2) - FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

                   int = -0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*wq
                   F2(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                   F2(i2) = F2(i2) - FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                END DO

                ! -tauU/h < [ q* ] , [ LM (B^n+1).d ] >
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

                ! tauU/h < [ q* ] , [ LM (B^n).d ] >
                DO id = 1, nDim
                  DO jd= 1, nDim
                     i1 = ii1 + (id-1)*nNodesPerElement
                     i2 = ii2 + (id-1)*nNodesPerElement
                     j1 = j   + (jd-1)*nNodesPerElement

                     int = -0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
                     F1(i1) = F1(i1) + FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                     F1(i2) = F1(i2) - FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                     int = -0.5d0*tauU/h*SF_i1M*SF_i2M*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
                     F1(i1) = F1(i1) + FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                     F1(i2) = F1(i2) - FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                     int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jM*(1.d0/lambdaM)*d_qM(jd)*wq
                     F2(i1) = F2(i1) + FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                     F2(i2) = F2(i2) - FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

                     int = 0.5d0*tauU/h*SF_i1P*SF_i2P*SF_jP*(1.d0/lambdaP)*d_qP(jd)*wq
                     F2(i2) = F2(i2) - FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                     F2(i1) = F2(i1) + FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                  END DO
                END DO


                ! -tauU/h*0.5 < [ q* ],LM [d^t D(B^n+1).d] >
                DO id =1, nDim
                  DO jd =1, nDim
                    Do kd = 1,nDim
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

                ! tauU/h*0.5 < [ q* ],LM [d^t D(B^n).d] >
                DO id =1, nDim
                  DO jd =1, nDim
                    Do kd = 1,nDim
                    i1 = ii1 + (id-1)*nNodesPerElement
                    i2 = ii2 + (id-1)*nNodesPerElement
                    j1 = j   + (jd-1)*nNodesPerElement

                     res=(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
                     int=tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
                     F1(i1) = F1(i1) +  FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                     F1(i2) = F1(i2) -  FM_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                     res=(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
                     int= -tauU/h*0.25d0*SF_i1M*SF_i2M*res*wq
                     F1(i1) = F1(i1) +  FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                     F1(i2) = F1(i2) -  FP_prec(j1)*(1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

                     res=(1.d0/lambdaM)*GSF_jM(kd)*d_qM(kd)*d_qM(jd)
                     int= -tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
                     F2(i1) = F2(i1) +  FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                     F2(i2) = F2(i2) -  FM_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int

                     res=(1.d0/lambdaP)*GSF_jP(kd)*d_qP(kd)*d_qP(jd)
                     int= tauU/h*0.25d0*SF_i1P*SF_i2P*res*wq
                     F2(i1) = F2(i1) +  FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
                     F2(i2) = F2(i2) -  FP_prec(j1)*(1.d0/lambdaP)*(XeP(ii1,id) - XeP(ii2,id))*int
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

             ! - < {q*} , j2_bis(n.~n) >
             DO id= 1, nDim
                i1 = ii1 + (id-1)*nNodesPerElement
                i2 = ii2 + (id-1)*nNodesPerElement

                int=-0.25d0*SF_i1M*SF_i2M*jn2_bis*nnt*wq
                F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
                int=-0.25d0*SF_i1P*SF_i2P*jn2_bis*nnt*wq
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


            ! a < [q*] , jt_bis >
            DO id= 1, nDim
               i1 = ii1 + (id-1)*nNodesPerElement
               i2 = ii2 + (id-1)*nNodesPerElement

               int=-0.5d0*tauU/h*SF_i1M*SF_i2M*jt_bis*wq
               F1(i1) = F1(i1) + (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int
               F1(i2) = F1(i2) - (1.d0/lambdaM)*(XeM(ii1,id) - XeM(ii2,id))*int

               int=0.5d0*tauU/h*SF_i1P*SF_i2P*jt_bis*wq
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

                ! tauU/h < [ q* ] , [ (T^n+1)* ] >
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

                ! -tauU/h < [ q* ] , [ (T^n)* ] >
                DO id = 1, nDim
                   DO jd = 1, nDim
                      i1 = ii1 + (id-1)*nNodesPerElement
                      i2 = ii2 + (id-1)*nNodesPerElement
                      j1 = jj1 + (jd-1)*nNodesPerElement
                      j2 = jj2 + (jd-1)*nNodesPerElement

                      int = -0.25d0*(1.d0/(lambdaM**2))*tauU/h*SF_i1M*SF_i2M*SF_j1M*SF_j2M*wq
                      F1(i1) = F1(i1)+FM_prec(j1)*(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                      F1(i1) = F1(i1)-FM_prec(j2)*(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                      F1(i2) = F1(i2)-FM_prec(j1)*(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                      F1(i2) = F1(i2)+FM_prec(j2)*(XeM(ii1,id)-XeM(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

                      int = 0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1M*SF_i2M*SF_j1P*SF_j2P*wq
                      F1(i1) = F1(i1) + FP_prec(j1)*(XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                      F1(i1) = F1(i1) - FP_prec(j2)*(XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                      F1(i2) = F1(i2) - FP_prec(j1)*(XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                      F1(i2) = F1(i2) + FP_prec(j2)*(XeM(ii1,id)-XeM(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int

                      int = 0.25d0*(1.d0/(lambdaM*lambdaP))*tauU/h*SF_i1P*SF_i2P*SF_j1M*SF_j2M*wq
                      F2(i1) = F2(i1) + FM_prec(j1)*(XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                      F2(i1) = F2(i1) - FM_prec(j2)*(XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                      F2(i2) = F2(i2) - FM_prec(j1)*(XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int
                      F2(i2) = F2(i2) + FM_prec(j2)*(XeP(ii1,id)-XeP(ii2,id))*(XeM(jj1,jd)-XeM(jj2,jd))*int

                      int = -0.25d0*(1.d0/(lambdaP**2))*tauU/h*SF_i1P*SF_i2P*SF_j1P*SF_j2P*wq
                      F2(i1) = F2(i1) + FP_prec(j1)*(XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                      F2(i1) = F2(i1) - FP_prec(j2)*(XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                      F2(i2) = F2(i2) - FP_prec(j1)*(XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                      F2(i2) = F2(i2) + FP_prec(j2)*(XeP(ii1,id)-XeP(ii2,id))*(XeP(jj1,jd)-XeP(jj2,jd))*int
                   END DO
                END DO


             END DO
          END DO
       END IF

    END DO

    Amm=0.5d0*Amm ; Amp=0.5d0*Amp ; Apm = 0.5d0*Apm ; App = 0.5d0*App
    F1 = 0.5d0*F1 ; F2 = 0.5d0*F2

    DEALLOCATE ( Np, Nm, Xq_M, Xq_P, XqM, XqP, X_M, X_P, Exa_p, Exa_m)
    DEALLOCATE ( XeP, XeM )
    DEALLOCATE (GSF_iM, GSF_jM, GSF_iP, GSF_jP, GSF_kM, GSF_kP)

  END SUBROUTINE ComputePMMatrices_CG_Embedded_Nico
  !======================================================================!

END MODULE EdgeIntegral_CG_Icing_Nicolson
