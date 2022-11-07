MODULE EdgeIntegral_CG_Icing_Explicit

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
  USE ElementIntegral_Mixed_Icing_Explicit

  IMPLICIT NONE

CONTAINS


  !=====================================================================================!
  SUBROUTINE ComputePMMatrices_CG_Embedded_Explicit(mesh, edg, edgM, edgP, eleM, eleP, F1, F2,&
  Floc_prec1,Floc_prec2, Data, iLS, Sol_prec, Sol)
    !=====================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)     :: Edg
    type(Edge), DIMENSION(:), INTENT(IN)     :: edgM, edgP
    type(element)           , INTENT(IN)     :: eleM, eleP
    REAL*8, DIMENSION(:)    , INTENT(OUT)    :: F1,F2
    REAL*8, DIMENSION(:)    , INTENT(IN)  	 :: Floc_prec1, Floc_prec2
    type(DataStructure)     , INTENT(IN)     :: Data
    INTEGER                 , INTENT(IN)     :: iLS
    type(MeshType)          , INTENT(IN)     :: Mesh
    type(SolStructure) 		, INTENT(IN)     :: Sol_prec,Sol

    !----------------------------------------------------------
    type(edge)    							 :: edgIntM, edgIntP
    type(element)							 :: ELE_inter
    !----------------------------------------------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE    :: Np, Nm, Xq_M, Xq_P, Exa_p, Exa_m
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE    :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE    :: XeM, XeP
    REAL*8, DIMENSION(:)    , ALLOCATABLE    :: XxqP, XxqM
    REAL*8, DIMENSION(data%nDim)   		     :: d_qP, n_qP, d_qM, n_qM
    REAL*8, DIMENSION(data%nDim-1) 		     :: tntM, tntP
    REAL*8, DIMENSION(data%nDim,data%nDim-1) :: t_qP, t_qM
    REAL*8, DIMENSION(data%nDim+1) 			 :: exa1_q, exa2_q, exa1_q_prec, exa2_q_prec
    REAL*8, DIMENSION(data%nDim)   			 :: x_qP, x_qM
    REAL*8, DIMENSION(:)    , ALLOCATABLE	 :: GSF_iM, GSF_jM,GSF_iP, GSF_jP
    REAL*8, DIMENSION(:)    , ALLOCATABLE    :: GSF_kM, GSF_kp, Exa_p_prec, Exa_M_prec
	LOGICAL								  	 :: testjumpT
    !-------------------------------------------------------------
    REAL*8  								 :: Nrm, wq, jt, pD_M, pD_P, nnt,jt_bis, jn2_bis, Dt, hhP, hhM 
    REAL*8 									 :: SF_iM, SF_iP, SF_jM, SF_jP, a, ttM, ttP,res, mu_qM, mu_qP,tauPM,tauPP 
    REAL*8  								 :: SF_i1M, SF_i1P, SF_i2P, SF_i2M, SF_j1M, SF_j2P, SF_j1P,SF_j2M,alphaM,alphaP
    REAL*8  								 :: tauU, h, LS_qP, LS_qM, lambdaM, lambdaP,res1,pD_M_prec, pD_P_prec, a2M, a2P 
    REAL*8  								 :: lambdaMax,dBx,dBy
    REAL*8  								 :: Sigma_bP, Sigma_bM, jn2, L, int, SF_kM, SF_kP,res2, jump,Jumpfx,JumpFy
    REAL*8  								 :: Gam
    !------------------------------------------------------------
    INTEGER 								 :: iM, iP, i, jM, jP,iq, i1, j1, j, id, jd, ii
    INTEGER 								 :: ii1,ii2, jj,jj1,jj2,k,i2,k1, kd ,ld,j2,Ns_1, Ns_2
    INTEGER 								 :: nDim, nNodesPerElement, nQuad, nEdgesPerElement, nVar
    !--------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nVar = nDim + 1  ; nEdgesPerElement = eleM%nEdgesPerElement ; Dt = Data%Dt 
    hhM = eleM%A**(1.d0/nDim) ; hhP = eleP%A**(1.d0/nDim)

    ALLOCATE ( Np(nDim),Nm(nDim), Xq_M(nDim), Xq_P(nDim))
    ALLOCATE (  XxqP(nDim), XxqM(nDim) ) ; ALLOCATE ( Exa_p(nVar), Exa_M(nVar)  )
    ALLOCATE ( Exa_p_prec(nVar), Exa_M_prec(nVar)  )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( XeM(nNodesPerElement,nDim), XeP(nNodesPerElement,nDim) )
    ALLOCATE ( GSF_iM(nDim), GSF_jM(nDim), GSF_iP(nDim), GSF_jP(nDim) )
    ALLOCATE ( GSF_kP(nDim), GSF_kM(nDim))
	
    
	Ns_1 = edg%Vertex(1) ; Ns_2 = edg%Vertex(2)
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
	Np = Np/Nrm  ;  Nm = -Np
	
	!value of the characteristic length 
	h = (eleM%A + eleP%A)/(2.d0*L)
	
	IF(lambdaM<lambdaP) THEN 
		lambdaMax = lambdaP
	ELSE 
		lambdaMax = lambdaM
	END IF 
	tauU = alphaIP1*Nrm*lambdaMax  !stab from the Nitsche's method 	
	Gam = 0.d0 !1.d0/(lambdaMax*h)*0.5d0 !stab Nitsche for the the condition on DT = sigma 
	
	CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim)
	CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)
	!xqM,xqP coordinates on the reference triangle 
	
	!Initialization 
	F1 = 0.d0 ; F2 = 0.d0 ; nnt = 0.d0
	
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
     
		CALL computeNodalDistance(x_qP, d_qP, LS_qP, Data, iLS)
		CALL getNormalAndTangents(d_qP, n_qP, t_qP, nP, nnt, tntP, nDim)

		IF (nnt>1 .OR. nnt<0) THEN
			PRINT*,"Produit des deux normales", eleP%num,eleP%ref,eleM%Num, eleM%ref,"value", nnt
			PRINT*,"Problème sur le calcul des normales" 
			!Stop
		END IF

		n_qM = -n_qP
		d_qM = -d_qP
		nM = - nP 
		
		!condition sur le flux, condition de stefan 
		!IF(icas==001) THEN 
			!pas de de déplacement en y donc pas de normal sur la condition 
			!jn2 = Icing_rho*Icing_Lm*Icing_chi*sqrt(Icing_alpha_l/(Data%t-Data%DT)) ! with the analytical solution
			!CALL JumpPrevNewInter(Data, Mesh, Sol, Sol_prec,JumpFx,JumpFy) 	
			!jn2 = JumpFx
			!jump at time n on the surrogate at time n+1 
			CALL ExactSolJumpt(eleP%L, exa1_q_prec, x_qP + d_qP ,Data%T-Data%Dt, eleP, nDim) !solution a the previous time 
			CALL ExactSolJumpt(eleM%L, exa2_q_prec, x_qM - d_qM ,Data%T-Data%Dt, eleM, nDim)

			jn2 = DOT_PRODUCT(exa1_q_prec(1:nDim) - exa2_q_prec(1:nDim), n_qP) ! the condition is defined with the normal vector 
			!PRINT*,"Value of b jump",jn2,JumpFx

		!ELSE
		!	PRINT*,"Error in ComputePMMatrices_CG_Embedded_Explicit"
		!	STOP
		!END IF
	
		 
		!probleme ici sur le calcul du saut au temps précédent 
		!on a bien le saut de temperature nul pour le cas 001
		!IF(icas==001) THEN 
			!computation for time n on the interface n+1 
			!jump = 0.d0 
			!CALL ValJumpT(Data, Mesh, Sol_prec, edgIntP, jump, d_qP, iq,lambdaM,lambdaP) 
			!jt = jump
			
			CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T-Data%Dt, eleP, nDim)
			CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T-Data%Dt, eleM, nDim)
			jt =  exa1_q(nDim+1) -  exa2_q(nDim+1)
			!PRINT*,"Value of jt",jT,jump
        !ELSE
		!	PRINT*,"Error in ComputePMMatrices_CG_Embedded_Explicit"
		!	STOP
		!END IF
				
		CALL WeightQuad(edg%eleType, wq, iq)
		wq = wq*L
		
		! terme a l'interface
		DO i = 1, nNodesPerElement
			CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL GSFEval(GSF_iM, xq_M, i,eleM)
			CALL SFEval(eleP%eleType, SF_iP, xq_P, i) ; CALL GSFEval(GSF_iP, xq_P, i,eleP)

			DO j = 1, nNodesPerElement
				CALL SFEval(eleM%eleType, SF_jM, xq_M, j) ; CALL GSFEval(GSF_jM, xq_M, j,eleM)
				CALL SFEval(eleP%eleType, SF_jP, xq_P, j) ; CALL GSFEval(GSF_jP, xq_P, j,eleP)


				!  - < { q } , [ lambda DT^n ].n(n.~n) >
				!DO id = 1, nDim !for the use of the gradient and the normal vector 
				!	i1 = i 
				!	j1 = j 
				!	
				!	F1(i1) = F1(i1) - 0.5d0*SF_iM*GSF_jM(id)*n_qM(id)*nnt*Floc_prec1(j1)*lambdaM*wq
				!	F1(i1) = F1(i1) - 0.5d0*SF_iM*GSF_jP(id)*n_qP(id)*nnt*Floc_prec2(j1)*lambdaP*wq
				!	
				!	F2(i1) = F2(i1) - 0.5d0*SF_iP*GSF_jM(id)*n_qM(id)*nnt*Floc_prec1(j1)*lambdaM*wq
				!	F2(i1) = F2(i1) - 0.5d0*SF_iP*GSF_jP(id)*n_qP(id)*nnt*Floc_prec2(j1)*lambdaP*wq
	!
	!			END DO
				
				!   < [ q ] , { lambda/pc DT^n }.n~ >
				DO id = 1, nDim !for the use of the gradient and the normal vector 
	
					F1(i) = F1(i) + 0.5d0*SF_iM*GSF_jM(id)*nM(id)*Floc_prec1(j)*lambdaM*wq
					F1(i) = F1(i) - 0.5d0*SF_iM*GSF_jP(id)*nP(id)*Floc_prec2(j)*lambdaP*wq
					
					F2(i) = F2(i) - 0.5d0*SF_iP*GSF_jM(id)*nM(id)*Floc_prec1(j)*lambdaM*wq 
					F2(i) = F2(i) + 0.5d0*SF_iP*GSF_jP(id)*nP(id)*Floc_prec2(j)*lambdaP*wq
	
				END DO
				
				!   < { q } , [ lambda/pc DT^n ].n~ >
				DO id = 1, nDim 
	
					F1(i) = F1(i) + lambdaM*SF_iM*0.5d0*GSF_jM(id)*nM(id)*Floc_prec1(j)*wq
					F1(i) = F1(i) + lambdaP*SF_iM*0.5d0*GSF_jP(id)*nP(id)*Floc_prec2(j)*wq
					
					F2(i) = F2(i) + lambdaM*SF_iP*0.5d0*GSF_jM(id)*nM(id)*Floc_prec1(j)*wq
					F2(i) = F2(i) + lambdaP*SF_iP*0.5d0*GSF_jP(id)*nP(id)*Floc_prec2(j)*wq
	
				END DO
	
					! -  < { lambda/pc  Dq }.n~ , [ T^n ] >
				!DO id = 1, nDim 
	!
	!				F1(i) = F1(i) + lambdaM*GSF_iM(id)*nM(id)*0.5d0*SF_jM*Floc_prec1(j)*wq
	!				F1(i) = F1(i) - lambdaM*GSF_iM(id)*nM(id)*0.5d0*SF_jP*Floc_prec2(j)*wq
	!				
	!				F2(i) = F2(i) - lambdaP*GSF_iP(id)*nP(id)*0.5d0*SF_jM*Floc_prec1(j)*wq
	!				F2(i) = F2(i) + lambdaP*GSF_iP(id)*nP(id)*0.5d0*SF_jP*Floc_prec2(j)*wq
	!
	!			END DO
							

				!Nitsche's stab for the second jump condition
			
				!- Gam < [L Dq.n^~ + D(lambda Dq)n^~] , [L DT.n^~ + D(lambda DT)n^~]  >
				!only one term for a P1 resolution 
				DO id = 1, nDim 
					DO jd = 1, nDim
					
						F1(i) = F1(i) - lambdaM**2*Gam*GSF_iM(id)*nM(id)*GSF_jM(jd)*nM(jd)*Floc_prec1(j)*wq
						F1(i) = F1(i) - lambdaM*lambdaP*Gam*GSF_iM(id)*nM(id)*GSF_jP(jd)*nP(jd)*Floc_prec2(j)*wq
						
						F2(i) = F2(i) - lambdaP*lambdaM*Gam*GSF_iP(id)*nP(id)*GSF_jM(jd)*nM(jd)*Floc_prec1(j)*wq
						F2(i) = F2(i) - lambdaP**2*Gam*GSF_iP(id)*nP(id)*GSF_jP(jd)*nP(jd)*Floc_prec2(j)*wq
					
					END DO 
				END DO
				
			
			
				!   < { q } , [ DB^n d ].n(n.~n) > !should be 0 too because of the P1 formulation 
			!	DO id = 1, nDim
			!		DO jd=1, nDim
			!			i1 = i + nDim*nNodesPerElement
			!			j1 = j + (id-1)*nNodesPerElement
			!			
			!			F1(i1) = F1(i1) - 0.5d0*SF_iM*N_qm(id)*GSF_jM(jd)*d_qM(jd)*nnt*Floc_prec1(j1)*wq
			!			F1(i1) = F1(i1) + 0.5d0*SF_iM*N_qp(id)*GSF_jP(jd)*d_qP(jd)*nnt*Floc_prec2(j1)*wq
			!		
			!			F2(i1) = F2(i1) - 0.5d0*SF_iP*N_qm(id)*GSF_jM(jd)*d_qM(jd)*nnt*Floc_prec1(j1)*wq
			!			F2(i1) = F2(i1) + 0.5d0*SF_iP*N_qp(id)*GSF_jP(jd)*d_qP(jd)*nnt*Floc_prec2(j1)*wq
			!		
			!		END DO
			!	END DO

				!ajout d'un terme de stabilisation sur la temperature
				!- a <[q + Dq . d],[T^n + DT^n.d]>
				
				!- a <[q],[T^n]>

				F1(i) = F1(i) - tauU/h*SF_iM*SF_jM*Floc_prec1(j)*wq
				F1(i) = F1(i) + tauU/h*SF_iM*SF_jP*Floc_prec2(j)*wq
				
				F2(i) = F2(i) + tauU/h*SF_iP*SF_jM*Floc_prec1(j)*wq
				F2(i) = F2(i) - tauU/h*SF_iP*SF_jP*Floc_prec2(j)*wq

				
				! - a < [q] ,[ DT^n .d] >
				DO jd = 1, nDim

					int = tauU/h
					F1(i) = F1(i) + SF_iM*GSF_jM(jd)*d_qM(jd)*int*Floc_prec1(j)*wq
					
					int = tauU/h
					F1(i) = F1(i) + SF_iM*GSF_jP(jd)*d_qP(jd)*int*Floc_prec2(j)*wq

					int = tauU/h
					F2(i) = F2(i) - SF_iP*GSF_jM(jd)*d_qM(jd)*int*Floc_prec1(j)*wq

					int = tauU/h
					F2(i) = F2(i) - SF_iP*GSF_jP(jd)*d_qP(jd)*int*Floc_prec2(j)*wq
					
				END DO

				! - a < [  Dq .d] ,[T^n] >
				DO id = 1, nDim

					int = tauU/h
					F1(i) = F1(i) + GSF_iM(id)*d_qM(id)*SF_jM*int*Floc_prec1(j)*wq

					int = tauU/h
					F1(i) = F1(i) - GSF_iM(id)*d_qM(id)*SF_jP*int*Floc_prec2(j)*wq

					int = tauU/h
					F2(i) = F2(i) + GSF_iP(id)*d_qP(id)*SF_jM*int*Floc_prec1(j)*wq

					int = tauU/h
					F2(i) = F2(i) - GSF_iP(id)*d_qP(id)*SF_jP*int*Floc_prec2(j)*wq
				END DO

				! - a < [ Dq.d], [  DT^n .D ] >
				DO id = 1,nDim
					DO jd = 1,nDim

						int = tauU/h
						F1(i) = F1(i) - GSF_iM(id)*d_qM(id)*GSF_jM(jd)*d_qM(jd)*int*Floc_prec1(j)*wq

						int = tauU/h
						F1(i) = F1(i) - GSF_iM(id)*d_qM(id)*gSF_jP(jd)*d_qP(jd)*int*Floc_prec2(j)*wq

						int = tauU/h
						F2(i) = F2(i) - GSF_iP(id)*d_qP(id)*gSF_jM(jd)*d_qM(jd)*int*Floc_prec1(j)*wq

						int = tauU/h
						F2(i) = F2(i) - GSF_iP(id)*d_qP(id)*gSF_jP(jd)*d_qP(jd)*int*Floc_prec2(j)*wq
					END DO
				END DO
          
			END DO

			! - < {q} , j2^n (n.~n)>

			!F1(i) = F1(i) - 0.5d0*SF_iM*jn2*nnt*wq
			!F2(i) = F2(i) - 0.5d0*SF_iP*jn2*nnt*wq


			!ajout du terme sur le second membre issu de la stabilisation
			!a <[q+Dq.d],jT^n>

			!a  <[q],jT>
			F1(i) = F1(i) - tauU/h*SF_iM*jt*wq
			F2(i) = F2(i) + tauU/h*SF_iP*jt*wq

			!terme pour le second membre
            ! a < [ Dq.d] ,jT^n >
            DO id=1, nDim

				int= tauU/h*jT
				F1(i) = F1(i) + GSF_iM(id)*d_qM(id)*int*wq

				int= tauU/h*jT
				F2(i) = F2(i) + GSF_iP(id)*d_qP(id)*int*wq
            END DO
            
            !terme pour le second membre
            ! - Gam < [ L Dq . n^~] ,j2^n >
            DO id=1, nDim
           
				int = Gam
				F1(i) = F1(i) - int*lambdaM*GSF_iM(id)*nM(id)*jn2*wq
				
				int = Gam
				F2(i) = F2(i) - int*lambdaP*GSF_iP(id)*nP(id)*jn2*wq
            END DO
            

		END DO

    END DO
	
    DEALLOCATE ( Np, Nm, Xq_M, Xq_P, XqM, XqP, X_M, X_P, Exa_p, Exa_m)
    DEALLOCATE ( XeP, XeM )
    DEALLOCATE (GSF_iM, GSF_jM, GSF_iP, GSF_jP, GSF_kM, GSF_kP)

  END SUBROUTINE ComputePMMatrices_CG_Embedded_Explicit
  !======================================================================!


  END MODULE EdgeIntegral_CG_Icing_Explicit
