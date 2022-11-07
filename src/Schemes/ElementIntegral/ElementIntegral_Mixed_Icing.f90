MODULE ElementIntegral_Mixed_Icing

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE SourceTerm

  IMPLICIT NONE

CONTAINS

  !==============================================================!
  SUBROUTINE Compute_ElementIntegralIcing_Mixed(ele, Aloc, Floc, Floc_prec, Floc_2prec, Data)
  !==============================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: Floc
    REAL*8, DIMENSION(:)  , INTENT(IN)  :: Floc_prec, Floc_2prec 
    type(DataStructure)   , INTENT(IN)  :: Data
    !----------------------------------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, Xe, L_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_i, GSF_j, GSF_k
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_i1, GSF_i2, GSF_j1, GSF_j2
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Xq, X_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: phi_q, Gphi,phi_q_prec
    !------------------------------------------------------------
    REAL*8 :: SF_i, SF_j, h, SF_k
    REAL*8 :: SF_i1, SF_i2, SF_j1, SF_j2
    REAL*8 :: wq, mu_q, lambda, int1, int2
    REAL*8 :: tauU, tauP, int, Dt, T, a2, conv, alpha, alphaS 
    !---------------------------------
    INTEGER :: nNodesPerElement, nDim, nQuad
    INTEGER :: i, iq, j, id, jd, i1, i2, j1, j2, kd, md, nd
    INTEGER :: nVar, ii, k, k1, ld, nEdgesPerElement
    INTEGER :: ii1, ii2, jj, jj1, jj2
    !-----------------------------------------------

    nDim = Data%nDim ; nVar = nDim + 1 ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; Dt = Data%Dt ; T = Data%T
    int1 = 0.d0 ; int2 = 0.d0 ; ALLOCATE(Gphi(nDim))
    
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Lm1_q(nDim,nDim), L_q(nDim,nDim) )
    ALLOCATE ( GSF_i(nDim), GSF_j(nDim), GSF_k(nDim), Xq(nDim), X_q(nDim) )
    ALLOCATE ( phi_q(nVar), phi_q_prec(nVar), Xe(nNodesPerElement,nDim) )
    ALLOCATE ( GSF_i1(nDim), GSF_i2(nDim), GSF_j1(nDim), GSF_j2(nDim) )
    Aloc = 0.d0 ; Floc = 0.d0

	!paramètre physique pour le cas test 001
	IF(icas == 001 .OR. Data%ResisTance) THEN ! we also do the definition with the resistance test 
		IF(ele%ref == 1) THEN 
			conv = Icing_rho*Icing_c1
		ELSE IF (ele%ref == 2) THEN 
			conv = Icing_rho*Icing_c2
		ELSE 
			conv = Icing_c3*Icing_rhoAirfoil ! rho is not the same in the solid airfoil 
		END IF 
	ELSE 
		conv = 1.0 ! not defined for the case where it not stefan 
	END IF 
    
	!on definit le rho*c avec le terme de discrétisation temporelle 
    SELECT CASE (TRIM(ADJUSTL(Data%TScheme)))
    CASE ("EULER")
       a2 = 1.d0*conv
    CASE ("NICOLSON")
       a2 = 2.d0*conv
    CASE ("BDF2") 
	   a2 = 1.0*conv
	CASE DEFAULT 
		PRINT*,"Error in the choice of the time scheme"
		STOP
    END SELECT

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, Data%icas)
    END DO

    Xe = ele%Coor 
    lambda = ele%lambda
	h = ele%A**(1.d0/nDim)
	
    tauU = 0.5d0  ! momemtum stabilisation
    alphaS = Data%alphaS
    tauP = alphaS*h*h/lambda*0.5d0 ! coefficient for the DIV-DIV STAB  

    CALL GetNQuad(ele%eleType, nQuad) ! gives 3 in 2D number of quadrature node 

    DO iq = 1, nQuad

		CALL CoorQuad(ele%eleType, Xq, iq)
		CALL WeightQuad(ele%eleType, wq, iq)
		wq = wq*ele%A

		X_q = 0.d0
		DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, Xq, i)
			DO id = 1, nDim
				X_q(id) = X_q(id) + ele%Coor(i,id)*SF_i
			END DO
		END DO

		CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, Data%icas)
		CALL SetNodalMatrixPermeability(L_q, ele, data%icas)

        DO i = 1, nNodesPerElement
            CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
            DO j = 1, nNodesPerElement
                CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

                ! ( w , Lm1 B^(n+1) )
                DO id = 1, nDim
                    i1 = i + (id-1)*nNodesPerElement
                    j1 = j + (id-1)*nNodesPerElement
                    Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*SF_i*SF_j*wq
                END DO

                ! - ( div w , T^n+1 ) !terme de la formulation d'origine sans stab
                DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) - GSF_i(id)*SF_j*wq
                END DO

				! ( q , div B^n+1 )
				DO id = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (id-1)*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) + SF_i*GSF_j(id)*wq
				END DO

				IF(Data%TScheme=="EULER" .OR. Data%Tscheme == "NICOLSON") THEN 
					!ajout du terme avec la dependance en temps
					! a2/Dt ( q , T^n+1 )
					i1 = i + nDim*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) + (a2/Dt)*SF_i*SF_j*wq

					! a2/dT ( q, T^n )
					i1 = i + nDim*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Floc(i1) = Floc(i1) + (a2/Dt)*SF_i*SF_j*Floc_prec(j1)*wq
					
				ELSE IF(Data%Tscheme=="BDF2") THEN 
					!ajout du terme avec la dependance en temps
					! 3/2 Dt *rho c ( q , T^n+1 )
					i1 = i + nDim*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) + (a2/Dt)*3.d0/2.d0*SF_i*SF_j*wq

					! 2/dT *rho c ( q, T^n )
					i1 = i + nDim*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Floc(i1) = Floc(i1) + 2.d0*(a2/Dt)*SF_i*SF_j*Floc_prec(j1)*wq
					
					! -1/2 dT *rho c ( q, T^n-1 )
					i1 = i + nDim*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Floc(i1) = Floc(i1) - 0.5*(a2/Dt)*SF_i*SF_j*Floc_2prec(j1)*wq
				ELSE 
					PRINT*,"Error in the choice of the time SCHEME"
					STOP
				ENDIF 
				
			END DO

			! ( q , phi^n+1 )
			CALL Compute_SourceTerm(phi_q, X_q, T, ele, Data%icas, nDim)
			IF(DATA%TSCHEME == "NICOLSON") THEN
				CALL Compute_SourceTerm(phi_q_prec, X_q, T-Data%Dt, ele, Data%icas, nDim)
				phi_q(nDim+1) = 0.5d0*(phi_q(nDim+1)+phi_q_prec(nDim+1))
			END IF

			i1 = i + nDim*nNodesPerElement
			Floc(i1) = Floc(i1) + SF_i*phi_q(nDim+1)*wq

		END DO

		IF ( data%TestFunctionEnrichment ) THEN
			DO ii = 1, nEdgesPerElement
				ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
				CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)
				! ( q* , phi^n+1 )
				DO id = 1, nDim
					i1 = ii1 + (id-1)*nNodesPerElement
					i2 = ii2 + (id-1)*nNodesPerElement
					int = 0.5d0*SF_i1*SF_i2*phi_q(nDim+1)*wq
					Floc(i1) = Floc(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
					Floc(i2) = Floc(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
				END DO

			END DO
		END IF

		!ajout des termes de stabilisation
		DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
			DO j = 1, nNodesPerElement
				CALL SFEVal(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

				! - tauU ( Lm1 w , B^n+1 )
				DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + (id-1)*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) - tauU*(1.d0/lambda)*SF_i*SF_j*wq
				END DO

				! - tauU ( Lm1 w , L GT^n+1 )
				DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) - tauU*SF_i*GSF_j(id)*wq
				END DO

				! tauU ( Gq , L GT^n+1 )
				DO id = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) + tauU*GSF_i(id)*lambda*GSF_j(id)*wq
				END DO

				! tauU ( Gq , B^n+1 )
				DO id = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (id-1)*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) + tauU*GSF_i(id)*SF_j*wq
				END DO

				! tauP ( div w , div B^n+1 )
				DO id = 1, nDim
					DO jd = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + (jd-1)*nNodesPerElement
						Aloc(i1,j1) = Aloc(i1,j1) + tauP*GSF_i(id)*GSF_j(jd)*wq
					END DO
				END DO

				IF (DATA%Tscheme=="EULER" .OR. Data%Tscheme=="NICOLSON") THEN 
					!terme de stabilisation div-div
					! tauP a2/dt ( div w , T^n+1 )
					DO id = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement
						Aloc(i1,j1) = Aloc(i1,j1) +(a2/Dt)*tauP*GSF_i(id)*SF_j*wq
					END DO

					! tauP  a2/dT ( div w , T^n )
					DO id = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement
						Floc(i1) = Floc(i1) + tauP*(a2/Dt)*GSF_i(id)*SF_j*Floc_prec(j1)*wq
					END DO
				
				ELSE IF (Data%TSCHEME == "BDF2") THEN 
					!ajout du terme avec la dependance en temps
					! tauP a2/dt 3.d0/2.d0 ( div w , T^n+1 )
					DO id = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement
						Aloc(i1,j1) = Aloc(i1,j1) + 3.d0/2.d0*(a2/Dt)*tauP*GSF_i(id)*SF_j*wq
					END DO

					! tauP  2.d0*a2/dT ( div w , T^n )
					DO id = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement
						Floc(i1) = Floc(i1) + 2.d0*tauP*(a2/Dt)*GSF_i(id)*SF_j*Floc_prec(j1)*wq
					END DO
					
					! tauP  0.5d0*a2/dT ( div w , T^n-1 )
					DO id = 1, nDim
						i1 = i + (id-1)*nNodesPerElement
						j1 = j + nDim*nNodesPerElement
						Floc(i1) = Floc(i1) - 0.5d0*tauP*(a2/Dt)*GSF_i(id)*SF_j*Floc_2prec(j1)*wq
					END DO
				END IF

			END DO

			! tauP ( div w , phi^n+1 ) !associé au terme de stabilisation div-div
			DO id = 1, nDim
				i1 = i + (id-1)*nNodesPerElement
				Floc(i1) = Floc(i1) + tauP*GSF_i(id)*phi_q(nDim+1)*wq
			END DO

		END DO

		!pas de modif faite sur l'enrichissement pour l'instant
        ! Pressure enrichment components
        IF ( Data%PressureEnrichment ) THEN
            DO i = 1, nNodesPerElement
				CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
				DO ii = 1, nEdgesPerElement
					j = Permut(ii,1) ; k = Permut(ii,2)
					CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
					CALL SFEval(ele%eleType, SF_k, xq, k) ; CALL GSFEval(GSF_k, xq, k, ele)


					!  - ( div w , T^n+1* ) !terme de la formulation d'origine
					DO id = 1, nDim
					   DO jd = 1, nDim
						 i1 = i + (id-1)*nNodesPerElement
						 j1 = j + (jd-1)*nNodesPerElement
						 k1 = k + (jd-1)*nNodesPerElement
						 int = -0.5d0*GSF_i(id)*SF_j*SF_k*wq
						 Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
						 Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
					   END DO
					END DO

					! tauU ( - Lm1 w , L GT^(n+1)* )
					DO id = 1, nDim
						DO jd = 1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement
							k1 = k + (jd-1)*nNodesPerElement
							int = -tauU*0.5d0*SF_i*(SF_j*GSF_k(id) + SF_k*GSF_j(id))*wq
							Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
						END DO
					END DO


					! tauU ( Gq , L GT^n+1* )
					DO id = 1, nDim
						DO jd = 1, nDim
							i1 = i + nDim*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement
							k1 = k + (jd-1)*nNodesPerElement
							int = 0.5d0*tauU*GSF_i(id)*(SF_k*GSF_j(id) + SF_j*GSF_k(id))*wq
							Aloc(i1,j1) = Aloc(i1,j1) + (Xe(j,jd) - Xe(k,jd))*int
							Aloc(i1,k1) = Aloc(i1,k1) - (Xe(j,jd) - Xe(k,jd))*int
						END DO
					END DO

					IF(Data%TSCHEME == "NICOLSON" .OR. Data%TSCHEME== "EULER") THEN 
						!tauP a2/Dt ( div w , T^n+1* ) terme de stabilisation
						DO id = 1, nDim
						   DO jd = 1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement
								k1 = k + (jd-1)*nNodesPerElement
								int = (a2/Dt)*tauP*0.5d0*GSF_i(id)*SF_j*SF_k*wq
								Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
								Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							END DO
						END DO

						!tauP a2/Dt ( div w , T^n* ) terme de stabilisation
						DO id = 1, nDim
							DO jd = 1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement
								k1 = k + (jd-1)*nNodesPerElement
								int = (a2/Dt)*tauP*0.5d0*GSF_i(id)*SF_j*SF_k*wq
								Floc(i1) = Floc(i1) + (Floc_prec(j1)/lambda)*(Xe(j,jd) - Xe(k,jd))*int
								Floc(i1) = Floc(i1) - (Floc_prec(k1)/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							END DO
						END DO
					ELSE IF (DATA%TSCHEME=="BDF2") THEN 
						!tauP a2/Dt ( div w , T^n+1* ) terme de stabilisation
						DO id = 1, nDim
						    DO jd = 1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement
								k1 = k + (jd-1)*nNodesPerElement
								int = (a2/Dt)*tauP*0.5d0*GSF_i(id)*SF_j*SF_k*wq
								Aloc(i1,j1) = Aloc(i1,j1) + 3.d0/2.d0*(1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
								Aloc(i1,k1) = Aloc(i1,k1) - 3.d0/2.d0*(1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							END DO
						END DO

						!tauP a2/Dt ( div w , T^n* ) terme de stabilisation
						DO id = 1, nDim
							DO jd = 1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement
								k1 = k + (jd-1)*nNodesPerElement
								int = (a2/Dt)*tauP*0.5d0*GSF_i(id)*SF_j*SF_k*wq
								Floc(i1) = Floc(i1) + 2.d0*(Floc_prec(j1)/lambda)*(Xe(j,jd) - Xe(k,jd))*int
								Floc(i1) = Floc(i1) - 2.d0*(Floc_prec(k1)/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							END DO
						END DO
						
						!tauP a2/Dt ( div w , T^n-1* ) terme de stabilisation
						DO id = 1, nDim
							DO jd = 1, nDim
								i1 = i + (id-1)*nNodesPerElement
								j1 = j + (jd-1)*nNodesPerElement
								k1 = k + (jd-1)*nNodesPerElement
								int = (a2/Dt)*tauP*0.5d0*GSF_i(id)*SF_j*SF_k*wq
								Floc(i1) = Floc(i1) - 0.5d0*(Floc_2prec(j1)/lambda)*(Xe(j,jd) - Xe(k,jd))*int
								Floc(i1) = Floc(i1) + 0.5d0*(Floc_2prec(k1)/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							END DO
						END DO
					
					END IF 


					IF(Data%Tscheme=="NICOLSON" .OR. Data%TSCHEME=="EULER") THEN 
						! a2/Dt ( q , T^(n+1)* )
						DO id = 1, nDim
							i1 = i + nDim*nNodesPerElement
							j1 = j + (id-1)*nNodesPerElement
							k1 = k + (id-1)*nNodesPerElement
							int = 0.5d0*(a2/Dt)*SF_i*SF_j*SF_k*wq
							Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,id) - Xe(k,id))*int
							Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,id) - Xe(k,id))*int
						END DO

						! a2/Dt ( q , T^n* )
						DO id = 1, nDim
							i1 = i + nDim*nNodesPerElement
							j1 = j + (id-1)*nNodesPerElement
							k1 = k + (id-1)*nNodesPerElement
							int = (1.d0/lambda)*0.5d0*(a2/Dt)*SF_i*SF_j*SF_k*wq
							Floc(i1) = Floc(i1) + Floc_prec(j1)*(Xe(j,id) - Xe(k,id))*int
							Floc(i1) = Floc(i1) - Floc_prec(k1)*(Xe(j,id) - Xe(k,id))*int
						END DO
					ELSE IF(Data%TSCHEME == "BDF2") THEN 
						! a2/Dt ( q , T^(n+1)* )
						DO id = 1, nDim
							i1 = i + nDim*nNodesPerElement
							j1 = j + (id-1)*nNodesPerElement
							k1 = k + (id-1)*nNodesPerElement
							int = 0.5d0*(a2/Dt)*SF_i*SF_j*SF_k*wq
							Aloc(i1,j1) = Aloc(i1,j1) + 3.d0/2.d0*(1.d0/lambda)*(Xe(j,id) - Xe(k,id))*int
							Aloc(i1,k1) = Aloc(i1,k1) - 3.d0/2.d0*(1.d0/lambda)*(Xe(j,id) - Xe(k,id))*int
						END DO

						! a2/Dt ( q , T^n* )
						DO id = 1, nDim
							i1 = i + nDim*nNodesPerElement
							j1 = j + (id-1)*nNodesPerElement
							k1 = k + (id-1)*nNodesPerElement
							int = (1.d0/lambda)*0.5d0*(a2/Dt)*SF_i*SF_j*SF_k*wq
							Floc(i1) = Floc(i1) + 2.d0*Floc_prec(j1)*(Xe(j,id) - Xe(k,id))*int
							Floc(i1) = Floc(i1) - 2.d0*Floc_prec(k1)*(Xe(j,id) - Xe(k,id))*int
						END DO
					
					
						! a2/Dt ( q , T^n-1* )
						DO id = 1, nDim
							i1 = i + nDim*nNodesPerElement
							j1 = j + (id-1)*nNodesPerElement
							k1 = k + (id-1)*nNodesPerElement
							int = (1.d0/lambda)*0.5d0*(a2/Dt)*SF_i*SF_j*SF_k*wq
							Floc(i1) = Floc(i1) - 0.5d0*Floc_2prec(j1)*(Xe(j,id) - Xe(k,id))*int
							Floc(i1) = Floc(i1) + 0.5d0*Floc_2prec(k1)*(Xe(j,id) - Xe(k,id))*int
						END DO
					
					END IF 
				END DO
			END DO
		END IF

       ! Test function enrichment
       IF ( Data%TestFunctionEnrichment ) THEN
          DO ii = 1, nEdgesPerElement
            ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
            CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL GSFEval(GSF_i1, xq, ii1, ele)
            CALL SFEval(ele%eleType, SF_i2, xq, ii2) ; CALL GSFEval(GSF_i2, xq, ii2, ele)
               DO j = 1, nNodesPerElement
                  CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

                  ! ( q* , div B^n+1 )
                  DO id = 1, nDim
                     DO jd = 1, nDim
                        i1 = ii1 + (id-1)*nNodesPerElement
                        i2 = ii2 + (id-1)*nNodesPerElement
                        j1 = j   + (jd-1)*nNodesPerElement
                        int = 0.5d0*SF_i1*SF_i2*GSF_j(jd)*wq
                        Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                        Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                      END DO
                  END DO

				IF(DATA%TSCHEME == "NICOLSON" .OR. DATA%TSCHEME =="EULER") THEN 
					! a2/Dt ( q* , T^(n+1) )
					DO id = 1, nDim
						i1 = ii1 + (id-1)*nNodesPerElement
						i2 = ii2 + (id-1)*nNodesPerElement
						j1 = j   + nDim*nNodesPerElement

						int = 0.5d0*(a2/Dt)*SF_i1*SF_i2*SF_j*wq
						Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
					END DO

                    ! a2/DT ( q* , T^n)
                    DO id = 1, nDim
						i1 = ii1 + (id-1)*nNodesPerElement
						i2 = ii2 + (id-1)*nNodesPerElement
						j1 = j   + nDim*nNodesPerElement

						int = 0.5d0*(a2/Dt)*SF_i1*SF_i2*SF_j*Floc_prec(j1)*wq
						Floc(i1) = Floc(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						Floc(i2) = Floc(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                    END DO
				ELSE IF (DATA%TSCHEME=="BDF2") THEN 
				
					! a2/Dt * 3.d0/2.d0*( q* , T^(n+1) )
					DO id = 1, nDim
						i1 = ii1 + (id-1)*nNodesPerElement
						i2 = ii2 + (id-1)*nNodesPerElement
						j1 = j   + nDim*nNodesPerElement

						int = 0.5d0*(a2/Dt)*SF_i1*SF_i2*SF_j*wq
						Aloc(i1,j1) = Aloc(i1,j1) + 3.d0/2.d0*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						Aloc(i2,j1) = Aloc(i2,j1) - 3.d0/2.d0*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
					END DO

                    ! a2/DT * 2* ( q* , T^n)
                    DO id = 1, nDim
						i1 = ii1 + (id-1)*nNodesPerElement
						i2 = ii2 + (id-1)*nNodesPerElement
						j1 = j   + nDim*nNodesPerElement

						int = 0.5d0*(a2/Dt)*SF_i1*SF_i2*SF_j*Floc_prec(j1)*wq
						Floc(i1) = Floc(i1) + 2.d0*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						Floc(i2) = Floc(i2) - 2.d0*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                    END DO
                    
                    ! -0.5d0*a2/DT ( q* , T^n-1)
                    DO id = 1, nDim
						i1 = ii1 + (id-1)*nNodesPerElement
						i2 = ii2 + (id-1)*nNodesPerElement
						j1 = j   + nDim*nNodesPerElement

						int = 0.5d0*(a2/Dt)*SF_i1*SF_i2*SF_j*Floc_2prec(j1)*wq
						Floc(i1) = Floc(i1) - 0.5d0*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						Floc(i2) = Floc(i2) + 0.5d0*(1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                    END DO
				
				END IF 

                 !terme de stabilisation enrichie

                 ! tauU ( Gq* , B^n+1 )
                 DO id = 1, nDim
                   DO kd = 1, nDim
                    i1 = ii1 + (id-1)*nNodesPerElement
                    i2 = ii2 + (id-1)*nNodesPerElement
                    j1 = j   + (kd-1)*nNodesPerElement
                    int = tauU*0.5d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*SF_j*wq
                    Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                    Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
                  END DO
                 END DO

                 ! tauU ( Gq* , L GT^n+1)
                 DO id = 1, nDim
                    DO jd = 1, nDim
                       i1 = ii1 + (id-1)*nNodesPerElement
                       i2 = ii2 + (id-1)*nNodesPerElement
                       j1 = j   + nDim*nNodesPerElement
                       int = tauU*0.5d0*(SF_i1*GSF_i2(jd) + GSF_i1(jd)*SF_i2)*GSF_j(jd)*wq
                       Aloc(i1,j1) = Aloc(i1,j1) + (Xe(ii1,id) - Xe(ii2,id))*int
                       Aloc(i2,j1) = Aloc(i2,j1) - (Xe(ii1,id) - Xe(ii2,id))*int
                    END DO
                 END DO

                END DO
              END DO
            END IF

        IF ( Data%TestFunctionEnrichment .AND. Data%PressureEnrichment ) THEN
            DO ii = 1, nEdgesPerElement
                ii1 = Permut(ii,1) ; ii2 = Permut(ii,2)
                CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL GSFEval(GSF_i1, xq, ii1, ele)
                CALL SFEval(ele%eleType, SF_i2, xq, ii2) ; CALL GSFEval(GSF_i2, xq, ii2, ele)
                DO jj = 1, nEdgesPerElement
					jj1 = Permut(jj,1) ; jj2 = Permut(jj,2)
					CALL SFEval(ele%eleType, SF_j1, xq, jj1) ; CALL GSFEval(GSF_j1, xq, jj1, ele)
					CALL SFEval(ele%eleType, SF_j2, xq, jj2) ; CALL GSFEval(GSF_j2, xq, jj2, ele)


					! tauU( Gq* , L GT* )
					DO id = 1, nDim
						DO jd = 1, nDim
							DO kd = 1, nDim
								i1 = ii1 + (id-1)*nNodesPerElement
								i2 = ii2 + (id-1)*nNodesPerElement
								j1 = jj1 + (jd-1)*nNodesPerElement
								j2 = jj2 + (jd-1)*nNodesPerElement

								int = tauU*0.25d0*(SF_i1*GSF_i2(kd) + GSF_i1(kd)*SF_i2)*&
									(SF_j1*GSF_j2(kd) + GSF_j1(kd)*SF_j2)*wq

								Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
								(Xe(jj1,jd) - Xe(jj2,jd))*int
								Aloc(i2,j1) = Aloc(i2,j1) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
								(Xe(jj1,jd) - Xe(jj2,jd))*int
								Aloc(i1,j2) = Aloc(i1,j2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
								(Xe(jj1,jd) - Xe(jj2,jd))*int
								Aloc(i2,j2) = Aloc(i2,j2) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*&
								(Xe(jj1,jd) - Xe(jj2,jd))*int
							END DO
						END DO
					END DO

					IF(Data%TSCHEME == "NICOLSON" .OR. Data%TSCHEME == "EULER") THEN 
						! a2/Dt ( q* , T^(n+1)* )
						DO id = 1, nDim
							DO jd = 1, nDim
							   i1 = ii1 + (id-1)*nNodesPerElement
							   i2 = ii2 + (id-1)*nNodesPerElement
							   j1 = jj1 + (jd-1)*nNodesPerElement
							   j2 = jj2 + (jd-1)*nNodesPerElement

							   int = (1.d0/(lambda**2))*0.25d0*(a2/Dt)*SF_i1*SF_i2*SF_j1*SF_j2*wq
							   Aloc(i1,j1) = Aloc(i1,j1) + (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Aloc(i1,j2) = Aloc(i1,j2) - (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Aloc(i2,j1) = Aloc(i2,j1) - (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Aloc(i2,j2) = Aloc(i2,j2) + (Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							END DO
						END DO

						! a2/Dt ( q* , T^n* )
						DO id = 1, nDim
							DO jd = 1, nDim
							   i1 = ii1 + (id-1)*nNodesPerElement
							   i2 = ii2 + (id-1)*nNodesPerElement
							   j1 = jj1 + (jd-1)*nNodesPerElement
							   j2 = jj2 + (jd-1)*nNodesPerElement

							   int = (1.d0/(lambda**2))*0.25d0*(a2/Dt)*SF_i1*SF_i2*SF_j1*SF_j2*wq
							   Floc(i1) = Floc(i1) + Floc_prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i1) = Floc(i1) - Floc_prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i2) = Floc(i2) - Floc_prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i2) = Floc(i2) + Floc_prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							END DO
						END DO
					ELSE IF (DATA%TSCHEME=="BDF2") THEN 
						! a2/Dt 3/2 ( q* , T^(n+1)* )
						DO id = 1, nDim
							DO jd = 1, nDim
							   i1 = ii1 + (id-1)*nNodesPerElement
							   i2 = ii2 + (id-1)*nNodesPerElement
							   j1 = jj1 + (jd-1)*nNodesPerElement
							   j2 = jj2 + (jd-1)*nNodesPerElement

							   int = (1.d0/(lambda**2))*0.25d0*(a2/Dt)*SF_i1*SF_i2*SF_j1*SF_j2*wq
							   Aloc(i1,j1) = Aloc(i1,j1) + 3.d0/2.d0*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Aloc(i1,j2) = Aloc(i1,j2) - 3.d0/2.d0*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Aloc(i2,j1) = Aloc(i2,j1) - 3.d0/2.d0*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Aloc(i2,j2) = Aloc(i2,j2) + 3.d0/2.d0*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							END DO
						END DO

						! a2/Dt 2*( q* , T^n* )
						DO id = 1, nDim
							DO jd = 1, nDim
							   i1 = ii1 + (id-1)*nNodesPerElement
							   i2 = ii2 + (id-1)*nNodesPerElement
							   j1 = jj1 + (jd-1)*nNodesPerElement
							   j2 = jj2 + (jd-1)*nNodesPerElement

							   int = (1.d0/(lambda**2))*0.25d0*(a2/Dt)*SF_i1*SF_i2*SF_j1*SF_j2*wq
							   Floc(i1) = Floc(i1) + 2.d0*Floc_prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i1) = Floc(i1) - 2.d0*Floc_prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i2) = Floc(i2) - 2.d0*Floc_prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i2) = Floc(i2) + 2.d0*Floc_prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							END DO
						END DO
					
					
						! a2/Dt 2*( q* , T^n-1* )
						DO id = 1, nDim
							DO jd = 1, nDim
							   i1 = ii1 + (id-1)*nNodesPerElement
							   i2 = ii2 + (id-1)*nNodesPerElement
							   j1 = jj1 + (jd-1)*nNodesPerElement
							   j2 = jj2 + (jd-1)*nNodesPerElement

							   int = (1.d0/(lambda**2))*0.25d0*(a2/Dt)*SF_i1*SF_i2*SF_j1*SF_j2*wq
							   Floc(i1) = Floc(i1) - 0.5d0*Floc_2prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i1) = Floc(i1) + 0.5d0*Floc_2prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i2) = Floc(i2) + 0.5d0*Floc_2prec(j1)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							   Floc(i2) = Floc(i2) - 0.5d0*Floc_2prec(j2)*(Xe(ii1,id)-Xe(ii2,id))*(Xe(jj1,jd)-Xe(jj2,jd))*int
							END DO
						END DO
					
					END IF 

				END DO
			END DO
        END IF
    END DO
    
    DEALLOCATE ( Lm1_q, GSF_i, GSF_j, GSF_k, Xe, Xq, X_q, phi_q, L_q )
    DEALLOCATE ( GSF_i1, GSF_i2, GSF_j1, GSF_j2 )

  END SUBROUTINE Compute_ElementIntegralIcing_Mixed
  !================================================================!

END MODULE ElementIntegral_Mixed_Icing
