MODULE NeumannConformal_Mixed

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactSolutiont
  USE ExactFunctionJumpt


  IMPLICIT NONE

CONTAINS

  !====================================================================!
  SUBROUTINE NeumannConf_Mixed(Aloc, F, ele, edgList, Data, iFace, iBC)
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
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, L_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, xq, x_q
    !-------------------------------------------------------------------------
    REAL*8 :: L, wq, x_iq, y_iq, z_iq, SF_i1, SF_i2, SF_k
    REAL*8 :: SF_i, SF_j, Nrm, h, hN_q, a, int, lambda, hN_q_bis
    !-----------------------------------------------------
    INTEGER :: i, j, iq, id, jd, nDim, i1, j1, ii1, ii2, k,jj,j2
    INTEGER :: nNodesPerElement, nQuad, nVar, ii, i2, k1, kd
    INTEGER :: nEdgesPerElement
    LOGICAL :: NonAbiatic
    !------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement
    nVar = nDim + 1 ; Aloc = 0.d0 ; F = 0.d0
    	
    ALLOCATE ( Lm1(nNodesPerElement, nDim, nDim) )
    ALLOCATE ( Lm1_q(nDim,nDim), L_q(nDim,nDim), Xe(nNodesPerElement, nDim) )
    ALLOCATE ( Exa_q(nVar), N(nDim), GSF_i(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )

    Xe = ele%Coor
    N = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; N = N/Nrm ; L = Nrm
    h = ele%A/L
	lambda = ele%lambda !depends only of the element 
	
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


		IF(icas==001) THEN 
			hN_q = Data%BCVal(iBC) !valable sur le bord abiatic 2 et 4 
			IF(hN_q==0) THEN 
				!on ne se trouve pas sur le bord 3 ou la condition est special
				NonAbiatic = .FALSE.
			ELSE !voir cette condition 
				hN_q = -Icing_Tr*Icing_hte
				NonAbiatic = .TRUE.
			END IF 
		ELSE IF ( icas /= 0) THEN
			CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
			hN_q = DOT_PRODUCT(Exa_q(1:nDim),N)
			IF(Data%TSCHEME=="NICOLSON") THEN
				Exa_q = 0.d0
				CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T-Data%Dt, ele, nDim)
				hN_q_bis = DOT_PRODUCT(Exa_q(1:nDim),N)
				hn_q =(hn_q+hN_q_bis)*0.5d0
			END IF
		ELSE
			hN_q = Data%BCVal(iBC)
		END IF

       DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
			DO j = 1, nNodesPerElement
				CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
				!condition sur un bord ou on a pas de condition
				! de Dirichlet
				! < w.n , T >
				DO id = 1, nDim
					i1 = i + (id-1)*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) + SF_i*SF_j*N(id)*wq
				END DO

				!terme de Nitsche
				! - < q , B.n >
				DO id = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (id-1)*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) - SF_i*SF_j*N(id)*wq
				END DO

			END DO
	
			! -< q , hN >
			i1 = i + nDim*nNodesPerElement
			F(i1) = F(i1) - SF_i*hN_q*wq

			IF ( data%pressureEnrichment ) THEN
				DO ii = 1, nEdgesPerElement
					j = permut(ii,1) ; k = permut(ii,2)
					CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)

					! < w.nt , T* >
					DO id = 1, nDim
						DO jd = 1, nDim
							i1 = i + (id-1)*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement
							k1 = k + (jd-1)*nNodesPerElement
							int = 0.5d0*SF_i*SF_j*SF_k*N(id)*wq
							Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							
							   
						END DO
					END DO
					
				END DO
			END IF
			
		END DO
		
		IF ( data%TestFunctionEnrichment ) THEN
			DO ii = 1, nEdgesPerElement
				ii1 = permut(ii,1) ; ii2 = permut(ii,2)
				CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)
              ! - < q* , hN >
				DO id=1,nDim
					i1 = ii1 + (id-1)*nNodesPerElement
					i2 = ii2 + (id-1)*nNodesPerElement
					int = -0.5d0*SF_i1*SF_i2*hN_q*wq
					F(i1) = F(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
					F(i2) = F(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
				END DO
				
				
				DO j = 1, nNodesPerElement
					CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
					! - < q* , B.n >
					DO id = 1, nDim
						DO jd=1,nDim
							 i1 = ii1 + (id-1)*nNodesPerElement
							 i2 = ii2 + (id-1)*nNodesPerElement
							 j1 = j   + (jd-1)*nNodesPerElement
							 int = -0.5d0*SF_i1*SF_i2*SF_j*N(jd)*wq
							 Aloc(i1,j1) = Aloc(i1,j1) + (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
							 Aloc(i2,j1) = Aloc(i2,j1) - (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						END DO
                    END DO
                END DO
            END DO
        END IF


		IF(icas==001 .AND. NonAbiatic) THEN 
		    DO i = 1, nNodesPerElement
			    CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
			    DO j = 1, nNodesPerElement
					CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
			
					!  hte < q , T >
					i1 = i + nDim*nNodesPerElement
					j1 = j + nDim*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) + SF_i*SF_j*Icing_hte*wq
					  
				END DO

				IF ( data%pressureEnrichment ) THEN
					DO ii = 1, nEdgesPerElement
						j = permut(ii,1) ; k = permut(ii,2)
						CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)

						!  hte < q , T* >
						DO jd = 1, nDim
							i1 = i + nDim*nNodesPerElement
							j1 = j + (jd-1)*nNodesPerElement
							k1 = k + (jd-1)*nNodesPerElement
							int = Icing_hte*0.5d0*SF_i*SF_j*SF_k*wq
							Aloc(i1,j1) = Aloc(i1,j1) + (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
							Aloc(i1,k1) = Aloc(i1,k1) - (1.d0/lambda)*(Xe(j,jd) - Xe(k,jd))*int
						END DO			
					END DO
				END IF
			END DO
			
		    IF ( data%TestFunctionEnrichment ) THEN
				DO ii = 1, nEdgesPerElement
					ii1 = permut(ii,1) ; ii2 = permut(ii,2)
					CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)
					DO j = 1, nNodesPerElement
						CALL SFEval(ele%eleType, SF_j, xq, j)    ; CALL GSFEval(GSF_j, xq, j, ele)
					 
						!  hte < q* , T >
						DO id = 1, nDim
							i1 = ii1 + (id-1)*nNodesPerElement
							i2 = ii2 + (id-1)*nNodesPerElement
							j1 = j   + nDim*nNodesPerElement
							int = Icing_hte*0.5d0*SF_i1*SF_i2*SF_j*wq
							Aloc(i1,j1) = Aloc(i1,j1) + (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
							Aloc(i2,j1) = Aloc(i2,j1) - (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						END DO
					END DO
				  
				  
					DO jj = 1, nEdgesPerElement
						j = permut(jj,1) ; k = permut(jj,2)
						CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL SFEval(ele%eleType, SF_k, xq, k)

				  
						!  hte < q* , T* >
						DO id = 1, nDim
							DO jd=1,nDim
								i1 = ii1 + (id-1)*nNodesPerElement
								i2 = ii2 + (id-1)*nNodesPerElement
								j1 = j   + (jd-1)*nNodesPerElement
								j2 = k   + (jd-1)*nNodesPerElement
						 
								int = Icing_hte*(1.d0/(lambda**2))*0.25d0*SF_i1*SF_i2*SF_j*SF_k*wq

								Aloc(i1,j1) = Aloc(i1,j1) + (Xe(ii1,id)-Xe(ii2,id))*(Xe(j,jd)-Xe(k,jd))*int
								Aloc(i1,j2) = Aloc(i1,j2) - (Xe(ii1,id)-Xe(ii2,id))*(Xe(j,jd)-Xe(k,jd))*int
								Aloc(i2,j1) = Aloc(i2,j1) - (Xe(ii1,id)-Xe(ii2,id))*(Xe(j,jd)-Xe(k,jd))*int
								Aloc(i2,j2) = Aloc(i2,j2) + (Xe(ii1,id)-Xe(ii2,id))*(Xe(j,jd)-Xe(k,jd))*int
							END DO
						END DO
					END DO 
				END DO
			END IF
		END IF
    END DO

    DEALLOCATE ( Lm1_q, L_q, Exa_q, N, GSF_i, GSF_j, xq, x_q )

  END SUBROUTINE NeumannConf_Mixed
  !=============================================================!
  
  
  !====================================================================!
  SUBROUTINE NeumannHeatFluxJump(Aloc, F, ele, edgList, Data, iFace, iBC)
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
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, L_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, xq, x_q
    !-------------------------------------------------------------------------
    REAL*8 :: L, wq, x_iq, y_iq, z_iq, SF_i1, SF_i2, SF_k
    REAL*8 :: SF_i, SF_j, Nrm, h, hN_q, a, int, lambda, hN_q_bis
    !-----------------------------------------------------
    INTEGER :: i, j, iq, id, jd, nDim, i1, j1, ii1, ii2, k,jj,j2
    INTEGER :: nNodesPerElement, nQuad, nVar, ii, i2, k1, kd
    INTEGER :: nEdgesPerElement
    LOGICAL :: NonAbiatic
    !------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement
    nVar = nDim + 1 ; Aloc = 0.d0 ; F = 0.d0
    	
    ALLOCATE ( Lm1(nNodesPerElement, nDim, nDim) )
    ALLOCATE ( Lm1_q(nDim,nDim), L_q(nDim,nDim), Xe(nNodesPerElement, nDim) )
    ALLOCATE ( Exa_q(nVar), N(nDim), GSF_i(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )

    Xe = ele%Coor
    N = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; N = N/Nrm ; L = Nrm
    h = ele%A/L
	lambda = ele%lambda !depends only of the element 
	
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


		IF ( icas /= 0) THEN
			CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T, ele, nDim)
			hN_q = DOT_PRODUCT(Exa_q(1:nDim),N)
		ELSE
			hN_q = Data%BCVal(iBC)
		END IF

        DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
			DO j = 1, nNodesPerElement
				CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

				!terme de Nitsche
				! - < q , B.n >
				DO id = 1, nDim
					i1 = i + nDim*nNodesPerElement
					j1 = j + (id-1)*nNodesPerElement
					Aloc(i1,j1) = Aloc(i1,j1) - SF_i*SF_j*N(id)*wq
				END DO

			END DO
	
			! -< q , hN >
			i1 = i + nDim*nNodesPerElement
			F(i1) = F(i1) - SF_i*hN_q*wq
			
		END DO
		
		IF ( data%TestFunctionEnrichment ) THEN
			DO ii = 1, nEdgesPerElement
				ii1 = permut(ii,1) ; ii2 = permut(ii,2)
				CALL SFEval(ele%eleType, SF_i1, xq, ii1) ; CALL SFEval(ele%eleType, SF_i2, xq, ii2)
              ! - < q* , hN >
				DO id=1,nDim
					i1 = ii1 + (id-1)*nNodesPerElement
					i2 = ii2 + (id-1)*nNodesPerElement
					int = -0.5d0*SF_i1*SF_i2*hN_q*wq
					F(i1) = F(i1) + (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
					F(i2) = F(i2) - (1.d0/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
				END DO
				
				
				DO j = 1, nNodesPerElement
					CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
					! - < q* , B.n >
					DO id = 1, nDim
						DO jd=1,nDim
							 i1 = ii1 + (id-1)*nNodesPerElement
							 i2 = ii2 + (id-1)*nNodesPerElement
							 j1 = j   + (jd-1)*nNodesPerElement
							 int = -0.5d0*SF_i1*SF_i2*SF_j*N(jd)*wq
							 Aloc(i1,j1) = Aloc(i1,j1) + (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
							 Aloc(i2,j1) = Aloc(i2,j1) - (1/lambda)*(Xe(ii1,id) - Xe(ii2,id))*int
						END DO
                    END DO
                END DO
            END DO
        END IF

    END DO

    DEALLOCATE ( Lm1_q, L_q, Exa_q, N, GSF_i, GSF_j, xq, x_q )

  END SUBROUTINE NeumannHeatFluxJump
  !=============================================================!
  
END MODULE NeumannConformal_Mixed
