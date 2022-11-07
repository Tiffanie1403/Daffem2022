MODULE NeumannConformal_Mixed_Explicit

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactSolutiont
  USE ExactFunctionJumpt


  IMPLICIT NONE

CONTAINS

  !====================================================================!
  SUBROUTINE NeumannConf_Mixed_Explicit(F, Floc_prec, ele, edgList, Data, iFace, iBC)
  !====================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    REAL*8, DIMENSION(:)    , INTENT(IN)  	:: Floc_prec
    type(element)           , INTENT(IN)    :: ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !-------------------------------------------------------------------------
    type(Edge) :: Edg
    !-------------------------------------------------------------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, L_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, xq, x_q
    !-------------------------------------------------------------------------
    REAL*8 								  :: L, wq, x_iq, y_iq, z_iq, SF_i1, SF_i2, SF_k, Dt, hh, tauP 
    REAL*8 								  :: SF_i, SF_j, Nrm, h, hN_q, a, int, lambda, hN_q_bis, a2
    !-------------------------------------------------------------------------
    INTEGER 							  :: i, j, iq, id, jd, nDim
    INTEGER 							  :: nNodesPerElement, nQuad, nVar, ii, i2, k1, kd
    INTEGER 							  :: nEdgesPerElement
    LOGICAL 							  :: NonAbiatic
    !-------------------------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; Dt = Data%Dt
    nVar = nDim + 1 ; F = 0.d0 ; hh = ele%A**(1.d0/nDim)

    	
    ALLOCATE ( Lm1(nNodesPerElement, nDim, nDim),Lm1_q(nDim,nDim), L_q(nDim,nDim), Xe(nNodesPerElement, nDim) )
    ALLOCATE ( Exa_q(nVar), N(nDim), GSF_i(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )

    Xe  = ele%Coor
    N   = -ele%N(iFace,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; N = N/Nrm ; L = Nrm
    h = ele%A/L ; lambda=ele%lambda 
	
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
			
			CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T-Data%Dt, ele, nDim)
			hN_q = DOT_PRODUCT(Exa_q(1:nDim),N)
			
		ELSE
			hN_q = Data%BCVal(iBC)
		END IF
		
		
        DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
			DO j = 1, nNodesPerElement
				CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
				
				!terme de Nitsche
				!  - < q , lambda DT^n .n >
			!	DO id = 1, nDim
			!		F(i) = F(i) - lambda*SF_i*GSF_j(id)*N(id)*Floc_prec(j)*wq
			!	END DO
			END DO
			!   - < q , hN >
			F(i) = F(i) - SF_i*hN_q*wq
       END DO  
    END DO

    !Robin condition 
    IF(icas==001 .AND. NonAbiatic) THEN 
		DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
			DO j = 1, nNodesPerElement
				CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)
		
				!  -hte  < q , T^n >

				F(i) = F(i) - SF_i*SF_j*Floc_prec(j)*Icing_hte*wq
			
			END DO
		END DO 
	END IF 
    	
    DEALLOCATE ( Lm1_q, L_q, Exa_q, N, GSF_i, GSF_j, xq, x_q )

  END SUBROUTINE NeumannConf_Mixed_Explicit
  !=============================================================!
  
END MODULE NeumannConformal_Mixed_Explicit
