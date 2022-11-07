MODULE DirichletIcingConformal_Mixed_Explicit

  USE Types
  USE PublicVar
  USE Permeability
  USE Quadrature
  USE ExactFunctionJumpt


  IMPLICIT NONE

CONTAINS

  !======================================================================!
  SUBROUTINE DirichletIcingConf_Mixed_Explicit(F, Floc_prec, ele, edgList, Data, iFace, iBC)
  !======================================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:)    , INTENT(INOUT) :: F
    REAL*8, DIMENSION(:)    , INTENT(IN)  	:: Floc_prec
    type(element)           , INTENT(IN)    :: ele
    type(edge), DIMENSION(:), INTENT(IN)    :: edgList
    type(DataStructure)     , INTENT(IN)    :: Data
    INTEGER                 , INTENT(IN)    :: iBC, iFace
    !-----------------------------------------------------
    type(Edge) :: Edg
    !------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, L_q, Xe
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_q, N, GSF_i, GSF_j, xq, x_q, GSF_i1, GSF_i2
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_j1, GSF_j2, GSF_k
    !-----------------------------------------------------------------------
    REAL*8 :: L, wq, x_iq, y_iq, z_iq, int, lambda,hn_q,hN_q_bis, pD_q_bis, alpha, Dt, hh 
    REAL*8 :: SF_i, SF_j, Nrm, h, pD_q, a, SF_k, SF_i1, SF_j1, SF_i2, SF_j2, a2, conv, tauP 
    !----------------------------------------------------------------------
    INTEGER :: i, j, iq, id, jd, nDim
    INTEGER :: jj, nNodesPerElement, nQuad, nVar, nEdgesPerElement
    !----------------------------------------------------------

    nDim = Data%nDim; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; nVar = nDim + 1 ; Dt = Data%Dt 

	lambda=ele%lambda
	hh = ele%A**(1.d0/nDim)

    ALLOCATE ( Lm1_q(nDim,nDim), L_q(nDim,nDim), Lm1(nNodesPerElement,nDim,nDim) )
    ALLOCATE ( Exa_q(nVar), N(nDim), GSF_i(nDim), GSF_i1(nDim), GSF_k(nDim),  &
    GSF_i2(nDim),GSF_j1(nDim), GSF_j2(nDim), GSF_j(nDim), xq(nDim), x_q(nDim) )
    ALLOCATE ( Xe(nNodesPerElement,nDim) )

	
    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, icas)
    END DO
    
    Xe = ele%Coor !coordonnes des points definissant l'element de reference
    F = 0.d0 ! Initialization 
    
    N = - ele%N(iFace,:) ! recuperation de la normale
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + N(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L = Nrm ; h = ele%A/L
    N = N/Nrm ! normalisation
    
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
				!passage a l'element courant car calcul fait sur l'element de reference
				x_q(id) = x_q(id) + SF_i*ele%Coor(i,id) 
			END DO
		END DO

		!coefficient for the stabilization
        a = alphaIP1*lambda/h**2 
       
		CALL ExactSolJumpt(ele%L, Exa_q, x_q, Data%T-Data%Dt, ele, nDim)
		pD_q = Exa_q(nDim+1)

		DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, xq, i)


			!stabilization term for the dirichlet condition
			DO j = 1, nNodesPerElement
				CALL SFEval(ele%eleType, SF_j, xq, j)
            
				! - a < q , T^n >
				F(i) = F(i) - a*SF_i*SF_j*Floc_prec(j)*wq
				
				!!  + < lambda DT .n, q >
				Do id =1, nDim 
					F(i) = F(i) + lambda*SF_i*GSF_j(id)*N(id)*Floc_prec(j)*wq
				END DO 
				

			END DO

			! a  < q , pD >
			F(i) = F(i) + a*SF_i*pD_q*wq

		END DO
    END DO
     
    DEALLOCATE ( Lm1_q, L_q, Exa_q, N, GSF_i, GSF_j, xq, x_q )

  END SUBROUTINE DirichletIcingConf_Mixed_Explicit
  !=============================================================!

END MODULE DirichleticingConformal_Mixed_Explicit
