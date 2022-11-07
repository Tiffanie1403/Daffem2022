MODULE ElementIntegral_Mixed_Icing_Explicit

  USE Types
  USE PublicVar
  USE Quadrature
  USE Permeability
  USE SourceTerm
  USE Algebra

  IMPLICIT NONE

CONTAINS

  !==============================================================!
  SUBROUTINE Compute_ElementIntegralIcing_Mixed_Explicit(ele, Aloc, Floc, Floc_prec, Data)
  !==============================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: Aloc
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: Floc
    REAL*8, DIMENSION(:)  , INTENT(IN)  :: Floc_prec
    type(DataStructure)   , INTENT(IN)  :: Data
    !----------------------------------------------
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Lm1
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: Lm1_q, Xe, L_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: GSF_i, GSF_j
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Xq, X_q
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: phi_q,phi_q_prec
    !------------------------------------------------------------
    REAL*8 :: SF_i, SF_j, h
    REAL*8 :: wq, mu_q, lambda, int1, int2
    REAL*8 :: tauU, tauP, int, Dt, T,alpha
    !---------------------------------
    INTEGER :: nNodesPerElement, nDim, nQuad
    INTEGER :: i, iq, j, id, jd
    INTEGER :: nVar, k, k1, ld, nEdgesPerElement
    !-----------------------------------------------
		
    nDim = Data%nDim ; nVar = nDim + 1 ; nNodesPerElement = ele%nNodesPerElement
    nEdgesPerElement = ele%nEdgesPerElement ; Dt = Data%Dt ; T=Data%T
    int1=0.d0 ; int2=0.d0 
    
    ALLOCATE ( Lm1(nNodesPerElement,nDim,nDim), Lm1_q(nDim,nDim), L_q(nDim,nDim) )
    ALLOCATE ( GSF_i(nDim), GSF_j(nDim), Xq(nDim), X_q(nDim) )
    ALLOCATE ( phi_q(nVar),phi_q_prec(nVar), Xe(nNodesPerElement,nDim))
    
    h = ele%A**(1.d0/nDim)
    tauU = 0.5d0  ! momemtum stabilisation
		
    Aloc = 0.d0 ; Floc = 0.d0

    DO i = 1, nNodesPerElement
       CALL SetNodalInverseMatrixPermeability(Lm1(i,:,:), ele, Data%icas)
    END DO

    Xe = ele%Coor 
    lambda=ele%lambda
	

    CALL GetNQuad(ele%eleType, nQuad) !number of quadrature node 

    DO iq = 1, nQuad

		CALL CoorQuad(ele%eleType, Xq, iq) !coordinate of the quadrature node 
		CALL WeightQuad(ele%eleType, wq, iq)
		wq = wq*ele%A

		X_q = 0.d0
		DO i = 1, nNodesPerElement
			CALL SFEval(ele%eleType, SF_i, Xq, i)
			DO id = 1, nDim
				!envoie les coordonées du noeuds sur l'element de référence 
				X_q(id) = X_q(id) + ele%Coor(i,id)*SF_i ! used for the value of the second member f 
			END DO
		END DO

		CALL SetNodalInverseMatrixPermeability(Lm1_q, ele, Data%icas)
		CALL SetNodalMatrixPermeability(L_q, ele, data%icas)

        DO i = 1, nNodesPerElement
            CALL SFEval(ele%eleType, SF_i, xq, i) ; CALL GSFEval(GSF_i, xq, i, ele)
            DO j = 1, nNodesPerElement
                CALL SFEval(ele%eleType, SF_j, xq, j) ; CALL GSFEval(GSF_j, xq, j, ele)

                ! the test function is only P1 so the second derivative should be nulle 
				!  lambda Dt/pc ( q , Grad(DT^n) )
				!DO id = 1, nDim
				!	i1 = i + nDim*nNodesPerElement
				!	j1 = j + (id-1)*nNodesPerElement
				!	Floc(i1) = Floc(i1) + lambda*SF_i*GSF_j(id)*Floc_prec(j1)*wq 
				!END DO

				! - lambda ( Dq , DT^n )
				DO id = 1, nDim
					Floc(i) = Floc(i) - lambda*GSF_i(id)*GSF_j(id)*Floc_prec(j)*wq 
				END DO

			END DO ! end on j 

			! ( q , phi^n )
			CALL Compute_SourceTerm(phi_q, X_q, Data%T-Data%Dt, ele, Data%icas, nDim) !second member at time n 
			Floc(i) = Floc(i) + SF_i*phi_q(nDim+1)*wq
		END DO
		    
    END DO
    DEALLOCATE ( Lm1_q, GSF_i, GSF_j, Xe, Xq, X_q, phi_q, L_q )

  END SUBROUTINE Compute_ElementIntegralIcing_Mixed_Explicit
  
  !================================================================!
  
  
  
  !================================================================!
  SUBROUTINE ReconstructDiv(ele,Floc_prec,nDim,nNodesPerElement,dBx,dBy)
	
	IMPLICIT NONE
    type(element)          , INTENT(IN)     :: ele
    REAL*8, DIMENSION(:)   , INTENT(IN)     :: Floc_prec
	REAL*8			       , INTENT(OUT)    :: dBx,dBy 
	INTEGER				   , INTENT(IN)		:: nDim, nNodesPerElement
	!------------------------------------------------
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: A, Am1
	REAL*8, DIMENSION(:)  , ALLOCATABLE :: Bx, By,Solx,Soly
	!------------------------------------------------
	REAL*8								:: xi, yi
	INTEGER								:: i,j 
	

	ALLOCATE(A(nNodesPerElement,nDim+1),Am1(nNodesPerElement,nDim+1))
	ALLOCATE(Bx(nNodesPerElement),By(nNodesPerElement))
	ALLOCATE(Solx(nNodesPerElement), Soly(nNodesPerElement))

	A = 0.d0 ; Bx = 0.d0 ; By = 0.d0 ; Solx = 0.d0 ; Soly = 0.d0 
	
	DO i =1, nNodesPerElement
		
		xi = ele%Coor(i,1) ; yi = ele%Coor(i,2)
		 
		A(i,1) = 1.d0 ; A(i,2) = xi ; A(i,3) = yi
		Bx(i) = Floc_prec(1+(i-1)*3)
		By(i) = Floc_prec(2+(i-1)*3)
		
		 
	END DO
	
	CALL Inverse3x3Matrix(A, Am1)
	
	DO i =1,nNodesPerElement
		DO j=1, nNodesPerElement
	
		Solx(i) = Solx(i) + Am1(i,j)*Bx(j)  
		Soly(i) = Soly(i) + Am1(i,j)*By(j)
		
		END DO 
	END DO 
	
	dBx = Solx(2) !dx Bx
	dBy = Soly(3) !dy By  
	
	DEALLOCATE(A,AM1,Bx,By,Solx,soly)
	
	 
  END SUBROUTINE ReconstructDiv
  !================================================================!


END MODULE ElementIntegral_Mixed_Icing_Explicit
