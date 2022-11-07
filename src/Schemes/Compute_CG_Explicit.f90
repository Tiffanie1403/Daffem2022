MODULE Compute_CG_Explicit

  USE Types
  USE PublicVar
  USE SparseMatrix_CG
  USE BCSelection
  USE ReadMesh_mod
  USE Affichage_M
  USE ElementIntegral_Mixed_Icing_Explicit
  USE EdgeIntegral_CG_Icing_Explicit
  USE DirichletIcingConformal_Mixed_Explicit 
  USE NeumannConformal_Mixed_Explicit
  USE Compute_CG 

  IMPLICIT NONE

CONTAINS

	!==================================================!
	SUBROUTINE Compute_LHS_RHS_CG_Explicit(Data, Mesh, Sol, Sol_prec, Mat)
	!==================================================!

	IMPLICIT NONE
    type(DataStructure), INTENT(INOUT) :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SolStructure) , INTENT(IN)    :: Sol_prec
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    !-------------------------------------------
    type(element) :: ele, eleM, eleP, eleBC, Adj
    type(Edge)    :: Edg
    !--------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Aloc
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Floc,Floc1,Floc2
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Floc_prec, Floc_prec1, Floc_prec2

    !--------------------------------------------------
    CHARACTER(len=20) :: BCType
    !--------------------------
    INTEGER :: AdjNum, nDim, NeActive, iLS, Nu_i 
    INTEGER :: Nv, Ne, Nb, nVar, Ns, ib, iFace, iFaceM, iFaceP
    INTEGER :: ie, i, j, k, iBC, ii, nNodesPerElement, Ndof, Nv_bis
    !--------------------------------------------------
    LOGICAL :: eleSolve, eleID, esM, esP
    !--------------------------

    IF (Data%Icing .eqv. .TRUE.) THEN
       Ns=Mesh%NarEmb + 1 + mesh%PtEmb ! nombre de noeud sur l'interface
    ELSE
       Ns=0
    END IF

    nDim = Data%nDim
    Ne = Mesh%Ne 
    Nv = Mesh%NvActive + Ns 
    Nv_bis = Mesh%NvActive + Data%MaxRowEmbedded
    Nb = Mesh%Nb 
    nVar = Mat%NVar 
    Ndof =  Mesh%NvActive + Data%MaxRowEmbedded
    Mat%Coeff = 0.d0 
    Sol%RHS = 0.d0 
    NeActive = Mesh%NeActive
	Sol%DualArea = 0.d0 
	
    !boucle sur les elements pour les integrales sur Omega
    DO ie = 1, Ne
		ele = Mesh%ele(ie)		
		nNodesPerElement = ele%nNodesPerElement
		
		!filling of the vector for the dual area associated to each node of the mesh 
		DO i = 1, nNodesPerElement 
			Nu_i = ele%SurVertex(i)
			IF(ele%ref==1) THEN
				IF(icas==001) THEN 
					Sol%DualArea(Nu_i) = Sol%DualArea(Nu_i) + ele%A/3.d0*Icing_rho*Icing_c1
				ELSE 
					Sol%DualArea(Nu_i) = Sol%DualArea(Nu_i) + ele%A/3.d0
				END IF 
			ELSE
				IF(icas==001) THEN 
					Sol%DualArea(Nu_i) = Sol%DualArea(Nu_i) + ele%A/3.d0*Icing_rho*Icing_c2
				ELSE 
					Sol%DualArea(Nu_i) = Sol%DualArea(Nu_i) + ele%A/3.d0
				END IF 
			END IF 
		END DO 
		
		
		ALLOCATE ( Aloc(nVar*nNodesPerElement,nVar*nNodesPerElement) ) !allocation de la matrice locale sur l'element
		ALLOCATE ( Floc(nVar*nNodesPerElement), Floc_prec(nVar*nNodesPerElement) )
		Aloc = 0.d0 ; Floc = 0.d0 ; Floc_prec = 0.d0 

		CALL Remplissage_FlocPrec(ele, Floc_prec, Sol_prec, nDim,Data)
       
		IF ( ele%solve ) THEN
			CALL Compute_ElementIntegralIcing_Mixed_Explicit(ele, Aloc, Floc, Floc_prec, Data) 
			CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv_bis, Ne)
			CALL AssembleRHSCG(Sol%RHS, Floc, ele, Ndof, Ne, nVar)           
		END IF
		DEALLOCATE ( Aloc, Floc, Floc_prec)
    END DO

    !Conformal BC
    !boucle sur les aretes de bords interface comprise
    
    DO ib = 1, mesh%nb
		edg = mesh%boundary(ib) !numerotation de la ib ieme arete
		ele = mesh%ele(edg%tri(1)) ! on se place sur le premier element possédant cette arete de bord, il peut y avoir deux elements si on est sur l'interface
		nNodesPerElement = ele%nNodesPerElement
		IF ( ele%Solve ) THEN
		
			nNodesPerElement = ele%nNodesPerElement
			ALLOCATE ( Aloc(nVar*nNodesPerElement,nVar*nNodesPerElement) )
			ALLOCATE ( Floc(nVar*nNodesPerElement) , Floc_prec(nVar*nNodesPerElement))
			Aloc = 0.d0 ; Floc = 0.d0 ; Floc_prec = 0.d0 

			CALL SelectBC2(edg, BCType, Data, iBC) 
			
			CALL GetFace(ele, edg, iFace) 
			CALL Remplissage_FlocPrec(ele,Floc_prec, Sol_prec,nDim,Data)

			SELECT CASE ( TRIM(ADJUSTL(BCType)) )

			CASE ( "Dirichlet", "dirichlet" )

				CALL DirichletIcingConf_Mixed_Explicit(Floc , Floc_prec, ele, mesh%edg(ele%edg), Data, iFace, iBC)
				CALL AssembleRHSCG(Sol%RHS, Floc, ele, Ndof, Ne, nVar)

			CASE ( "Neumann", "neumann" )
	
				CALL NeumannConf_Mixed_Explicit(Floc ,Floc_prec, ele, mesh%edg(ele%edg), Data, iFace, iBC)
				CALL AssembleRHSCG(Sol%RHS, Floc, ele, Ndof, Ne, nVar)

			CASE ( "NeumannJump","neumannJump","neumannjump","Neumannjump" )
				PRINT*,"I shouldn't use this part here" 
				!STOP 
			CASE DEFAULT
				IF(Data%Icing .AND. Data%Embedded) THEN
					IF(mesh%boundary(ib)%ref == Data%N_BC+ 1) THEN
					!do nothing
					!permet de laisser les aretes de l'interface dans le fichier mesh
					END IF
				ELSE
				   PRINT*, " ** ERROR in BC type... (Compute_CG_Explicit.f90)"
				   STOP
				END IF
			END SELECT
			DEALLOCATE ( Aloc, Floc, Floc_prec )
		END IF
    END DO


    !Embedded BC
    !Condition sur les aretes de la surrogate
    IF ( Data%Embedded .AND. Data%Icing ) THEN
        DO ib = 1, mesh%nSurB !mesh%nSurB est le nombre d'arete sur la surrogate
			iLS= mesh%EmbBoundary(ib)%EmbLS
			!on se place sur une des aretes de la surrogate
			edg = mesh%embBoundary(ib) !mesh%embBoundary enregistre les aretes sur la surrogate
			
			IF(edg%tri(2)/=-1 .AND. edg%tri(1) /=-1) THEN ! on verifie que l'on se trouve pas sur une arete de bord 
				nNodesPerElement = nDim + 1
				ALLOCATE ( Aloc(nVar*nNodesPerElement, nVar*nNodesPerElement) )
				ALLOCATE ( Floc(nVar*nNodesPerElement) )

				SELECT CASE ( TRIM(ADJUSTL(mesh%embBoundary(ib)%BC_Type)) )
				!on peut d'abord ce concentrer sur le cas NeumannJump car c'est celui qui nous intéresse dans notre cas
				CASE ( "Dirichlet","dirichlet","Neumann","neumann" )
					PRINT*,"I shouldn't use this part here / Compute_LHS_RHS_CG_Explicit" 
					STOP 
				CASE ( "NeumannJump","neumannJump","neumannjump","Neumannjump" )

					IF (mesh%ele(edg%tri(2))%ref == IcingZone2_ref) THEN
						eleP = mesh%ele(edg%tri(1))
						eleM = mesh%ele(edg%tri(2))
					ELSE
						eleP = mesh%ele(edg%tri(2))
						eleM = mesh%ele(edg%tri(1))
					END IF

					nNodesPerElement = eleM%nNodesPerElement

					ALLOCATE ( Floc1(nVar*nNodesPerElement), Floc2(nVar*nNodesPerElement) )
					ALLOCATE ( Floc_prec1(nVar*nNodesPerElement), Floc_prec2(nVar*nNodesPerElement) )
					Floc_prec1 = 0.d0 ; Floc_prec2 = 0.d0 
					Floc1 = 0.d0      ; Floc2 = 0.d0 
					
					CALL Remplissage_FlocPrec(eleM, Floc_prec1, Sol_prec, nDim,Data)
					CALL Remplissage_FlocPrec(eleP, Floc_prec2, Sol_prec, nDim,Data)

					CALL ComputePMMatrices_CG_Embedded_Explicit(Mesh, Edg,  Mesh%Edg(eleM%Edg), Mesh%Edg(eleP%Edg), &
							EleM, eleP, Floc1, Floc2, Floc_prec1, Floc_prec2, Data, iLS, Sol_prec, Sol)
					
					CALL AssembleRHSCG(Sol%RHS, Floc1, eleM, Ndof, Ne, nVar)
					CALL AssembleRHSCG(Sol%RHS, Floc2, eleP, Ndof, Ne, nVar)

					DEALLOCATE ( Floc1, Floc2 )
					DEALLOCATE ( Floc_prec1, Floc_prec2)

				CASE DEFAULT
					PRINT*, " ** ERROR in embedded BC type... (Compute_CG_Explicit.f90)"
					STOP
				END SELECT
				DEALLOCATE(Aloc, Floc)

			END IF 
        END DO

	END IF !fin du cas si embedded
	Sol%RHSTmp = Sol%RHS !voir a quoi sert cette variable

	END SUBROUTINE Compute_LHS_RHS_CG_Explicit
	!=============================================!

END MODULE Compute_CG_Explicit
