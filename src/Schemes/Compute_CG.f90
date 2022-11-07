MODULE Compute_CG

  USE Types
  USE PublicVar
  USE ElementIntegral_Mixed
  USE SparseMatrix_CG
  USE BCSelection
  USE DirichletConformal_Mixed
  USE NeumannConformal_Mixed
  USE DirichletEmbedded_Mixed
  USE NeumannEmbedded_Mixed
  USE ReadMesh_mod
  USE SparseMatrix_DG
  USE EdgeIntegral_DG
  USE NeumannJumpConformal_Mixte
  USE NeumannIcingConformal_Mixed
  USE DirichletIcingConformal_Mixed
  USE EdgeIntegral_CG_Icing
  USE ElementIntegral_Mixed_Icing
  USE Affichage_M
  USE ElementNewtonIntegral_Mixed_Icing
  
  IMPLICIT NONE

CONTAINS

	!==================================================!
	SUBROUTINE Compute_LHS_RHS_CG(Data, Mesh, Sol, Sol_prec, Sol_2prec, Mat, PhysicalInterface )
	!==================================================!

	IMPLICIT NONE
    type(DataStructure), INTENT(INOUT) :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SolStructure) , INTENT(IN)    :: Sol_prec, Sol_2prec 
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    type(VectorInterface), INTENT(IN)  :: PhysicalInterface 

    !-------------------------------------------
    type(element) :: ele, eleM, eleP, eleBC, Adj
    type(Edge)    :: Edg
    !--------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: Aloc, Amm, Amp, Apm, App
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Floc, Floc1, Floc2, Floc_prec, Floc_prec2nd, Floc_2prec
    REAL*8, DIMENSION(:)  , ALLOCATABLE :: Floc_prec1, Floc_prec2

    !--------------------------------------------------
    CHARACTER(len=20) :: BCType
    !--------------------------
    INTEGER :: AdjNum, nDim, NeActive, iLS
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
    Ne = Mesh%Ne ; Nv = Mesh%NvActive + Ns ; Nv_bis = Mesh%NvActive + Data%MaxRowEmbedded
    Nb = Mesh%Nb ; nVar = Mat%NVar ; Ndof =  Mesh%NvActive + Data%MaxRowEmbedded
    Mat%Coeff = 0.d0 ; Sol%RHS = 0.d0
    NeActive = Mesh%NeActive

    !boucle sur les elements pour les integrales sur Omega
    DO ie = 1, Ne
    
		ele = Mesh%ele(ie)
		nNodesPerElement = ele%nNodesPerElement
		ALLOCATE ( Aloc(nVar*nNodesPerElement,nVar*nNodesPerElement) ) !allocation de la matrice locale sur l'element
		ALLOCATE ( Floc(nVar*nNodesPerElement),  Floc_prec(nVar*nNodesPerElement),  Floc_2prec(nVar*nNodesPerElement) )
			
		Aloc = 0.d0 ; Floc = 0.d0 ; Floc_prec = 0.d0

		CALL Remplissage_FlocPrec(ele, Floc_prec, Sol_prec, nDim, Data) ! for the derivative in time on T 
		CALL Remplissage_FlocPrec(ele, Floc_2prec, Sol_2prec, nDim, Data) ! for the derivative in time on T 

		IF ( ele%solve ) THEN
			IF( .NOT. Data%Icing ) THEN
				CALL Compute_ElementIntegral_Mixed(ele, Aloc, Floc, Data) !pas encore modifier pour prendre en compte Floc_prec
			ELSE	
				CALL Compute_ElementIntegralIcing_Mixed(ele, Aloc, Floc, Floc_prec,Floc_2prec, Data)
			END IF
          
			CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv_bis, Ne)
			CALL AssembleRHSCG(Sol%RHS, Floc, ele, Ndof, Ne, nVar)

		END IF
		
		!we check if the element is in the solid airfoil 
		IF (Data%Resistance) THEN 
			IF(ele%ref == IcingZone3_ref) THEN 
				!need to identify if an element has an edge which is also a part of the resistance 
				!then we need a loop on the edges of the element 
				! then to identify which edge
				! then to do the correct addition to the matrix 
				!edg = ?? ! need to define wich edge we send 
				
				!need to modify this function to take into account the new condition 
				CALL SelectBC2(edg, BCType, Data, iBC) ! Resort le type de condition définit sur l'arete
				CALL GetFace(ele, edg, iFace) ! need to see if i have to do modification to here 
				CALL NeumannHeatFluxJump(Aloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
				
				CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv_bis, Ne)
				CALL AssembleRHSCG(Sol%RHS, Floc, ele, Ndof, Ne, nVar)
				
			END IF 
			!then add the part to add the condition on the flux for the term defining the resistance 
		
		END IF 
		DEALLOCATE ( Aloc, Floc, Floc_prec, Floc_2prec)
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
			Aloc = 0.d0 ; Floc = 0.d0

			CALL SelectBC2(edg, BCType, Data, iBC) ! Resort le type de condition définit sur l'arete
			! iBC est la reference de la condition associé à l'arête 
			CALL GetFace(ele, edg, iFace) !retourne le numero de l'arete dans Iface associé à la numerotation de l'arte dans la lsite des aretes de bord de l'element

			CALL Remplissage_FlocPrec(ele,Floc_prec, Sol_prec,nDim,Data)

			SELECT CASE ( TRIM(ADJUSTL(BCType)) )

			CASE ( "Dirichlet", "dirichlet" )

				IF( .NOT. Data%Icing ) THEN
					CALL DirichletConf_Mixed(Aloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
				ELSE
					CALL DirichletIcingConf_Mixed(Aloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
				END IF
				

				CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv_bis, Ne)
				CALL AssembleRHSCG(Sol%RHS, Floc, ele, Ndof, Ne, nVar)

			CASE ( "Neumann", "neumann" )
				IF ( .NOT. Data%Icing ) THEN
					CALL NeumannConf_Mixed(Aloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
				ELSE
					CALL NeumannConf_Mixed(Aloc, Floc, ele, mesh%edg(ele%edg), Data, iFace, iBC)
				END IF
				

				CALL AssembleSparseMatrixCG(Mat, Aloc, ele, Nv_bis, Ne)
				CALL AssembleRHSCG(Sol%RHS, Floc, ele, Ndof, Ne, nVar)

			CASE ( "NeumannJump","neumannJump","neumannjump","Neumannjump" )
				eleM = mesh%ele(edg%tri(2)) ; eleP = mesh%ele(edg%tri(1))
				nNodesPerElement = eleM%nNodesPerElement
				ALLOCATE ( Amm(nVar*nNodesPerElement,nVar*nNodesPerElement) )
				ALLOCATE ( Amp(nVar*nNodesPerElement,nVar*nNodesPerElement) )
				ALLOCATE ( Apm(nVar*nNodesPerElement,nVar*nNodesPerElement) )
				ALLOCATE ( App(nVar*nNodesPerElement,nVar*nNodesPerElement) )
				
				Amm=0.d0 ; Amp=0.d0 ; Apm=0.d0 ; App=0.d0
				ALLOCATE ( Floc1(nVar*nNodesPerElement), Floc2(nVar*nNodesPerElement) )
				Floc1 = 0.d0 ; Floc2=0.d0

				!pas de condition de saut ailleurs qu'à l'interface pour le moment
				CALL ComputePMMatrices_CG(Mesh, iBC, Edg,  Mesh%Edg(eleM%Edg), Mesh%Edg(eleP%Edg), &
				EleM, eleP, Amm, Amp, Apm, App, Floc1, Floc2, Data)
				
				CALL AssembleSparseMatrixInterface(Mat, Amm, Amp, Apm, App, EleM, EleP, Nv_bis, NeActive)
				CALL AssembleRHSCG(Sol%RHS, Floc1, eleM, Ndof, Ne, nVar)
				CALL AssembleRHSCG(Sol%RHS, Floc2, eleP, Ndof, Ne, nVar)

				DEALLOCATE ( Amm, Amp, Apm, App, Floc1, Floc2 )
			CASE DEFAULT
				IF(Data%Icing .AND. Data%Embedded) THEN
					IF(mesh%boundary(ib)%ref == Data%N_BC+ 1) THEN
					!do nothing
					!permet de laisser les aretes de l'interface dans le fichier mesh
					END IF
				ELSE
				   PRINT*, ie, TRIM(ADJUSTL(BCType))
				   PRINT*, " ** ERROR in BC type... (Compute_CG.f90)"
				   PRINT*, "    Available :"
				   PRINT*, "      + Dirichlet/dirichlet"
				   PRINT*, "      + Neumann/neumann"
				   PRINT*, "      + NeumannJump/neumannJump/neumannjump/Neumannjump"
				   STOP
				END IF
			END SELECT
			DEALLOCATE ( Aloc, Floc,Floc_prec )
		END IF
    END DO

    ! Embedded BC
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
				CASE ( "Dirichlet","dirichlet" )
					CALL DirichletEmb_Mixed(ALoc, Floc, ele, Mesh%edg(ele%edg), &
					Sol%dist(1,ele%Vertex,:), data, iFace, mesh%embBoundary(ib)%embLS)
				CASE ( "Neumann","neumann" )
					CALL NeumannEmb_Mixed(Aloc, Floc, ele, mesh%edg(ele%edg), &
					Sol%dist(1,ele%vertex,:), data, iFace, mesh%embBoundary(ib)%embLS)
				CASE ( "NeumannJump","neumannJump","neumannjump","Neumannjump" )
					!on recupere les deux elements associe a l'arete courante
					IF (mesh%ele(edg%tri(2))%ref == IcingZone2_ref) THEN !we define eleP as the left element at the surrogate 
						eleP = mesh%ele(edg%tri(1))
						eleM = mesh%ele(edg%tri(2))
					ELSE
						eleP = mesh%ele(edg%tri(2))
						eleM = mesh%ele(edg%tri(1))
					END IF

					nNodesPerElement = eleM%nNodesPerElement
					ALLOCATE ( Amm(nVar*nNodesPerElement,nVar*nNodesPerElement) )
					ALLOCATE ( Amp(nVar*nNodesPerElement,nVar*nNodesPerElement) )
					ALLOCATE ( Apm(nVar*nNodesPerElement,nVar*nNodesPerElement) )
					ALLOCATE ( App(nVar*nNodesPerElement,nVar*nNodesPerElement) )
					ALLOCATE ( Floc1(nVar*nNodesPerElement), Floc2(nVar*nNodesPerElement) )

					IF (Data%SetEquation == 1) THEN 
					CALL ComputePMMatrices_CG_Embedded(PhysicalInterface, Mesh, Edg,  Mesh%Edg(eleM%Edg), Mesh%Edg(eleP%Edg), &
							EleM, EleP, Amm, Amp, Apm, App, Floc1, Floc2, Data, iLS, Sol_prec, Sol)
					ELSE IF (Data%SetEquation == 2) THEN 
					CALL ComputePMMatrices_CG_Embedded2(PhysicalInterface, Mesh, Edg,  Mesh%Edg(eleM%Edg), Mesh%Edg(eleP%Edg), &
							EleM, EleP, Amm, Amp, Apm, App, Floc1, Floc2, Data, iLS, Sol_prec, Sol)
					ELSE
						PRINT*,"Error in the choice of the interface conditions"
						STOP 
					END IF 				
					CALL AssembleSparseMatrixInterface(Mat, Amm, Amp, Apm, App, EleM, EleP, Nv_bis, NeActive)
					CALL AssembleRHSCG(Sol%RHS, Floc1, eleM, Ndof, Ne, nVar)
					CALL AssembleRHSCG(Sol%RHS, Floc2, eleP, Ndof, Ne, nVar)

					DEALLOCATE ( Amm, Amp, Apm, App, Floc1, Floc2 )

				CASE DEFAULT
					PRINT*, " ** ERROR in embedded BC type... (Compute_CG.f90)"
					PRINT*, " Available :"
					PRINT*, "   + NeumannJump/neumannJump/neumannjump/Neumannjump"
					STOP
				END SELECT
				DEALLOCATE(Aloc, Floc)

			END IF 
       END DO
	END IF !fin du cas si embedded

	END SUBROUTINE Compute_LHS_RHS_CG
	!=============================================!

	!==================================================!
	SUBROUTINE Remplissage_FlocPrec(ele,Floc_prec,Sol_prec,nDim,Data)
    !==================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)    :: ele
    type(SolStructure)    , INTENT(IN)    :: Sol_prec
    REAL*8, DIMENSION(:)  , INTENT(INOUT) :: Floc_prec
    INTEGER               , INTENT(IN)    :: nDim
    type(DataStructure)   , INTENT(IN) 	  :: Data
    !--------------------------------------------------
    INTEGER :: i, nNodesPerElement
    !--------------------------------------------------

    nNodesPerElement = ele%nNodesPerElement
    IF (Data%Space_Scheme=="CG") THEN !  1+nDim variable 
		DO i = 1 ,nNodesPerElement
			Floc_prec(i)   = Sol_prec%Bcg(ele%SurVertex(i),1)
			Floc_prec(i+(nDim-1)*nNodesPerElement) = Sol_prec%Bcg(ele%SurVertex(i),2)
			Floc_prec(i+nDim*nNodesPerElement) = Sol_prec%Pcg(ele%SurVertex(i))
		END DO
	ELSE IF (Data%Space_Scheme=="CG-Primal") THEN 
		DO i = 1 ,nNodesPerElement
			Floc_prec(i) = Sol_prec%Pcg(ele%SurVertex(i)) !only one variable 
		END DO
	ELSE 
		PRINT*,"Error in Remplissage_Floc_Prec"
		STOP
	
	END IF 

	END SUBROUTINE Remplissage_FlocPrec
	!=============================================!
  
   !==================================================!
   SUBROUTINE Remplissage_FlocNewtown(Data,ele,Floc_inter,Sol,nDim)
   !==================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)    :: ele
    type(SolStructure)    , INTENT(IN)    :: Sol
    REAL*8, DIMENSION(:)  , INTENT(INOUT) :: Floc_inter
    INTEGER               , INTENT(IN)    :: nDim
    type(DataStructure)   , INTENT(IN) 	  :: Data

    !--------------------------------------------------
    INTEGER :: i, nNodesPerElement
    !--------------------------------------------------
    
	
    nNodesPerElement = ele%nNodesPerElement

    DO i = 1 ,nNodesPerElement
		Floc_inter(i)   						= Sol%Bcg_Newton(ele%SurVertex(i),1)
		Floc_inter(i+(nDim-1)*nNodesPerElement) = Sol%Bcg_Newton(ele%SurVertex(i),2)
		Floc_inter(i+nDim*nNodesPerElement)		= 0 ! the value of T isn't used so this part is to keep the same construction than for the other local vector at time n 
    END DO


	END SUBROUTINE Remplissage_FlocNewtown
	!=============================================!


END MODULE Compute_CG
