MODULE PostProcess

  USE Types
  USE PublicVar
  USE Data_mod
  !USE FluxReconstruction
  !USE Compute_MSEG
  USE GradientReconstruction
  USE Initialisation
  USE PolyReconstruction
  USE SparseMatrix_CG
  USE Distance
  USE Algebra
  USE Data_mod
  USE EmbeddedVariableInitialization
  USE ExactSolutiont
  USE ExactFunctionJumpt

  IMPLICIT NONE

CONTAINS

  !===================================================!
  SUBROUTINE PostProcessSolution(Data, Sol, Sol_prec, Sol_2prec, Mesh, Mat)
  !===================================================!

  !*************************************************************************
  !  Post Precessing after resolution of the sytem
  !  Deals with flux reconstruction depending of type of approximation
  !  i.e. continuous galerking or discontinuous galerking
  !
  !  Parameters:
  !
  !    Input,Output, type (SparseMatrix) Mat, Structure for the global matrix
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input,Output, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input,Output, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(INOUT)  :: Data
    type(SolStructure) , INTENT(INOUT) 	:: Sol
    type(SolStructure) , INTENT(INOUT) 	:: Sol_prec
    type(SolStructure) , INTENT(INOUT) 	:: Sol_2prec
    type(MeshType)     , INTENT(INOUT) 	:: Mesh
    type(SparseMatrix) , INTENT(IN) 	:: Mat
    !-------------------------------------------

    CALL MUMPS2Sol(Data, Mesh, Sol,Sol_prec, Sol%RHS)
    ! Sol%RHS est en faite U

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG-Primal" )
       PRINT*, "   ==> Primal formulation --> Flux reconstruction "
       !separe le flux de la temperature en deux variable Sol%Bcg et Sol%pcg et met à jour si nicolson par changement de variable 
		IF(Data%Space_Type == "EXPLICIT") THEN
			!reconstruction of B 
			CALL GreenGaussReconstruction_Primal(Sol, Mesh, Data) !reconstruction of the flux B 
		ELSE
			PRINT*,"Implicit Resolution - Problem in data input file for Primal formulation"
			STOP 
		END IF
    CASE ( "DG-Primal" )
       PRINT*, "   ==> Primal formulation --> Flux reconstruction "
       !CALL FluxReconstruction_DG(Sol, Mesh, Data)
    CASE ( "CG", "DG" )
       ! Do nothing
    !CASE ( "MSEG", "MSEG-Stab" )
      ! CALL CG2EG(Data, Mesh, Sol)
      ! CALL PostProcessMS(Data, Mesh, Sol)
    CASE DEFAULT
       CALL PrintError("PostProcess")
    END SELECT
    
    !changement de variable comme calcul intermediaire si Nicolson
    IF(DATA%TSCHEME=="NICOLSON" .AND. Data%Space_Type =="IMPLICIT" ) THEN
		Sol%Pcg = 2*Sol%Pcg - Sol_prec%Pcg
		Sol%Bcg = 2*Sol%Bcg - Sol_prec%Bcg
	END IF

	!Mise à jour du temps précédent au temps courant
	CALL MAJSolTmp(Data, Sol, Sol_prec, Sol_2prec, Mesh)!ajouter le tableau de complémentaire 
		
  END SUBROUTINE PostProcessSolution
  !===================================================!
  
  
  !===================================================!
  SUBROUTINE PostProcessSolutionNewton(Data, Sol, Sol_prec, Sol_Nicolson, Mesh, Mat)
  !===================================================!

  !*************************************************************************
  !  Post Precessing after resolution of the sytem
  !  Deals with flux reconstruction depending of type of approximation
  !  i.e. continuous galerking or discontinuous galerking
  !
  !  Parameters:
  !
  !    Input,Output, type (SparseMatrix) Mat, Structure for the global matrix
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input,Output, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input,Output, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(INOUT)  :: Data
    type(SolStructure) , INTENT(INOUT) 	:: Sol,Sol_Nicolson
    type(SolStructure) , INTENT(INOUT) 	:: Sol_prec
    type(MeshType)     , INTENT(INOUT) 	:: Mesh
    type(SparseMatrix) , INTENT(IN) 	:: Mat
    !-------------------------------------------

    CALL MUMPS2Sol(Data, Mesh, Sol,Sol_prec, Sol%RHS)
    ! Sol%RHS est en faite U

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG-Primal" )
       PRINT*, "   ==> Primal formulation --> Flux reconstruction "
       !separe le flux de la temperature en deux variable Sol%Bcg et Sol%pcg et met à jour si nicolson par changement de variable 
		IF(Data%Space_Type == "EXPLICIT") THEN
			!reconstruction of B 
			CALL GreenGaussReconstruction_Primal(Sol, Mesh, Data) !reconstruction of the flux B 
		ELSE
			PRINT*,"Implicit Resolution - Problem in data input file for Primal formulation"
			STOP 
		END IF
    CASE ( "DG-Primal" )
       PRINT*, "   ==> Primal formulation --> Flux reconstruction "
       !CALL FluxReconstruction_DG(Sol, Mesh, Data)
    CASE ( "CG", "DG" )
       ! Do nothing
    !CASE ( "MSEG", "MSEG-Stab" )
      ! CALL CG2EG(Data, Mesh, Sol)
      ! CALL PostProcessMS(Data, Mesh, Sol)
    CASE DEFAULT
       CALL PrintError("PostProcess")
    END SELECT
    
      !changement de variable comme calcul intermediaire si Nicolson
    IF(DATA%TSCHEME=="NICOLSON" .AND. Data%Space_Type =="IMPLICIT" ) THEN
		Sol_Nicolson%Pcg = 2*Sol%pcg - Sol_prec%Pcg
		Sol_Nicolson%Bcg = 2*Sol%Bcg - Sol_prec%Bcg
	END IF
	
	!Mise à jour du temps précédent au temps courant
	CALL MAJSolTmp(Data, Sol, Sol_prec, Sol_prec, Mesh)!ajouter le tableau de complémentaire 
		
  END SUBROUTINE PostProcessSolutionNewton
  !===================================================!
  

  !===================================================!
  SUBROUTINE MUMPS2Sol(Data, Mesh, Sol,Sol_prec, U) ! U correspond à Sol%RHS
  !===================================================!

  !*************************************************************************
  !  Solves the system with MUMPS
  !
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(:)  U, solution vector after resolution
  !    Input,Output, type(SolStructur!Mesh%NarEmb_prec = Mesh%NarEmb !number of nodes at the interface from the previous simulation 
	!Mesh%Tab_NoeudComp_prec = 0 ;  Mesh%Tab_NoeudComp_prec = mesh%Tab_NoeudComp
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure)  , INTENT(IN)    :: Data
    type(MeshType)       , INTENT(IN)    :: Mesh
    type(SolStructure)   , INTENT(INOUT) :: Sol
    type(SolStructure)   , INTENT(IN) 	 :: Sol_prec
    REAL*8, DIMENSION(:) , INTENT(INOUT) :: U
    !-----------------------------------------------------------------------------------
    type(element) 						 :: ele
    !-----------------------------------------------------------------------------------
    INTEGER, DIMENSION(Mesh%Nv)          :: Compt
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: Tab_i
    !-----------------------------------------------------------------------------------
    INTEGER 							 :: Ne, Nv, NDof, NU_i, Ndim, Ned
    INTEGER								 :: ie, i, Pos_i, id, nVar, surNU_i, surPos_i
    INTEGER 							 :: nNodesPerElement,Ns,cpt,j
    !-----------------------------------------------------------------------------------

    IF (Data%Icing .eqv. .TRUE.) THEN
      Ns=Mesh%NB_VertInterf
    ELSE
      Ns=0
    END IF

    Nv = Mesh%Nv ; Ne = Mesh%Ne ; Ned = Mesh%Ned ; Ndim = Data%Ndim
    nVar = Sol%nVar ; cpt=1

    ALLOCATE(Tab_i(Ns,2))

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ("DG")
       NDof = (nDim+1)*mesh%NeActive
       DO ie = 1, Ne
          ele = Mesh%Ele(ie)
          IF ( ele%solve ) THEN
             nNodesPerElement = ele%nNodesPerElement
             DO i = 1, nNodesPerElement
                NU_i = ele%Vertex(i) ; surNU_i = ele%surVertex(i)
                Pos_i = ele%Pos(i) ; surPos_i = ele%surPos(i)
                DO id = 1, nDim
                   Sol%Bdg(Pos_i,id) = U(surPos_i + (id-1)*nDof)
                END DO
                Sol%pdg(Pos_i) = U(surPos_i + nDim*NDof)
             END DO
          END IF
       END DO
    CASE ("DG-Primal")
       DO ie = 1, Ne
          ele = mesh%ele(ie)
          IF ( ele%solve ) THEN
             nNodesPerElement = mesh%ele(ie)%nNodesPerElement
             DO i = 1, nNodesPerElement
                Pos_i = ele%Pos(i) ; surPos_i = ele%surPos(i)
                Sol%pdg(Pos_i) = Sol%pdg(Pos_i) + U(surPos_i)
             END DO
          END IF
       END DO
    CASE ("CG-Primal")

        Sol%Pcg = 0.d0
		DO i = 1, Nv+Ns
			Sol%Pcg(i) =  Sol_prec%Pcg(i) + Data%Dt*U(i)/Sol%DualArea(i)
        END DO
        U = 0.d0 !reset of the vector U for next simulation 
       
    CASE ("CG", "MSEG", "MSEG-Stab")
       Ndof = Nv + Data%MaxRowEmbedded
       Sol%Bcg = 0.d0   ;  Sol%Pcg = 0.d0

       DO i = 1, Nv
          IF ( Mesh%Vertex(i)%Active ) THEN
             DO id = 1, Data%nDim
                Sol%Bcg(i,id) = U(i + (id-1)*NDof)
             END DO
             Sol%Pcg(i) = U(i + nDim*NDof)
          END IF
       END DO

       DO i = Nv+1, Nv+Ns
         DO id = 1, Data%nDim
          Sol%Bcg(i,id) = U(i+(id-1)*NDof)
         END DO
         Sol%Pcg(i) =  U(i+nDim*NDof)
      END DO
      U = 0.d0 !reset of the vector U for next simulation 

    CASE DEFAULT
       CALL PrintError("MUMPS2Sol in PostProcess")
    END SELECT

  END SUBROUTINE MUMPS2Sol
  !=================================================!

  !========================================!
  SUBROUTINE PostProcessMS(Data, Mesh, Sol)
  !========================================!

  !*************************************************************************
  !  Post processing for MS case
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType), INTENT(IN)         :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    !------------------------------------------
    type(element) :: ele
    !---------------------
    INTEGER :: Nv, Ne, nDof, NU_i, nDim
    INTEGER :: ie, i, Pos_i, id
    INTEGER :: nNodesPerElement, nDofPerElement
    !--------------------------------------------

    Nv = Mesh%Nv ; Ne = Mesh%Ne ; Ndim = Data%Ndim
    Sol%Bdg = 0. ; Sol%pdg = 0.
    nDim = Data%nDim

    DO ie = 1, Ne
       ele = Mesh%Ele(ie)
       nNodesPerElement = ele%nNodesPerElement
       nDofPerElement = nNodesPerElement + 1
       DO i = 1, nNodesPerElement
          NU_i = ele%Vertex(i)
          Pos_i = ele%Pos(i)
          DO id = 1, nDim
             Sol%Bdg(Pos_i,id) =  Sol%Bcg(NU_i,id) + Sol%Beg(ie,id)
          END DO
          Sol%pdg(Pos_i) =  Sol%pcg(NU_i) + Sol%peg(ie)
       END DO
    END DO

  END SUBROUTINE PostProcessMS
  !========================================!

  !========================================!
  SUBROUTINE MAJSolTmp(Data, Sol, Sol_prec, Sol_2prec, Mesh)
  !========================================!

  !*************************************************************************
  !  Update the previous time with the currenPostProcessSolutiont time after the resolution
  !  for the next time step
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol_prec, Structure for the solution vector, previous time
  !    Input, type(SolStructure) Sol, Structure for the solution vector, current time
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !GreenGaussReconstruction_Primal
  !*************************************************************************

    IMPLICIT NONE
    type(DataStructure)   , INTENT(INOUT) :: Data 
    type(MeshType)        , INTENT(INOUT) :: Mesh 
    type(SolStructure)    , INTENT(IN)    :: Sol
    type(SolStructure)    , INTENT(INOUT) :: Sol_prec
    type(SolStructure)    , INTENT(INOUT) :: Sol_2prec

    !----------------------------------------
    type(element) :: ele
    !------------------------------------------
    Integer :: ie, i, NU_i
    INTEGER :: Nv, id, nDim, nNodesPerElement,Ns
    !----------------------------------------------

    ! Il faut faire le calcul du nouveau tableau de complementaire
    IF (Data%Icing .eqv. .TRUE.) THEN
      Ns = Mesh%NB_VertInterf
    ELSE
      Ns = 0
    END IF
    
	! enregistrement du nombre de noeud avant deplacement de l'Interface
	Mesh%NarEmb_prec = Mesh%NarEmb !number of nodes at the interface from the previous simulation 
	Mesh%Tab_NoeudComp_prec = 0.d0 ;  Mesh%Tab_NoeudComp_prec = mesh%Tab_NoeudComp

    Nv = Mesh%Nv ; nDim = Data%nDim
	Sol_prec%dist = Sol%dist ! on conserve les informations sur la distance 
	Sol_prec%Nv = Sol%Nv
	
    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG", "CG-Primal" )
      !Mise à jour du temps precedent par le temps courant pour la prochaine itération en temps
      Sol_2prec%Pcg = 0.d0 ; Sol_2prec%Bcg = 0.d0
      
      Sol_2prec%Bcg(1:Ns+Nv,1) = Sol_prec%Bcg(1:Ns+Nv,1)
      Sol_2prec%Bcg(1:Ns+Nv,2) = Sol_prec%Bcg(1:Ns+Nv,2)
      Sol_2prec%Pcg(1:Ns+Nv)   = Sol_prec%Pcg(1:Ns+Nv)
      
      Sol_prec%Pcg = 0.d0 ; Sol_prec%Bcg = 0.d0
      
      Sol_prec%Bcg(1:Ns+Nv,1) = Sol%Bcg(1:Ns+Nv,1)
      Sol_prec%Bcg(1:Ns+Nv,2) = Sol%Bcg(1:Ns+Nv,2)
      Sol_prec%Pcg(1:Ns+Nv)   = Sol%Pcg(1:Ns+Nv)
    CASE DEFAULT
       PRINT*,"Only CG case available"
       CALL PrintError("Error in MAJSolTmp  ")
    END SELECT
    
    
  END SUBROUTINE MAJSolTmp
  !========================================!

END MODULE PostProcess
