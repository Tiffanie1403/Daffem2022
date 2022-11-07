MODULE Initialisation

  USE Types
  USE PublicVar
  USE ReadMesh_Mod
  USE SparseMatrix_CG
  USE SparseMatrix_DG
  USE Distance
  USE Algebra
  USE Data_mod
  USE EmbeddedVariableInitialization
  USE ExactSolutiont
  USE ExactFunctionJumpt
  USE PolyReconstruction
  USE IcingTest_mod 
  USE Velocity

  IMPLICIT NONE

CONTAINS

  !========================================================!
  SUBROUTINE Simulation_Tinit(Sol_prec, Sol_2prec, Sol, Mesh, Data, Mat, PhysicalInterface)
  !========================================================!

    IMPLICIT NONE
    
    type(DataStructure)  , INTENT(INOUT)  :: Data
    type(SolStructure)   , INTENT(OUT)    :: Sol_prec, Sol, Sol_2prec 
    type(MeshType)       , INTENT(OUT)    :: Mesh
    type(SparseMatrix)   , INTENT(INOUT)  :: Mat
    type(VectorInterface), INTENT(OUT)    :: PhysicalInterface 
    !-----------------------------------------------------------------
    INTEGER 							  :: i, ie, Nu_i
    REAL*8							  	  :: x0, eps, distNode
    !-----------------------------------------------------------------

    
    !Definition of the value Icing_chi 
	IF(icas==001 ) THEN 
		!Definition of the variable chi used in the definition of the interface position
		x0 = 1E-5! initial guess
		eps= 1E-15
		CALL Newton_Interface(x0,eps)
		Icing_chi = x0
		!PRINT*,"Value of chi",Icing_chi 
	END IF 
	    
	!specificities of the 001 case 
	IF(icas==001) THEN 
		!association of a initial position to an initial time 
		!Value of the initial time knowing the initial position x 
		!IF(.NOT. Data%DefInterfaceNode) THEN 
		Data%T = Icing_L1**2/(4*Icing_chi**2*Icing_alpha_l)  ! T_init 
		!ELSE 
		!	Data%T     = 0.d0
		!	PRINT*,"We start the simulation with one phase only" 
		!END IF 
		Data%Tinit = Data%T 
		
		IF(Data%T==0.d0 .AND. .NOT. Data%DefInterfaceNode) THEN 
			PRINT*,"Warning, the definition of Dt needs to be update T is zero"
			STOP
		ELSE IF (DATA%TSCHEME=="NICOLSON" .OR. DATA%TSCHEME == "BDF2") THEN 
			!on choisi la valeur du T initial comme valeur de t
			Data%Dt = Data%h*sqrt(Data%T)/(Icing_chi*sqrt(Icing_alpha_l))*0.5d0
			PRINT*,"Value of Dt", Data%dt 	
		ELSE IF (DATA%TSCHEME=="EULER") THEN 
			!h**2 to keep a second order accuracy
			!Data%Dt = Data%h**2*sqrt(Data%T)/(2*Icing_chi*sqrt(Icing_alpha_l))
			Data%Dt = Data%h*sqrt(Data%T)/(Icing_chi*sqrt(Icing_alpha_l))*0.5d0
		END IF 
		Data%Tmax = 10.2741*0.5d0 !10*Data%dt + Data%T !initialization of the final time to get 10 time steps !10.2741
	ELSE
		!definition of the initial time and final time 
		!Data%T     = Data%IntitialTime
		Data%Tinit = Data%T 
		Data%Tmax  = Data%FinalTime !10*Data%dt + Data%T !initialization of the final time to get 10 time steps !10.2741
		! defined in GetData Data%Dt 
	END IF
		
	!on ne peut pas definir de tour max si dt n'est pas constant 
	Data%tourMax = aint((Data%Tmax-Data%T)/Data%Dt) !on defnit le nombre d'itérations en temps
	PRINT*,"Data%tourMax : ",Data%tourMax-1,Data%Dt

	IF(Data%tourMax==0) THEN 
		PRINT*,"Value of Tmax not big enough compared to the value of DT"
		STOP
	END IF 
	
	! on va modifier Data%Tmax pour correspondre au temp final actuel de la resolution 
	Data%Tmax = Data%tourMax*Data%Dt !to prevent any lost of precision from the definition of tourMax 
    PRINT*,"Value of Tmax", Data%Tmax
    Data%tourMax = Data%tourMax - 1 ! we lost one loop to initialize two values !!voir cette ligne !!!! 
    ! Read the mesh
    CALL ReadMesh(Data, Mesh)
 
    IF (Data%Icing) THEN
      RefZoneAChanger = IcingZone2_ref
    END IF

    Sol_prec%Nv  = Mesh%Nv  ; Sol_prec%Ne = Mesh%Ne
    Sol_prec%Ned = Mesh%Ned ; Sol_prec%Nb = Mesh%Nb    
    
    ! This is modified later if embedded
    Mesh%NvActive = Mesh%Nv           ; Mesh%NeActive = Mesh%Ne
    Sol_prec%NeActive = Mesh%NeActive ; Sol_prec%NvActive = Mesh%NvActive

    ! initialisation de variable sur la structure solution courant
    Sol%Nv = Mesh%Nv ; Sol%Ne = Mesh%Ne             ; Sol%Ned = Mesh%Ned
    Sol%Nb = Mesh%Nb ; Sol%NeActive = Mesh%NeActive ; Sol%NvActive = Mesh%NvActive

    ! Define boundary conditions
    CALL SetBC(Mesh)  
    CALL AllocEmbeddedVariables_Tinit(Data,Sol,Mesh) 
    CALL AllocEmbeddedVariables_Tinit(Data,Sol_prec,Mesh)
    !initialisation de Sol_2prec par copie de Sol_prec 
    Sol_2prec = Sol_prec ! if needed by the choice of the time scheme, but at least the structure exists 
    
    ! initialisation avant modification dans la routine de tag
    DO ie = 1, mesh%Ne
        DO i = 1, mesh%ele(ie)%nNodesPerElement
			Mesh%ele(ie)%surVertex(i)       = Mesh%ele(ie)%vertex(i)
			Mesh%ele(ie)%surPos(i)          = Mesh%ele(ie)%Pos(i)
			Mesh%ele(ie)%SurVertex_Sauvg(i) = 0
        END DO
    END DO
	
	Data%Tour = 0 
	
	IF(Data%DefInterfaceNode) THEN !definition of the interface by a set of node 
		!Icing_L3 height of the domaine, means only take into account horizontal displacement 
		! initialization of the vector 
		ALLOCATE(PhysicalInterface%Coor(Icing_NodePhysicalInterface,2))       ; PhysicalInterface%Coor(:,:)       = 0.d0 
		ALLOCATE(PhysicalInterface%Coor_prec(Icing_NodePhysicalInterface,2))  ; PhysicalInterface%Coor_prec(:,:)  = 0.d0 
		ALLOCATE(PhysicalInterface%Coor_2prec(Icing_NodePhysicalInterface,2)) ; PhysicalInterface%Coor_2prec(:,:) = 0.d0 
		
		distNode = Icing_L3/(Icing_NodePhysicalInterface-1) ! on definit les point de maniere equidistant sur l'axe des ordonnées 
		
		IF (icas==001) THEN 
			PhysicalInterface%Coor(:,1) = 2.d0*Icing_chi*sqrt(Icing_alpha_l*Data%T) ! x coordinate
		ELSE 
			PhysicalInterface%Coor(:,1) = Icing_L1 ! starting with one phase so it should work 
		END IF 
		DO i = 1, Icing_NodePhysicalInterface
			PhysicalInterface%Coor(i,2) = (i-1)*distNode  ! y coordinate 
		END DO 
	ELSE 
		IF (icas==001) THEN 
			!modification of the interface position for the initial time to get the initial position with the correct configuration 
			CALL Interface_deplacement(Data)
		END IF 
	END IF 
    
    ! If embedded, define geometry and tags
    IF ( Data%Embedded .OR. Data%Immersed ) THEN
		CALL ComputeDistance(PhysicalInterface, Data, Mesh, Sol_2prec)
		IF ( Data%Embedded .AND. .NOT. Data%Icing ) CALL tagElement2Solve(Data, Mesh)
		IF ( Data%Immersed .OR. Data%Icing ) CALL tagElement2SolveImmersed(Data, Mesh)
		Sol_2prec%NeActive = Mesh%NeActive ; Sol_2prec%NvActive = Mesh%NvActive
		Mesh%NarEmb= mesh%NSurB !nombre d'arete sur la surrogate
    END IF
	! copy for the next loop
	ALLOCATE (Mesh%Tab_NoeudComp_prec(Data%MaxRowEmbedded,2))
	Mesh%Tab_NoeudComp_prec = 0.d0 
	Mesh%Tab_NoeudComp_prec = Mesh%Tab_NoeudComp  
	
	! Enregistrement des surVertex
	IF (Data%Icing) THEN
	    DO ie = 1, mesh%Ne
			DO i = 1, mesh%ele(ie)%nNodesPerElement
				Nu_i=mesh%ele(ie)%vertex(i)
				IF (mesh%vertex(Nu_i)%OnInterface .AND. .NOT. Data%Embedded) THEN !cas conforme
					IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
						mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
						mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
					END IF
				ELSE IF (mesh%vertex(Nu_i)%OnSurrogate .AND.  Data%Embedded) THEN ! cas embedded
					IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
						mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
						mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
					END IF
			    END IF    
		    END DO
	    END DO
	END IF
    mesh%vertex(:)%active = .TRUE.

	CALL AllocSol_Init(Data, Sol_prec, Sol_2prec, Mesh, 1) ;  CALL InitSol0(Mesh, Sol_2prec, Data)
	
	IF(Data%TSCHEME=="BDF2") THEN ! we also need another vector in the initialization procedure
	
		Data%T = Data%T+Data%Dt ! Modification of the current time
		PRINT*,"initialisation of the n vector at time", Data%T  
		IF(Data%DefInterfaceNode) THEN !definition of the interface by a set of node 
			!Icing_L3 height of the domaine, means only take into account horizontal displacement 		
			PhysicalInterface%Coor_2prec = PhysicalInterface%Coor 
			!part to change when the points are not following the same behaviour 
			IF (icas==001) THEN 
				PhysicalInterface%Coor(:,1) = 0.d0 ! only the x coordinate is put on hold 
				PhysicalInterface%Coor(:,1) = 2.d0*Icing_chi*sqrt(Icing_alpha_l*Data%T) ! x coordinate
			ELSE IF(Data%ResisTance) THEN 
				PRINT*,"Need to put something here" ! I need to update every point separatly 
				STOP 
			ELSE 
				DO i = 1, Icing_NodePhysicalInterface
					PhysicalInterface%Coor(i,1) = PhysicalInterface%Coor(i,1) + Data%Speed*Data%Dt 	!only a displacement of the x coordinate 	
				END DO
			END IF
		ELSE 
			!modification of the interface position - same function for any case 
			Icing_L1_prec = Icing_L1 ! keep the value of the interface at time n-1 
			CALL Interface_deplacement(Data) ! deplacement de la position de la véritable interface
		END IF 
		
		!Remise à O
		!Every data concerning the previous solution on the mesh for the complementary nodes 
		Mesh%Tab_NoeudComp_prec = 0 ;  Mesh%Tab_NoeudComp_prec = Mesh%Tab_NoeudComp !voir cette ligne !!!!
		Mesh%NarEmb_prec        = Mesh%NarEmb
		Mesh%Tab_NoeudComp      = 0 ; Mesh%NarEmb = 0.d0
		
		! remise à 0 aprés le calcul fait dans la boucle précédente
		DO i = 1, mesh%Nv
		   mesh%vertex(i)%Complementaire = 0
		END DO

		DO ie = 1, mesh%Ne
			DO i = 1, mesh%ele(ie)%nNodesPerElement
				Mesh%ele(ie)%surVertex(i) = Mesh%ele(ie)%vertex(i)
				Mesh%ele(ie)%surPos(i) = Mesh%ele(ie)%Pos(i)
				Mesh%ele(ie)%SurVertex_Sauvg(i) = 0
			END DO
		END DO

		Data%Tour=1
		! If embedded, define geometry and tags
		! plus define number of active elements
		IF ( Data%Embedded .OR. Data%Immersed ) THEN
			CALL computeDistance(PhysicalInterface, Data, Mesh, Sol_prec)
			IF ( Data%Embedded .AND. .NOT. Data%Icing ) CALL tagElement2Solve(Data, Mesh)
			IF ( Data%Immersed .OR. Data%Icing ) CALL tagElement2SolveImmersed(Data, Mesh)
			Sol_prec%NeActive = Mesh%NeActive ; Sol_prec%NvActive = Mesh%NvActive
			Mesh%NarEmb= mesh%NSurB !nombre d'arete sur la surrogate
		END IF
		Data%Tour = 0

		IF (Data%Icing .eqv. .TRUE.) THEN
			! Enregistrement des surVertex
			DO ie = 1, mesh%Ne
				DO i = 1, mesh%ele(ie)%nNodesPerElement
					Nu_i=mesh%ele(ie)%vertex(i)
					IF (mesh%vertex(Nu_i)%OnInterface .AND. .NOT. Data%Embedded) THEN
						IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
						   mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
						   mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
						END IF
					ELSE IF (mesh%vertex(Nu_i)%OnSurrogate .AND.  Data%Embedded) THEN
						IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
						   mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
						   mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
						END IF
					END IF
				END DO
			END DO
		END IF
		
		mesh%vertex(:)%active = .TRUE.
		CALL InitSol0(Mesh, Sol_prec, Data)
		
		PhysicalInterface%Coor_prec = PhysicalInterface%Coor 
		 
		Sol_2prec%Dist = Sol_prec%Dist !the vector has been modified for the new time so the distance vector needs to be updated too 
		CALL AllocSol_PrecModif(Data, Sol_2prec, Mesh) !maj du vecteur temps precedent
		
		Mesh%Tab_NoeudComp_prec = 0 ;  Mesh%Tab_NoeudComp_prec = Mesh%Tab_NoeudComp !voir cette ligne !!!!
		Mesh%NarEmb_prec = Mesh%NarEmb
	ELSE
		Sol_prec = Sol_2prec ! because we do not use a second order scheme 
		PhysicalInterface%Coor_prec = PhysicalInterface%Coor 
		Mesh%Tab_NoeudComp_prec = 0 ;  Mesh%Tab_NoeudComp_prec = Mesh%Tab_NoeudComp
		Mesh%NarEmb_prec = Mesh%NarEmb
	END IF 
		
	!Initialisation Sol 
    CALL AllocSol(Data, Sol, Mat, Mesh) ; CALL InitVar(Data, Sol, Mat) 
    
    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    
    CASE ( 'CG-Primal', 'CG', 'MSEG', 'MSEG-Stab' )
        CALL FillIndCG(Mat, Mesh)
    CASE ( 'DG', 'DG-Primal' )
        CALL FillIndDG(Mat, Mesh)
    CASE DEFAULT
        CALL PrintError("SimulationInitialisation")
        
    END SELECT
            
  END SUBROUTINE Simulation_Tinit
  !========================================================!
  
  !========================================================!
  SUBROUTINE SimulationInitialisation(PhysicalInterface, Data, Sol, Sol_prec, Sol_2prec,&
   Mesh, Mat, SpeedJump, JumpF, SpeedJump2, JumpF2)
    !========================================================!

    IMPLICIT NONE
    type(DataStructure)  , INTENT(INOUT) 		:: Data
    type(VectorInterface), INTENT(INOUT)        :: PhysicalInterface 
    type(SolStructure)   , INTENT(INOUT)		:: Sol, Sol_prec, Sol_2prec
    type(MeshType)       , INTENT(INOUT) 		:: Mesh
    type(SparseMatrix)   , INTENT(INOUT) 		:: Mat
    REAL*8			     , INTENT(IN)    		:: SpeedJump, JumpF, SpeedJump2, JumpF2  !value of the velocity and jump to impose ro move the interface 
    !-----------------------------------------------------------------------------------------------------------------
    INTEGER          				   		    :: i, ie, Nu_i, Ns, k ,nb, nbMaxITER
    REAL*8           				   		    :: decel_coef, JumpF_prev, Gam,v_prec,v2_prec,COMP_VAL,truePos
    REAL*8 							  		    :: JumpFx,JumpFy,JumpNewGEo,TolSousIter
    LOGICAL        					  		    :: test, Ttest
    REAL*8 							   		    :: TestVelocity,ModifiedJump, MeanT, L1_inter
    REAL*8									    :: SpeedJumpx, SpeedJumpy
    !-----------------------------------------------------------------------------------------------------------------
    REAL*8   , DIMENSION(:,:), ALLOCATABLE 	    :: SolJump_Bcg
    REAL*8   , DIMENSION(:)  , ALLOCATABLE 	    :: SolJump_Pcg,X_interface
    !-----------------------------------------------------------------------------------------------------------------
	
	
	Data%ErrorT  = 0.d0  ;  Data%ErrorBx = 0.d0  ;  Data%ErrorBy = 0.d0
	ALLOCATE(X_interface(2)) ; X_interface = 0.d0
	!The modification has to be done on the update of the interface 
	!IF(Data%DefInterfaceNode) THEN !definition of the interface by a set of node 
		IF(icas==001) THEN
			! i need to define a function to get the correct value to move every point with the good speed 
			IF(Data%InterfaceScheme == "EulerExp") THEN 
				DO i = 1 , Icing_NodePhysicalInterface
					!PRINT*,"iter",i,"/",icing_NodePhysicalInterface
					!SpeedJump need to be a vector of dimension 2
					!SpeedJumpx = 0.d0 ; SpeedJumpy = 0.d0 ; JumpFx = 0.d0 ; JumpFy = 0.d0
					PhysicalInterface%Coor_2prec(i,:) = PhysicalInterface%Coor_prec(i,:)
					PhysicalInterface%Coor_prec(i,:)  = PhysicalInterface%Coor(i,:)
					!CALL LocalMeanFluxReconstruct(Data, Mesh, Sol_prec, PhysicalInterface%Coor(i,:), SpeedJumpx, SpeedJumpy, JumpFx, JumpFy)
					PhysicalInterface%Coor(i,1) = PhysicalInterface%Coor(i,1) + SpeedJump*Data%Dt !2.d0*icing_chi*sqrt(icing_alpha_l*Data%T) !PhysicalInterface%Coor(i,1) + SpeedJumpx*Data%Dt 
					!PhysicalInterface%Coor(i,2) = PhysicalInterface%Coor(i,2) !+ SpeedJumpy*Data%Dt 
				END DO  
				!Data%LS_L(:,1) = 2*PhysicalInterface%Coor(1,1)
			!ELSE IF (Data%InterfaceScheme == "BDF2Ext") THEN 
			!STOP 
			!see to implement the bdf2 scheme for this part 
			!	L1_inter  = Icing_L1 !intermediate step to init icing_L1_prec 		
			!	Icing_L1 = 2.d0/3.d0*Data%Dt*(2.d0*SpeedJump - SpeedJump2)+ 4.d0/3.d0*Icing_L1 - 1.d0/3.d0*Icing_L1_prec ! a changer avec le BDF2		
			!	Data%LS_L(:,1) = 2.d0*Icing_L1
			!	WRITE(111,*)Data%T,Icing_L1,2.d0*Icing_chi*sqrt(Icing_alpha_l*Data%T)
			!	Icing_L1_prec = L1_inter 
			!ELSE 
			!	PRINT*,"Error in the choice of the interface SCHEME"
			!	STOP 
			END IF 	
		END IF 
	!	ELSE IF(Data%SpeedCst) THEN 
	!		!independent of the system constant velocity 
	!		!for this case only the x coordinate is changed 
	!		PhysicalInterface%Coor(:,1) = PhysicalInterface%Coor(:,1) + Data%Speed*Data%Dt 	!only a displacement of the x coordinate 	
	!	ELSE ! we need to reconstruct a value for each point defininf the true interface
	!!		IF(Data%Resistance) THEN 
	!			IF(Data%InterfaceScheme == "EulerExp") THEN 
	!				DO i =2, Icing_NodePhysicalInterface-1 ! do not touch the points at the extremities of the domainto be sure that we always have one connected to the exterior boundary 
	!				! maybe see if i can still move following the x coordinate but not the y one 
	!					X_interface(1) = PhysicalInterface%Coor(i,1)
	!					X_interface(2) = PhysicalInterface%Coor(i,2)
	!					IF (Data%TOUR == 1) THEN !not displacement 
	!						SpeedJumpx = 0.d0 ; SpeedJumpy = 0.d0 
	!					ELSE 
	!						CALL MeanFluxReconstructNodeVal(Data, Mesh, Sol_prec, SpeedJumpx, SpeedJumpy, X_interface) !NO NEED OF THE NORMAL HERE 
	!						IF (SpeedJumpx==0) THEN 
	!							 PRINT*,"There is maybe a problem, need to check"
	!						END IF 
	!					END IF 
	!					PhysicalInterface%Coor(i,1) = PhysicalInterface%Coor(i,1) + SpeedJumpx*Data%Dt
	!					PhysicalInterface%Coor(i,2) = PhysicalInterface%Coor(i,2) + SpeedJumpy*Data%Dt
	!				END DO 
				!ELSE 
				!	PRINT*,"Not Implemented yet"
				!	STOP
	!			END IF 
			!ELSE 
			!	PRINT*,"What am I doing here, pb in the choice of the test case"
			!	END IF 
	!		END IF 
	!	END IF 
!	ELSE 
!		IF(icas==001) THEN
!			!True simulation 
!			!Displacement of the front after the different itteration with caculated solution 
!			IF(Data%InterfaceScheme == "EulerExp") THEN 
!				Icing_L1 = Icing_L1 + SpeedJump*Data%Dt ! a changer avec le BDF2
!				Data%LS_L(:,1) = 2.d0*Icing_L1	
!			ELSE IF (Data%InterfaceScheme == "BDF2Ext") THEN 
!				L1_inter  = Icing_L1 !intermediate step to init icing_L1_prec 
!				
!				Icing_L1 = 2.d0/3.d0*Data%Dt*(2.d0*SpeedJump - SpeedJump2)+ 4.d0/3.d0*Icing_L1 - 1.d0/3.d0*Icing_L1_prec ! a changer avec le BDF2		
!				Data%LS_L(:,1) = 2.d0*Icing_L1
!				WRITE(111,*)Data%T,Icing_L1,2.d0*Icing_chi*sqrt(Icing_alpha_l*Data%T)
!				Icing_L1_prec = L1_inter 
!			ELSE 
!				PRINT*,"Error in the choice of the interface SCHEME"
!				STOP 
!			END IF 
!			PRINT*,"New interface position at time",Data%T, Icing_L1,"true pos",2.d0*Icing_chi*sqrt(Icing_alpha_l*Data%T)
!		ELSE
!			!cas test polynome
!			CALL Interface_deplacement(Data) ! deplacement de la position de la véritable interface
!		END IF
	!END IF 	
	
	!Remise à O
	!Every data concerning the previous solution on the mesh for the complementary nodes 
	Mesh%Tab_NoeudComp_prec = 0 ;  Mesh%Tab_NoeudComp_prec = Mesh%Tab_NoeudComp !voir cette ligne !!!!
	Mesh%NarEmb_prec = Mesh%NarEmb
	Mesh%Tab_NoeudComp = 0 ; Mesh%NarEmb = 0.d0
	
	! remise à 0 aprés le calcul fait dans la boucle précédente
	DO i = 1, mesh%Nv
	   mesh%vertex(i)%Complementaire = 0
	END DO

	DO ie = 1, mesh%Ne
		DO i = 1, mesh%ele(ie)%nNodesPerElement
			Mesh%ele(ie)%surVertex(i) = Mesh%ele(ie)%vertex(i)
			Mesh%ele(ie)%surPos(i) = Mesh%ele(ie)%Pos(i)
			Mesh%ele(ie)%SurVertex_Sauvg(i) = 0
		END DO
	END DO

	! If embedded, define geometry and tags
	! plus define number of active elements
	IF ( Data%Embedded .OR. Data%Immersed ) THEN
		CALL computeDistance(PhysicalInterface, Data, Mesh, Sol)
		IF ( Data%Embedded .AND. .NOT. Data%Icing ) CALL tagElement2Solve(Data, Mesh)
		IF ( Data%Immersed .OR. Data%Icing ) CALL tagElement2SolveImmersed(Data, Mesh)
		Sol%NeActive = Mesh%NeActive ; Sol%NvActive = Mesh%NvActive
		Mesh%NarEmb= mesh%NSurB !nombre d'arete sur la surrogate
	END IF

	IF (Data%Icing .eqv. .TRUE.) THEN
		! Enregistrement des surVertex
		DO ie = 1, mesh%Ne
			DO i = 1, mesh%ele(ie)%nNodesPerElement
				Nu_i=mesh%ele(ie)%vertex(i)
				IF (mesh%vertex(Nu_i)%OnInterface .AND. .NOT. Data%Embedded) THEN
					IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
					   mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
					   mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
					END IF
				ELSE IF (mesh%vertex(Nu_i)%OnSurrogate .AND.  Data%Embedded) THEN
					IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
					   mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
					   mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
					END IF
				END IF
			END DO
		END DO
	END IF
	
	mesh%vertex(:)%active = .TRUE.
	Sol_prec%Dist = Sol%Dist !the vector has been modified for the new time so the disctance vector needs to be updated too 
	CALL AllocSol_PrecModif(Data, Sol_prec, Mesh) !maj du vecteur temps precedent

	Sol_2prec%Dist = Sol%Dist !the vector has been modified for the new time so the disctance vector needs to be updated too 
	CALL AllocSol_PrecModif(Data, Sol_2prec, Mesh) !maj du vecteur temps precedent
	

	IF (icas==001) THEN
		Data%Jump = JumpF
	ENDIF 

	IF(DATA%TSCHEME=="NICOLSON" ) THEN
		JumpFx = 0.d0 ; JumpFy = 0.d0
		CALL JumpPrevNewInter(Data, Mesh, Sol, Sol_prec,JumpFx,JumpFy) 	
		Data%PrevJumpSurrogate   = JumpFx !Definition of the flux jump value for the previous time on the new geometry  
		Data%PrevJumpSurrogate_x = JumpFx !Definition of the flux jump value for the previous time on the new geometry  
		Data%PrevJumpSurrogate_y = JumpFy !Definition of the flux jump value for the previous time on the new geometry  
		!CALL MeanT_Reconstruct22(Data, Mesh, Sol_prec, MeanT)
		CALL MeanTemp(Data, Mesh, Sol_prec, MeanT) 
		Data%MeanT = MeanT !mean Temperature from the previous time on the new geometry 
		PRINT*,"Value of MeanTemp",Data%MeanT
	END IF 

    CALL InitVar(Data, Sol, Mat)
    CALL FillIndCG(Mat, Mesh) ! depends on the Mesh and the complementary nodes 

  END SUBROUTINE SimulationInitialisation
  !========================================================!

 !========================================================!
  SUBROUTINE InitialisationNewton(Data_Newton,Sol_prec,Sol_prec_Newton,Mesh,Mesh_Newton)
    !========================================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(INOUT)  :: Data_Newton
    type(SolStructure) , INTENT(INOUT)  :: Sol_prec_Newton
    type(SolStructure) , INTENT(IN)     :: Sol_prec
    type(MeshType)     , INTENT(INOUT)  :: Mesh_Newton 
    type(MeshType)     , INTENT(IN)     :: Mesh
    !----------------------------------------
	
	
	Data_Newton%Tour = Data_Newton%Tour + 1
	Data_Newton%T = Data_Newton%T + Data_Newton%Dt ! Copy for Newton iteration 
	Data_Newton%NewtonProcess =  .TRUE.
	 
	IF(Data_Newton%tour==1) THEN 
		ALLOCATE(Sol_prec_Newton%Bcg_Newton(Sol_prec%Nv + Data_Newton%MaxRowEmbedded,2))
	END IF 
	
	Sol_prec_Newton%Bcg_Newton = 0.d0 
		
	Data_Newton%iterJump = -1 ! first iteration
	!we stay on a constant value of Icing_Tau 
	!we stay on a constant value of Icing_Tau 
	!Icing_Tau = Icing_Tau_init ! reinitialization of the value for the new time step
	!Icing_PrevFixPosition = Icing_L1 !saving of the previous interface position 
	
	Mesh_Newton = Mesh 
	Sol_prec_Newton = Sol_prec
	Sol_prec_Newton%Bcg_Newton = Sol_prec%Bcg ! before the first pseudo Newton iteration 


  END SUBROUTINE InitialisationNewton
  !========================================================!

  !=================================!
  SUBROUTINE InitVar(Data, Sol, Mat)
    !=================================!

    !*************************************************************************
    !  Initialization of variables in the structures Sol and Mat
    !  Parameters:
    !
    !    Input, type(DataStructure) Data, Data Structure from dom.data
    !    Input,Output, type(SolStructure) Sol, Sol Structure
    !    Input,Output, type(SolStructure) Mat, Structure of the matrix DISCRETIZATION
    !
    !*************************************************************************

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    !--------------------------------------------

  !  PRINT*, '  --  Initialisation Solution Vector  --'

    Mat%Ind = 0  
	Mat%Coeff = 0.d0 
	Sol%RHS = 0.d0


    IF ( Data%Embedded ) THEN
       Sol%EleStatus = 1 ; Sol%NodeStatus = 1
    END IF

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG", "CG-Primal" )
       Sol%pcg = 0.d0 !; Sol%pcg_inter = 0.d0! solution vector for the t    
       Sol%Bcg = 0.d0 !; Sol%Bcg_inter = 0.d0 ! solution vector for the flux
       
    CASE ( "DG", "DG-Primal" )
       Sol%pdg = 0.d0
       Sol%Bdg = 0.d0
    END SELECT

  !  PRINT*, '  --  Initialisation Solution Vector done --  '
  !  PRINT*, ' '

  END SUBROUTINE InitVar
  !=================================!

  !========================================!
  SUBROUTINE InitSol0(Mesh, Sol, Data)
    !========================================!

    !*************************************************************************
    !   Initialisation Solution at initial time Tinit
    !  Parameters:
    !    Input, type(DataStructure) Data, Data Structure from dom.data
    !    Input, type(MeshType) Mesh, Mesh Structure
    !    Output, type(SolStructure) Sol, Sol Structure
    !
    !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    !----------------------------------------------
    type(element) :: ele
    !----------------------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: Exa
    !----------------------------------------------
    REAL*8 :: x, y, z
    !----------------------------------------------
    INTEGER :: ie, i, NU_i, Nu_i2, Ne, Nv 
    INTEGER :: id, nDim, nNodesPerElement,Ns   

    !----------------------------------------------
    LOGICAL :: EleSolve
    !----------------------------------------------

    PRINT*, '  --  Initialisation Solution at initial time   --', Data%T

    IF (Data%Icing .eqv. .TRUE.) THEN
       Ns=Mesh%NB_VertInterf
    ELSE
       Ns=0
    END IF

    IF ( Data%Embedded ) THEN
       Sol%EleStatus = 1 ; Sol%NodeStatus = 1
    END IF

    nDim = Data%nDim ; Ne = Mesh%Ne ; Nv = mesh%nv 
    ALLOCATE ( Exa(nDim+1) ) ; Exa=0.d0

	IF(Data%Tinit == 0.d0 ) THEN ! means we start with one phase only 
		SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )

		CASE ( "CG" )
			Sol%Pcg(:)   = Icing_Tinit ! we consider that the temperature is uniform before switching on the resistance 
			Sol%Bcg(:,:) = 0.d0 ! no heat flux yet 
		CASE DEFAULT
		   CALL PrintError("Error in InitSol0")
		   CALL PrintError("Only CG scheme available in Primal and Mixed form")
		END SELECT
	ELSE 
		SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )

		CASE ( "CG" )
		Sol%Pcg = 0.d0 ; Sol%Bcg = 0.d0 ; Sol%RHS=0.d0
		DO ie = 1, Ne
			ele = Mesh%Ele(ie) ; eleSolve = .TRUE.
			nNodesPerElement = ele%nNodesPerElement
			IF (.NOT. Data%Embedded .AND. Data%Icing) THEN
				DO i = 1, Data%N_LS
					IF ( ele%Status(i) /= -1 ) eleSolve = .FALSE.
				END DO
			END IF
			IF ( eleSolve ) THEN
				DO i = 1, nNodesPerElement
					NU_i = ele%Vertex(i) 
					NU_i2= Ele%SurVertex(i)!on reste que sur du icing 
					!PRINT*,Nu_i,Nu_i2
					CALL ExactSolJumpt(ele%L, Exa(:), mesh%Vertex(NU_i)%X ,Data%T, ele, nDim)
						
					DO id = 1, nDim
						Sol%Bcg(NU_i2,id) = Exa(id)
					END DO
						
					Sol%Pcg(NU_i2) = Exa(nDim+1)
						
				 END DO
			END IF
		END DO	
		CASE ("CG-Primal") !initialization only for the T variable 
		Sol%pcg = 0.d0 ; Sol%RHS=0.d0
			DO ie = 1, Ne
				ele = Mesh%Ele(ie) ; eleSolve = .TRUE.
				nNodesPerElement = ele%nNodesPerElement
				IF (.NOT. Data%Embedded .AND. Data%Icing) THEN
					DO i = 1, Data%N_LS
						IF ( ele%Status(i) /= -1 ) eleSolve = .FALSE.
					END DO
				END IF

				IF ( eleSolve ) THEN
					DO i = 1, nNodesPerElement
						NU_i = ele%Vertex(i) 
						NU_i2=Ele%SurVertex(i)!on reste que sur du icing 
				
						CALL ExactSolJumpt(ele%L, Exa(:), mesh%Vertex(NU_i)%X ,Data%T, ele, nDim)
												
						Sol%Pcg(NU_i2) = Exa(nDim+1)
							
					 END DO
				END IF
			END DO
			
		CASE DEFAULT
		   CALL PrintError("Error in InitSol0")
		   CALL PrintError("Only CG scheme available in Primal and Mixed form")
		END SELECT
	END IF 


    DEALLOCATE (Exa) ; DEALLOCATE (Sol%LS) 
    DEALLOCATE (Sol%NodeStatus) ; DEALLOCATE (Sol%EleStatus) 


    PRINT*, '  --  Initialisation Solution at initial time ', Data%T, ' done -- '
    PRINT*, ' '

  END SUBROUTINE InitSol0
  !========================================!

  !==================================!
  SUBROUTINE AllocSol(Data, Sol, Mat, Mesh)
    !==================================!

    !*************************************************************************
    !  Allocation of different variables in the different structures
    !  Parameters:
    !
    !    Input, type(DataStructure) Data, Data Structure from dom.data
    !    Input,Output, type(SolStructure) Sol, Sol Structure
    !    Input,Output, type(SparseMatrix) Mat, Structure of the matrix DISCRETIZATION
    !    Input,Output, type(SolStructure) Mesh, Mesh Structure from dom.mesh
    !
    !*************************************************************************

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    type(MeshType)     , INTENT(INOUT) :: Mesh

    !----------------------------------------
    INTEGER :: i, NDof
    INTEGER :: NMax, NMaxTri, NMAxExtra, Ntot
    INTEGER :: Nn, Ntri, nDim, nVar, NvActive, NeActive
    INTEGER :: Nv, Ne, Ned, Ns
    INTEGER :: nNodesPerElement
    !-----------------------------------------

    !  PRINT*, '  --  Allocation  --  '

    Nv = Sol%Nv ; Ne = Sol%Ne ; Ned = Sol%Ned ; nDim = Data%nDim
    nNodesPerElement = nDim + 1
    NeActive = Sol%NeActive ; NvActive = Sol%NvActive

    IF (Data%Icing .eqv. .TRUE.) THEN
       Ns=Mesh%NB_VertInterf
    ELSE
       Ns=0
    END IF


    ALLOCATE ( Mesh%MappingReal2Sur(Nv + Data%MaxRowEmbedded) )
    ALLOCATE ( Mesh%MappingSur2Real(Nv + Data%MaxRowEmbedded) )
    DO i = 1, Nv + Data%MaxRowEmbedded  

       Mesh%mappingreal2Sur(i) = i
       Mesh%mappingSur2Real(i) = i
    END DO

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG" )
       NMax = 70 ; Mat%Nmax = Nmax
       NVar = nDim + 1 ; Mat%NVar = NVar
       ALLOCATE ( Sol%Pcg( Nv + Data%MaxRowEmbedded ) )
       ALLOCATE ( Sol%Bcg( Nv + Data%MaxRowEmbedded , nDim ) )
       ALLOCATE ( Mat%Ind( Nv + Data%MaxRowEmbedded , NMax ) )
       ALLOCATE ( Sol%RHS( NVar*(Nv + Data%MaxRowEmbedded)) )
       ALLOCATE ( Mat%Coeff( NVar*(NvActive + Data%MaxRowEmbedded ) , NVar*Nmax ) )
    CASE ( "CG-Primal" )
       Nmax = 70 ; Mat%Nmax = Nmax
       Nvar = 1 ; Mat%nVar = nVar
       ALLOCATE ( Sol%pcg( Nv + Data%MaxRowEmbedded ) )
       ALLOCATE ( Sol%Bcg( Nv + Data%MaxRowEmbedded , nDim ) )
       ALLOCATE ( Sol%RHS( Nv + Data%MaxRowEmbedded ) )
       ALLOCATE ( Sol%DualArea( Nv + Data%MaxRowEmbedded ) )
       !ALLOCATE ( Sol%rhsN( NvActive + Ns ) ) ! peut etre modifier la dimension
       !ALLOCATE ( Sol%rhs0( NvActive + Ns) ) ! peut etre modifier la dimension
       !ALLOCATE ( Sol%rhsTmp( NvActive + Ns ) ) ! peut etre modifier la dimension
       ALLOCATE ( Mat%Ind( Nv + Data%MaxRowEmbedded , NMax) )
       ALLOCATE ( Mat%Coeff( NvActive + Data%MaxRowEmbedded , NMax ) )
    CASE ( "DG" )
       Nmax = 70 ; Mat%Nmax = Nmax
       nVar = nDim + 1 ; Mat%NVar = nVar
       ALLOCATE ( Sol%Bdg( nNodesPerElement*(Ne+Ns) , nDim ) )
       ALLOCATE ( Sol%pdg(nNodesPerElement*(Ne+Ns)) )
       ALLOCATE ( Sol%RHS(nVar*nNodesPerElement*(NeActive+Ns)) )
       ALLOCATE ( Mat%Ind(nNodesPerElement*(NeActive+Ns) ,NMax) ) !modif faite sur la taille du tableau
       ALLOCATE ( Mat%Coeff(nNodesPerElement*nVar*(NeActive+Ns), nVar*NMax) )
    CASE ( "DG-Primal" )
       Nmax = 70
       Mat%Nmax = Nmax
       Mat%NVar = 1 ; nVar = 1
       ALLOCATE ( Sol%pdg(nVar*Ne*nNodesPerElement) )
       ALLOCATE ( Sol%RHS(nVar*NeActive*nNodesPerElement) )
       ALLOCATE ( Sol%rhsN(nVar*NeActive*nNodesPerElement ) )
       ALLOCATE ( Sol%rhs0(nVar*NeActive*nNodesPerElement) )
       ALLOCATE ( Sol%Bdg(nVar*Ne*nNodesPerElement,nDim) )
       ALLOCATE ( Mat%Ind(nNodesPerElement*NeActive ,NMax) )
       ALLOCATE ( Mat%Coeff(nNodesPerElement*NeActive, NMax) )
    CASE ( "MSEG", "MSEG-Stab" )
       nMax = 70
       mat%nMAx = nMax
       Mat%nVar = nDim + 1
       ALLOCATE ( Sol%Bcg(sol%nv, nDim) )
       ALLOCATE ( Sol%Beg(Sol%Nv + Ne,nDim) )
       ALLOCATE ( Sol%pcg(Sol%Nv ) )
       ALLOCATE ( Sol%peg(Sol%Nv + Ne) )
       ALLOCATE ( Sol%RHS((nDim+1)*Sol%Nv) )
       ALLOCATE ( Sol%Bdg((nDim+1)*Ne,nDim) )
       ALLOCATE ( Sol%pdg((nDim+1)*Ne) )
       ALLOCATE ( Mat%Ind(Nv,NMax) )
       ALLOCATE ( Mat%Coeff((nDim+1)*Nv, (nDim+1)*NMax) )
    CASE DEFAULT
       CALL PrintError("AllocSol")
    END SELECT
    Sol%nVar = Mat%nVar
    !  PRINT*, '  --  Allocation done  --  '

  END SUBROUTINE AllocSol
  !========================================!


  !==================================!
  SUBROUTINE AllocSol_init(Data, Sol, Sol_prec, Mesh, k)
    !==================================!

    !*************************************************************************
    !  Allocation of the second solution structure to keep the previous time
    !  Parameters:
    !
    !    Input, type(DataStructure) Data, Data Structure from dom.data
    !    Input,Output, type(SolStructure) Sol, Sol Structure
    !    Input, type(SolStructure) Mesh, Mesh Structure from dom.mesh
    !
    !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(SolStructure) , INTENT(INOUT) :: Sol,Sol_prec
    type(MeshType)     , INTENT(IN)    :: Mesh
    INTEGER			   , INTENT(IN)    :: k 

    !----------------------------------------
    INTEGER :: nDim, nVar, MaxRowEmbedded
    INTEGER :: Nv, Ns, i, NMax, nNodesPerElement
    !-----------------------------------------

    PRINT*, '  --  Allocation  --  '

    Nv = Mesh%Nv ; nDim = Mesh%nDim ; nNodesPerElement = nDim + 1

    IF (Data%Icing) THEN
       Ns=Mesh%NB_VertInterf
    ELSE
       Ns=0
    END IF

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG" )
        NMax = 70 ; NVar = nDim + 1 ; MaxRowEmbedded = Data%MaxRowEmbedded
        ALLOCATE ( Sol%Bcg( Nv + MaxRowEmbedded , nDim ) ) ! vecteur solution flux à nDim composantes
        ALLOCATE ( Sol%Pcg( Nv + MaxRowEmbedded ) ) ! vecteur solution température
        ALLOCATE ( Sol%RHS( NVar*(Nv + MaxRowEmbedded)) ) ! vecteur solution avant decoupage par variables not necessary at n time
        IF (k == 1) THEN 
			ALLOCATE ( Sol_prec%Bcg( Nv + MaxRowEmbedded , nDim ) ) ! vecteur solution flux à nDim composantes
			ALLOCATE ( Sol_prec%Pcg( Nv + MaxRowEmbedded ) ) ! vecteur solution température
			ALLOCATE ( Sol_prec%RHS( NVar*(Nv + MaxRowEmbedded)) ) ! vecteur solution avant decoupage par variables not necessary at n time
        END IF 
    CASE ( "CG-Primal" )   
       Nmax = 70 ; Nvar = 1 ; MaxRowEmbedded = Data%MaxRowEmbedded
       ALLOCATE ( Sol%Pcg( Nv + MaxRowEmbedded ) )
       ALLOCATE ( Sol%Bcg( Nv + MaxRowEmbedded , nDim ) )
       ALLOCATE ( Sol%RHS( Nv + MaxRowEmbedded ) )
       
    CASE DEFAULT
       PRINT*,"Error in Scheme choice "
       PRINT*,"Only CG case is available in Primal and Mixed Form"
       CALL PrintError("Error in AllocSol_init")
    END SELECT
    Sol%nVar = nVar

    PRINT*, '  --  Allocation done for T init --  '

  END SUBROUTINE AllocSol_init
  !========================================!

  !==================================!
  SUBROUTINE AllocSol_PrecModif(Data, Sol, Mesh)
  !==================================!
    ! modification du vecteur solution au temps n pour etre de la taille du vecteur au temps n+1
    IMPLICIT NONE
    type(DataStructure), INTENT(INOUT) 		:: Data
    type(SolStructure) , INTENT(INOUT) 		:: Sol
    type(MeshType)     , INTENT(INOUT)      :: Mesh

    !----------------------------------------
    INTEGER             	 :: nDim, nVar, Ns_prec, Nie , k
    INTEGER             	 :: Nv, Ns, i, NMax, nNodesPerElement, j, Nu_i, Nu_j, Nu_k
    INTEGER             	 :: vertcs, vertcs_comp, nb_comTab, choix, typeCond 
    LOGICAL             	 :: test, InIt, InIt1, InIt2,VertcsOnBound
    INTEGER, DIMENSION(2,2)	 :: typeV ! on va regarder le type des deux aretes associé 
    REAL*8, DIMENSION(3) 	 :: Exa 
    type(Element)            :: ele_original  !element etudie
    CHARACTER(len=20)        :: BCType1, BCType2
    REAL*8              	 :: val1,val2,val3
    !-----------------------------------------
    REAL*8   , DIMENSION(:,:), ALLOCATABLE :: Sol_inter_Bcg
    REAL*8   , DIMENSION(:)  , ALLOCATABLE :: Sol_inter_Pcg
    !-----------------------------------------

    !  PRINT*, '  --  Allocation  --  '

    Nv = Sol%Nv ; nDim = Data%nDim ; nNodesPerElement = nDim + 1
    Mesh%Vertex(:)%EmptyComp = .TRUE. ; Test = .FALSE. ! pour savoir si une valeur n'a pas encore était reconstruite sur le noeud 
	typeV = 0 
    IF (Data%Icing) THEN
		!ici ce sont des variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb           !nombre de noeud sur l'interface au temps n+1
		Ns_prec=Mesh%NarEmb_prec + 1 + mesh%PtEmb 
    ELSE
       Ns=0 ; Ns_prec=0
    END IF

    CALL VERIF_Comp(Mesh,test) !on verifie que la surrogate n'est pas la meme qu'au temps precedent
		! il suffit de comparer les tableaux avec complementaires
		IF (.NOT. Test) THEN  

			PRINT*,"Displacement of the interface, Geometry Modification "

			SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
			CASE ( 'CG','CG-Primal' )
				! vecteurs intermediaires pour la partie noeuds complementaire 
				ALLOCATE(Sol_inter_Bcg(Mesh%NB_VertInterf,2)) ; ALLOCATE(Sol_inter_Pcg(Mesh%NB_VertInterf))
				! Mesh%NB_VertInterf REPRESENTE SEULEMENT LE nombre de noeuds avec un complementaire 
				Sol_inter_Bcg = 0.d0 ; Sol_inter_Pcg = 0.d0
				PRINT*,"**** Deplacement of the shifted interface ****"
				!il faut récupérer la valeur du tableau de complémentaire précédent
				DO j= 1, Ns_prec
					vertcs = mesh%Tab_NoeudComp_prec(j,1) !ref du noeud original
					vertcs_comp = mesh%Tab_NoeudComp_prec(j,2) ! ref associe au noeud comp
					! on verifie d'abord que ce noeud n'est pas dans le nouveau tableau de complementaires
					! et qu'il est donc toujour sur la surrogate 
					choix = 0
					!on déplace les valeurs déja présente dans la partie complementaire au temps precedents
					CALL LookForComp(Mesh,InIt,vertcs,nb_comTab,choix) ! on regarde si vertcs a un noeud dans la nvelle table
					!le cas vertcs_comp /= vertcs est present pour traiter les noeuds presents à l'interface sans complemetaire
					! cas des noeuds de bords gauche si on commence à t=0
					IF (InIt .AND. vertcs_comp /= vertcs) THEN ! cas ou le complementaire est toujours présent, !! il peut etre a un emplacement different
						!ce qu'on veut ces mettre les valeurs complementaires au bon endroit avec leur nouvelle numérotation
						Sol_inter_Bcg(nb_comTab-Nv,1) = Sol%Bcg(vertcs_comp,1)
						Sol_inter_Bcg(nb_comTab-Nv,2) = Sol%Bcg(vertcs_comp,2)
						Sol_inter_Pcg(nb_comTab-Nv)   = Sol%Pcg(vertcs_comp)
					END IF
				END DO

				!Now we look at the current complementaries table
				j = 0
				DO i=1, Ns
				
					j=j+1
					Nu_i = mesh%Tab_NoeudComp(j,1) ; Nu_j = mesh%Tab_NoeudComp(j,2)
					Init= .FALSE.
					choix = 1 ! on regarde si il est dans la table precedente
					CALL LookForComp(Mesh,InIt,Nu_i,nb_comTab,choix)
					IF (.NOT. INIT) THEN ! valeurs non modifies
						IF(Data%IterJump==0) THEN 
							nbChangInterface = nbChangInterface+1 ! variable globale pour savoir combien de fois on a bougé l'interface 
						END IF 
						!On copie la valeur initale dans celle du complementaire car le noeud complementare est associé
						! à la valeur à droite de l'interface ou les élements non pas changer de référence
						!par contre les noeuds de gauche eux sont associé à des elements qui ont changé de references avec le mouvement de l'interface
						Sol_inter_Bcg(Nu_j-nv,1) = Sol%Bcg(Nu_i,1)
						Sol_inter_Bcg(Nu_j-Nv,2) = Sol%Bcg(Nu_i,2)
						Sol_inter_Pcg(Nu_j-Nv)   = Sol%Pcg(Nu_i)
					 
				    	! il faut d'abord faire un test pou savoir si le noeud que l'on regarde se trouve sur le bord 
						VertcsOnBound = .FALSE.
						CALL VertcsPos(Mesh, Nu_i, VertcsOnBound, typeV)
						IF(.NOT. VertcsOnBound) THEN !cas ou le noeud ne se trouve pas sur le bord  
							!car un noeud sur le bord doit d'abord verifié la condition de bord 
							!modification de la valeur du noeud initial
							! la valeur que l'on reconstruit c'est la valeur de gauche dans la zone référencé 1 
							CALL PolyConstrQR(Data,Mesh,Sol,Nu_i,val1,val2,val3)
							
							Sol%Pcg(Nu_i)   = val1 !valeur reconstruite pour la temperature
							Sol%Bcg(Nu_i,1) = val2 !valeur reconstruite pour le flux en x
							Sol%Bcg(Nu_i,2) = val3 !valeur reconstruite pour le flux en y
							
						ELSE ! deplacement d'un noeud de bord, la valeur de la condition de bord prédomine 

							!PRINT*,"Reconstruction of a node at the boundary, ref",Nu_i
							
							BCType1 = TRIM(ADJUSTL(Data%Bc(TypeV(1,1)))) ! condition associé à la reference de bord 
							BCType2 = TRIM(ADJUSTL(Data%Bc(TypeV(2,1)))) 
							! en fonction de type de condition de bord on va se 
							IF(BCType1 == BCType2) THEN 
							
								!pour reconstruire les valeurs non données par les conditions de bords 
								CALL PolyConstrQR(Data,Mesh,Sol,Nu_i,val1,val2,val3)
			
								!IF (icas==001) THEN !cas specific stefan 
								!
								!	!on recupere la condition sur le fichier de donnée 
								!	SELECT CASE ( TRIM(ADJUSTL(BCType1)) )
							!
							!		CASE ( "Dirichlet", "dirichlet" )
							!										
							!		IF(	TypeV(1,1)==1) THEN					
							!			Sol%Pcg(Nu_i)   = Icing_Tw  
							!			Sol%Bcg(Nu_i,1) = val2 
							!			Sol%Bcg(Nu_i,2) = val3 
							!		ELSE 
							!			PRINT*,"error for boundary reconstruction test 001"
							!			STOP
							!		END IF 
							!	 
							!		CASE ( "Neumann", "neumann" )
							!	
							!		SELECT CASE ( TypeV(1,1) )
							!																		
							!			CASE(2) !bord de reference 2
							!			! la condition de Neumann donne -By 
!
!											Sol%Pcg(Nu_i)   = val1 
!											Sol%Bcg(Nu_i,1) = val2 
!											Sol%Bcg(Nu_i,2) = 0 
!										
!										CASE(4) !bord de reference 4
!										! la condition de Neumann donne By 
!
!											Sol%Pcg(Nu_i)   = val1  
!											Sol%Bcg(Nu_i,1) = val2 
!											Sol%Bcg(Nu_i,2) = 0 
!										
!										CASE DEFAULT
!										   PRINT*, " ** ERROR inAllocSol_PrecModif case 001"
!										   STOP
!										END SELECT
!												
!								CASE DEFAULT
!								   PRINT*, " ** ERROR inAllocSol_PrecModif case 001"
!								   STOP
!								END SELECT	
						   
!								ELSE !cas general 
							
									!on ne se trouve pas sur un coin
									! besoin d'un element de référence pour aller chercher la valeur de condition de bord 		
									DO Nie= 1, mesh%Ne ! boucle sur le nombre d'élément
										DO k = 1, mesh%ele(Nie)%nNodesPerElement ! boucle sur le nombre de noeud associé à l'élement
											Nu_k=mesh%ele(Nie)%vertex(k) ! noeud sur l'element
											IF(Nu_k == Nu_i .AND. mesh%ele(Nie)%ref == 1 ) THEN
												ele_original = mesh%ele(Nie) ! on veut l'enregistrement d'un element associé au noeud ou on a la valeur manquante pour l'utiliser pour le calcul de la solution exacte 
											END IF
										END DO 
									END DO 
											
									CALL ExactSolJumpt(ele_original%L, Exa,Mesh%Vertex(Nu_i)%X,Data%T-data%Dt, ele_original, nDim)
									
									SELECT CASE ( TRIM(ADJUSTL(BCType1)) )
									
										CASE ( "Dirichlet", "dirichlet" )
																									
											Sol%Pcg(Nu_i)   = Exa(3)  
											Sol%Bcg(Nu_i,1) = val2 
											Sol%Bcg(Nu_i,2) = val3 
										 
										CASE ( "Neumann", "neumann" )
										
											SELECT CASE ( TypeV(1,1) )
											
												CASE(1) !bord 1
												! la condition de Neumann donne -Bx
																													
													Sol%Pcg(Nu_i)   = val1  
													Sol%Bcg(Nu_i,1) = Exa(1) 
													Sol%Bcg(Nu_i,2) = val3 								

												
												CASE(2) !bord de reference 2
												! la condition de Neumann donne -By 

													Sol%Pcg(Nu_i)   = val1 
													Sol%Bcg(Nu_i,1) = val2 
													Sol%Bcg(Nu_i,2) = Exa(2) 
												
												CASE(3) !bord de reference 3 
												! la condition de Neumann donne Bx 
												
													Sol%Pcg(Nu_i)   = val1 
													Sol%Bcg(Nu_i,1) = Exa(1) 
													Sol%Bcg(Nu_i,2) = val3 
												
												CASE(4) !bord de reference 4
												! la condition de Neumann donne By 

													Sol%Pcg(Nu_i)   = val1  
													Sol%Bcg(Nu_i,1) = val2 
													Sol%Bcg(Nu_i,2) = Exa(2) 
												
												CASE DEFAULT
												   PRINT*, " ** ERROR inAllocSol_PrecModif"
												   STOP
												END SELECT
														
										CASE DEFAULT
										   PRINT*, " ** ERROR inAllocSol_PrecModif"
										   STOP
								   END SELECT	
						   
								!END IF 
							ELSE 
							! noeud sur un coin pour un domaine rectangulaire , !! a modifier continuité T 
							!! 
							! a definir 
							!! 
							END IF ! fin de la reconstruction ou il manquait une valeur 
						END IF 
					END IF  
					! cas avec interface sur le bord 
					IF (INIT .AND. Nu_i/= Nu_j .AND. Nu_i == nb_comTab ) THEN !i;e que l'on avait un noeud de bord sans complementaire 
						!on recupére la condition associé aux noeUD 
								! à remplir 
					END IF 	 ! fin (INIT .AND. Nu_i/= Nu_j .AND. Nu_i == nb_comTab ) 
			  END DO
			  
			  
			  !L'interface ne peut passer par plus d'un élément à la fois ce qui veut dire que tous les noeuds possède
			  ! à un moment une valeur complémentaire
			  !Cependant il existe une exception pour les noeuds de bords gauche quand on considère une interface confondus avec le bord 
			  ! il faut ajouter que le noeud que l'on regarde était bien sur le bord 
		   !   j=0 
			!   DO i=1, Ns_prec
			!	  j=j+1
			!	  Nu_i = mesh%Tab_NoeudComp_prec(j,1) ; Nu_j = mesh%Tab_NoeudComp_prec(j,2)
			!	  
			!	  Init1= .FALSE.
			!	  choix = 1 ! on regarde si il est dans la table precedente
			!	  CALL LookForComp(Mesh,InIt1,Nu_i,nb_comTab,choix)
			 ! 
			  !	  InIt2 = .FALSE. 
				!  choix = 0 ! on regarde si il est dans la table courante
				!  CALL LookForComp(Mesh,InIt2,Nu_i,k,choix)
				 ! 
				!  IF (INIT1 .AND. Nu_i == nb_comTab .AND. .NOT. InIt2 ) THEN
				!	! ca veut dire que le noeud Nu_i n'a jamais eu de complementaire on doit donc reconstruire le noeud pour le définir dans la zone 1 
				!	CALL PolyConstr7(Data,Mesh,Sol,Nu_i,val1,val2,val3)
				!	
				!	! besoin d'un element de référence pour la solution exacte 
				!	DO Nie= 1, mesh%Ne ! boucle sur le nombre d'élément
				!		DO k = 1, mesh%ele(Nie)%nNodesPerElement ! boucle sur le nombre de noeud associé à l'élement
				!			Nu_k=mesh%ele(Nie)%vertex(k) ! noeud sur l'element
				!			IF(Nu_k == Nu_i .AND. mesh%ele(Nie)%ref == 1 ) THEN
				!				ele_original = mesh%ele(Nie) ! on veut l'enregistrement d'un element associé au noeud ou on a la valeur manquante pour l'utiliser pour le calcul de la solution exacte 
				!			END IF
				!		END DO 
				!	END DO 
				!			
				!	CALL ExactSolJumpt(ele_original%L, Exa,Mesh%Vertex(Nu_i)%X, &
				!		Data%T-data%Dt, ele_original, nDim)
				!					
				!	IF(DATA%ContinuousTmp) THEN ! on n'a pas a reconstruit la temperature 
				!		!Sol%Pcg(Nu_i)   = Exa(3) ! valeur de bord 
				!		Sol%Bcg(Nu_i,1) = val2 !valeur reconstruite pour le flux en x
				!		Sol%Bcg(Nu_i,2) = val3 !valeur reconstruite pour le flux en y
				!	ELSE 
				!		Sol%Pcg(Nu_i)   = Exa(3) ! valeur de bord 
				!		Sol%Bcg(Nu_i,1) = val2 !valeur reconstruite pour le flux en x
				!		Sol%Bcg(Nu_i,2) = val3 !valeur reconstruite pour le flux en y
				!	END IF 		     
				 !  END IF 
				!END DO 

			   !la partie Sol%Bcg(Nv,:) a deja etait modifié de meme pour Sol%Pcg(Nv)
			   Sol%Bcg(Nv+1:Nv+Mesh%NB_VertInterf,1)= Sol_inter_Bcg(:,1)
			   Sol%Bcg(Nv+1:Nv+Mesh%NB_VertInterf,2)= Sol_inter_Bcg(:,2)
			   Sol%Pcg(Nv+1:Nv+Mesh%NB_VertInterf)= Sol_inter_Pcg(:)
			   
			   !remise a 0 de la partie non utilisé
			   Sol%Bcg(Nv + Mesh%NB_VertInterf + 1 : Nv + data%MaxRowEmbedded ,:) = 0.d0
			   Sol%Pcg(Nv + Mesh%NB_VertInterf + 1 : Nv + data%MaxRowEmbedded ) = 0.d0
			   

							
			   DEALLOCATE(Sol_inter_Bcg) ; DEALLOCATE(Sol_inter_Pcg)
		CASE DEFAULT
			CALL PrintError("AllocSol_PrecModif")
		END SELECT
    END IF
			
			
      PRINT*, '  -- Structural Modification after moving the surrogate done for Sol at time n --  '
	END SUBROUTINE AllocSol_PrecModif
	!========================================!
	
	
  !==================================!
  SUBROUTINE AllocSol_PrecModif_Newton(Data, Sol, Mesh)
  !==================================!
    ! modification du vecteur solution au temps n pour etre de la taille du vecteur au temps n+1
    IMPLICIT NONE
    type(DataStructure), INTENT(INOUT) 		:: Data
    type(SolStructure) , INTENT(INOUT) 		:: Sol
    type(MeshType)     , INTENT(INOUT)      :: Mesh

    !----------------------------------------
    INTEGER             	 				:: nDim, nVar, Ns_prec, Nie , k
    INTEGER             	 				:: Nv, Ns, i, NMax, nNodesPerElement, j, Nu_i, Nu_j, Nu_k
    INTEGER             					:: vertcs, vertcs_comp, nb_comTab, choix, typeCond 
    LOGICAL             	 				:: test, InIt, InIt1, InIt2,VertcsOnBound
    INTEGER, DIMENSION(2,2)	 				:: typeV ! on va regarder le type des deux aretes associé 
    REAL*8, DIMENSION(3) 					:: Exa 
    type(Element)            				:: ele_original  !element etudie
    CHARACTER(len=20)       				:: BCType1, BCType2
    REAL*8              	 				:: val1,val2,val3
    !-----------------------------------------
    REAL*8   , DIMENSION(:,:), ALLOCATABLE  :: Sol_inter_Bcg
    !-----------------------------------------

    !  PRINT*, '  --  Allocation  --  '

    Nv = Sol%Nv ; nDim = Data%nDim ; nNodesPerElement = nDim + 1
    Mesh%Vertex(:)%EmptyComp = .TRUE. ; Test = .FALSE. ! pour savoir si une valeur n'a pas encore était reconstruite sur le noeud 
	typeV = 0 
	
    IF (Data%Icing) THEN
		!ici ce sont des variables pour se déplacer dans le tableau 
		Ns = Mesh%NarEmb + 1 + Mesh%PtEmb           
		Ns_prec = Mesh%NarEmb_prec_inter + 1 + Mesh%PtEmb 
    ELSE
        Ns=0 ; Ns_prec=0
    END IF

    CALL VERIF_Comp_Newton(Mesh,test) !on verifie que la surrogate n'est pas la meme qu'au temps precedent
		! il suffit de comparer les tableaux avec complementaires
		IF (.NOT. Test) THEN  
			SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
			CASE ( 'CG','CG-Primal' )
				! vecteurs intermediaires pour la partie noeuds complementaire 
				ALLOCATE(Sol_inter_Bcg(Mesh%NB_VertInterf,2)) ; Sol_inter_Bcg = 0.d0 
				PRINT*,"**** Deplacement of the shifted interface ****"
				!il faut récupérer la valeur du tableau de complémentaire précédent
				DO j= 1, Ns_prec
					vertcs      = mesh%Tab_NoeudComp_prec_inter(j,1) !ref du noeud original
					vertcs_comp = mesh%Tab_NoeudComp_prec_inter(j,2) !ref associe au noeud comp
					
					choix = 0
					CALL LookForComp(Mesh,InIt,vertcs,nb_comTab,choix) ! on regarde si vertcs a un noeud dans la nouvelle table

					IF (InIt .AND. vertcs_comp /= vertcs) THEN ! cas ou le complementaire est toujours présent, !! il peut etre a un emplacement different
						!ce qu'on veut ces mettre les valeurs complementaires au bon endroit avec leur nouvelle numérotation
						Sol_inter_Bcg(nb_comTab-Nv,1) = Sol%Bcg_Newton(vertcs_comp,1)
						Sol_inter_Bcg(nb_comTab-Nv,2) = Sol%Bcg_Newton(vertcs_comp,2)
					END IF
				END DO

				!Now we look at the current complementaries table
				j = 0
				DO i = 1, Ns
				
					j=j+1
					Nu_i = mesh%Tab_NoeudComp(j,1) ; Nu_j = mesh%Tab_NoeudComp(j,2)
					Init= .FALSE.
					choix = 2 ! on regarde si il est dans la table precedente mais sur les iter en k 
					CALL LookForComp(Mesh,InIt,Nu_i,nb_comTab,choix)
					IF (.NOT. INIT) THEN ! valeurs non modifies
						IF(Data%IterJump==0) THEN 
							nbChangInterface = nbChangInterface+1 ! variable globale pour savoir combien de fois on a bougé l'interface 
						END IF 
						!On copie la valeur initale dans celle du complementaire car le noeud complementare est associé
						!par contre les noeuds de gauche eux sont associé à des elements qui ont changé de references avec le mouvement de l'interface
						Sol_inter_Bcg(Nu_j-nv,1) = Sol%Bcg_Newton(Nu_i,1)
						Sol_inter_Bcg(Nu_j-Nv,2) = Sol%Bcg_Newton(Nu_i,2)
					 
					 
						! il faut d'abord faire un test pou savoir si le noeud que l'on regarde se trouve sur le bord 
						VertcsOnBound = .FALSE.
						CALL VertcsPos(Mesh, Nu_i, VertcsOnBound, typeV)
						IF(.NOT. VertcsOnBound) THEN  

							CALL PolyConstrQR(Data,Mesh,Sol,Nu_i,val1,val2,val3)

							Sol%Bcg_Newton(Nu_i,1) = val2 !valeur reconstruite pour le flux en x
							Sol%Bcg_Newton(Nu_i,2) = val3 !valeur reconstruite pour le flux en y
							
						ELSE ! deplacement d'un noeud de bord, la valeur de la condition de bord prédomine 

							
							BCType1 = TRIM(ADJUSTL(Data%Bc(TypeV(1,1)))) ; BCType2 = TRIM(ADJUSTL(Data%Bc(TypeV(2,1)))) 

							IF(BCType1 == BCType2) THEN 
							
								CALL PolyConstrQR(Data,Mesh,Sol,Nu_i,val1,val2,val3)
			
								IF (icas==001) THEN !cas specific stefan 
								
									SELECT CASE ( TRIM(ADJUSTL(BCType1)) )
							
									CASE ( "Dirichlet", "dirichlet" )
																	
									IF(	TypeV(1,1)==1) THEN					
										Sol%Bcg_Newton(Nu_i,1) = val2 
										Sol%Bcg_Newton(Nu_i,2) = val3 
									ELSE 
										PRINT*,"error for boundary reconstruction test 001"
										STOP
									END IF 
								 
									CASE ( "Neumann", "neumann" )
								
									SELECT CASE ( TypeV(1,1) )
																									
										CASE(2) !bord de reference 2
										! la condition de Neumann donne -By 

											Sol%Bcg_Newton(Nu_i,1) = val2 
											Sol%Bcg_Newton(Nu_i,2) = 0 
										
										CASE(4) !bord de reference 4
										! la condition de Neumann donne By 

											Sol%Bcg_Newton(Nu_i,1) = val2 
											Sol%Bcg_Newton(Nu_i,2) = 0 
										
										CASE DEFAULT
										   PRINT*, " ** ERROR inAllocSol_PrecModif case 001"
										   STOP
										END SELECT
												
								CASE DEFAULT
								   PRINT*, " ** ERROR inAllocSol_PrecModif case 001"
								   STOP
								END SELECT	
						   
								ELSE !cas general 
							
									!on ne se trouve pas sur un coin
									! besoin d'un element de référence pour aller chercher la valeur de condition de bord 		
									DO Nie= 1, mesh%Ne ! boucle sur le nombre d'élément
										DO k = 1, mesh%ele(Nie)%nNodesPerElement ! boucle sur le nombre de noeud associé à l'élement
											Nu_k=mesh%ele(Nie)%vertex(k) ! noeud sur l'element
											IF(Nu_k == Nu_i .AND. mesh%ele(Nie)%ref == 1 ) THEN
												ele_original = mesh%ele(Nie) ! on veut l'enregistrement d'un element associé au noeud ou on a la valeur manquante pour l'utiliser pour le calcul de la solution exacte 
											END IF
										END DO 
									END DO 
											
									CALL ExactSolJumpt(ele_original%L, Exa, Mesh%Vertex(Nu_i)%X, Data%T-data%Dt, ele_original, nDim)
									
									SELECT CASE ( TRIM(ADJUSTL(BCType1)) )
									
										CASE ( "Dirichlet", "dirichlet" )
																									
											Sol%Bcg_Newton(Nu_i,1) = val2 
											Sol%Bcg_Newton(Nu_i,2) = val3 
										 
										CASE ( "Neumann", "neumann" )
										
											SELECT CASE ( TypeV(1,1) )
											
												CASE(1) !bord 1
												! la condition de Neumann donne -Bx
																													
													Sol%Bcg_Newton(Nu_i,1) = Exa(1) 
													Sol%Bcg_Newton(Nu_i,2) = val3 								

												
												CASE(2) !bord de reference 2
												! la condition de Neumann donne -By 

													Sol%Bcg_Newton(Nu_i,1) = val2 
													Sol%Bcg_Newton(Nu_i,2) = Exa(2) 
												
												CASE(3) !bord de reference 3 
												! la condition de Neumann donne Bx 
												
													Sol%Bcg_Newton(Nu_i,1) = Exa(1) 
													Sol%Bcg_Newton(Nu_i,2) = val3 
												
												CASE(4) !bord de reference 4
												! la condition de Neumann donne By 

													Sol%Bcg_Newton(Nu_i,1) = val2 
													Sol%Bcg_Newton(Nu_i,2) = Exa(2) 
												
												CASE DEFAULT
												   PRINT*, " ** ERROR inAllocSol_PrecModif"
												   STOP
												END SELECT
														
										CASE DEFAULT
										   PRINT*, " ** ERROR inAllocSol_PrecModif"
										   STOP
								   END SELECT	
						   
								END IF 
							ELSE 
								PRINT*,"Error in AllocSol_PrecModif_Newton"
								STOP
							END IF ! fin de la reconstruction ou il manquait une valeur 
						END IF 
					END IF  
					! cas avec interface sur le bord 
					IF (INIT .AND. Nu_i/= Nu_j .AND. Nu_i == nb_comTab ) THEN !i;e que l'on avait un noeud de bord sans complementaire 
						PRINT*,"Error in AllocSol_PrecModif_Newton"
						STOP
					END IF 	 ! fin (INIT .AND. Nu_i/= Nu_j .AND. Nu_i == nb_comTab ) 
			  END DO
	

			   !la partie Sol%Bcg(Nv,:) a deja etait modifié de meme pour Sol%Pcg(Nv)
			   Sol%Bcg_Newton(Nv+1:Nv+Mesh%NB_VertInterf,1)= Sol_inter_Bcg(:,1)
			   Sol%Bcg_Newton(Nv+1:Nv+Mesh%NB_VertInterf,2)= Sol_inter_Bcg(:,2)
			   
			   !remise a 0 de la partie non utilisé
			   Sol%Bcg_Newton(Nv + Mesh%NB_VertInterf + 1 : Nv + data%MaxRowEmbedded ,:) = 0.d0
			   
							
			   DEALLOCATE(Sol_inter_Bcg)
		CASE DEFAULT
			CALL PrintError("AllocSol_PrecModif")
		END SELECT
    END IF
			
			
      PRINT*, '  -- Structural Modification after moving the surrogate done for Sol at time n --  '
	END SUBROUTINE AllocSol_PrecModif_Newton
	!========================================!
	

END MODULE Initialisation
