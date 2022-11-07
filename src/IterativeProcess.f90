MODULE IterativeProcess

	USE Types
	USE PublicVar
	USE Data_mod
	USE Initialisation
	USE LHS_RHS_SchemeSelection
	USE SystemResolution
	USE PostProcess
    USE SparseMatrix_CG
    USE Distance
    USE EmbeddedVariableInitialization
    USE IcingTest_mod 
    USE Velocity
  	USE L2Norm
	USE WriteSol


CONTAINS


  !========================================================!
  SUBROUTINE CopyValForN_Sol(Sol_Newton,Sol,Mesh_Newton,Mesh,Data_Newton,Data)
  !========================================================!


    IMPLICIT NONE
    
    !--------------------------------------------------------------------

    type(SolStructure) , INTENT(INOUT)    :: Sol_Newton
    type(SolStructure) , INTENT(IN) 	  :: Sol
    
	type(MeshType)     , INTENT(INOUT) 	  :: Mesh
    type(MeshType)     , INTENT(INOUT)	  :: Mesh_Newton
    
    type(DataStructure), INTENT(INOUT)    :: Data_Newton
    type(DataStructure), INTENT(INOUT) 	  :: Data

    !--------------------------------------------------------------------
  
    !Every data concerning the previous solution on the mesh for the complementary nodes 
    Mesh_Newton%NarEmb_prec        = Mesh%NarEmb_prec
	Mesh_Newton%Tab_NoeudComp_prec = 0.d0 
	Mesh_Newton%Tab_NoeudComp_prec = Mesh%Tab_NoeudComp_prec
	Data_Newton%Dt = Data%Dt 

	!values saved to modified the jump value after the displacement of the interface 
	!Complementary node's info from the pseudo Newton part 
	!Icing_Tau = Data%h/Data%Dt*1.d0/0.6d0*1.d0
	!PRINT*,"Here is the value of Dtau",Icing_Tau 
		
	IF(Data_newton%IterJump==-1) THEN 
		Mesh_Newton%NarEmb_prec_inter        = Mesh%NarEmb_prec  
		Mesh_Newton%Tab_NoeudComp_prec_inter = mesh%Tab_NoeudComp_prec
		Sol_Newton%Bcg_Newton 		         = Sol%Bcg !the value from the previous pseudo newton iter is kept 

	ELSE 
		Mesh_Newton%NarEmb_prec_inter        = Mesh_Newton%NarEmb  
		Mesh_Newton%Tab_NoeudComp_prec_inter = Mesh_Newton%Tab_NoeudComp
		Sol_Newton%Bcg_Newton 		         = Sol_Newton%Bcg !the value from the previous pseudo newton iter is kept 
	END IF
	
	
	!from the prevois time n  
	!Solution not modified by the geometry -- starting point of the algorithm 
    Sol_Newton%Bcg  = Sol%Bcg
    Sol_Newton%Pcg  = Sol%Pcg
    Sol_Newton%Dist = Sol%Dist

       
  END SUBROUTINE CopyValForN_Sol
  !========================================================!
 	
 	
!  !========================================================!
!  SUBROUTINE UpdateFrontPseudoNewton(Sol_prec,Sol,Mesh,Data,Mat,nbMaxITER,TolSousIter,SpeedJump, JumpF) 
!  !========================================================!
!
!
 !   IMPLICIT NONE
  !  
!    type(SolStructure) , INTENT(INOUT) :: Sol_prec, Sol
!    type(SparseMatrix) , INTENT(INOUT) :: Mat
!    INTEGER		 	   , INTENT(IN)    :: nbMaxITER 
!    REAL*8       	   , INTENT(IN)    :: TolSousIter
!    REAL*8       	   , INTENT(OUT)   :: SpeedJump, JumpF
!    type(DataStructure), INTENT(INOUT) :: Data
!    type(MeshType)     , INTENT(INOUT) :: Mesh

!    !--------------------------------------------------------------------
!    REAL*8 				:: truePos , trueVelo, trueVelo_prev, JumpFx,JumpFy,TestVelocity
!    REAL*8 				:: v_prec, v2_prec, test_Velocity,testPos
!    REAL*8 				:: SpeedJump_diff, JumpF_diff,SpeedJump_test, JumpF_test
!    REAL*8     			:: SpeedJumpNew, JumpNew, interJ, TolSousIterPos,Speed_prec,testJump_r
!    INTEGER 			:: Ns_prec
!    !--------------------------------------------------------------------

!    LOGICAL     		:: test
 !   INTEGER 			:: i, ie , Nu_i, NarEmb_prec_inter
  !  type(MeshType)      :: Mesh_inter
	!type(DataStructure) :: Data_inter
!	type(SolStructure)  :: Sol_prec_inter, Sol_Nicolson
!	CHARACTER(len=99)   :: charTmp1, charTmp2,charTmp3
    !--------------------------------------------------------------------

    !Analytical solutions 
 !   truePos       = 2.d0*Icing_chi*sqrt(Icing_alpha_l*Data%T)
  !  trueVelo      = Icing_chi*sqrt(Icing_alpha_l/Data%T)
!	trueVelo_prev = Icing_chi*sqrt(Icing_alpha_l/(Data%T-Data%DT))!
!	Ns_prec = Mesh%NarEmb_prec_inter + 1 + Mesh%PtEmb 
!	SpeedJump = 0.d0 ; JumpF = 0.d0 
!	TolSousIterPos = 1.E-10
	
!	IF(Data%iterJump == -1) THEN !Initialisation for the first iteration 
!		Data%ResNewt = 0.d0 ; Data%ResNewtD = 0.d0 
!	END IF 
	
!	IF(icas==001) THEN
!		PRINT*,"---- Calculation of the approximative Jump value ----"
!		! the mean value for the jump and the velocity
!		IF(Data%iterJump == -1 ) THEN 
!			
!			CALL MeanFluxReconstruct22(Data, Mesh, Sol_prec, SpeedJump, JumpF) 
!						
!			Sol%Bcg_inter = Sol_prec%Bcg ; Sol%Pcg_inter = Sol_prec%Pcg  ! for pseudo newton 
!
!			Data%InitResNewt = 0.d0 
!			DO i = 1, Mesh%NarEmb + 1 + mesh%PtEmb 
!				Data%InitResNewt = Data%InitResNewt + Sol_prec%Bcg(i,1)**2 
!			END DO 	
!					
!		ELSE 
!			IF(DATA%TSCHEME=="EULER") THEN 
!				CALL MeanFluxReconstruct22(Data, Mesh, Sol, SpeedJump, JumpF) 
!				Sol%Bcg_inter = Sol%Bcg ; Sol%Pcg_inter = Sol%Pcg    ! for pseudo newton 
!			ELSE IF (Data%TSCHEME=="NICOLSON") THEN 
!			!we need the change of variable fo r the value of the real jump 
!				CALL MeanFluxReconstruct22Nicolson(Data, Mesh, Sol, SpeedJump, JumpF) 
!				Sol%Bcg_inter = Sol%Bcg ; Sol%Pcg_inter = Sol%Pcg    ! for pseudo newton 
!			END IF 
!			!JumpF = Icing_rho*Icing_Lm*Icing_chi*sqrt(Icing_alpha_l/Data%T)
!			!SpeedJump = Icing_chi*sqrt(Icing_alpha_l/Data%T)
!
!		END IF 
!		
!		PRINT*,"--- END Calculation of the jump value ---" 
!		
!		IF(Data%Space_Type =="IMPLICIT") THEN 
!		!	IF(Data%IterJump==-1) THEN !initialization of the method 
!	   !
!		!		Data%Tab_Prev2Jump(Data%tour,1) = SpeedJump
!		!		Data%Tab_Prev2Jump(Data%tour,2) = trueVelo_prev !only to compare 
!!
!!				IF(Data%Tour==1) THEN 
!!					! we need a value from before the initialization 
!!					v2_prec = Icing_chi*sqrt(Icing_alpha_l/(Data%T-(2*Data%DT))) 
!!					v_prec  = Data%Tab_Prev2Jump(Data%tour,1)
!!				ELSE
!!					v2_prec = Data%Tab_Prev2Jump(Data%tour-1,1) ! value from the previous time 
!!					v_prec  = Data%Tab_Prev2Jump(Data%tour,1)
!!				END IF 
!!				! polynomial reconstruction of order 1 for initialization of the velocity 
!!				!SpeedJump = 1.d0/Data%Dt*((Data%T-Data%Dt)*v2_prec - (Data%T-2*Data%Dt)*v_prec+Data%T*(v_prec-v2_prec))
!!				!modifier l'initialisation pour avoir quelques chose de plus proche de la vraie solution en prenant une approximation par le bas
!!				!JumpF = Icing_rho*Icing_Lm*SpeedJump
!!			END IF
!!			
	!		!interJ = 1.d0 
	!		!interJ = 0.95+(Data%iTerJump+1)*0.001
	!		!If(interJ>1.d0) interJ = 1.d0 
	!		!PRINT*,"Valeur of interJ", interJ 
	!		!JumpF = JumpF*interJ 
	!		!SpeedJump = JumpF/(Icing_rho*Icing_Lm) 
	!		
	!		!IF(Data%h<0.0001) THEN
			!	interJ = 1.d0 !0.99 !(1.01-Data%h*10**3)
	!		!	!PRINT*,"Value of interJ",interJ
	!		!	!IF(interJ>1.0) THEN 
	!		!	!	!interJ = 1.0d0
	!		!	!	PRINT*,"Error in the definition of the jump"
	!		!	!	STOP
	!		!	!	
	!		!	!ELSE
	!		!	!	IF (interJ>0.98) interJ = 0.97
	!		!		JumpF = JumpF*interJ ; SpeedJump = JumpF/(Icing_rho*Icing_Lm) 
	!		!	!END IF			
	!		!ELSE
	!		!	PRINT*,"The value of Data%h is too big for the simulation",Data%h 
	!		!	STOP 
	!		!END IF 
	!		
	!		!Only done after one iteration 
	!		!test on the jump value 
	!		IF (Data%IterJump /= -1 ) testJump_r = abs((Data%Icing_Prev_Velocity-JumpF)/Data%Icing_Prev_Velocity)
	!		
	!		!test on the velocity of the interface 
	!		!Speed_prec = (Icing_L1 - 4.d0/3.d0*Icing_L1_prec+1.d0/3.d0*icing_L1_2prec)*3.d0/2.d0/Data%Dt	
	!		
	!		IF (Data%IterJump /= -1 ) TestVelocity = abs((Data%Speed_prec-SpeedJump)/Data%Speed_prec)		! Test on the velocity
	!
	!		
	!		!IF (Data%IterJump /= -1 ) TestVelocity = abs((SpeedJump-Speed_prec)/Speed_prec)
	!		!test on the interface position 
	!		IF (Data%IterJump /= -1 ) TestPos = abs((Icing_L1-(Icing_PrevFixPosition+SpeedJump*Data%Dt))/Icing_L1)  ! test on the front position 
!
!
!			!IF (Data%IterJump /= -1 ) TestPos = abs((Icing_L1-(4.d0/3.d0*Icing_PrevFixPosition-1.d0/3.d0*Icing_L1_2prec+&
!			!2.d0/3.d0*SpeedJump*Data%Dt))/Icing_L1)  ! test on the front position
!		
!		END IF 
!		! For the explicit case only 1 iteration 
!		IF(Data%Space_Type =="EXPLICIT") THEN !the interface isn't moving anymore, explicit resolution only 1 iteration 
!				PRINT*,"Explicit Resolution, not the one taking into account"
!				STOP 
!				!Data%iterJump  = 0 ! we got back to zero and we can move the interface for real 
!				!Icing_L1       = Icing_PrevFixPosition ! initialization to the initial position 
!				!Data%LS_L(:,1) = 2.d0*Icing_L1
!		ELSE IF (Data%IterJump == -1 .OR. Data%iterJump <=1 ) THEN !.OR. interJ < 1.d0  ) THEN !it means first or second iter 
!			! Modif interface position 
!			
!			Icing_L1 = Icing_PrevFixPosition+SpeedJump*Data%Dt
!			
!			!Icing_L1 = 4.d0/3.d0*Icing_PrevFixPosition-1.d0/3.d0*Icing_L1_2prec+2.d0/3.d0*SpeedJump*Data%Dt 
!			Data%LS_L(:,1) = 2.d0*Icing_L1
!			
!			IF(Data%iterJump==-1) THEN 
!				Data%IterJump = 1 
!			ELSE 
!				Data%IterJump = Data%IterJump + 1 
!			END IF 
!		ELSE 	! we have at least done two iterations to compare 	
!			!IF(( TestPos<=TolSousIter .OR. Data%iterJump>=nbMaxITER)) THEN !the interface isn't moving anymore 
!			!IF( (Data%ResNewt<=TolSousIter .AND. TestPos<=TolSousIterPos) .OR. Data%iterJump>=nbMaxITER ) THEN !the interface isn't moving anymore 
!			!IF( (Data%ResNewt<=0.00001 .AND. TestPos<=TolSousIterPos) .OR. Data%iterJump>=nbMaxITER ) THEN !the interface isn't moving anymore 
!			!IF( TestPos<=TolSousIterPos .OR. Data%iterJump>=nbMaxITER ) THEN !the interface isn't moving anymore 
!			IF( testJump_r<=TolSousIterPos .OR. Data%iterJump>=nbMaxITER ) THEN !the interface isn't moving anymore 
!		
!				Data%iterJump=0 ! we got back to zero and we can move the interface for real 
!				
!				!go back to the initial value 
!				Icing_L1 = Icing_PrevFixPosition ! initialization to the initial position 
!				Data%LS_L(:,1) = 2.d0*Icing_L1
!				
!			!ELSE IF(TestPos>TolSousIter) THEN 
!			ELSE !IF(Data%ResNewt>TolSousIter) THEN 
!			
!				Data%iterJump = Data%iterJump+1 
!				Icing_L1 = Icing_PrevFixPosition+SpeedJump*Data%Dt
!				!Icing_L1 = 4.d0/3.d0*Icing_PrevFixPosition-1.d0/3.d0*Icing_L1_2prec+2.d0/3.d0*SpeedJump*Data%Dt 
!				Data%LS_L(:,1) = 2.d0*Icing_L1 
!			END IF
!		END IF 
!			
!	ELSE
!	    PRINT*,"Error in the implementation2 IterativeProcess/UpdateFrontPseudoNewton"
!		STOP
!	END IF
!	IF(Data%Space_Type =="IMPLICIT") THEN
!
!		IF (Data%iterJump/=1) THEN 	
!			PRINT*,"relative error / Velocity : ",TestVelocity,"Position : ",TestPos,"Jump B:",testJump_r,"error Bk", Data%ResNewt
!		END IF 
!		IF(Data%iterJump/=0) THEN 
!			PRINT*,"Iteration tour ", Data%tour,"/ iterJump Newton ",Data%iterJump
!			PRINT*,"Position of the new interface",Icing_L1, "expected : ", truePos,"prev :",Icing_PrevFixPosition
!			PRINT*,"Flux Jump imposed in the scheme", JumpF, "expected : ", Icing_rho*Icing_Lm*trueVelo
!			PRINT*,"Velocity", SpeedJump,"expected : ", trueVelo
!		END IF 
!	END IF 
!
!	!IF iterJump ==0 then we have found a correct guess so no need to go further 
!	IF(Data%iterJump/=0 .AND. Data%Space_Type =="IMPLICIT") THEN !we need to do a pseudo newton iteration 
!
!	!Step 1 : Definition of the new geometry 
!		!Remise à O
!		Mesh%Tab_NoeudComp = 0 ; Mesh%NarEmb = 0.d0 
!		
!		DO i = 1, mesh%Nv
!		   mesh%vertex(i)%Complementaire = 0
!		END DO
!
!		DO ie = 1, mesh%Ne
!			DO i = 1, mesh%ele(ie)%nNodesPerElement
!				Mesh%ele(ie)%surVertex(i) = Mesh%ele(ie)%vertex(i)
!				Mesh%ele(ie)%surPos(i) = Mesh%ele(ie)%Pos(i)
!				Mesh%ele(ie)%SurVertex_Sauvg(i) = 0
!			END DO
!		END DO
!
!		! If embedded, define geometry and tagsInitResNewt
!		! + define number of active elements
!		IF ( Data%Embedded .OR. Data%Immersed ) THEN
!			CALL computeDistance(PhysicalInterface, Data, Mesh, Sol)
!			IF ( Data%Embedded .AND. .NOT. Data%Icing ) CALL tagElement2Solve(Data, Mesh)
!			IF ( Data%Immersed .OR. Data%Icing ) CALL tagElement2SolveImmersed(Data, Mesh)
!			Sol%NeActive = Mesh%NeActive ; Sol%NvActive = Mesh%NvActive
!			Mesh%NarEmb = mesh%NSurB !nombre d'arete sur la surrogate
!		END IF
!
!		IF (Data%Icing .eqv. .TRUE.) THEN
!			! Enregistrement des surVertex
!			DO ie = 1, mesh%Ne
!				DO i = 1, mesh%ele(ie)%nNodesPerElement
!					Nu_i=mesh%ele(ie)%vertex(i)
!					IF (mesh%vertex(Nu_i)%OnInterface .AND. .NOT. Data%Embedded) THEN
!						IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
!						   mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
!						   mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
!						END IF
!					ELSE IF (mesh%vertex(Nu_i)%OnSurrogate .AND.  Data%Embedded) THEN
!						IF(mesh%ele(ie)%Ref==RefZoneAChanger) THEN
!						   mesh%ele(ie)%SurVertex_Sauvg(i) = mesh%ele(ie)%surVertex(i)
!						   mesh%ele(ie)%surVertex(i) = mesh%vertex(Nu_i)%Complementaire
!						END IF
!					END IF
!				END DO
!			END DO
!		END IF
!		mesh%vertex(:)%active = .TRUE.
!		!maj du vecteur temps precedent sur la nouvelle geometry pour T^n
!		CALL AllocSol_PrecModif(Data, Sol_prec, Mesh)
!		!Only to modif the B^k vector only to modify Sol%Bcg_Newton
!		CALL AllocSol_PrecModif_Newton(Data, Sol_prec, Mesh) 
!		! Modification of the jump for the new geometry, but not at the initialisation of the process 
!		Data%Jump = JumpF ! Jump of B at time n+1 on the Surrogate n+1 
!		
!		IF( Data%iterJump/=1 ) THEN
!			CALL JumpModifInter(Data, Sol, Mesh, SpeedJumpNew, JumpNew)
!			IF(JumpNew/=0) THEN !it means the interface didn't move 
!				Data%Jump = JumpNew  ! Jump of B at time n+1 on the Surrogate n+1 
!				PRINT*,"New value of the Jump",Data%Jump 
!		!	ELSE
!		!		PRINT*,"Error JumpNew is null"
!		!		STOP 
!			END IF 
!		END IF 
!
!		IF(DATA%TSCHEME=="NICOLSON" ) THEN
!			JumpFx = 0.d0 ; JumpFy = 0.d0
!			CALL JumpPrevNewInter(Data, Mesh, Sol, Sol_prec,JumpFx,JumpFy) 	
!			Data%PrevJumpSurrogate = JumpFx !Definition of the flux jump value for the previous time on the new geometry  
!		END IF 
!		
!		CALL InitVar(Data, Sol, Mat)
!		CALL FillIndCG(Mat, Mesh)
!
!		CALL RHS_LHS_SchemeSelection(Data, Mesh, Sol, Sol_prec, Sol_prec, Mat)
!		CALL SolveSystem(Data, Sol, Mesh, Mat, "Resolution")
!
!		IF(DATA%TSCHEME=="EULER") THEN 
!			PRINT*,"Need to take into account sol_2prec"
!			STOP 
!			CALL PostProcessSolution(Data, Sol, Sol_prec,Sol_prec, Mesh, Mat) !Separation of the solution into two vectors in Pcg And Bcg 
!		ELSE IF (Data%TSCHEME == "NICOLSON") THEN 
!			Sol_Nicolson = Sol
!			CALL PostProcessSolutionNewton(Data, Sol,Sol_prec,Sol_Nicolson, Mesh, Mat)
!			Sol%BcgN = Sol_Nicolson%Bcg ; Sol%PcgN = Sol_Nicolson%Pcg
!		END IF 
!
!		Data%ResNewt = 0.d0 ; Data%ResNewtD = 0.d0 
!		DO i = 1, Mesh%NarEmb + 1 + mesh%PtEmb 
!			Data%ResNewt  = Data%Resnewt  + (Sol%Bcg(i,1) - Sol_prec%Bcg_Newton(i,1))**2 
!			Data%ResNewtD = Data%ResNewtD +  Sol_prec%Bcg_Newton(i,1)**2 
!		END DO 
!
!		!Data%ResNewt = sqrt(Data%ResNewt/Data%ResNewtD) ! We want this term to go to zero 
!		Data%ResNewt = sqrt(Data%ResNewt/Data%InitResNewt) ! We want this term to go to zero 
!		!Data%ResNewt = sqrt(Data%ResNewt) ! We want this term to go to zero 
!		PRINT*,"Value of the residu of the Newton iter",Data%ResNewt 
!		!We need to check the residu of the term added for the convergence of the pseudo iteration 
!		
!		Data%Icing_Prev_Velocity =  Data%Jump ! for the norm of the Jump at the next iter 
!		Data%Speed_prec =  SpeedJump 
!	END IF 
!

 ! END SUBROUTINE UpdateFrontPseudoNewton
  !========================================================!

!==================================!
  SUBROUTINE JumpModifInter(DataR, SolR, MeshR, SpeedJumpNew, JumpNew)
  !==================================!
    ! modification du vecteur solution au temps n pour etre de la taille du vecteur au temps n+1
    IMPLICIT NONE
    type(DataStructure), INTENT(IN) 	:: DataR
    type(SolStructure) , INTENT(IN) 	:: SolR
    type(MeshType)     , INTENT(IN)     :: MeshR
	REAL*8			   , INTENT(OUT)	:: SpeedJumpNew, JumpNew
    !----------------------------------------
    type(DataStructure)			 :: Data
    type(SolStructure) 			 :: Sol
    type(MeshType)    	     	 :: Mesh
    INTEGER             		 :: nDim, nVar, Ns_prec, Nie , k
    INTEGER             		 :: Nv, Ns, i, NMax, nNodesPerElement, j, Nu_i, Nu_j, Nu_k
    INTEGER             		 :: vertcs, vertcs_comp, nb_comTab, choix, typeCond 
    LOGICAL             	 	 :: test, InIt, InIt1, InIt2,VertcsOnBound
    INTEGER, DIMENSION(2,2)	 	 :: typeV ! on va regarder le type des deux aretes associé 
    REAL*8, DIMENSION(3) 	 	 :: Exa 
    type(Element)            	 :: ele_original  !element etudie

    CHARACTER(len=20)        	 :: BCType1, BCType2
    REAL*8              		 :: val1,val2,val3
    !-----------------------------------------
    REAL*8   , DIMENSION(:,:), ALLOCATABLE :: Sol_inter_Bcg
    REAL*8   , DIMENSION(:)  , ALLOCATABLE :: Sol_inter_Pcg
    !-----------------------------------------

    !  PRINT*, '  --  Allocation  --  '
    Data = DataR ; Sol = SolR ; Mesh = MeshR
    Sol%Bcg = 0.d0 ; Sol%Pcg = 0.d0 
    Sol%Bcg = Sol%Bcg_inter ; Sol%Pcg = Sol%Pcg_inter 
    

	JumpNew = 0.d0
    Nv = Sol%Nv ; nDim = Data%nDim ; nNodesPerElement = nDim + 1
    Mesh%Vertex(:)%EmptyComp = .TRUE. ! pour savoir si une valeur n'a pas encore était reconstruite sur le noeud 
	typeV = 0 

    IF (Data%Icing) THEN
		!ici ce sont des variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb           !nombre de noeud sur l'interface au temps n+1
		Ns_prec=mesh%NarEmb_prec_inter + 1 + mesh%PtEmb 
    ELSE
       Ns=0 ; Ns_prec=0
    END IF
	
    CALL VERIF_Comp(Mesh,test) !on verifie que la surrogate n'est pas la meme qu'au temps precedent
	
	! il suffit de comparer les tableaux avec complementaires
	IF (.NOT. Test) THEN  
		! vecteurs intermediaires pour la partie noeuds complementaire 
		ALLOCATE(Sol_inter_Bcg(Mesh%NB_VertInterf,2)) ; ALLOCATE(Sol_inter_Pcg(Mesh%NB_VertInterf))
		! Mesh%NB_VertInterf REPRESENTE SEULEMENT LE nombre de noeuds avec un complementaire 
		Sol_inter_Bcg = 0.d0 ; Sol_inter_Pcg = 0.d0
		
		!il faut récupérer la valeur du tableau de complémentaire précédent
		DO j= 1, Ns_prec
			vertcs      = Mesh%Tab_NoeudComp_prec_inter(j,1) !ref du noeud original
			vertcs_comp = Mesh%Tab_NoeudComp_prec_inter(j,2) ! ref associe au noeud comp
			choix = 0
			CALL LookForComp(Mesh,InIt,vertcs,nb_comTab,choix) ! on regarde si vertcs a un noeud dans la nvelle table
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
			choix = 2 ! on regarde si il est dans la table precedente
			CALL LookForComp(Mesh,InIt,Nu_i,nb_comTab,choix)
			IF (.NOT. INIT) THEN ! valeurs non modifies

				Sol_inter_Bcg(Nu_j-nv,1) = Sol%Bcg(Nu_i,1)
				Sol_inter_Bcg(Nu_j-Nv,2) = Sol%Bcg(Nu_i,2)
				Sol_inter_Pcg(Nu_j-Nv)   = Sol%Pcg(Nu_i)
				 
				 
				VertcsOnBound = .FALSE.
				CALL VertcsPos(Mesh, Nu_i, VertcsOnBound, typeV)
				IF(.NOT. VertcsOnBound) THEN !cas ou le noeud ne se trouve pas sur le bord  
				
					CALL PolyConstrQR(Data,Mesh,Sol,Nu_i,val1,val2,val3)
					
					Sol%Pcg(Nu_i)   = val1 !valeur reconstruite pour la temperature
					Sol%Bcg(Nu_i,1) = val2 !valeur reconstruite pour le flux en x
					Sol%Bcg(Nu_i,2) = val3 !valeur reconstruite pour le flux en y
					
				ELSE ! deplacement d'un noeud de bord, la valeur de la condition de bord prédomine 

					BCType1 = TRIM(ADJUSTL(Data%Bc(TypeV(1,1)))) ! condition associé à la reference de bord 
					BCType2 = TRIM(ADJUSTL(Data%Bc(TypeV(2,1)))) 
					! en fonction de type de condition de bord on va se 
					IF(BCType1 == BCType2) THEN 
					
						!pour reconstruire les valeurs non données par les conditions de bords 
						CALL PolyConstrQR(Data,Mesh,Sol,Nu_i,val1,val2,val3)
						
						IF (icas==001) THEN !cas specific stefan 
						
							!on recupere la condition sur le fichier de donnée 
							SELECT CASE ( TRIM(ADJUSTL(BCType1)) )
					
							CASE ( "Dirichlet", "dirichlet" )
																
								IF(	TypeV(1,1)==1) THEN					
									Sol%Pcg(Nu_i)   = Icing_Tw  
									Sol%Bcg(Nu_i,1) = val2 
									Sol%Bcg(Nu_i,2) = val3 
								ELSE 
									PRINT*,"error for boundary reconstruction test 001"
									STOP
								END IF 
							 
								CASE ( "Neumann", "neumann" )
							
								SELECT CASE ( TypeV(1,1) )
																								
									CASE(2) !bord de reference 2
									! la condition de Neumann donne -By 

										Sol%Pcg(Nu_i)   = val1 
										Sol%Bcg(Nu_i,1) = val2 
										Sol%Bcg(Nu_i,2) = 0 
									
									CASE(4) !bord de reference 4
									! la condition de Neumann donne By 

										Sol%Pcg(Nu_i)   = val1  
										Sol%Bcg(Nu_i,1) = val2 
										Sol%Bcg(Nu_i,2) = 0 
									
									CASE DEFAULT
									   PRINT*, " ** ERROR inAllocSol_PrecModif case 001"
									   STOP
									END SELECT
											
							CASE DEFAULT
							   PRINT*, " ** ERROR inAllocSol_PrecModif case 001"
							   STOP
							END SELECT	
					   
							
						ELSE 
							PRINT*,"Error JumpModifInter"
							STOP 
						END IF 
					ELSE 
						PRINT*,"Error JumpModifInter"
						STOP 
					END IF ! fin de la reconstruction ou il manquait une valeur 
				END IF 
			END IF  
			! cas avec interface sur le bord 
			IF (INIT .AND. Nu_i/= Nu_j .AND. Nu_i == nb_comTab ) THEN !i;e que l'on avait un noeud de bord sans complementaire 
				PRINT*,"Error JumpModifInter"
				STOP 
		    END IF 	 ! fin (INIT .AND. Nu_i/= Nu_j .AND. Nu_i == nb_comTab ) 
	    END DO
	      
		!la partie Sol%Bcg(Nv,:) a deja etait modifié de meme pour Sol%Pcg(Nv)
		Sol%Bcg(Nv+1:Nv+Mesh%NB_VertInterf,1)= Sol_inter_Bcg(:,1)
		Sol%Bcg(Nv+1:Nv+Mesh%NB_VertInterf,2)= Sol_inter_Bcg(:,2)
		Sol%Pcg(Nv+1:Nv+Mesh%NB_VertInterf)= Sol_inter_Pcg(:)
		  
		!remise a 0 de la partie non utilisé
		Sol%Bcg(Nv + Mesh%NB_VertInterf + 1 : Nv + data%MaxRowEmbedded ,:) = 0.d0
		Sol%Pcg(Nv + Mesh%NB_VertInterf + 1 : Nv + data%MaxRowEmbedded )   = 0.d0
		   
		CALL MeanFluxReconstruct22(Data, Mesh, Sol, SpeedJumpNew, JumpNew) 
					
		DEALLOCATE(Sol_inter_Bcg) ; DEALLOCATE(Sol_inter_Pcg)
    END IF
			
  END SUBROUTINE JumpModifInter
  !========================================!

END MODULE IterativeProcess
