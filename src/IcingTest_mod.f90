MODULE IcingTest_mod

  USE Types
  USE PublicVar
  USE Locate

  IMPLICIT NONE

CONTAINS



  !====================================!
  SUBROUTINE H_enthalpy(ele, temp, H)
  !====================================!

  !*************************************************************************
  !   Define the value of the enthalphy depending of the mesh ref
  !*************************************************************************

	!attention voir pour la temperature au noeud souhaité 
    IMPLICIT NONE
    REAL*8       , INTENT(OUT)  :: h
    type(Element), INTENT(IN)   :: ele
    REAL*8       , INTENT(IN)   :: temp ! temperature sur le noeud considéré 

    !----------------------------------------

    
    IF(ele%ref ==1) THEN 
    
      H = icing_c2*(temp - Icing_Tm)

    ELSE ! ele%ref==2 
		H = Icing_c1*(temp - Icing_Tm) + Icing_Lm 
    END IF 


  END SUBROUTINE H_enthalpy
  !====================================!


  !====================================!
  SUBROUTINE L_function(ele, temp, lambda_icing )
  !====================================!

  !*************************************************************************
  !   Lambda function
  !*************************************************************************

	
    IMPLICIT NONE
    REAL*8       , INTENT(OUT)  :: lambda_icing 
    type(Element), INTENT(IN)   :: ele
    REAL*8       , INTENT(IN)   :: temp ! temperature sur le noeud considéré 

    !----------------------------------------
    REAL*8                      :: h

    
   CALL H_enthalpy(ele, temp, H)
   
   IF ( H<0) THEN 
		Lambda_icing = Icing_Lambda2
   ELSE IF (H> Icing_Lm) THEN 
		Lambda_icing = Icing_Lambda1
   ELSE 
		Lambda_icing = (H/icing_Lm)/Icing_Lambda1 + (1 - (H/Icing_Lm))/Icing_Lambda2
   END IF 
   

  END SUBROUTINE L_function
  !====================================!
  
  
  !====================================!
  SUBROUTINE Interface_deplacement(Data)
  !====================================!

  !*************************************************************************
  !  Modify the date values for the new position of the interface before the new time resolution 
  !*************************************************************************
  
   IMPLICIT NONE
    type(DataStructure), INTENT(INOUT) :: Data
    !----------------------------------------

	INTEGER  :: i 
	REAL*8   :: chi , t , alpha_l, pos_init
  
    !----------------------------------------
    
    t = Data%T ; alpha_l = Icing_Lambda1/(Icing_rho*Icing_c1)
    

	IF( icas/=001) THEN ! i.E qu'on prend les exemple précédents avec polynome comme fonction analytique 
	   IF (Data%SpeedCst) THEN ! si on a une interface mobile suivant une constante
			 SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(1,1))) )
			 CASE ( "Circle","circle","Sphere","sphere" )
				Icing_L3 = Icing_L3 + Data%Speed*Data%Dt  ! modification du rayon, fichier Icing
				DO i = 1, Data%N_LS
				   Data%LS_L(i,1) = Data%LS_L(i,1) + 2.d0*Data%Speed*Data%Dt  ! modification du rayon, fichier Embedded
				END DO
			 CASE ( "Box","box" )
				Icing_L1 = Icing_L1 + Data%Speed*Data%Dt  ! modification de la position en abscisse, fichier Icing
				DO i = 1, Data%N_LS
				   Data%LS_L(i,1) = Data%LS_L(i,1) + 2.d0*Data%Speed*Data%Dt  ! modification de la position en abscisse, fichier Embedded
				END DO
				PRINT*,"New Interface position at time", data%T," : ", Icing_L1
			 END SELECT
		END IF
	ELSE ! i.e. que l'on fait le cas test avec partie glace et eau 	
	!used only at the initialization step to define the initial interface position 
		SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(1,1))) )
		CASE ( "Circle","circle","Sphere","sphere" )
			PRINT*,"Error test case not available on a circle"
			STOP
		CASE ( "Box","box" )
			Icing_L1 = 2.d0*Icing_chi*sqrt(alpha_l*t) 
			Data%LS_L(:,1) = 2.d0*Icing_L1
		END SELECT
	END IF 
      
  END SUBROUTINE Interface_deplacement 
  !====================================!
  
  
	!====================================!
	SUBROUTINE Interface_Function(x,res1,res2)
	!====================================!

	!*************************************************************************
	!  Equation from the one needed to solve to get the coeffcient defined for the interface position 
	!*************************************************************************
  
	IMPLICIT NONE
	
	REAL*8  ,INTENT(IN) :: x
	REAL*8  ,INTENT(OUT):: res1,res2
    !---------------------------------------
	INTEGER  :: i   
	REAL*8   :: nu 
    !----------------------------------------
    nu = sqrt(Icing_alpha_l/Icing_alpha_s)
    
    !calcul de la fonction
    res1  = Stw/(exp(x**2)*erf(x)) - Sti/(nu*exp(nu**2*x**2)*erfc(nu*x)) -x*sqrt(PI) 
    
    res2 = -2*Stw*((x*exp(x**2)*erf(x) +(1.d0/sqrt(PI)))/(exp(x**2)**2*erf(x)**2))-sqrt(PI) &
    + 2.d0*Sti*((nu*x*exp(nu**2*x**2)*erfc(nu*x)-(1.d0/sqrt(PI)))/(exp(nu**2*x**2)**2*erfc(nu*x)**2))
	!calcul de sa dérivée 
	
      
	END SUBROUTINE Interface_Function
	!====================================!
  
  
  !======================================
  SUBROUTINE Newton_Interface(x0,eps) 
  !======================================
  
  !********************************
  !newthod method to find the coeffcient chi in the interface definition 
  !**********************************************************************
  
  IMPLICIT NONE 
  
	REAL*8,		INTENT(INOUT)	:: x0
	REAL*8,		INTENT(IN)		:: eps 
	!-------------------------------------------------------
	REAL*8 						:: res1,res2
	LOGICAL 					:: ARRET 
	!--------------------------------------------------------
	
	res1 = 0.d0 ; res2 = 0.d0 
	ARRET = .FALSE. 
	
	DO WHILE (.NOT. ARRET ) 
		
		res1 = 0.d0 !value of the function
		res2 = 0.d0 !value of the derivative 
		CALL Interface_Function(x0,res1,res2)
		IF(res1<eps) THEN !we got the right x0
			ARRET = .TRUE. 
		ELSE 
			!we verify that the derivative isn't null 
			IF(abs(res2)<= 1E-12) THEN 
				PRINT*,"Error in Newton_Interface the derivative is null"
				STOP 
			ElSE 
				x0 = x0 -res1/res2
			END IF 
		END IF 
	END DO 
	
    
	END SUBROUTINE Newton_Interface
	!==============================
  
  
  !===================================
  SUBROUTINE Coeff_hte(Data) 
  
  IMPLICIT NONE 
      type(DataStructure), INTENT(IN) :: Data
      !--------------------------------------------------
      
      REAL*8	:: res 



	res = Icing_Lambda2/erfc(Icing_L2/(2.d0*sqrt(Icing_alpha_s*Data%t)))&
	*1.d0/sqrt(PI*Icing_alpha_s*Data%T)*exp(-Icing_L2**2/(4.d0*Icing_alpha_s*Data%t))
  
	Icing_hte = res 
	PRINT*,"Voici la valeur de hte calculée",res 
	
  
  END SUBROUTINE Coeff_hte 
  !===================================
  
  !===================================
	SUBROUTINE MeanFlux(Data,Mesh, Sol, JumpF,SpeedJump) 
  
	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: JumpF,SpeedJump

    !-------------------------------------------------- 
    INTEGER :: Ne, Nv, nDim, nNodesPerElement, i
    INTEGER :: Nu_i, Nu_j, nbJumpf, Ns     
    REAL*8	:: analy_velocity,std_deviation
    CHARACTER(len=20) :: h_value
    !-------------------------------------------------- 
      
      
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim ; std_deviation = 0.d0 
    nNodesPerElement = nDim + 1 ; JumpF = 0.d0 ; nbJumpf = 0
	analy_velocity = Icing_chi*sqrt(Icing_Lambda1/(Icing_rho*Icing_c1))/sqrt(Data%T-Data%Dt) !analytical velocity 
	WRITE (h_value,'(F7.5)') Data%h
		
	IF (Data%Icing) THEN
		!ici ce sont des variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb           !nombre de noeud sur l'interface au temps n+1
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    
    !Condition sur les aretes de la surrogate
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ! ref mesh
		Nu_j = mesh%Tab_NoeudComp(i,2) ! ref complementary

		IF(icas ==001) THEN 
			!2d test but with a 1d displacement for the two points of the edge 
			IF((Sol%Bcg(Nu_i,1) - Sol%Bcg(Nu_j,1))<0) THEN 
				PRINT*,"Problem of sign in IcingTest_mod\MeanFlux function"
				PRINT*,"Means that the previous approximation is bad" 
			ELSE			
				JumpF = JumpF + (Sol%Bcg(Nu_i,1) - Sol%Bcg(Nu_j,1)) !we only look at the difference in x since 1d case 
				nbJumpf = nbJumpf + 1
			END IF 

		ELSE IF(nDim == 2) THEN 
			!!write something 
			PRINT*,"Not implemented yet look at MeanFlux"
			STOP 
		ELSE 
			PRINT*,"3D case not available in MeanFlux"
			STOP 
		END IF 

    END DO

	IF (nbJumpf==0) THEN 
		PRINT*,"Error the number of jump at the interface is null"
	ELSE 
		JumpF = JumpF/nbJumpf !mean value 
		SpeedJump = JumpF/(Icing_rho*Icing_Lm) ! mean velocity
	END IF 
		
  
  END SUBROUTINE MeanFlux
  !===================================
  
    
  !===================================
	SUBROUTINE ValJumpT(Data, Mesh, Sol, edg, Jump, MeanT, dis, k, lambdaM, lambdaP) 
  
	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 			:: Data
	type(MeshType)     , INTENT(IN) 			:: Mesh
	type(SolStructure) , INTENT(IN) 			:: Sol
	REAL*8 			   , INTENT(OUT)			:: Jump,MeanT
    type(Edge)         , INTENT(IN)    			:: Edg
    REAL*8, DIMENSION(data%nDim), INTENT(IN)    :: dis
	INTEGER 		   ,INTENT(IN)				:: k
	REAL*8			   ,INTENT(IN) 			    :: lambdaM,lambdaP

    !-------------------------------------------------- 
    INTEGER :: nDim, nNodesPerElement, i,j 
    INTEGER :: Nu_j1, Nu_j2, nbJumpf, Ns ,Nu_k1,NU_k2   
    LOGICAL :: FOUND
    INTEGER, DIMENSION(:), ALLOCATABLE :: NU 
    REAL*8, DIMENSION(2) :: dis2 
    REAL*8  :: Temp_Left, Temp_Right, Bx_Left, Bx_Right, By_Left, By_Right
    !-------------------------------------------------- 
      
    dis2 = dis
    nDim = Data%nDim ; nNodesPerElement = nDim + 1 ; Jump = 0.d0 ; MeanT = 0.d0 
    
    ALLOCATE(NU(2))
	Nu(1) = edg%vertex(1) ; Nu(2) = edg%vertex(2)
		
	IF (Data%Icing) THEN
		!ici ce sont des variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb           !nombre de noeud sur l'interface au temps n+1
    ELSE
       PRINT*,"We have a problem here Not icing case, Icing_testMod/ValJumpT"
       STOP
    END IF
    
    IF(k==1) THEN 
    
		!on veut trouver ça position dans le tableau de complémentaire 
		FOUND = .FALSE. 
		j= 0
		DO WHILE (.NOT. FOUND)
		j = j+1
			IF(j>Ns) THEN 		
				PRINT*,"ERROR in Icing_testMod/ValJumpT "
				STOP 
			END IF 
			IF(mesh%Tab_NoeudComp(j,1) == Nu(1)) THEN 
				FOUND = .TRUE. 
				Nu_k1 = mesh%Tab_NoeudComp(j,1) ! ref mesh
				Nu_k2 = mesh%Tab_NoeudComp(j,2) ! ref complementary
			END IF 
		END DO  
			 
		IF(nDim == 2) THEN 
			MeanT  = Sol%Pcg(Nu_k1) -1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis2(1) -1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis2(2)
			MeanT = MeanT + Sol%Pcg(Nu_k2) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis2(1) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis2(2)
			MeanT = MeanT*0.5d0
		
			Jump = (Sol%Pcg(Nu_k1) - Sol%Pcg(Nu_k2)) 
			Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis2(1) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis2(1))
			Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis2(2) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis2(2))
		ELSE 
			PRINT*,"3D case not available in Icing_testMod/ValJumpT"
			STOP 
		END IF 

	ELSE IF(k==2) THEN 
	
	! computation of the middle value 
		
	!	FOUND = .FALSE. 
	!	j= 0
	!	DO WHILE (.NOT. FOUND)
	!	j = j+1
	!	IF(j>Ns) THEN 		
	!			PRINT*,"ERROR in Icing_testMod/ValJumpT "
	!			STOP 
	!		END IF 
	!		IF(mesh%Tab_NoeudComp(j,1) == Nu(1)) THEN 
	!			FOUND = .TRUE. 
	!			Nu_k1 = mesh%Tab_NoeudComp(j,1) ! ref mesh
	!			Nu_k2 = mesh%Tab_NoeudComp(j,2) ! ref complementary
	!		END IF 
	!	END DO  
	!		 
	!	IF(nDim == 2) THEN 
	!	
	!		MeanT  = Sol%Pcg(Nu_k1) -1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis(1) -1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis(2)
	!		MeanT = MeanT + Sol%Pcg(Nu_k2) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis(1) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis(2)
	!		MeanT = MeanT*0.5d0
	!		
	!		Jump = (Sol%Pcg(Nu_k1) - Sol%Pcg(Nu_k2))
	!		Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis(1) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis(1))
	!		Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis(2) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis(2))
	!	ELSE 
	!		PRINT*,"3D case not available in Icing_testMod/ValJumpT"
	!		STOP 
	!	END IF 
	!	
	!	FOUND = .FALSE. 
	!	j= 0
	!	DO WHILE (.NOT. FOUND)
	!	j = j+1
	!		IF(j>Ns) THEN 		
	!			PRINT*,"ERROR in Icing_testMod/ValJumpT "
	!			STOP 
	!		END IF 
	!		IF(mesh%Tab_NoeudComp(j,1) == Nu(2)) THEN 
	!			FOUND = .TRUE. 
	!			Nu_k1 = mesh%Tab_NoeudComp(j,1) ! ref mesh
	!			Nu_k2 = mesh%Tab_NoeudComp(j,2) ! ref complementary
	!		END IF 
	!	END DO  
	!		 
	!	IF(nDim == 2) THEN 
	!	
	!		MeanT = MeanT + Sol%Pcg(Nu_k1) - 1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis(1) - 1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis(2)
	!		MeanT = MeanT + Sol%Pcg(Nu_k2) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis(1) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis(2)
	!		MeanT = MeanT*0.5d0
	!		
	!		Jump = Jump + (Sol%Pcg(Nu_k1) - Sol%Pcg(Nu_k2)) 
	!		Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis(1) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis(1))
	!		Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis(2) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis(2))
	!	ELSE 
	!		PRINT*,"3D case not available in Icing_testMod/ValJumpT"
	!		STOP 
	!	END IF 
		
	!	Jump  = Jump*0.5d0 !we take the mean of the value at the vertices 
	!	MeanT = MeanT*0.5d0
	
	! computation of the middle value 
		
		FOUND = .FALSE. 
		j= 0
		DO WHILE (.NOT. FOUND)
		j = j+1
		IF(j>Ns) THEN 		
				PRINT*,"ERROR in Icing_testMod/ValJumpT "
				STOP 
			END IF 
			IF(mesh%Tab_NoeudComp(j,1) == Nu(1)) THEN 
				FOUND = .TRUE. 
				Nu_k1 = mesh%Tab_NoeudComp(j,1) ! ref mesh
				Nu_k2 = mesh%Tab_NoeudComp(j,2) ! ref complementary
			END IF 
		END DO  
		
		FOUND = .FALSE. 
		j= 0
		DO WHILE (.NOT. FOUND)
		j = j+1
			IF(j>Ns) THEN 		
				PRINT*,"ERROR in Icing_testMod/ValJumpT "
				STOP 
			END IF 
			IF(mesh%Tab_NoeudComp(j,1) == Nu(2)) THEN 
				FOUND = .TRUE. 
				Nu_j1 = mesh%Tab_NoeudComp(j,1) ! ref mesh
				Nu_j2 = mesh%Tab_NoeudComp(j,2) ! ref complementary
			END IF 
		END DO  

			 
		IF(nDim == 2) THEN 
		
			Temp_Left = (Sol%Pcg(Nu_k1)+Sol%Pcg(Nu_j1))*0.5
			Bx_Left   = (Sol%Bcg(Nu_k1,1)+Sol%Bcg(Nu_j1,1))*0.5
			By_Left   = (Sol%Bcg(Nu_k1,2)+Sol%Bcg(Nu_j1,2))*0.5
			
			Temp_Right = (Sol%Pcg(Nu_k2)+Sol%Pcg(Nu_j2))*0.5
			Bx_Right   = (Sol%Bcg(Nu_k2,1)+Sol%Bcg(Nu_j2,1))*0.5
			By_Right   = (Sol%Bcg(Nu_k2,2)+Sol%Bcg(Nu_j2,2))*0.5
			

			MeanT  = (Temp_Left+Temp_Right)*0.5 - (1.d0/lambdaP*Bx_Left+1.d0/lambdaM*Bx_Right)*dis2(1)*0.5 &
			- (1.d0/lambdaP*By_Left+1.d0/lambdaM*By_Right)*0.5*dis2(2)
			
			Jump  = (Temp_Left - Temp_Right) - (1.d0/lambdaP*Bx_Left - 1.d0/lambdaM*Bx_Right)*dis2(1) &
			- (1.d0/lambdaP*By_Left - 1.d0/lambdaM*By_Right)*dis2(2)

		ELSE 
			PRINT*,"3D case not available in Icing_testMod/ValJumpT"
			STOP 
		END IF 
			
		
	ELSE 
	
		FOUND = .FALSE. 
		j= 0
		DO WHILE (.NOT. FOUND)
		j = j+1
			IF(j>Ns) THEN 		
				PRINT*,"ERROR in Icing_testMod/ValJumpT "
				STOP 
			END IF 
			IF(mesh%Tab_NoeudComp(j,1) == Nu(2)) THEN 
				FOUND = .TRUE. 
				Nu_k1 = mesh%Tab_NoeudComp(j,1) ! ref mesh
				Nu_k2 = mesh%Tab_NoeudComp(j,2) ! ref complementary
			END IF 
		END DO  
			 
		IF(nDim == 2) THEN 
		
			MeanT  = Sol%Pcg(Nu_k1) -1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis2(1) -1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis2(2)
			MeanT = MeanT + Sol%Pcg(Nu_k2) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis2(1) - 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis2(2)
			MeanT = MeanT*0.5d0
			
			Jump = (Sol%Pcg(Nu_k1) - Sol%Pcg(Nu_k2)) 
			Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,1)*dis2(1) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,1)*dis2(1))
			Jump = Jump +  (-1.d0/lambdaP*Sol%Bcg(Nu_k1,2)*dis2(2) + 1.d0/lambdaM*Sol%Bcg(Nu_k2,2)*dis2(2))
		ELSE 
			PRINT*,"3D case not available in Icing_testMod/ValJumpT"
			STOP 
		END IF 
	
	END IF
  
  END SUBROUTINE ValJumpT
  !===================================
  
      
  !===================================
	SUBROUTINE MeanTemp(Data, Mesh, Sol, MeanT) 
  
	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 			:: Data
	type(MeshType)     , INTENT(IN) 			:: Mesh
	type(SolStructure) , INTENT(IN) 			:: Sol
	REAL*8 			   , INTENT(OUT)			:: MeanT
    !-------------------------------------------------- 
    INTEGER :: nDim, i
    INTEGER :: Ns ,Nu_k1,NU_k2   
    REAL*8, DIMENSION(2) :: dis 
    !-------------------------------------------------- 
      
    
    nDim = Data%nDim ; MeanT = 0.d0 
    
	IF (Data%Icing) THEN
		!ici ce sont des variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb           !nombre de noeud sur l'interface au temps n+1
    ELSE
       PRINT*,"We have a problem here Not icing case, Icing_testMod/ValJumpT"
       STOP
    END IF
    
    DO i=1, Ns 
    
	Nu_k1 = mesh%Tab_NoeudComp(i,1) ! ref mesh
	Nu_k2 = mesh%Tab_NoeudComp(i,2) ! ref complementary
	
	dis(1) = Sol%dist(1,Nu_k1,1) 
	dis(2) = Sol%dist(1,Nu_k1,2)  
	PRINT*,"dis",dis(1),dis(2)
		 
	MeanT = MeanT + (Sol%Pcg(Nu_k1) - 1.d0/Icing_Lambda1*Sol%Bcg(Nu_k1,1)*dis(1) - 1.d0/Icing_Lambda1*Sol%Bcg(Nu_k1,2)*dis(2))*0.5d0
	MeanT = MeanT + (Sol%Pcg(Nu_k2) - 1.d0/Icing_Lambda2*Sol%Bcg(Nu_k2,1)*dis(1) - 1.d0/Icing_Lambda2*Sol%Bcg(Nu_k2,2)*dis(2))*0.5d0

	END DO 
	
	MeanT = MeanT/Ns 
	PRINT*,"Value of the Mean Temperature at the interface",MeanT
  
  END SUBROUTINE MeanTemp
  !===================================

END MODULE IcingTest_mod
