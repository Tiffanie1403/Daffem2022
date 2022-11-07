 MODULE Velocity

  USE Types
  USE PublicVar
  USE Locate
  USE Interpolation 
  
  USE Algebra
  USE Algebra2
  USE Data_mod
  USE ReadMesh_Mod
  USE GMRES 
  USE IcingTest_mod 
  
  
  IMPLICIT NONE

CONTAINS

	!===================================
	SUBROUTINE MeanFluxReconstruct(Data, Mesh, Sol, JumpF) 
  
	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: JumpF

    !-------------------------------------------------- 
    INTEGER 			  :: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			  :: Nu_i, Nu_j, nbJumpf, Ns, NU_k 
    REAL*8				  :: analy_velocity,interpSol, Beta_RightSide, Beta_LeftSide,JumpReconstruct
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               :: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL               :: found 
    REAL*8, ALLOCATABLE, DIMENSION(:)  :: SolTriangle, SolTriangley, n_interface
	REAL*8				 :: previous_JumpF, distance,interpSoly,norm,Beta_Left_True, Beta_Right_True
	CHARACTER(len=20) 	 :: h_value
	REAL*8				 :: JumpReconstruct_test1, JumpReconstruct_test2 ,alpha_l,alpha_s		

    !-------------------------------------------------- 
     
    ALLOCATE(X_cur(2),SolTriangle(3),SolTriangley(3),n_interface(2))
    SolTriangle =0.d0 ; previous_JumpF = 0.d0 ; SolTriangley =0.d0
    n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 
    Beta_RightSide = 0.d0 ; Beta_LeftSide = 0.d0 ; distance = 0.d0 
    JumpF = 0.d0 ; nbJumpf = 0 ; interpSol =0.d0
	alpha_l = Icing_Lambda1/(Icing_rho*Icing_c1)
	alpha_s = Icing_Lambda2/(Icing_rho*Icing_c2)


	WRITE (h_value,'(F7.5)') Data%h

	 
    !analytical velocity 
	analy_velocity = Icing_chi*sqrt(Icing_Lambda1/(Icing_rho*Icing_c1))/sqrt(Data%T-Data%Dt) 
		
	IF (Data%Icing) THEN
		!variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    
    !Condition sur les aretes de la surrogate
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 

		!Calculus of \beta^- by interpolation
		
		!First :Defining the point on the true interface
		
		X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! position of the node on the true interface 
		
		!Second :Finding the triangle where the node is 
		Nele = mesh%vertex(Nu_i)%ele1 !element of departure for the research 
		CALL LocateNode(X_cur, mesh, Nele, found)
		
		!Third: Calculus of \beta^- by interpolation
		! warning : we use the complementary value for the interpolation , interpolation only for the x coordinate 
		! Filling of the SolTRiangle Vector 
		SolTriangle = 0.d0 !; SolTriangley = 0.d0
		DO ik=1,3 
		!we need to check if the node is on the surrogate to take the compementary value in the case it is
			Nu_k = Mesh%ele(Nele)%Vertex(ik)
			IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
				NU_k = Mesh%Vertex(Nu_k)%Complementaire
			END IF 
			SolTriangle(ik) = Sol%Bcg(Nu_k,1)
		END DO 
		
		!interpolation of the left side for \beta^-
		interpSol = 0.d0
		CALL InterpolateSolution(X_cur, interpSol, SolTriangle, Mesh%ele(nele))
		Beta_RightSide =  interpSol

		!The next step is to compute the left side \beta^+ by extrapolation, only the x coordinate is considered
		CALL ConstrQR_Beta_LeftSide(Data, Mesh, Sol, X_cur, Beta_LeftSide)
		
		JumpReconstruct = Beta_LeftSide - Beta_RightSide

		
		IF(icas ==001) THEN 
			IF(JumpReconstruct <0) THEN 
				PRINT*,"Problem of sign in Velocity\MeanFluxReconstruct function"
			ELSE			
				JumpF = JumpF + JumpReconstruct
				nbJumpf = nbJumpf + 1
			END IF 

		ELSE IF(nDim == 2) THEN 
			!!write something 
			PRINT*,"Not implemented yet look at MeanFluxReconstruct"
			STOP 
		ELSE 
			PRINT*,"3D case not available in MeanFluxReconstruct"
			STOP 
		END IF 

    END DO

	IF (nbJumpf==0) THEN 
		PRINT*,"Error the number of jump at the interface is null"
	ELSE 
		JumpF = JumpF/nbJumpf !mean value 
		JumpF = JumpF/(Icing_rho*Icing_Lm) ! mean velocity		
		!PRINT*,"Numerical velocity",JumpF,"Analytical velocity",analy_velocity!,"time",Data%T

		!pour arreter le code a t environ 3 secondes
		!IF(Data%T-Data%Dt<3 .AND. Data%T > 3) THEN			
		!	DO i = 1, Ns 
         ! 	 
			!	Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2)

			!	X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! position of the node on the true interface 
				
			!	distance =sqrt((mesh%Vertex(Nu_i)%X(1)-X_cur(1))**2+(mesh%Vertex(Nu_i)%X(2)-X_cur(2))**2)
				
			!	Nele = mesh%vertex(Nu_i)%ele1 !element of departure for the research 
			!	CALL LocateNode(X_cur, mesh, Nele, found)
				
			!	SolTriangle =0.d0 
			!	DO ik=1,3 
			!		Nu_k = Mesh%ele(Nele)%Vertex(ik)
			!		IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
			!			NU_k = Mesh%Vertex(Nu_k)%Complementaire
			!		END IF 
			!		SolTriangle(ik)  = Sol%Bcg(Nu_k,1) !x coordinate
			!	END DO 
				
			!	interpSol = 0.d0
			!	CALL InterpolateSolution(X_cur, interpSol, SolTriangle, Mesh%ele(nele))
			!	Beta_RightSide =  interpSol
				
			!	CALL ConstrQR_Beta_LeftSide(Data, Mesh, Sol, X_cur, Beta_LeftSide)
			!	!real reconstruction
			!	JumpReconstruct = Beta_LeftSide - Beta_RightSide
			!	!real left value

			!	
			!	IF(JumpReconstruct<0) THEN 
			!	PRINT*,"Problem of sign in IcingTest_mod\MeanFlux function"
			!	PRINT*,"Means that the previous approximation is bad" 
			!	ELSE
			!		OPEN(unit = 77,file='speed_compare_defnew'//TRIM(h_value)//'.data',action='write',position='append')
			!		WRITE(77,*)Nu_i, mesh%Vertex(Nu_i)%X(1), mesh%Vertex(Nu_i)%X(2), Sol%dist(1,Nu_i,1), &
			!		JumpReconstruct/(Icing_rho*Icing_Lm), JumpF, analy_velocity
			!		CLOSE(77)
!
!				END IF 
!
!			END DO
!			PRINT*,"Writen"
!			STOP
!		END IF

		!OPEN(unit = 177,file='output_velocity_position_enriched_'//TRIM(h_value)//'.data',action='write',position='append')
		!WRITE(177,*)Data%T-Data%Dt,JumpF,analy_velocity,Icing_L1+JumpF,2*Icing_chi*sqrt(Icing_alpha_l*(Data%T-Data%Dt))
		!CLOSE(177)
	
		JumpF = JumpF*Data%Dt ! we already add Data%Dt to update the position of the front 
	END IF 
	
  
  END SUBROUTINE MeanFluxReconstruct
  !===================================
  
  
  !send the jump value and the velocity in two different variable 
  !===================================
	SUBROUTINE MeanFluxReconstruct22(Data, Mesh, Sol,SpeedJump,JumpF)
  	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: JumpF,SpeedJump

    !-------------------------------------------------- 
    INTEGER 			  				:: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			 				:: Nu_i, Nu_j, nbJumpf, Ns, Nu_k 
    REAL*8				  				:: analy_velocity, Beta_RightSide, Beta_LeftSide,JumpReconstruct
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               				:: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL              				:: found 
    REAL*8, ALLOCATABLE, DIMENSION(:) 	:: SolTriangle, n_interface
	REAL*8				 				:: alpha_l,alpha_s,res1,res2 
	REAL*8, DIMENSION(2) 				:: nP 
	CHARACTER(len=99) 					:: charTmp1, charTmp2,charTmp3
    !-------------------------------------------------- 
     
    WRITE(charTmp2,*) Data%iterJump
    WRITE(charTmp1,*) Data%Tour
    WRITE(charTmp3,*) Data%TSCHEME

    !open(unit = 331,file='MeanValue'//TRIM(ADJUSTL(charTmp1))//TRIM(ADJUSTL(charTmp2))//&
    !TRIM(ADJUSTL(charTmp3))//'.data', action="write",position="append")

    ALLOCATE(X_cur(2),SolTriangle(3),n_interface(2))
    SolTriangle =0.d0 ;  n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 ; Beta_RightSide = 0.d0 ; Beta_LeftSide = 0.d0 
    JumpF = 0.d0 ; nbJumpf = 0 ; SpeedJump = 0.d0 ; JumpReconstruct = 0.d0 
	res1 = (Icing_Tw-Icing_Tm)/erf(Icing_chi) 
	res2 = (Icing_Tm-Icing_Tinit)/erfc(Icing_chi*sqrt(Icing_alpha_l/Icing_alpha_s))
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    DO i = 1, Ns 
          	 
		Nu_i  = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
		X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! Sol_prec est definis via Sol en post processing 
		np    = Sol%dist(1,Nu_i,:)/sqrt(Sol%dist(1,Nu_i,1)**2+Sol%dist(1,Nu_i,2)**2) ! normal vector at the interface, nout used for one direction displacement 

		Nele  = mesh%vertex(Nu_i)%ele1 
		CALL LocateNode(X_cur, mesh, Nele, found)

		SolTriangle = 0.d0 
		DO ik=1,3 
			Nu_k = Mesh%ele(Nele)%Vertex(ik)
			IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
				NU_k = Mesh%Vertex(Nu_k)%Complementaire
			END IF 
			SolTriangle(ik) = Sol%Bcg(Nu_k,1)
		END DO 
		
		CALL InterpolateSolution(X_cur, Beta_RightSide, SolTriangle, Mesh%ele(nele))		
	
		CALL ConstrQR_Beta_LeftSide(Data, Mesh, Sol, X_cur, Beta_LeftSide)
				
		JumpReconstruct = Beta_LeftSide - Beta_RightSide

		!WRITE(331,*)X_cur(1),X_cur(2),JumpReconstruct

		IF(icas ==001) THEN 
			IF(JumpReconstruct <0.d0) THEN 
				PRINT*,"Problem of sign in Velocity\MeanFluxReconstruct function"
			ELSE			
				JumpF = JumpF + JumpReconstruct ! jump on the true interface 
				nbJumpf = nbJumpf + 1
			END IF 

		ELSE IF(nDim == 2) THEN 
			PRINT*,"Not implemented yet look at MeanFluxReconstruct22"
			STOP 
		ELSE 
			PRINT*,"3D case not available in MeanFluxReconstruct22"
			STOP 
		END IF 

    END DO

	IF (nbJumpf==0) THEN 
		PRINT*,"Error the number of jump at the interface is null"
	ELSE 
		JumpF  = JumpF/nbJumpf !mean jump value 
		SpeedJump = JumpF/(Icing_rho*Icing_Lm) ! value of the velocity with the value of the gradient to 
		PRINT*,"JumpF", JumpF, "SpeedJump", SpeedJump 
	END IF 
	
	!WRITE(331,*)"Mean Value",JumpF,"True Value",Icing_rho*Icing_Lm*Icing_chi*sqrt(Icing_alpha_l/Data%t)

	!CLOSE(331)
	
  END SUBROUTINE MeanFluxReconstruct22
  !============================================
  
  
   !========================================================================
	SUBROUTINE MeanFluxReconstructNodeVal(Data, Mesh, Sol, SpeedJumpx, SpeedJumpy, X)
  	IMPLICIT NONE 
  
	type(DataStructure)			   , INTENT(IN)  :: Data
	type(MeshType)    			   , INTENT(IN)  :: Mesh
	type(SolStructure) 			   , INTENT(IN)  :: Sol
	REAL*8  		  			   , INTENT(OUT) :: SpeedJumpx, SpeedJumpy
	REAL*8 , DIMENSION(2) 		   , INTENT(IN)	 :: X
    !-------------------------------------------------- 
    INTEGER 			  				:: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			 				:: Nu_i, Nu_j, nbJumpf, Ns, Nu_k 
    REAL*8				  				:: analy_velocity, Beta_RightSide_x,Beta_RightSide_y 
    REAL*8								:: Beta_LeftSide_x, Beta_LeftSide_y, JumpReconstructX, JumpReconstructY
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               				:: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL              				:: found 
    REAL*8, ALLOCATABLE, DIMENSION(:) 	:: SolTriangleX, SolTriangleY, n_interface
	REAL*8				 				:: alpha_l,alpha_s,res1,res2 
	REAL*8, DIMENSION(2) 				:: nP 
	CHARACTER(len=99) 					:: charTmp1, charTmp2,charTmp3
	REAL*8								:: JumpF,SpeedJump, JumpReconstruct 
    !-------------------------------------------------- 
     
    WRITE(charTmp2,*) Data%iterJump
    WRITE(charTmp1,*) Data%Tour
    WRITE(charTmp3,*) Data%TSCHEME

    ALLOCATE(X_cur(2), SolTriangleX(3),SolTriangleY(3), n_interface(2))
    SolTriangleX =0.d0 ; SolTriangleY =0.d0 ;  n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 ; Beta_RightSide_y = 0.d0 ; Beta_LeftSide_x = 0.d0 
    Beta_RightSide_x = 0.d0 ; Beta_LeftSide_y = 0.d0 
    JumpF = 0.d0 ; nbJumpf = 0  ; SpeedJump = 0.d0      ; JumpReconstruct = 0.d0 
	res1 = (Icing_Tw-Icing_Tm)/erf(Icing_chi) 
	res2 = (Icing_Tm-Icing_Tinit)/erfc(Icing_chi*sqrt(Icing_alpha_l/Icing_alpha_s))
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
          	 
	X_cur = X ! 
	Nele  = mesh%vertex(1)%ele1 ! point of departure to find the element where the node is located 
	CALL LocateNode(X_cur, mesh, Nele, found)

	SolTriangleX = 0.d0 ; 	SolTriangleY = 0.d0 

	DO ik = 1,3 
		Nu_k = Mesh%ele(Nele)%Vertex(ik)
		IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
			NU_k = Mesh%Vertex(Nu_k)%Complementaire
		END IF 
		SolTriangleX(ik) = Sol%Bcg(Nu_k,1)
		SolTriangleY(ik) = Sol%Bcg(Nu_k,2)
	END DO 
		
	CALL InterpolateSolution(X_cur, Beta_RightSide_x, SolTriangleX, Mesh%ele(nele))		! x coordinate for B 
	CALL InterpolateSolution(X_cur, Beta_RightSide_Y, SolTriangleY, Mesh%ele(nele))		! y coordinate for B 
	
	CALL ConstrQR_Beta_LeftSide(Data, Mesh, Sol, X_cur, Beta_LeftSide_x) ! need to change tha also to get y 
	CALL ConstrQR_Beta_LeftSide_y(Data, Mesh, Sol, X_cur, Beta_LeftSide_y) ! need to change tha also to get y 
			
	SpeedJumpx = Beta_LeftSide_x - Beta_RightSide_x
	SpeedJumpy = Beta_LeftSide_y - Beta_RightSide_y

	
  END SUBROUTINE MeanFluxReconstructNodeVal
  !====================================================================================================
  
  
  !====================================================================================================
  SUBROUTINE LocalMeanFluxReconstruct(Data, Mesh, Sol, X_cur, SpeedJumpx, SpeedJumpy, JumpFx, JumpFy)
  !in this case we do not use the mean of the jump but a value for each of the points defining the interface 
	
  	IMPLICIT NONE 
  
	type(DataStructure)  , INTENT(IN) 	:: Data
	type(MeshType)       , INTENT(IN) 	:: Mesh
	type(SolStructure)   , INTENT(IN) 	:: Sol
	REAL*8  		     , INTENT(OUT)	:: JumpFx, JumpFy, SpeedJumpx, SpeedJumpy
	REAL*8, DIMENSION (:), INTENT(IN)   :: X_cur 
    !---------------------------------------------------------------------------------------- 
    INTEGER 			  				:: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			 				:: Nu_i, Nu_j, nbJumpf, Ns, Nu_k 
    REAL*8				  				:: Beta_RightSide, Beta_LeftSide, JumpReconstruct
    REAL*8								:: JumpReconstruct_x, JumpReconstruct_y 
    INTEGER               				:: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL              				:: found 
    REAL*8, ALLOCATABLE, DIMENSION(:) 	:: SolTriangle, n_interface
    !----------------------------------------------------------------------------------------
     
    ALLOCATE(SolTriangle(3),n_interface(2))
    SolTriangle =0.d0 ;  n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 ; Beta_RightSide = 0.d0 ; Beta_LeftSide = 0.d0 

    !analytical velocity 
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    Nu_i = mesh%Tab_NoeudComp(1,1) 
	Nele = mesh%vertex(Nu_i)%ele1 

	CALL LocateNode(X_cur, Mesh, Nele, found)
	
	SolTriangle = 0.d0 
	DO ik=1,3 
		Nu_k = Mesh%ele(Nele)%Vertex(ik)
		IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
			NU_k = Mesh%Vertex(Nu_k)%Complementaire
		END IF 
		SolTriangle(ik) = Sol%Bcg(Nu_k,1)
	END DO 
		
	CALL InterpolateSolution(X_cur, Beta_RightSide, SolTriangle, Mesh%ele(nele))		
	CALL ConstrQR_Beta_LeftSide(Data, Mesh, Sol, X_cur, Beta_LeftSide)
				
	JumpReconstruct_x = Beta_LeftSide - Beta_RightSide
	!PRINT*,JumpReconstruct_x,icing_rho*icing_Lm*Icing_chi*sqrt(Icing_alpha_l/(Data%T-Data%Dt))

	SolTriangle = 0.d0 
	DO ik=1,3 
		Nu_k = Mesh%ele(Nele)%Vertex(ik)
		IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
			NU_k = Mesh%Vertex(Nu_k)%Complementaire
		END IF 
		SolTriangle(ik) = Sol%Bcg(Nu_k,2)
	END DO 
		
	CALL InterpolateSolution(X_cur, Beta_RightSide, SolTriangle, Mesh%ele(nele))	
		
	CALL ConstrQR_Beta_LeftSide_y(Data, Mesh, Sol, X_cur, Beta_LeftSide)
				
	JumpReconstruct_y = Beta_LeftSide - Beta_RightSide

	IF(icas==001) THEN !from the stefan condition 
		SpeedJumpx = JumpReconstruct_x/(Icing_Lm*Icing_rho)
		SpeedJumpy = JumpReconstruct_y/(Icing_Lm*Icing_rho)
	ELSE 
		PRINT*,"Not implemented yet"
		STOP 
	END IF 
	
	IF(JumpReconstruct <0.d0) THEN 
		PRINT*,"Problem of sign in Velocity\MeanFluxReconstruct function"
		STOP 
	END IF 
	JumpFx = JumpReconstruct_x
	JumpFy = JumpReconstruct_y
   
  END SUBROUTINE LocalMeanFluxReconstruct
  !===================================
  
  
   !send the jump value and the velocity in two different variable 
  !===================================
	SUBROUTINE MeanFluxReconstruct22Nicolson(Data, Mesh, Sol,SpeedJump,JumpF)
  	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: JumpF,SpeedJump

    !-------------------------------------------------- 
    INTEGER 			  				:: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			 				:: Nu_i, Nu_j, nbJumpf, Ns, Nu_k 
    REAL*8				  				:: analy_velocity, Beta_RightSide, Beta_LeftSide,JumpReconstruct
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               				:: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL              				:: found 
    REAL*8, ALLOCATABLE, DIMENSION(:) 	:: SolTriangle, n_interface
	REAL*8				 				:: alpha_l,alpha_s,res1,res2 
	REAL*8, DIMENSION(2) 				:: nP 
	type(SolStructure) 					:: Sol_copy
	CHARACTER(len=99) 					:: charTmp1, charTmp2,charTmp3
    !-------------------------------------------------- 
    
    WRITE(charTmp2,*) Data%iterJump
    WRITE(charTmp1,*) Data%Tour
    WRITE(charTmp3,*) Data%TSCHEME

    open(unit = 331,file='MeanValue'//TRIM(ADJUSTL(charTmp1))//TRIM(ADJUSTL(charTmp2))//&
    TRIM(ADJUSTL(charTmp3))//'.data', action="write",position="append")


     Sol_copy = Sol 
     Sol_copy%Bcg = Sol_copy%BcgN ; Sol_copy%Pcg = Sol_copy%PcgN 
     
    ALLOCATE(X_cur(2),SolTriangle(3),n_interface(2))
    SolTriangle =0.d0 ;  n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 ; Beta_RightSide = 0.d0 ; Beta_LeftSide = 0.d0 
    JumpF = 0.d0 ; nbJumpf = 0 ; SpeedJump = 0.d0 ; JumpReconstruct = 0.d0 
	res1 = (Icing_Tw-Icing_Tm)/erf(Icing_chi) 
	res2 = (Icing_Tm-Icing_Tinit)/erfc(Icing_chi*sqrt(Icing_alpha_l/Icing_alpha_s))

    !analytical velocity 
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
		X_cur = mesh%Vertex(Nu_i)%X + Sol_copy%dist(1,Nu_i,:) ! Sol_prec est definis via Sol en post processing 
		np = Sol_copy%dist(1,Nu_i,:)/sqrt(Sol_copy%dist(1,Nu_i,1)**2+Sol_copy%dist(1,Nu_i,2)**2) ! normal vector at the interface, nout used for one direction displacement 
		
		Nele = mesh%vertex(Nu_i)%ele1 
		CALL LocateNode(X_cur, mesh, Nele, found)
		
		SolTriangle = 0.d0 
		DO ik=1,3 
			Nu_k = Mesh%ele(Nele)%Vertex(ik)
			IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
				NU_k = Mesh%Vertex(Nu_k)%Complementaire
			END IF 
			SolTriangle(ik) = Sol_copy%Bcg(Nu_k,1)
		END DO 
		
		CALL InterpolateSolution(X_cur, Beta_RightSide, SolTriangle, Mesh%ele(nele))			
		CALL ConstrQR_Beta_LeftSide(Data, Mesh, Sol_copy, X_cur, Beta_LeftSide)
				
		JumpReconstruct = Beta_LeftSide - Beta_RightSide

		WRITE(331,*)X_cur(1),X_cur(2),JumpReconstruct


		IF(icas ==001) THEN 
			IF(JumpReconstruct <0.d0) THEN 
				!PRINT*,"Problem of sign in Velocity\MeanFluxReconstruct function"
			ELSE			
				JumpF = JumpF + JumpReconstruct ! jump on the true interface 
				nbJumpf = nbJumpf + 1
			END IF 

		ELSE IF(nDim == 2) THEN 
			PRINT*,"Not implemented yet look at MeanFluxReconstruct22"
			STOP 
		ELSE 
			PRINT*,"3D case not available in MeanFluxReconstruct22"
			STOP 
		END IF 

    END DO

	IF (nbJumpf==0 .OR. JumpF == 0.d0 ) THEN 
		PRINT*,"Error the number of jump at the interface is null"
		STOP 
	ELSE 
		JumpF  = JumpF/nbJumpf !mean jump value 
		SpeedJump = JumpF/(Icing_rho*Icing_Lm) ! value of the velocity with the value of the gradient to 
	END IF 
	
	WRITE(331,*)"Mean Value",JumpF,"True Value",Icing_rho*Icing_Lm*Icing_chi*sqrt(Icing_alpha_l/Data%t)

	close(331)
	
  END SUBROUTINE MeanFluxReconstruct22Nicolson
  !===================================
  
  
   
  !===================================
	SUBROUTINE MeanFluxReconstruct_LowerVers(Data, Mesh, Sol,SpeedJump,JumpF)
  	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: JumpF,SpeedJump

    !-------------------------------------------------- 
    INTEGER 			  :: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			  :: Nu_i, Nu_j, nbJumpf, Ns, NU_k 
    REAL*8				  :: analy_velocity,interpSol, Beta_RightSide, Beta_LeftSide,JumpReconstruct
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               :: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL               :: found 
    REAL*8, ALLOCATABLE, DIMENSION(:)  :: SolTriangle, n_interface
	REAL*8				 :: previous_JumpF, distance,interpSoly,norm,Beta_Left_True, Beta_Right_True
	CHARACTER(len=20) 	 :: h_value
	REAL*8				 :: JumpReconstruct_test1, JumpReconstruct_test2 ,alpha_l,alpha_s		

    !-------------------------------------------------- 
     
    ALLOCATE(X_cur(2),SolTriangle(3),n_interface(2))
    SolTriangle =0.d0 ; previous_JumpF = 0.d0 ; n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 ; Beta_RightSide = 0.d0 ; Beta_LeftSide = 0.d0 ; distance = 0.d0 
    JumpF = -10000000000.d0 ; nbJumpf = 0 ; interpSol =0.d0 ; SpeedJump = 0.d0


    !analytical velocity 
	analy_velocity = Icing_chi*sqrt(Icing_Lambda1/(Icing_rho*Icing_c1))/sqrt(Data%T-Data%Dt) 
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
		X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! Sol_prec est definis via Sol en post processing 
		Nele = mesh%vertex(Nu_i)%ele1 
		CALL LocateNode(X_cur, mesh, Nele, found)
		
		SolTriangle = 0.d0 
		DO ik=1,3 
			Nu_k = Mesh%ele(Nele)%Vertex(ik)
			IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
				NU_k = Mesh%Vertex(Nu_k)%Complementaire
			END IF 
			SolTriangle(ik) = Sol%Bcg(Nu_k,1)
		END DO 
		
		interpSol = 0.d0
		CALL InterpolateSolution(X_cur, interpSol, SolTriangle, Mesh%ele(nele))
		Beta_RightSide =  interpSol

		CALL ConstrQR_Beta_LeftSide(Data, Mesh, Sol, X_cur, Beta_LeftSide)
		JumpReconstruct = Beta_LeftSide - Beta_RightSide

		IF(icas ==001) THEN 
			!PRINT*,"jump comp here", JumpReconstruct, JumpF 
			IF(JumpReconstruct >JumpF .AND. JumpReconstruct > 0.d0) THEN 		
				JumpF = JumpReconstruct
				nbJumpf = 1
			ELSE IF (JumpReconstruct < 0.d0) THEN 
				PRINT*,"PB here Velocity\MeanFluxReconstruct_LowerVers"
				STOP
			END IF 

		ELSE IF(nDim == 2) THEN 
			PRINT*,"Not implemented yet look at MeanFluxReconstruct22"
			STOP 
		ELSE 
			PRINT*,"3D case not available in MeanFluxReconstruct22"
			STOP 
		END IF 

    END DO

	IF (nbJumpf==0) THEN 
		PRINT*,"Error the number of jump at the interface is null in MeanFluxReconstruct_LowerVers"
		STOP
	ELSE 
		JumpF = JumpF/nbJumpf !mean jump value 
		SpeedJump = JumpF/(Icing_rho*Icing_Lm) ! value of the velocity 
	END IF 
	
  END SUBROUTINE MeanFluxReconstruct_LowerVers
  !===================================
  
	!==================================!
	SUBROUTINE ConstrQR_Beta_LeftSide(Data, Mesh, Sol, X_cur, Beta_LeftSide)
	!==================================!

	IMPLICIT NONE

	  type(DataStructure)	, INTENT(IN)  		:: Data
	  type(MeshType)    	, INTENT(IN)  		:: Mesh
	  type(SolStructure)	, INTENT(IN)  		:: Sol
	  REAL*8, DIMENSION(:)  , INTENT(IN) 		:: X_cur
	  REAL*8             	, INTENT(OUT)  		:: Beta_LeftSide

	  !----------------------------------------
	  INTEGER                               :: j, Ne, i,valcomp1, lengTab, k, leng_recherche,taille_stencil
	  INTEGER                               :: valcomp2, ie, nu_i, choix, Ndim, nu_j, ib, ils 
	  INTEGER                               :: Nu_k , Nv, leng_tabvtc
	  LOGICAL                               :: InIt1, InIt2, test, TEst_ele, AJOUT , test_edg
	  LOGICAL                               :: Test1 , Test2 ,IFoundYou 
	  type(Element)                         :: ele_original 
	  REAL*8                                :: xx, yy, lambda, max, test_val, residu 
	  REAL*8                                :: x_test, y_test, x_k, x_j, x_i,y_i , y_k , y_j, dist,mintab ,x,y
	  REAL*8                                :: residu1, residu2, residu3 
	  INTEGER                               :: test_tab, tour,lengdif, numb_min, chang, cmpt , Nie,degre_rsct, taille 
	  INTEGER, DIMENSION(:)  , ALLOCATABLE  :: Tab_stencil, Tab_vtc, TabRecherch, TabTridist
	  REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: A, Am1
	  REAL*8, DIMENSION(:)  , ALLOCATABLE   :: test_vertcs, B1, Xa1,B,Xa, B2, Xa2
	  type(Edge)    						:: Edg

	  !----------------------------------------

	Ne = Mesh%Ne  ;  Ndim = Data%Ndim  ;  Nv = Sol%Nv 
	
	lengdif=0 
	leng_recherche = 0
	dist = 0.d0
	degre_rsct = 6 
	taille_stencil = 12 !nombre de noeuds pour la résolution 
	lengTab = 0 
	tour = 0 
	Test_ele = .FALSE.
	xx = X_cur(1)  ;  yy =  X_cur(2)
  
  
	lambda = Icing_Lambda1 ! la reconstruction s'effectue dans la zone +

	ALLOCATE(TabRecherch(Nv)) ; TabRecherch = 0 ; leng_recherche = 0
	  
	!on enregistre tous les noeuds qui sont dans la zone 1 
	!car on reconstruit une valeur à gauche 
	DO ie= 1, Ne 
		DO i = 1, mesh%ele(ie)%nNodesPerElement 
			IF (mesh%ele(ie)%ref == IcingZone1_Ref) THEN ! si on se trouve dans la bonne zone et pas sur le noeud ou il manque une valeur
				Nu_i=mesh%ele(ie)%vertex(i) ;  AJOUT = .TRUE. 
				dist = sqrt((xx-Mesh%vertex(Nu_i)%X(1))**2+(yy-Mesh%vertex(Nu_i)%X(2))**2)
				IF(dist <= 10*Data%h) THEN ! on concentre la recherche sur une longueur maximum au noeud (cas de maillage trés fin!!)
					IF (leng_recherche/=0) THEN !on doit verifié que si le tableau n'est pas vide le noeud n'y est pas deja
						DO k=1, leng_recherche
							IF(Nu_i == TabRecherch(k)) THEN 
								AJOUT = .FALSE.
							END IF 
						END DO 
					END IF 
					IF (AJOUT) THEN 
						leng_Recherche = leng_Recherche +1 
						TabRecherch(leng_Recherche)= Nu_i
					END IF 
				END IF
			END IF 
		END DO 
	END DO 
	
	
	! Mainteant on va trier le tableau par distance par rapport a X_cur  
	! Le nouveau tableau sera TabTridist
	ALLOCATE(TabTridist(leng_Recherche)) ; TabTridist = 0 
	chang = 0 
	DO WHILE (chang /= leng_Recherche) !on doit faire autant de changement qu'il y a de valeur dans le tableau 
		mintab = 10000  ;  numb_min = 0
		DO k = 1, leng_Recherche
			IF(TabRecherch(k)/= 0) THEN !on met à 0 les valeurs que l'on a déja enlevé 
				x_k= mesh%Vertex(TabRecherch(k))%X(1) ; y_k = mesh%Vertex(TabRecherch(k))%X(2)
				IF ( mintab > sqrt((xx-x_k)**2+(yy-y_k)**2) ) THEN 
					mintab = sqrt((xx-x_k)**2+(yy-y_k)**2)
					numb_min = k ! position du noeud dansle tableau 
				END IF
			END IF 
		END DO  
		chang = chang + 1
		TabTridist(chang) = TabRecherch(numb_min)
		TabRecherch(numb_min) = 0 
	END DO 
	
	taille = 12 ! choix de la longueur du stencil 
	ALLOCATE(Tab_stencil(taille)) ; Tab_stencil = 0.d0 ; leng_tabvtc = 0
	k = 0

    DO WHILE (leng_tabvtc < taille .AND. k < leng_Recherche) 
		k = k+1 
		IF(TabTridist(k) /= 0 ) THEN 
			leng_tabvtc = leng_tabvtc + 1
			Tab_stencil(leng_tabvtc)=TabTridist(k) 
		END IF 
    END DO 
    
    taille_stencil = leng_tabvtc

	IF(taille_stencil < 6) THEN
		PRINT*, "probleme sur la definition du stencil pas assez de neoud",taille_stencil
		STOP
	END IF 
	
	DEALLOCATE(TabTridist,TabRecherch)
	
	!Initialisation of the systems matrices 
	ALLOCATE(A(taille_stencil,6)) !matrix
	A = 0.d0 
	ALLOCATE(B1(taille_stencil)) 

	B1 = 0.d0 

	 
	DO i = 1, taille_stencil 
		Nu_i = Tab_stencil(i)
		x = mesh%Vertex(Nu_i)%X(1)  ; y =  mesh%Vertex(Nu_i)%X(2)
			
		!matrice du systeme lineaire reconstruction pour T 
		A(i,1)= 1        	 ; A(i,2)= x    			; A(i,3)= y
		A(i,4)= x*y    	     ; A(i,5)= x**2 			; A(i,6)= y**2

		B1(i)=Sol%bcg(Nu_i,1) !the x coordinate 

	END DO 
    
    ALLOCATE(Xa1(3))
	Xa1 = 0  ! initial guess 
	CALL SolveSystemQR(A(:,1:3),B1,Xa1,taille_stencil,3)
	Beta_LeftSide = Xa1(1) + Xa1(2)*xx  + Xa1(3)*yy
	
	DEALLOCATE(Xa1,A,b1,Tab_stencil)
	

	END SUBROUTINE ConstrQR_Beta_LeftSide
	!========================================!


	!==================================!
	SUBROUTINE ConstrQR_Beta_LeftSide_y(Data, Mesh, Sol, X_cur, Beta_LeftSide)
	!==================================!

	IMPLICIT NONE

	  type(DataStructure)	, INTENT(IN)  		:: Data
	  type(MeshType)    	, INTENT(IN)  		:: Mesh
	  type(SolStructure)	, INTENT(IN)  		:: Sol
	  REAL*8, DIMENSION(:)  , INTENT(IN) 		:: X_cur
	  REAL*8             	, INTENT(OUT)  		:: Beta_LeftSide

	  !----------------------------------------
	  INTEGER                               :: j, Ne, i,valcomp1, lengTab, k, leng_recherche,taille_stencil
	  INTEGER                               :: valcomp2, ie, nu_i, choix, Ndim, nu_j, ib, ils 
	  INTEGER                               :: Nu_k , Nv, leng_tabvtc
	  LOGICAL                               :: InIt1, InIt2, test, TEst_ele, AJOUT , test_edg
	  LOGICAL                               :: Test1 , Test2 ,IFoundYou 
	  type(Element)                         :: ele_original 
	  REAL*8                                :: xx, yy, lambda, max, test_val, residu 
	  REAL*8                                :: x_test, y_test, x_k, x_j, x_i,y_i , y_k , y_j, dist,mintab ,x,y
	  REAL*8                                :: residu1, residu2, residu3 
	  INTEGER                               :: test_tab, tour,lengdif, numb_min, chang, cmpt , Nie,degre_rsct, taille 
	  INTEGER, DIMENSION(:)  , ALLOCATABLE  :: Tab_stencil, Tab_vtc, TabRecherch, TabTridist
	  REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: A, Am1
	  REAL*8, DIMENSION(:)  , ALLOCATABLE   :: test_vertcs, B1, Xa1,B,Xa, B2, Xa2
	  type(Edge)    						:: Edg

	  !----------------------------------------

	Ne = Mesh%Ne  ;  Ndim = Data%Ndim  ;  Nv = Sol%Nv 
	
	lengdif=0 
	leng_recherche = 0
	dist = 0.d0
	degre_rsct = 6 
	taille_stencil = 12 !nombre de noeuds pour la résolution 
	lengTab = 0 
	tour = 0 
	Test_ele = .FALSE.
	xx = X_cur(1)  ;  yy =  X_cur(2)
  
  
	lambda = Icing_Lambda1 ! la reconstruction s'effectue dans la zone +

	ALLOCATE(TabRecherch(Nv)) ; TabRecherch = 0 ; leng_recherche = 0
	  
	!on enregistre tous les noeuds qui sont dans la zone 1 
	!car on reconstruit une valeur à gauche 
	DO ie= 1, Ne 
		DO i = 1, mesh%ele(ie)%nNodesPerElement 
			IF (mesh%ele(ie)%ref == IcingZone1_Ref) THEN ! si on se trouve dans la bonne zone et pas sur le noeud ou il manque une valeur
				Nu_i=mesh%ele(ie)%vertex(i) ;  AJOUT = .TRUE. 
				dist = sqrt((xx-Mesh%vertex(Nu_i)%X(1))**2+(yy-Mesh%vertex(Nu_i)%X(2))**2)
				IF(dist <= 10*Data%h) THEN ! on concentre la recherche sur une longueur maximum au noeud (cas de maillage trés fin!!)
					IF (leng_recherche/=0) THEN !on doit verifié que si le tableau n'est pas vide le noeud n'y est pas deja
						DO k=1, leng_recherche
							IF(Nu_i == TabRecherch(k)) THEN 
								AJOUT = .FALSE.
							END IF 
						END DO 
					END IF 
					IF (AJOUT) THEN 
						leng_Recherche = leng_Recherche +1 
						TabRecherch(leng_Recherche)= Nu_i
					END IF 
				END IF
			END IF 
		END DO 
	END DO 
	
	
	! Mainteant on va trier le tableau par distance par rapport a X_cur  
	! Le nouveau tableau sera TabTridist
	ALLOCATE(TabTridist(leng_Recherche)) ; TabTridist = 0 
	chang = 0 
	DO WHILE (chang /= leng_Recherche) !on doit faire autant de changement qu'il y a de valeur dans le tableau 
		mintab = 10000  ;  numb_min = 0
		DO k = 1, leng_Recherche
			IF(TabRecherch(k)/= 0) THEN !on met à 0 les valeurs que l'on a déja enlevé 
				x_k= mesh%Vertex(TabRecherch(k))%X(1) ; y_k = mesh%Vertex(TabRecherch(k))%X(2)
				IF ( mintab > sqrt((xx-x_k)**2+(yy-y_k)**2) ) THEN 
					mintab = sqrt((xx-x_k)**2+(yy-y_k)**2)
					numb_min = k ! position du noeud dansle tableau 
				END IF
			END IF 
		END DO  
		chang = chang + 1
		TabTridist(chang) = TabRecherch(numb_min)
		TabRecherch(numb_min) = 0 
	END DO 
	
	taille = 12 ! choix de la longueur du stencil 
	ALLOCATE(Tab_stencil(taille)) ; Tab_stencil = 0.d0 ; leng_tabvtc = 0
	k = 0

    DO WHILE (leng_tabvtc < taille .AND. k < leng_Recherche) 
		k = k+1 
		IF(TabTridist(k) /= 0 ) THEN 
			leng_tabvtc = leng_tabvtc + 1
			Tab_stencil(leng_tabvtc)=TabTridist(k) 
		END IF 
    END DO 
    
    taille_stencil = leng_tabvtc

	IF(taille_stencil < 6) THEN
		PRINT*, "probleme sur la definition du stencil pas assez de neoud",taille_stencil
		STOP
	END IF 
	
	DEALLOCATE(TabTridist,TabRecherch)
	
	!Initialisation of the systems matrices 
	ALLOCATE(A(taille_stencil,6)) !matrix
	A = 0.d0 
	ALLOCATE(B1(taille_stencil)) 

	B1 = 0.d0 

	 
	DO i = 1, taille_stencil 
		Nu_i = Tab_stencil(i)
		x = mesh%Vertex(Nu_i)%X(1)  ; y =  mesh%Vertex(Nu_i)%X(2)
			
		!matrice du systeme lineaire reconstruction pour T 
		A(i,1)= 1        	 ; A(i,2)= x    			; A(i,3)= y
		A(i,4)= x*y    	     ; A(i,5)= x**2 			; A(i,6)= y**2

		B1(i)=Sol%bcg(Nu_i,2) !the x coordinate 

	END DO 
    
    ALLOCATE(Xa1(3))
	Xa1 = 0  ! initial guess 
	CALL SolveSystemQR(A(:,1:3),B1,Xa1,taille_stencil,3)
	Beta_LeftSide = Xa1(1) + Xa1(2)*xx  + Xa1(3)*yy
	
	DEALLOCATE(Xa1,A,b1,Tab_stencil)
	

	END SUBROUTINE ConstrQR_Beta_LeftSide_y
	!========================================!

	!==================================!
	SUBROUTINE ConstrQR_T_LeftSide(Data, Mesh, Sol, X_cur, T_LeftSide)
	!==================================!

	IMPLICIT NONE

	  type(DataStructure)	, INTENT(IN)  		:: Data
	  type(MeshType)    	, INTENT(IN)  		:: Mesh
	  type(SolStructure)	, INTENT(IN)  		:: Sol
	  REAL*8, DIMENSION(:)  , INTENT(IN) 		:: X_cur
	  REAL*8             	, INTENT(OUT)  		:: T_LeftSide

	  !----------------------------------------
	  INTEGER                               :: j, Ne, i,valcomp1, lengTab, k, leng_recherche,taille_stencil
	  INTEGER                               :: valcomp2, ie, nu_i, choix, Ndim, nu_j, ib, ils 
	  INTEGER                               :: Nu_k , Nv, leng_tabvtc
	  LOGICAL                               :: InIt1, InIt2, test, TEst_ele, AJOUT , test_edg
	  LOGICAL                               :: Test1 , Test2 ,IFoundYou 
	  type(Element)                         :: ele_original 
	  REAL*8                                :: xx, yy, lambda, max, test_val, residu 
	  REAL*8                                :: x_test, y_test, x_k, x_j, x_i,y_i , y_k , y_j, dist,mintab ,x,y
	  REAL*8                                :: residu1, residu2, residu3 
	  INTEGER                               :: test_tab, tour,lengdif, numb_min, chang, cmpt , Nie,degre_rsct, taille 
	  INTEGER, DIMENSION(:)  , ALLOCATABLE  :: Tab_stencil, Tab_vtc, TabRecherch, TabTridist
	  REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: A, Am1
	  REAL*8, DIMENSION(:)  , ALLOCATABLE   :: test_vertcs, B1, Xa1,B,Xa, B2, Xa2
	  type(Edge)    						:: Edg

	  !----------------------------------------

	Ne = Mesh%Ne  ;  Ndim = Data%Ndim  ;  Nv = Sol%Nv 
	
	lengdif=0 
	leng_recherche = 0
	dist = 0.d0
	degre_rsct = 6 
	taille_stencil = 12 !nombre de noeuds pour la résolution 
	lengTab = 0 
	tour = 0 
	Test_ele = .FALSE.
	xx = X_cur(1)  ;  yy =  X_cur(2)
  
  
	lambda = Icing_Lambda1 ! la reconstruction s'effectue dans la zone +

	ALLOCATE(TabRecherch(Nv)) ; TabRecherch = 0 ; leng_recherche = 0
	  
	!on enregistre tous les noeuds qui sont dans la zone 1 
	!car on reconstruit une valeur à gauche 
	DO ie= 1, Ne 
		DO i = 1, mesh%ele(ie)%nNodesPerElement 
			IF (mesh%ele(ie)%ref == IcingZone1_Ref) THEN ! si on se trouve dans la bonne zone et pas sur le noeud ou il manque une valeur
				Nu_i=mesh%ele(ie)%vertex(i) ;  AJOUT = .TRUE. 
				dist = sqrt((xx-Mesh%vertex(Nu_i)%X(1))**2+(yy-Mesh%vertex(Nu_i)%X(2))**2)
				IF(dist <= 10*Data%h) THEN ! on concentre la recherche sur une longueur maximum au noeud (cas de maillage trés fin!!)
					IF (leng_recherche/=0) THEN !on doit verifié que si le tableau n'est pas vide le noeud n'y est pas deja
						DO k=1, leng_recherche
							IF(Nu_i == TabRecherch(k)) THEN 
								AJOUT = .FALSE.
							END IF 
						END DO 
					END IF 
					IF (AJOUT) THEN 
						leng_Recherche = leng_Recherche +1 
						TabRecherch(leng_Recherche)= Nu_i
					END IF 
				END IF
			END IF 
		END DO 
	END DO 
	
	
	! Mainteant on va trier le tableau par distance par rapport a X_cur  
	! Le nouveau tableau sera TabTridist
	ALLOCATE(TabTridist(leng_Recherche)) ; TabTridist = 0 
	chang = 0 
	DO WHILE (chang /= leng_Recherche) !on doit faire autant de changement qu'il y a de valeur dans le tableau 
		mintab = 10000  ;  numb_min = 0
		DO k = 1, leng_Recherche
			IF(TabRecherch(k)/= 0) THEN !on met à 0 les valeurs que l'on a déja enlevé 
				x_k= mesh%Vertex(TabRecherch(k))%X(1) ; y_k = mesh%Vertex(TabRecherch(k))%X(2)
				IF ( mintab > sqrt((xx-x_k)**2+(yy-y_k)**2) ) THEN 
					mintab = sqrt((xx-x_k)**2+(yy-y_k)**2)
					numb_min = k ! position du noeud dansle tableau 
				END IF
			END IF 
		END DO  
		chang = chang + 1
		TabTridist(chang) = TabRecherch(numb_min)
		TabRecherch(numb_min) = 0 
	END DO 
	
	taille = 12 ! choix de la longueur du stencil 
	ALLOCATE(Tab_stencil(taille)) ; Tab_stencil = 0.d0 ; leng_tabvtc = 0
	k = 0

    DO WHILE (leng_tabvtc < taille .AND. k < leng_Recherche) 
		k = k+1 
		IF(TabTridist(k) /= 0 ) THEN 
			leng_tabvtc = leng_tabvtc + 1
			Tab_stencil(leng_tabvtc)=TabTridist(k) 
		END IF 
    END DO 
    
    taille_stencil = leng_tabvtc

	IF(taille_stencil < 6) THEN
		PRINT*, "probleme sur la definition du stencil pas assez de neoud",taille_stencil
		STOP
	END IF 
	
	DEALLOCATE(TabTridist,TabRecherch)
	
	!Initialisation of the systems matrices 
	ALLOCATE(A(taille_stencil,6)) !matrix
	A = 0.d0 
	ALLOCATE(B1(taille_stencil)) 

	B1 = 0.d0 

	 
	DO i = 1, taille_stencil 
		Nu_i = Tab_stencil(i)
		x = mesh%Vertex(Nu_i)%X(1)  ; y =  mesh%Vertex(Nu_i)%X(2)
			
		!matrice du systeme lineaire reconstruction pour T 
		A(i,1)= 1        	 ; A(i,2)= x    			; A(i,3)= y
		A(i,4)= x*y    	     ; A(i,5)= x**2 			; A(i,6)= y**2

		B1(i)=Sol%Pcg(Nu_i) !the x coordinate 

	END DO 
    
    ALLOCATE(Xa1(6))
	Xa1 = 0  ! initial guess 
	CALL SolveSystemQR(A,B1,Xa1,taille_stencil,6)
	T_LeftSide = Xa1(1) + Xa1(2)*xx  + Xa1(3)*yy + Xa1(4)*yy*xx+Xa1(5)*xx**2+Xa1(6)*yy**2
	
	DEALLOCATE(Xa1,A,b1,Tab_stencil)
	

	END SUBROUTINE ConstrQR_T_LeftSide
	!========================================!
	
	!===================================
	SUBROUTINE MeanT_Reconstruct22(Data, Mesh, Sol,MeanT)
  	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: MeanT
    !-------------------------------------------------- 
    INTEGER 			  				:: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			 				:: Nu_i, Nu_j, nbJumpf, Ns, Nu_k 
    REAL*8				  				:: T_RightSide, T_LeftSide, TReconstruct
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               				:: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL              				:: found 
    REAL*8, ALLOCATABLE, DIMENSION(:) 	:: SolTriangle, n_interface
	REAL*8, DIMENSION(2) 				:: nP 
    !-------------------------------------------------- 
     

    ALLOCATE(X_cur(2),SolTriangle(3),n_interface(2))
    SolTriangle =0.d0 ;  n_interface = 0.d0 ; TReconstruct = 0.d0  ; MeanT = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 ; T_RightSide = 0.d0 ; T_LeftSide = 0.d0 

    !analytical velocity 
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
		X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! Sol_prec est definis via Sol en post processing 
		np = Sol%dist(1,Nu_i,:)/sqrt(Sol%dist(1,Nu_i,1)**2+Sol%dist(1,Nu_i,2)**2) ! normal vector at the interface, nout used for one direction displacement 

		Nele = mesh%vertex(Nu_i)%ele1 
		CALL LocateNode(X_cur, mesh, Nele, found)

		SolTriangle = 0.d0 
		DO ik=1,3 
			Nu_k = Mesh%ele(Nele)%Vertex(ik)
			IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
				NU_k = Mesh%Vertex(Nu_k)%Complementaire
			END IF 
			SolTriangle(ik) = Sol%Pcg(Nu_k)
		END DO 
		
		CALL InterpolateSolution(X_cur, T_RightSide, SolTriangle, Mesh%ele(nele))		
	
		CALL ConstrQR_T_LeftSide(Data, Mesh, Sol, X_cur, T_LeftSide) !à changer 
		!PRINT*,"C",T_LeftSide,T_RightSide 
		TReconstruct = (T_LeftSide + T_RightSide)*0.5d0
		MeanT = MeanT + TReconstruct ! jump on the true interface 
		
    END DO

	MeanT  = MeanT/Ns  !mean jump value 
	!PRINT*,"Value before exit",meanT
  END SUBROUTINE MeanT_Reconstruct22
  !===================================

	!===================================
	SUBROUTINE JumpPrevNewInter (Data, Mesh, Sol, Sol_prec, JumpFx,JumpFy) 
  
	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol_prec, Sol
	REAL*8  		   , INTENT(OUT)	:: JumpFx, JumpFy

    !-------------------------------------------------- 
    INTEGER 			  :: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			  :: Nu_i, Nu_j, nbJumpf, Ns, NU_k 
    REAL*8				  :: analy_velocity,interpSol, Beta_RightSide_x,Beta_RightSide_y, Beta_LeftSide_x,Beta_LeftSide_y
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               :: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL               :: found 
    REAL*8, ALLOCATABLE, DIMENSION(:)  :: SolTriangle_x, SolTriangle_y, n_interface
	REAL*8				 :: distance,interpSoly,norm,Beta_Left_True, Beta_Right_True
	REAL*8				 :: JumpReconstruct_x, JumpReconstruct_y	

    !-------------------------------------------------- 


    ALLOCATE(X_cur(2),SolTriangle_x(3),SolTriangle_y(3))
    SolTriangle_x = 0.d0 ;SolTriangle_y = 0.d0

    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 
    Beta_RightSide_x = 0.d0 ; Beta_LeftSide_x = 0.d0 ; distance = 0.d0 
    JumpFx = 0.d0 ; nbJumpf = 0 ; interpSol =0.d0
	Beta_RightSide_y = 0.d0 ; Beta_LeftSide_y = 0.d0 ; JumpFy = 0.d0
	
	IF (Data%Icing) THEN
		!variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    !Condition sur les aretes de la surrogate
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 

		!Calculus of \beta^- by interpolation
		
		!First :Defining the point on the true interface
		
		X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! position of the node on the true interface 
		
		!Second :Finding the triangle where the node is 
		Nele = mesh%vertex(Nu_i)%ele1 !element of departure for the research 
		CALL LocateNode(X_cur, mesh, Nele, found)
		!Third: Calculus of \beta^- by interpolation
		SolTriangle_x = 0.d0 ; SolTriangle_y = 0.d0
		DO ik=1,3 
		!we need to check if the node is on the surrogate to take the compementary value in the case it is
			Nu_k = Mesh%ele(Nele)%Vertex(ik)
			IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
				NU_k = Mesh%Vertex(Nu_k)%Complementaire
			END IF 
			SolTriangle_x(ik) = Sol_prec%Bcg(Nu_k,1)
			SolTriangle_y(ik) = Sol_prec%Bcg(Nu_k,2)
		END DO 
		
		!interpolation of the left side for \beta^-
		interpSol = 0.d0
		CALL InterpolateSolution(X_cur, interpSol, SolTriangle_x, Mesh%ele(nele))
		Beta_RightSide_x =  interpSol

		interpSol = 0.d0
		CALL InterpolateSolution(X_cur, interpSol, SolTriangle_y, Mesh%ele(nele))
		Beta_RightSide_y =  interpSol

		!The next step is to compute the left side \beta^+ by extrapolation, only the x coordinate is considered
		CALL ConstrQR_Beta_LeftSide2(Data, Mesh, Sol_prec, X_cur, Beta_LeftSide_x, Beta_LeftSide_y)
		
		JumpReconstruct_x = Beta_LeftSide_x - Beta_RightSide_x
		JumpReconstruct_y = Beta_LeftSide_y - Beta_RightSide_y	
		
		JumpFx = JumpFx + JumpReconstruct_x
		JumpFy = JumpFy + JumpReconstruct_y

		nbJumpf = nbJumpf + 1

    END DO

	IF (nbJumpf==0) THEN 
		PRINT*,"Error the number of jump at the interface is null"
	ELSE 
		JumpFx = JumpFx/nbJumpf !mean value 
		JumpFy = JumpFy/nbJumpf !mean value 
	END IF 
	
	  
	END SUBROUTINE JumpPrevNewInter 
	!===================================
  
	!===================================
	SUBROUTINE JumpPrevNewInter_Specific (Data, Mesh, Sol, Sol_prec,Ns_1,Ns_2, JumpFx,JumpFy) 
  
		IMPLICIT NONE 
	  
		type(DataStructure), INTENT(IN)   :: Data
		type(MeshType)     , INTENT(IN)   :: Mesh
		type(SolStructure) , INTENT(IN)   :: Sol_prec, Sol
		REAL*8  		   , INTENT(OUT)  :: JumpFx, JumpFy
		INTEGER			   , INTENT(IN)   :: Ns_1,Ns_2 

		!-------------------------------------------------- 
		INTEGER 			  :: Ne, Nv, nDim, nNodesPerElement, i, ik
		INTEGER 			  :: Nu_i, Nu_j, nbJumpf, Ns, NU_k 
		REAL*8				  :: analy_velocity,interpSol, Beta_RightSide_x,Beta_RightSide_y, Beta_LeftSide_x,Beta_LeftSide_y
		REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
		INTEGER               :: nele ! numerotation of the element countaining the node on the true interface
		LOGICAL               :: found 
		REAL*8, ALLOCATABLE, DIMENSION(:)  :: SolTriangle_x, SolTriangle_y, n_interface
		REAL*8				 :: previous_JumpF, distance,interpSoly,norm,Beta_Left_True, Beta_Right_True
		REAL*8				 :: JumpReconstruct_x, JumpReconstruct_y	

		!-------------------------------------------------- 
		 
		ALLOCATE(X_cur(2),SolTriangle_x(3),SolTriangle_y(3))
		SolTriangle_x =0.d0 ;SolTriangle_y =0.d0; previous_JumpF = 0.d0
		n_interface = 0.d0
		 
		Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
		nNodesPerElement = nDim + 1 
		Beta_RightSide_x = 0.d0 ; Beta_LeftSide_x = 0.d0 ; distance = 0.d0 
		JumpFx = 0.d0 ; nbJumpf = 0 ; interpSol =0.d0
		Beta_RightSide_y = 0.d0 ; Beta_LeftSide_y = 0.d0	
		JumpFy = 0.d0
		
		IF (Data%Icing) THEN
			!variables pour se déplacer dans le tableau 
			Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
		ELSE
		   PRINT*,"We have a problem here Not icing case, MeanFlux function"
		   STOP
		END IF
		
		!Condition sur les aretes de la surrogate
		DO i = 1, Ns 
				 
			Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
			IF(Nu_i == Ns_1 .OR. Nu_i == Ns_2) THEN ! on prend seulement les deux noeuds de l'interface pour l'élément considéré 
				!Calculus of \beta^- by interpolation
				
				!First :Defining the point on the true interface
				
				X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! position of the node on the true interface 
				
				!Second :Finding the triangle where the node is 
				Nele = mesh%vertex(Nu_i)%ele1 !element of departure for the research 
				CALL LocateNode(X_cur, mesh, Nele, found)
				!Third: Calculus of \beta^- by interpolation
				SolTriangle_x = 0.d0 ; SolTriangle_y = 0.d0
				DO ik=1,3 
				!we need to check if the node is on the surrogate to take the compementary value in the case it is
					Nu_k = Mesh%ele(Nele)%Vertex(ik)
					IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
						NU_k = Mesh%Vertex(Nu_k)%Complementaire
					END IF 
					SolTriangle_x(ik) = Sol_prec%Bcg(Nu_k,1)
					SolTriangle_y(ik) = Sol_prec%Bcg(Nu_k,2)
				END DO 
				
				!interpolation of the left side for \beta^-
				interpSol = 0.d0
				CALL InterpolateSolution(X_cur, interpSol, SolTriangle_x, Mesh%ele(nele))
				Beta_RightSide_x =  interpSol

				interpSol = 0.d0
				CALL InterpolateSolution(X_cur, interpSol, SolTriangle_y, Mesh%ele(nele))
				Beta_RightSide_y =  interpSol

				!The next step is to compute the left side \beta^+ by extrapolation, only the x coordinate is considered
				CALL ConstrQR_Beta_LeftSide2(Data, Mesh, Sol_prec, X_cur, Beta_LeftSide_x,Beta_LeftSide_y)
				
				JumpReconstruct_x = Beta_LeftSide_x - Beta_RightSide_x
				JumpReconstruct_y = Beta_LeftSide_y - Beta_RightSide_y

				
				
				JumpFx = JumpFx + JumpReconstruct_x
				JumpFy = JumpFy + JumpReconstruct_y

				nbJumpf = nbJumpf + 1
			END IF 
		END DO

		IF (nbJumpf==0) THEN 
			PRINT*,"Error the number of jump at the interface is null"
		ELSE 
			JumpFx = JumpFx/nbJumpf !mean value 
			JumpFy = JumpFy/nbJumpf !mean value 

		END IF 
		


	  
	END SUBROUTINE JumpPrevNewInter_Specific
	!===================================
  
  
!==================================!
SUBROUTINE ConstrQR_Beta_LeftSide2(Data, Mesh, Sol, X_cur, Beta_LeftSide_x, Beta_LeftSide_y)
!==================================!

  IMPLICIT NONE

  type(DataStructure)	, INTENT(IN)  		:: Data
  type(MeshType)    	, INTENT(IN)  		:: Mesh
  type(SolStructure)	, INTENT(IN)  		:: Sol
  REAL*8, DIMENSION(:)  , INTENT(IN) 		:: X_cur
  REAL*8             	, INTENT(OUT)  		:: Beta_LeftSide_x, Beta_LeftSide_y

  !----------------------------------------
  INTEGER                               :: j, Ne, i,valcomp1, lengTab, k, leng_recherche,taille_stencil
  INTEGER                               :: valcomp2, ie, nu_i, choix, Ndim, nu_j, ib, ils 
  INTEGER                               :: Nu_k , Nv, leng_tabvtc
  LOGICAL                               :: InIt1, InIt2, test, TEst_ele, AJOUT , test_edg
  LOGICAL                               :: Test1 , Test2 ,IFoundYou 
  type(Element)                         :: ele_original 
  REAL*8                                :: xx, yy, lambda, max, test_val, residu 
  REAL*8                                :: x_test, y_test, x_k, x_j, x_i,y_i , y_k , y_j, dist,mintab ,x,y
  REAL*8                                :: residu1, residu2, residu3 
  INTEGER                               :: test_tab, tour,lengdif, numb_min, chang, cmpt , Nie,degre_rsct, taille 
  INTEGER, DIMENSION(:)  , ALLOCATABLE  :: Tab_stencil, Tab_vtc, TabRecherch, TabTridist
  REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: A, Am1
  REAL*8, DIMENSION(:)  , ALLOCATABLE   :: test_vertcs, B1, Xa1,B,Xa, B2, Xa2
  type(Edge)    						:: Edg
  !----------------------------------------

	Ne = Mesh%Ne  ;  Ndim = Data%Ndim  ;  Nv = Sol%Nv 
	
	lengdif=0 
	leng_recherche = 0
	dist = 0.d0
	degre_rsct = 6 
	taille_stencil = 12 !nombre de noeuds pour la résolution 
	lengTab = 0 
	tour = 0 
	Test_ele = .FALSE.
	xx = X_cur(1)  ;  yy =  X_cur(2)
  
  
	lambda = Icing_Lambda1 ! la reconstruction s'effectue dans la zone +

	ALLOCATE(TabRecherch(Nv)) ; TabRecherch = 0 ; leng_recherche = 0
	  
	!on enregistre tous les noeuds qui sont dans la zone 1 
	!car on reconstruit une valeur à gauche 
	DO ie= 1, Ne 
		DO i = 1, mesh%ele(ie)%nNodesPerElement 
			IF (mesh%ele(ie)%ref == IcingZone1_Ref) THEN ! si on se trouve dans la bonne zone et pas sur le noeud ou il manque une valeur
				Nu_i=mesh%ele(ie)%vertex(i) ;  AJOUT = .TRUE. 
				dist = sqrt((xx-Mesh%vertex(Nu_i)%X(1))**2+(yy-Mesh%vertex(Nu_i)%X(2))**2)
				IF(dist <= 10*Data%h) THEN ! on concentre la recherche sur une longueur maximum au noeud (cas de maillage trés fin!!)
					IF (leng_recherche/=0) THEN !on doit verifié que si le tableau n'est pas vide le noeud n'y est pas deja
						DO k=1, leng_recherche
							IF(Nu_i == TabRecherch(k)) THEN 
								AJOUT = .FALSE.
							END IF 
						END DO 
					END IF 
					IF (AJOUT) THEN 
						leng_Recherche = leng_Recherche +1 
						TabRecherch(leng_Recherche)= Nu_i
					END IF 
				END IF
			END IF 
		END DO 
	END DO 
	
	
	! Mainteant on va trier le tableau par distance par rapport a X_cur  
	! Le nouveau tableau sera TabTridist
	ALLOCATE(TabTridist(leng_Recherche)) ; TabTridist = 0 
	chang = 0 
	DO WHILE (chang /= leng_Recherche) !on doit faire autant de changement qu'il y a de valeur dans le tableau 
		mintab = 10000  ;  numb_min = 0
		DO k = 1, leng_Recherche
			IF(TabRecherch(k)/= 0) THEN !on met à 0 les valeurs que l'on a déja enlevé 
				x_k= mesh%Vertex(TabRecherch(k))%X(1) ; y_k = mesh%Vertex(TabRecherch(k))%X(2)
				IF ( mintab > sqrt((xx-x_k)**2+(yy-y_k)**2) ) THEN 
					mintab = sqrt((xx-x_k)**2+(yy-y_k)**2)
					numb_min = k ! position du noeud dansle tableau 
				END IF
			END IF 
		END DO  
		chang = chang + 1
		TabTridist(chang) = TabRecherch(numb_min)
		TabRecherch(numb_min) = 0 
	END DO 
	
	taille = 12 ! choix de la longueur du stencil 
	ALLOCATE(Tab_stencil(taille)) ; Tab_stencil = 0.d0 ; leng_tabvtc = 0
	k = 0

    DO WHILE (leng_tabvtc < taille .AND. k < leng_Recherche) 
		k = k+1 
		IF(TabTridist(k) /= 0 ) THEN 
			leng_tabvtc = leng_tabvtc + 1
			Tab_stencil(leng_tabvtc)=TabTridist(k) 
		END IF 
    END DO 
    
    taille_stencil = leng_tabvtc

	IF(taille_stencil < 6) THEN
		PRINT*, "probleme sur la definition du stencil pas assez de neoud",taille_stencil
		STOP
	END IF 
	
	DEALLOCATE(TabTridist,TabRecherch)
	
	!Initialisation of the systems matrices 
	ALLOCATE(A(taille_stencil,6)) ; A = 0.d0 
	ALLOCATE(B1(taille_stencil),B2(taille_stencil)) 
	B1 = 0.d0 ; B2 = 0.d0 

	 
	DO i = 1, taille_stencil 
		Nu_i = Tab_stencil(i)
		x = mesh%Vertex(Nu_i)%X(1)  ; y =  mesh%Vertex(Nu_i)%X(2)
			
		!matrice du systeme lineaire reconstruction pour T 
		A(i,1)= 1        	 ; A(i,2)= x    			; A(i,3)= y
		A(i,4)= x*y    	     ; A(i,5)= x**2 			; A(i,6)= y**2

		B1(i)=Sol%bcg(Nu_i,1) !the x coordinate 
		B2(i)=Sol%bcg(Nu_i,2) !the x coordinate 

	END DO 
    
    ALLOCATE(Xa1(3))
	Xa1 = 0  ! initial guess 
	CALL SolveSystemQR(A(:,1:3),B1,Xa1,taille_stencil,3)
	Beta_LeftSide_x = Xa1(1) + Xa1(2)*xx  + Xa1(3)*yy
	
	Xa1 = 0  ! initial guess 
	CALL SolveSystemQR(A(:,1:3),B2,Xa1,taille_stencil,3)
	Beta_LeftSide_y = Xa1(1) + Xa1(2)*xx  + Xa1(3)*yy

	DEALLOCATE(Xa1,A,B1,B2,Tab_stencil)
	

END SUBROUTINE ConstrQR_Beta_LeftSide2
!========================================!

 
  !===================================
  SUBROUTINE JumpPrevNewInter_Test (Data, Mesh, Sol, Sol_prec,JumpFx,JumpFy) 
  
	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN)   :: Data
	type(MeshType)     , INTENT(IN)   :: Mesh
	type(SolStructure) , INTENT(IN)   :: Sol_prec, Sol
	REAL*8  		   , INTENT(OUT)  :: JumpFx, JumpFy
    !-------------------------------------------------- 
    INTEGER 			  :: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			  :: Nu_i, Nu_j, nbJumpf, Ns, NU_k 
    REAL*8				  :: analy_velocity,interpSol, Beta_RightSide_x,Beta_RightSide_y, Beta_LeftSide_x,Beta_LeftSide_y
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               :: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL               :: found 
    REAL*8, ALLOCATABLE, DIMENSION(:)  :: SolTriangle_x, SolTriangle_y, n_interface
	REAL*8				 :: previous_JumpF, distance,interpSoly,norm,Beta_Left_True, Beta_Right_True
	REAL*8				 :: JumpReconstruct_x, JumpReconstruct_y	

    !-------------------------------------------------- 
     
    ALLOCATE(X_cur(2),SolTriangle_x(3),SolTriangle_y(3))
    SolTriangle_x =0.d0 ;SolTriangle_y =0.d0; previous_JumpF = 0.d0
    n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim 
    nNodesPerElement = nDim + 1 
    Beta_RightSide_x = 0.d0 ; Beta_LeftSide_x = 0.d0 ; distance = 0.d0 
    JumpFx = 0.d0 ; nbJumpf = 0 ; interpSol =0.d0
	Beta_RightSide_y = 0.d0 ; Beta_LeftSide_y = 0.d0	
	JumpFy = 0.d0
	
	IF (Data%Icing) THEN
		!variables pour se déplacer dans le tableau 
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    !Condition sur les aretes de la surrogate
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
			!Calculus of \beta^- by interpolation
			
			!First :Defining the point on the true interface
			
			X_cur = mesh%Vertex(Nu_i)%X + Sol%dist(1,Nu_i,:) ! position of the node on the true interface 
			
			!Second :Finding the triangle where the node is 
			Nele = mesh%vertex(Nu_i)%ele1 !element of departure for the research 
			CALL LocateNode(X_cur, mesh, Nele, found)
			!Third: Calculus of \beta^- by interpolation
			SolTriangle_x = 0.d0 ; SolTriangle_y = 0.d0
			DO ik=1,3 
			!we need to check if the node is on the surrogate to take the compementary value in the case it is
				Nu_k = Mesh%ele(Nele)%Vertex(ik)
				IF (Mesh%Vertex(Nu_k)%OnSurrogate) THEN 
					NU_k = Mesh%Vertex(Nu_k)%Complementaire
				END IF 
				SolTriangle_x(ik) = Sol_prec%Bcg(Nu_k,1)
				SolTriangle_y(ik) = Sol_prec%Bcg(Nu_k,2)
			END DO 
			
			!interpolation of the left side for \beta^-
			interpSol = 0.d0
			CALL InterpolateSolution(X_cur, interpSol, SolTriangle_x, Mesh%ele(nele))
			Beta_RightSide_x =  interpSol

			interpSol = 0.d0
			CALL InterpolateSolution(X_cur, interpSol, SolTriangle_y, Mesh%ele(nele))
			Beta_RightSide_y =  interpSol

			!The next step is to compute the left side \beta^+ by extrapolation, only the x coordinate is considered
			CALL ConstrQR_Beta_LeftSide2(Data, Mesh, Sol_prec, X_cur, Beta_LeftSide_x,Beta_LeftSide_y)
			
			JumpReconstruct_x = Beta_LeftSide_x - Beta_RightSide_x
			JumpReconstruct_y = Beta_LeftSide_y - Beta_RightSide_y

			
			
			JumpFx = JumpFx + JumpReconstruct_x
			JumpFy = JumpFy + JumpReconstruct_y

			nbJumpf = nbJumpf + 1
    END DO

	IF (nbJumpf==0) THEN 
		PRINT*,"Error the number of jump at the interface is null"
	ELSE 
		JumpFx = JumpFx/nbJumpf !mean value 
		JumpFy = JumpFy/nbJumpf !mean value 

	END IF 
	
  END SUBROUTINE JumpPrevNewInter_Test
  !===================================
  
   !send the jump value and the velocity in two different variable 
  !===================================
	SUBROUTINE MeanFluxReconstruct_Surrogate(Data, Mesh, Sol,SpeedJump,JumpF)
  	IMPLICIT NONE 
  
	type(DataStructure), INTENT(IN) 	:: Data
	type(MeshType)     , INTENT(IN) 	:: Mesh
	type(SolStructure) , INTENT(IN) 	:: Sol
	REAL*8  		   , INTENT(OUT)	:: JumpF,SpeedJump

    !-------------------------------------------------- 
    INTEGER 			  :: Ne, Nv, nDim, nNodesPerElement, i, ik
    INTEGER 			  :: Nu_i, Nu_j, nbJumpf, Ns, NU_k 
    REAL*8				  :: analy_velocity,interpSol, Beta_RightSide, Beta_LeftSide,JumpReconstruct
    REAL*8, ALLOCATABLE , DIMENSION (:) :: X_cur ! Current X 
    INTEGER               :: nele ! numerotation of the element countaining the node on the true interface
    LOGICAL               :: found 
    REAL*8, ALLOCATABLE, DIMENSION(:)  :: SolTriangle, n_interface
	REAL*8				 :: previous_JumpF, distance,interpSoly,norm,Beta_Left_True, Beta_Right_True
	CHARACTER(len=20) 	 :: h_value
	REAL*8				 :: JumpReconstruct_test1, JumpReconstruct_dx ,alpha_l,alpha_s,res1,res2 
	REAL*8				 :: Beta_RightSide_dx, Beta_LeftSide_dx,JumpReconstruct2, JumpF2
	REAL*8				 :: beta_RightSide2, Beta_LeftSide2, JumpDist, TrueJump, True_Left, True_Right 
	REAL*8, DIMENSION(2) :: nP 
    !-------------------------------------------------- 
     
    ALLOCATE(X_cur(2),SolTriangle(3),n_interface(2))
    SolTriangle =0.d0 ; previous_JumpF = 0.d0 ; n_interface = 0.d0
     
    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim ; JumpDist = 0.d0 
    nNodesPerElement = nDim + 1 ; Beta_RightSide = 0.d0 ; Beta_LeftSide = 0.d0 ; distance = 0.d0 
    JumpF = 0.d0 ; nbJumpf = 0 ; interpSol =0.d0 ; SpeedJump = 0.d0 ; JumpF2 = 0.d0 
	JumpReconstruct = 0.d0 ; JumpF2 = 0.d0 ; JumpReconstruct2 = 0.d0 ; JumpReconstruct_dx = 0.d0 
	res1 = (Icing_Tw-Icing_Tm)/erf(Icing_chi) 
	res2 = (Icing_Tm-Icing_Tinit)/erfc(Icing_chi*sqrt(Icing_alpha_l/Icing_alpha_s))


    !analytical velocity 
		
	IF (Data%Icing) THEN
		Ns=Mesh%NarEmb + 1 + mesh%PtEmb          
    ELSE
       PRINT*,"We have a problem here Not icing case, MeanFlux function"
       STOP
    END IF
    
    DO i = 1, Ns 
          	 
		Nu_i = mesh%Tab_NoeudComp(i,1) ; Nu_j = mesh%Tab_NoeudComp(i,2) 
		X_cur = mesh%Vertex(Nu_i)%X !+ Sol%dist(1,Nu_i,:) ! Sol_prec est definis via Sol en post processing 		
		
		Beta_RightSide =  Sol%Bcg(Nu_j,1) ;	Beta_LeftSide = Sol%Bcg(Nu_i,1) 	

	
		Beta_RightSide_dx = -X_cur(1)/(2.d0*Icing_alpha_s*Data%t)*Icing_lambda2*res2/sqrt(PI*Icing_alpha_s*Data%T)&
		*exp(-X_cur(1)**2/(4.d0*Icing_alpha_s*Data%t))*Sol%dist(1,Nu_i,1)

		Beta_LeftSide_dx = -X_cur(1)/(2.d0*Icing_alpha_l*Data%t)*Icing_lambda1*res1/sqrt(PI*Icing_alpha_l*Data%T)&
		*exp(-X_cur(1)**2/(4.d0*Icing_alpha_l*Data%t))*Sol%dist(1,Nu_i,1)
		

		True_Left        = -Icing_lambda1*(-res1/(sqrt(PI*Icing_alpha_l*Data%T))*exp(-X_cur(1)**2/(4.d0*Icing_alpha_l*Data%T)))
		   
		True_Right       = -Icing_lambda2*(-res2/(sqrt(PI*Icing_alpha_s*Data%T))*exp(-X_cur(1)**2/(4.d0*Icing_alpha_s*Data%T)))
				
		TrueJump         =  (True_Left - True_Right)+ (Beta_LeftSide_dx - Beta_RightSide_dx )
		JumpReconstruct  =  Beta_LeftSide - Beta_RightSide 
		JumpReconstruct2 =  Beta_LeftSide_dx - Beta_RightSide_dx 
		JumpDist         =  JumpReconstruct + JumpReconstruct2 
		
		!PRINT*,"Difference jump with dist", JumpReconstruct, JumpDist, TrueJump 


		IF(icas ==001) THEN 
			IF(JumpReconstruct <0.d0) THEN 
				PRINT*,"Problem of sign in Velocity\MeanFluxReconstruct_Surrogate"
				STOP
			ELSE			
				JumpF = JumpF + JumpReconstruct ! jump on the true interface 
				!JumpF2 = JumpF2 + TrueJump 
				JumpF2 = JumpF2 + TrueJump
				nbJumpf = nbJumpf + 1
			END IF 

		ELSE IF(nDim == 2) THEN 
			PRINT*,"Not implemented yet look at MeanFluxReconstruct_Surrogate"
			STOP 
		ELSE 
			PRINT*,"3D case not available in MeanFluxReconstruct_Surrogate"
			STOP 
		END IF 

    END DO

	IF (nbJumpf==0) THEN 
		PRINT*,"Error the number of jump at the interface is null"
	ELSE 
		JumpF  = JumpF/nbJumpf !mean jump value 
		JumpF2 = JumpF2/nbJumpf 
		PRINT*,"Difference final jump value", JumpF, JumpF2 
		SpeedJump = JumpF/(Icing_rho*Icing_Lm) ! value of the velocity with the value of the gradient to 
	END IF 
	Jumpf = Jumpf2 
	
  END SUBROUTINE MeanFluxReconstruct_Surrogate
  !===================================
 
END MODULE Velocity
