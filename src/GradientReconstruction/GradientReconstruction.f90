MODULE GradientReconstruction

  USE Types
  USE PublicVar
  USE Algebra

  IMPLICIT NONE

CONTAINS

  !==========================================================!
  SUBROUTINE GreenGaussReconstruction_Primal(Sol, Mesh, Data)
  !==========================================================!

  !*************************************************************************
  !   Builds flux in the case of a problem in primal form
  !	  Modified to reconstruct B with the correct lambda value and the SurVertex structure 
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************


    IMPLICIT NONE
    type(SolStructure) , INTENT(INOUT)  :: Sol
    type(MeshType)     , INTENT(IN)     :: Mesh
    type(DataStructure), INTENT(IN)     :: Data
    !-------------------------------------------------------------------------------
    type(Element)						:: ele
    !---------------------
    REAL*8, DIMENSION(:), ALLOCATABLE   :: denom
    REAL*8, DIMENSION(data%nDim) 		:: GradP1, Gphi
    !-------------------------------------------------------------------------------
    REAL*8								:: ui, lambda, Bx_Recons,By_Recons
    !-------------
    INTEGER 							:: Ne, Nv, nNodesPerElement, nDim, ir 
    INTEGER 							:: i, ie, NU_i, iv, Ns, Nu_j
    !-------------------------------------------------------------------------------


    IF (Data%Icing .eqv. .TRUE.) THEN 
      Ns=Mesh%NB_VertInterf
    ELSE
      Ns=0
    END IF
    
    Ne = mesh%Ne ; Nv = mesh%Nv ; nDim = Data%nDim
    
    ALLOCATE(denom(Nv+Ns)) ; denom = 0.d0
    Sol%Bcg = 0.d0 
    	
    DO ie = 1, Ne
        ele = mesh%ele(ie) 
        nNodesPerElement = ele%nNodesPerElement

        IF ( ele%Solve ) THEN
        
			lambda = ele%lambda !value depending on the element position 
			!to know where to perform the reconstruction 
			        
			gradP1 = 0.d0 ; ui = 0.d0 
				
			!reconstruction of the gradient in the element 
			DO i = 1, nNodesPerElement
				NU_i = ele%SurVertex(i) !with the complementary notation 
				ui = Sol%Pcg(NU_i)
				GradP1 = GradP1 + ui*ele%N(i,:)/(nDim*ele%A)
				
				denom(NU_i) = denom(NU_i) + ele%A ! aire du triangle 
			END DO !end definition of the gradient in the element 
			
			GradP1 = -lambda*GradP1 !definition of B from the value of D.T 
			
			DO i = 1, nNodesPerElement
				NU_i = ele%SurVertex(i)
				Sol%Bcg(NU_i,:) = Sol%Bcg(NU_i,:) + GradP1*ele%A
			END DO
							
		END IF
    END DO

	DO iv = 1, Nv+Ns
		IF ( denom(iv) < 1e-10 ) STOP
		Sol%Bcg(iv,:) = Sol%Bcg(iv,:)/denom(iv)
	END DO
    
    !IF(icas==001) THEN 
	!
	!	DO i = 1, nNodesPerElement
	!		NU_i = ele%SurVertex(i) ; 	NU_j = ele%Vertex(i)
	!		!check if a value for Nu_i already exists 
	!		Already_reconstruct = .FALSE. 
	!		
	!		DO ir =1, Nv+Ns 
	!			IF(ValReconstruct(ir)==Nu_i) Already_reconstruct = .TRUE. 
	!		END DO 
	!		
	!		IF(.NOT. Already_reconstruct) THEN 
	!			ValReconstruct(taille) = Nu_i 
	!			Bx_Recons = 0.d0 ; By_Recons = 0.d0 
	!			CALL PolyConstrQR_Primal(Data,Mesh,Sol,Nu_i,Nu_j,Bx_Recons,By_Recons,lambda,areaId)
	!			!CALL Reconstruction_Primal(Sol, Mesh, Data,ele,Nu_i,Bx_Recons,By_Recons)
	!			Sol%Bcg(NU_i,1) = Bx_Recons
	!			Sol%Bcg(NU_i,2) = By_Recons
	!		END IF 
	!		
	!	END DO 
	!	
	!ELSE 

  END SUBROUTINE GreenGaussReconstruction_Primal
  !==========================================================!
  
  !==================================!
  SUBROUTINE PolyConstrQR_Primal(Data,Mesh,Sol,vertcs,vertcs_nonComp,Bx_Recons,By_Recons,lambda,Zone_Ref)
  !==================================!

  IMPLICIT NONE

  type(MeshType)     , INTENT(IN)	  	:: Mesh
  type(DataStructure), INTENT(IN)  	    :: Data
  REAL*8             , INTENT(OUT)  	:: Bx_Recons, By_Recons 
  REAL*8             , INTENT(IN)  		:: lambda
  INTEGER            , INTENT(INOUT)    :: vertcs
  INTEGER            , INTENT(IN) 	    :: vertcs_nonComp
  INTEGER            , INTENT(IN) 		:: Zone_Ref !area of reconstruction
  type(SolStructure) , INTENT(IN)  		:: Sol

  !----------------------------------------
  INTEGER                               :: j, Ne, i,valcomp1, lengTab, k, leng_recherche,taille_stencil
  INTEGER                               :: valcomp2, ie, nu_i, choix, Ndim, nu_j, ib, ils 
  INTEGER                               :: Nu_k , Nv, leng_tabvtc, vertcs_original
  LOGICAL                               :: InIt1, InIt2, test, TEst_ele, AJOUT , test_edg
  LOGICAL                               :: Test1 , Test2 ,IFoundYou 
  type(Element)                         :: ele_original 
  REAL*8                                :: xx, yy, x, y, max, test_val, residu, exact_BX 
  REAL*8                                :: U1,U2,U3, eps,resu1,resu2,pente,h_ref, exact_T
  REAL*8                                :: x_test, y_test, x_k, x_j, x_i,y_i , y_k , y_j, dist,mintab 
  REAL*8                                :: residu1, residu2, residu3, Nrm, res1
  INTEGER                               :: test_tab, tour,lengdif, numb_min, chang, cmpt , Nie,degre_rsct, taille 
  INTEGER								:: id 
  INTEGER, DIMENSION(:,:)  , ALLOCATABLE  :: Tab_stencil
  REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: A, Am1, AB,invAB,TabRecherch,TabTridist
  REAL*8, DIMENSION(:)  , ALLOCATABLE   :: B, Bdef , Xa, Cdef, test_vertcs, Exa, Exa_prec, B1,B2,cdef2, cdef3 
  type(Edge)    						:: Edg
  type(element) 						:: eleM, eleP

  !----------------------------------------

	Ne = Mesh%Ne ; Ndim = Data%Ndim ; Nv = Sol%Nv ; lengdif=0 ; leng_recherche = 0
	vertcs_original = vertcs ; dist = 0.d0
	degre_rsct = 6 ; taille_stencil = 12 !nombre de noeuds pour la résolution 
	lengTab = 0 ; tour = 0 ; Test_ele = .FALSE.
	xx = mesh%Vertex(vertcs_nonComp)%X(1)  ; yy =  mesh%Vertex(vertcs_nonComp)%X(2)
  
	ALLOCATE(TabRecherch(Nv,2)) ; TabRecherch = 0 ; leng_recherche = 0
	
	
	
	DO ie= 1, Ne 
		DO i = 1, mesh%ele(ie)%nNodesPerElement 
			IF (mesh%ele(ie)%ref == Zone_Ref .AND. mesh%ele(ie)%SurVertex(i)/=vertcs_original) THEN ! si on se trouve dans la bonne zone et pas sur le noeud ou il manque une valeur				
				
				Nu_i = mesh%ele(ie)%SurVertex(i) ; Nu_j = mesh%ele(ie)%Vertex(i)  
				AJOUT = .TRUE. ; dist = 0.d0 
				dist = sqrt((xx-Mesh%vertex(Nu_j)%X(1))**2+(yy-Mesh%vertex(Nu_j)%X(2))**2)		
				
				Nrm = 0.d0
				DO id = 1, nDim
					Nrm = Nrm + (-mesh%ele(ie)%N(1,id))**2
				END DO
				Nrm = sqrt(Nrm)
				!value of the characteristic length 
				h_ref = mesh%ele(ie)%A/Nrm ! valeur caracteristique de l'élément 
	
				IF(dist <= 8*h_ref) THEN ! on concentre la recherche sur une longueur maximum au noeud (cas de maillage trés fin!!)
					IF (leng_recherche/=0) THEN !on doit verifié que si le tableau n'est pas vide le noeud n'y est pas deja
						DO k= 1, leng_recherche
							IF(Nu_i == TabRecherch(k,1)) THEN 
								AJOUT = .FALSE.
							END IF 
						END DO 
					END IF 
					IF (AJOUT) THEN 
						leng_Recherche = leng_Recherche +1 
						TabRecherch(leng_Recherche,1)= Nu_i
						TabRecherch(leng_Recherche,2)= Nu_j
					END IF 	
				END IF			
			END IF 	
		END DO 
	END DO 

	
	! Mainteant on va trier le tableau par distance par rapport a vertcs_original  
	! Le nouveau tableau sera TabTridist
	ALLOCATE(TabTridist(leng_Recherche,2)) ; TabTridist = 0 
	x= mesh%Vertex(vertcs_nonComp)%X(1) ; y = mesh%Vertex(vertcs_nonComp)%X(2)
	chang = 0 
	taille_stencil = 12 ! choix de la longueur du stencil 

	
	DO WHILE (chang /= taille_stencil+1) !on doit faire autant de changement qu'il y a de valeur dans le tableau 
		mintab = 10000  ;  numb_min = 0
		DO k = 1, leng_Recherche
			IF(TabRecherch(k,1)/= 0) THEN !on met à 0 les valeurs que l'on a déja enlevé 
				Nu_j = TabRecherch(k,2)
				x_k = mesh%Vertex(Nu_j)%X(1) ; y_k = mesh%Vertex(Nu_j)%X(2)
				IF ( mintab > sqrt((x-x_k)**2+(y-y_k)**2) ) THEN 
					mintab = sqrt((x-x_k)**2+(y-y_k)**2)
					numb_min = k ! position du noeud dans le tableau 
				END IF
			END IF 
		END DO  
		chang = chang + 1
		TabTridist(chang,1) = TabRecherch(numb_min,1)
		TabTridist(chang,2) = TabRecherch(numb_min,2)
		TabRecherch(numb_min,:) = 0.d0 
	END DO 	
	
    
	DEALLOCATE(TabRecherch)
	
	!Initialisation of the systems matrices 
	ALLOCATE(A(taille_stencil,6)) ; ALLOCATE(B(taille_stencil)) 
	A = 0.d0 ; B = 0.d0 

	!reconstruction quadratique pour la température, reconstruction linéaire pour le flux 
	 
	DO i = 1, taille_stencil 
		Nu_i = TabTridist(i,1) ; Nu_j = TabTridist(i,2)
		x = mesh%Vertex(Nu_j)%X(1)  ; y =  mesh%Vertex(Nu_j)%X(2)
			
		!matrice du systeme lineaire reconstruction pour T 
		A(i,1)= 1        	 ; A(i,2)= x    			; A(i,3)= y
		A(i,4)= x*y    	     ; A(i,5)= x**2 			; A(i,6)= y**2

		! second membre pour la temperature 
		IF(Zone_Ref==1) THEN 
			res1 = (Icing_Tw-Icing_Tm)/erf(Icing_chi)
			exact_T = Icing_Tw - res1*erf(x/(2.d0*sqrt(Icing_alpha_l*Data%t)))
		ELSE 	
			res1 = (Icing_Tm-Icing_Tinit)/erfc(Icing_chi*sqrt(Icing_alpha_l/Icing_alpha_s))
			exact_T = Icing_Tinit + res1*erfc(x/(2.d0*sqrt(icing_alpha_s*Data%t)))
		END IF 
		
		B(i)=Sol%pcg(Nu_i)
		!B(i)=exact_T
			
	END DO 
            
    ALLOCATE(Xa(6))
	Xa = 0 ! initial guess 
	!
	CALL SolveSystemQR(A,b,Xa,taille_stencil,6)
	xx = mesh%Vertex(vertcs_nonComp)%X(1)  ; yy =  mesh%Vertex(vertcs_nonComp)%X(2)

	Bx_Recons = -lambda*( Xa(2) + Xa(4)*yy  + 2*Xa(5)*xx)
	By_Recons = -lambda*( Xa(3) + Xa(4)*xx  + 2*Xa(6)*yy)
	
	
!	IF(Zone_Ref==1) THEN 
!		res1 = (Icing_Tw-Icing_Tm)/erf(Icing_chi)
!		exact_BX = -lambda*(-res1/(sqrt(PI*Icing_alpha_l*Data%t))*exp(-xx**2/(4.d0*Icing_alpha_l*Data%T))) !Bx
!		exact_T = Icing_Tw - res1*erf(xx/(2.d0*sqrt(Icing_alpha_l*Data%t)))
!	ELSE 	
!		res1 = (Icing_Tm-Icing_Tinit)/erfc(Icing_chi*sqrt(Icing_alpha_l/Icing_alpha_s))
!		exact_BX = -lambda*(-res1/(sqrt(PI*Icing_alpha_s*Data%T))*exp(-xx**2/(4.d0*Icing_alpha_s*Data%T)))
!		exact_T = Icing_Tinit + res1*erfc(xx/(2.d0*sqrt(icing_alpha_s*Data%t)))
!	END IF 
	
	!PRINT*,"T",	Sol%pcg(vertcs_original),"ex T",exact_T,Zone_Ref	
	!PRINT*,"B recons",Bx_Recons,By_Recons,"ex",exact_BX
!	IF(icas == 001) THEN 
!		By_Recons = 0.d0
!	END IF 
!	 Bx_Recons = exact_BX

!	DEALLOCATE(Xa,A,B)

END SUBROUTINE PolyConstrQR_Primal
!========================================!

  
   !==========================================================!
  SUBROUTINE Reconstruction_Primal(Sol, Mesh, Data,ele,vertcs,Bx,By)
  !==========================================================!

  !*************************************************************************
  !   Builds flux in the case of a problem in primal form
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************


    IMPLICIT NONE
    type(SolStructure) , INTENT(IN)		:: Sol
    type(MeshType)     , INTENT(IN)    	:: Mesh
    type(DataStructure), INTENT(IN)    	::Data
    type(Element)      , INTENT(IN)  	:: ele  !element etudie
    REAL*8 			   , INTENT(OUT)  	:: Bx, By 
    INTEGER			   , INTENT(IN) 	:: vertcs 
    !-----------------------------------------

    REAL*8, DIMENSION(mesh%Nv)   :: denom
    REAL*8, DIMENSION(data%nDim) :: GradP1, Gphi
    REAl*8 , DIMENSION(data%nDim):: Bcg
    !--------------------------------------------
    REAL*8 :: ui
    !-------------
    INTEGER :: Ne, Nv, nNodesPerElement, nDim
    INTEGER :: i, ie, NU_i, iv
    !-----------------------------------------

    Ne = mesh%Ne ; Nv = mesh%Nv ; nDim = data%nDim

    Bcg = 0.d0 ; denom = 0.d0

    nNodesPerElement = ele%nNodesPerElement
    gradP1 = 0.d0
    DO i = 1, nNodesPerElement
        NU_i = ele%Vertex(i)
        ui = sol%pcg(NU_i)
        GradP1 = GradP1 + ui*ele%N(i,:)/((nDim)*ele%A)
        denom(NU_i) = denom(NU_i) + ele%A
    END DO
    DO i = 1, nNodesPerElement
        NU_i = ele%Vertex(i)
        IF(nu_i == vertcs) THEN 
			Bcg = Bcg + GradP1*ele%A
		END IF 
    END DO


    IF ( denom(vertcs) > 1e-10 ) THEN
        Bcg = -Bcg/denom(vertcs)
    END IF
    Bx= Bcg(1) ; By= Bcg(2)


  END SUBROUTINE Reconstruction_Primal
  !==========================================================!
  
  
  !==================================!
  SUBROUTINE LookForComp2(Mesh,InIt,vertcs,i,test)
    !==================================!

    IMPLICIT NONE

    type(MeshType)     , INTENT(IN)   :: Mesh
    LOGICAL            ,INTENT(OUT) :: InIt
    Integer             , INTENT(IN)  :: vertcs
    INTEGER             , INTENT(OUT) :: i
    INTEGER             , INTENT(IN)  :: test
    !----------------------------------------
    INTEGER       :: j, Ns,k
    !----------------------------------------



    IF (test == 0) THEN ! on regarde dans le tableau courant
       Ns = Mesh%NarEmb + 1 + mesh%PtEmb
       InIt = .FALSE.
       j=0 ; i=0
       DO WHILE (.NOT. InIt .AND. j<Ns)
          j=j+1
          IF (mesh%Tab_NoeudComp(j,1)==vertcs) THEN
             InIt = .TRUE.
             i= mesh%Tab_NoeudComp(j,2) !la ref est assoie a la position dans le vecteur
          END IF
       END DO
    ELSE IF(test==1) THEN !on regarde dans le tableau au temps precedent
       Ns = Mesh%NarEmb_prec + 1 + mesh%PtEmb
       InIt = .FALSE.
       j=0 ; i=0
       DO WHILE (.NOT. InIt .AND. j<Ns)
          j=j+1
          IF (mesh%Tab_NoeudComp_prec(j,1)==vertcs) THEN
             InIt = .TRUE.
             i= mesh%Tab_NoeudComp_prec(j,2)
          END IF
       END DO
    ELSE IF (test == 2) THEN
		Ns = Mesh%NarEmb_prec_inter + 1 + mesh%PtEmb
		InIt = .FALSE.
		j=0 ; i=0
		DO WHILE (.NOT. InIt .AND. j<Ns)
			j=j+1
			IF (Mesh%Tab_NoeudComp_prec_inter(j,1)==vertcs) THEN
				InIt = .TRUE.
				i= Mesh%Tab_NoeudComp_prec_inter(j,2)
			END IF
		END DO
	ELSE
		PRINT*,"Error in LookForComp, choice ", test ,"not available"
		STOP
    END IF


END SUBROUTINE LookForComp2
!========================================!
  
  
END MODULE GradientReconstruction
