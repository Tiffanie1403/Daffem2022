MODULE PolyReconstruction
! routine de reconstruction de polynome à l'interface mobile
! definit aussi des routines utilisés dans la modification de la structure du vecteur solution au temps n pour le temps n+1
  USE Types
  USE Algebra
  USE Algebra2
  USE Data_mod
  USE ExactFunctionJumpt
  USE ReadMesh_Mod
  USE SVD_Mod
  USE GMRES 
  USE GradientReconstruction
  
  USE PublicVar
  USE Quadrature
  USE Distance
  USE ComputeNormalsAndTangents

  IMPLICIT NONE

CONTAINS

  !==================================!
  SUBROUTINE VERIF_Comp(Mesh,test)
    !==================================!

    IMPLICIT NONE

    type(MeshType)     , INTENT(IN)  :: Mesh
    LOGICAL            ,INTENT(OUT)  :: test
    !----------------------------------------
    INTEGER       :: i, Ns, Ns_prec
    !----------------------------------------

    !ici on utilise la veritable valeur de Ns pour définir si l'interface à changé 
    Ns = Mesh%NarEmb + 1 + mesh%PtEmb ! il s'agit de la taille du tableau de complementaire 
    Ns_prec = Mesh%NarEmb_prec + 1 + mesh%PtEmb

    test = .TRUE. ! on considére pour l'intialisation que les deux tableau sont identiques
    i=0
    IF(Ns == Ns_prec) THEN
       DO WHILE (test .AND. i < Ns)
          i=i+1
          IF (mesh%Tab_NoeudComp(i,1)/= mesh%Tab_NoeudComp_prec(i,1)) THEN
             test = .FALSE.
          END IF
       END DO
    ELSE !si les tableaux ne sont pas de la meme taille ils sont forcement differents
       test = .FALSE.
    END IF

  END SUBROUTINE VERIF_Comp
  !========================================!


!==================================!
  SUBROUTINE VERIF_Comp_Newton(Mesh,test)
    !==================================!

	!There is an intermediate interface with the pseudo newton process 
    IMPLICIT NONE

    type(MeshType)     , INTENT(IN)  :: Mesh
    LOGICAL            ,INTENT(OUT)  :: test
    !----------------------------------------
    INTEGER       :: i, Ns, Ns_prec_inter
    !----------------------------------------

    !ici on utilise la veritable valeur de Ns pour définir si l'interface à changé 
    Ns = Mesh%NarEmb + 1 + mesh%PtEmb ! il s'agit de la taille du tableau de complementaire 
    Ns_prec_inter = Mesh%NarEmb_prec_inter  + 1 + mesh%PtEmb

    test = .TRUE. ! on considére pour l'intialisation que les deux tableau sont identiques
    i=0
    IF(Ns == Ns_prec_inter) THEN
       DO WHILE (test .AND. i < Ns)
          i=i+1
          IF (mesh%Tab_NoeudComp(i,1) /= mesh%Tab_NoeudComp_prec_inter(i,1)) THEN
             test = .FALSE.
          END IF
       END DO
    ELSE !si les tableaux ne sont pas de la meme taille ils sont forcement differents
       test = .FALSE.
    END IF

  END SUBROUTINE VERIF_Comp_Newton
  !========================================!
  

  !==================================!
  SUBROUTINE LookForComp(Mesh,InIt,vertcs,i,test)
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


END SUBROUTINE LookForComp
!========================================!

!==================================!
SUBROUTINE DirectNeigh(Mesh,vertcs,Nu_i,test)
!==================================!
! Dans la recherche du stencil de reconstruction on ne vérifie que le noeud étudié n'est pas un noeud direct du noeud courant
! car si à l'étape précédente il n'a pas était choisi cela veut dire qu'il fait parti des nouveaux complémentaires et donc
! la valeur pour la reconstruction n'est pas connu en ce point
  IMPLICIT NONE

  type(MeshType)     , INTENT(IN)  :: Mesh
  INTEGER            , INTENT(IN)  :: vertcs, Nu_i
  logical            , INTENT(OUT) :: test
  !----------------------------------------
  INTEGER                              :: ie, i, Nu_j, j, Nu_k
  INTEGER                              :: Ne, Nu_p
  !----------------------------------------

 Ne = Mesh%Ne
 ie = 0 ; test = .TRUE.
  DO WHILE ( ie < Ne .AND. test)
    ie = ie + 1
     DO i = 1, mesh%ele(ie)%nNodesPerElement
        Nu_j=mesh%ele(ie)%vertex(i)
        IF (Nu_j == vertcs) THEN ! si on se trouve sur mon nouveau noeud d'interface
           DO j = 1, mesh%ele(ie)%nNodesPerElement
              IF (i /=j) THEN
                 Nu_k = mesh%ele(ie)%vertex(j)
                 IF (Nu_k == Nu_i) THEN
                   test = .FALSE.
                   ! les deux noeuds on un elements en commun
                 END IF
            END IF
          END DO
        END IF
        IF (Nu_j == Nu_i) THEN ! si on se trouve sur mon nouveau noeud d'interface
           DO j = 1, mesh%ele(ie)%nNodesPerElement
              IF (i /=j) THEN
                 Nu_k = mesh%ele(ie)%vertex(j)
                 IF (Nu_k == vertcs) THEN
                   test = .FALSE.
                   ! les deux noeuds on un elements en commun
                 END IF
            END IF
          END DO
        END IF
      END do
    END DO



END SUBROUTINE DirectNeigh
!========================================!

!==================================!
SUBROUTINE PolyConstr7(Data,Mesh,Sol,vertcs,val1,val2,val3)
!==================================!

  IMPLICIT NONE

  type(MeshType)     , INTENT(INOUT)  	:: Mesh
  type(DataStructure), INTENT(INOUT)  	:: Data
  REAL*8             ,INTENT(OUT)  		:: val1,val2,val3
  INTEGER            , INTENT(INOUT) 	:: vertcs
  type(SolStructure) , INTENT(IN)  		:: Sol

  !----------------------------------------
  INTEGER                               :: j, Ne, i,valcomp1, lengTab, k, leng_recherche,taille_stencil
  INTEGER                               :: valcomp2, ie, nu_i, choix, Ndim, nu_j, ib, ils 
  INTEGER                               :: Nu_k , Nv, leng_tabvtc, vertcs_original
  LOGICAL                               :: InIt1, InIt2, test, TEst_ele, AJOUT , test_edg
  LOGICAL                               :: Test1 , Test2 
  type(Element)                         :: ele_original 
  REAL*8                                :: xx, yy, x, y, lambda, max, test_val, residu 
  REAL*8                                :: U1,U2,U3, eps,resu1,resu2,pente
  REAL*8                                :: x_test, y_test, x_k, x_j, x_i,y_i , y_k , y_j, dist,mintab 
  REAL*8                                :: residu1, residu2, residu3 
  INTEGER                               :: test_tab, tour,lengdif, numb_min, chang, cmpt , Nie,degre_rsct, taille 
  INTEGER, DIMENSION(:)  , ALLOCATABLE  :: Tab_stencil, Tab_vtc, TabRecherch, TabTridist
  REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: A, Am1, AB,invAB
  REAL*8, DIMENSION(:)  , ALLOCATABLE   :: B, Bdef , Xa, Cdef, test_vertcs, Exa, Exa_prec, B1,B2, Xa1,Xa2,cdef2, cdef3 
  type(Edge)    						:: Edg
  type(element) 						:: eleM, eleP

  !----------------------------------------

  Ne = Mesh%Ne ; Ndim = Data%Ndim ; Nv = Sol%Nv ; lengdif=0 ; leng_recherche = 0
  vertcs_original = vertcs ; dist = 0.d0
  degre_rsct = 6 ; taille_stencil = 12!nombre de noeuds pour la résolution 
  lengTab = 0 ; tour = 0 ; Test_ele = .FALSE.
  xx = mesh%Vertex(vertcs_original)%X(1)  ; yy =  mesh%Vertex(vertcs_original)%X(2)
  
  
	lambda = Icing_Lambda1 ! la reconstruction s'effectue dans la zone 1 

	ALLOCATE(TabRecherch(Nv)) ; TabRecherch = 0 ; leng_recherche = 0
	  
	!on enregistre tous les noeuds qui sont dans la zone 1 
	DO ie= 1, Ne 
		DO i = 1, mesh%ele(ie)%nNodesPerElement 
			IF (mesh%ele(ie)%ref == IcingZone1_Ref) THEN ! si on se trouve dans la bonne zone
				Nu_i=mesh%ele(ie)%vertex(i) 
				AJOUT = .TRUE. 
				IF (leng_recherche/=0) THEN !on doit verifié que si le tableau n'est pas vide le noeud n'y ets pas deja
					DO k=1, leng_recherche
						IF(Nu_i == TabRecherch(k)) THEN 
							AJOUT = .FALSE.
						END IF 
					END DO 
				END IF 
				IF ( Nu_i==vertcs_original ) THEN ! noeud ou il manque une valeur 
					AJOUT = .FALSE. 
				END IF 
				IF (AJOUT) THEN 
					leng_Recherche = leng_Recherche +1 
					TabRecherch(leng_Recherche)= Nu_i
				END IF 
			END IF 
		END DO 
	END DO 
	
	
	! Mainteant on va trier le tableau par distance par rapport a vertcs_original  
	! Le nouveau tableau sera TabTridist
	ALLOCATE(TabTridist(leng_Recherche)) ; TabTridist = 0 
	x= mesh%Vertex(vertcs_original)%X(1) ; y = mesh%Vertex(vertcs_original)%X(2)
	chang = 0 
	DO WHILE (chang /= leng_Recherche) !on doit faire autant de changement qu'il y a de valeur dans le tableau 
		mintab = 10000  ;  numb_min = 0
		DO k = 1, leng_Recherche
			IF(TabRecherch(k)/= 0) THEN !on met à 0 les valeurs que l'on a déja enlevé 
				x_k= mesh%Vertex(TabRecherch(k))%X(1) ; y_k = mesh%Vertex(TabRecherch(k))%X(2)
				IF ( mintab > sqrt((x-x_k)**2+(y-y_k)**2) ) THEN 
					mintab = sqrt((x-x_k)**2+(y-y_k)**2)
					numb_min = k ! position du noeud dansle tableau 
				END IF
			END IF 
		END DO  
		chang = chang + 1
		TabTridist(chang) = TabRecherch(numb_min)
		TabRecherch(numb_min) = 0 
	END DO 
	
	
	!on va enlever de Tabtridist les noeuds ou il manque aussi une valeur 
	Do j = 1, leng_Recherche
		IF (TabTridist(j)/=0) THEN 
			Nu_j = 	TabTridist(j)	   
			
			choix = 1 ! table precedente 
			CALL LookForComp(Mesh,InIt1,Nu_j,valcomp1,choix) 
			choix = 0 ! table courante 
			CALL LookForComp(Mesh,InIt2,Nu_j,valcomp2,choix) 
			IF(InIt1 .AND. InIt2) THEN 
				IF(valcomp1==Nu_j .AND. valcomp2/=Nu_j) THEN !alors on a un nvx cmpt de bord
					IF(mesh%vertex(Nu_j)%EmptyComp) THEN ! la valeur n'a pas etait reconstruite
						TabTridist(j) = 0 
					END IF 				
				END IF
			ELSE IF (.NOT. InIt1 .AND. InIt2) THEN !nouveaux noeuds à l'interface 
				IF(mesh%vertex(Nu_j)%EmptyComp) THEN ! la valeur n'a pas etait reconstruite
					TabTridist(j) = 0 
				END IF 
			END IF
		END IF 
	  END DO 
	! maintenant on peut definir le stencil 
	

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

	IF(leng_tabvtc < 6) THEN
		PRINT*, "probleme sur la definition du stencil pas assez de neoud"
		STOP
	END IF 

	 ALLOCATE(A(6,6),AB(3,3)) ; ALLOCATE(B(6),B1(3),B2(3)) 
	 A = 0.d0 ; B = 0.d0 
	 B1= 0.d0  ; B2= 0.d0 ; AB = 0.d0 

	 !reconstruction quadratique pour la température, reconstruction linéaire pour le flux 
	 
	 DO i = 1, taille_stencil 
			Nu_i = Tab_stencil(i)
			x = mesh%Vertex(Nu_i)%X(1)  ; y =  mesh%Vertex(Nu_i)%X(2)
			
			!matrice du systeme lineaire reconstruction pour T 
			A(1,1)= A(1,1)+1        	 ; A(1,2)= A(1,2)+x    			; A(1,3)= A(1,3)+y
			A(1,4)= A(1,4)+x*y    	     ; A(1,5)= A(1,5)+x**2 			; A(1,6)= A(1,6)+y**2

			A(2,1)= A(2,1)+x       		 ; A(2,2)= A(2,2)+x**2 			; A(2,3)= A(2,3)+y*x
			A(2,4)= A(2,4)+(x**2)*y 	 ; A(2,5)= A(2,5)+x**3 			; A(2,6)= A(2,6)+x*(y**2)

			A(3,1)= A(3,1)+y        	 ; A(3,2)= A(3,2)+x*y           ; A(3,3)= A(3,3)+y**2
			A(3,4)= A(3,4)+x*(y**2) 	 ; A(3,5)= A(3,5)+(x**2)*y      ; A(3,6)= A(3,6)+y**3

			A(4,1)= A(4,1)+x*y           ; A(4,2)= A(4,2)+(x**2)*y 		; A(4,3)= A(4,3)+x*(y**2)
			A(4,4)= A(4,4)+(x**2)*(y**2) ; A(4,5)= A(4,5)+(x**3)*y 		; A(4,6)= A(4,6)+x*(y**3)

			A(5,1)= A(5,1)+x**2     	 ; A(5,2)= A(5,2)+x**3          ; A(5,3)= A(5,3)+(x**2)*y
			A(5,4)= A(5,4)+(x**3)*y 	 ; A(5,5)= A(5,5)+x**4          ; A(5,6)= A(5,6)+(x**2)*(y**2)
				
			A(6,1)= A(6,1)+y**2     	 ; A(6,2)= A(6,2)+(y**2)*x      ; A(6,3)= A(6,3)+y**3
			A(6,4)= A(6,4)+x*(y**3) 	 ; A(6,5)= A(6,5)+(x**2)*(y**2) ; A(6,6)= A(6,6)+y**4
 
			! second membre pour la temperature 
			B(1)=B(1)+Sol%pcg(Nu_i)
			B(2)=B(2)+Sol%pcg(Nu_i)*x
			B(3)=B(3)+Sol%pcg(Nu_i)*y
			B(4)=B(4)+Sol%pcg(Nu_i)*x*y
			B(5)=B(5)+Sol%pcg(Nu_i)*(x**2)
			B(6)=B(6)+Sol%pcg(Nu_i)*(y**2)
			
			!matrice du systeme lineaire reconstruction lineaire pour le flux 
			AB(1,:)=A(1,1:3);
			AB(2,:)=A(2,1:3);
			AB(3,:)=A(3,1:3);
			
			! second membre pour le flux en x 	
			B1(1)=B1(1)+Sol%bcg(Nu_i,1)
			B1(2)=B1(2)+Sol%bcg(Nu_i,1)*x
			B1(3)=B1(3)+Sol%bcg(Nu_i,1)*y
			
			! second membre pour le flux en y 
			B2(1)=B2(1)+Sol%bcg(Nu_i,2)
			B2(2)=B2(2)+Sol%bcg(Nu_i,2)*x
			B2(3)=B2(3)+Sol%bcg(Nu_i,2)*y
			
			
	END DO 
	
	  ALLOCATE(Am1(6,6),invAB(3,3)) 
      Am1=0.d0 ; invAB = 0.d0 
      CALL InverseMatrix_maxPivot(A,Am1)
      CALL InverseMatrix_maxPivot(AB,invAB)
            
            
      ALLOCATE(Xa(6),Xa1(3),Xa2(3))
	  Xa = 0 ; Xa1= 0 ; Xa2= 0 ! initial guess 

	 ! produit de l'inverse avec le second membre 
	 !pour la temperature on est sur un vecteur solution de taille 6 
      Do i = 1, 6
		Do k= 1, 6
			Xa(i)  = Xa(i) + Am1(i,k)*B(k)
			IF(i<=3) THEN 
				IF(k<=3) THEN 
					Xa1(i) = Xa1(i)+ invAB(i,k)*B1(k)
					Xa2(i) = Xa2(i)+ InvAB(i,k)*B2(k)
				END IF
			END IF
		END DO
	  END DO 

    !CAlcul du résidu 
    ALLOCATE(Cdef(6),Cdef2(3),Cdef3(3))
    Cdef = 0.d0 ; cdef2 = 0.d0 ; Cdef3 = 0.d0 
    
    Do i = 1, 6
		Do j = 1, 6
			Cdef(i)  = Cdef(i)  + Xa(j)*A(i,j)
			IF(i<=3) THEN 
				IF(j<=3) THEN 
					Cdef2(i) = Cdef2(i) + Xa1(j)*Ab(i,j)
					Cdef3(i) = Cdef3(i) + Xa2(j)*Ab(i,j)
				END IF 
			END IF 
		END DO 
	END DO 
	

	residu1 = 0.d0 ; residu2 = 0.d0 ; residu3 = 0.d0 
	Cdef = Cdef - B ; cdef2 = cdef2 - b1 ; cdef3 = cdef3 - b2 
	Do i = 1, 6
		residu1 = residu1 + Cdef(i)**2
		IF(i<=3) THEN 
			residu2 = residu2 + Cdef2(i)**2
			residu3 = residu3 + Cdef3(i)**2
		END IF 
	END DO  
	residu1= sqrt(residu1) ; residu2=sqrt(residu2) ; residu3 = sqrt(residu3) 
	
	!IF(residu1 > 0.0001 .OR. residu2 > 0.0001  .OR. residu3 > 0.0001 ) THEN 
		!PRINT*,"Le residu de la reconstruction est trop élévé"
	!	PRINT*,"residu temperature ", residu1
	!	PRINT*,"residu flux  x", residu2
	!	PRINT*,"residu flux y" , residu3
		!STOP
	!END IF 

	 val1 = Xa(1) +Xa(2)*xx + Xa(3)*yy + Xa(4)*xx*yy  + Xa(5)*(xx**2)  + Xa(6)*(yy**2)
	!val2 = -lambda*( Xa(2) + Xa(4)*yy  + 2*Xa(5)*xx)
	!val3 = -lambda*( Xa(3) + Xa(4)*xx  + 2*Xa(6)*yy)

	!reconstruction separé pour le flux independante de la reconstruction sur T 
	!plus de precision avec reconstruction séparé 
	val2 =  Xa1(1) + Xa1(2)*xx + Xa1(3)*yy
	val3 =  Xa2(1) + Xa2(2)*xx + Xa2(3)*yy

	! on a reconstruit la valeur, elle n'est plus manquante 
	Do i =1, mesh%nv
		IF(Mesh%vertex(i)%num == vertcs_original) THEN 
			Mesh%vertex(i)%EmptyComp = .FALSE. 
		END IF 
	END DO 

  !solution exacte
  ALLOCATE(Exa(3)) ; Exa=0
  
  !on a besoin de la reference d'un element pour calculer la solution exacte 
	DO Nie= 1, mesh%Ne ! boucle sur le nombre d'élément
		DO k = 1, mesh%ele(Nie)%nNodesPerElement ! boucle sur le nombre de noeud associé à l'élement
			Nu_k=mesh%ele(Nie)%vertex(k) ! noeud sur l'element
			IF(Nu_k == vertcs_original  .AND. mesh%ele(Nie)%ref == 1 ) THEN
					ele_original = mesh%ele(Nie) ! on veut l'enregistrement d'un element associé au noeud ou on a la valeur manquante pour l'utiliser pour le calcul de la solution exacte 
			END IF
		END DO 
	END DO 
	CALL ExactSolJumpt(ele_original%L, Exa,Mesh%Vertex(vertcs_original)%X, &
								Data%T-data%Dt, ele_original, nDim)	  
	
	!methode qui semble reconstruire correctement le flux 			
	CALL ReconstructionVariableT(mesh,Sol,vertcs_original,val1,val2,val3)
	!CALL ReconstructionVariableB(mesh,Sol,vertcs_original,val2,val3)
		
	Data%ErrorT  = Data%ErrorT   + abs(Exa(3)-val1)
	Data%ErrorBx = Data%ErrorBx  + abs(Exa(1)-val2)
	Data%ErrorBy = Data%ErrorBy  + abs(Exa(2)-val3)
	
  DEALLOCATE(Tab_stencil) ;   DEALLOCATE(A,Am1) ; DEALLOCATE(B,B1,B2) 
  DEALLOCATE(Xa,Xa1,Xa2) ; DEALLOCATE(TabTridist,TabRecherch)

END SUBROUTINE PolyConstr7
!========================================!

!==================================!
SUBROUTINE PolyConstrQR(Data,Mesh,Sol,vertcs,val1,val2,val3)
!==================================!

  IMPLICIT NONE

  type(MeshType)     , INTENT(INOUT)  	:: Mesh
  type(DataStructure), INTENT(INOUT)  	:: Data
  REAL*8             ,INTENT(OUT)  		:: val1,val2,val3
  INTEGER            , INTENT(INOUT) 	:: vertcs
  type(SolStructure) , INTENT(IN)  		:: Sol

  !----------------------------------------
  INTEGER                               :: j, Ne, i,valcomp1, lengTab, k, leng_recherche,taille_stencil
  INTEGER                               :: valcomp2, ie, nu_i, choix, Ndim, nu_j, ib, ils 
  INTEGER                               :: Nu_k , Nv, leng_tabvtc, vertcs_original
  LOGICAL                               :: InIt1, InIt2, test, TEst_ele, AJOUT , test_edg
  LOGICAL                               :: Test1 , Test2 ,IFoundYou 
  type(Element)                         :: ele_original 
  REAL*8                                :: xx, yy, x, y, lambda, max, test_val, residu 
  REAL*8                                :: U1,U2,U3, eps,resu1,resu2,pente
  REAL*8                                :: x_test, y_test, x_k, x_j, x_i,y_i , y_k , y_j, dist,mintab 
  REAL*8                                :: residu1, residu2, residu3 
  INTEGER                               :: test_tab, tour,lengdif, numb_min, chang, cmpt , Nie,degre_rsct, taille 
  INTEGER, DIMENSION(:)  , ALLOCATABLE  :: Tab_stencil, Tab_vtc, TabRecherch, TabTridist
  REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: A, Am1, AB,invAB
  REAL*8, DIMENSION(:)  , ALLOCATABLE   :: B, Bdef , Xa, Cdef, test_vertcs, Exa, Exa_prec, B1,B2, Xa1,Xa2,cdef2, cdef3 
  type(Edge)    						:: Edg
  type(element) 						:: eleM, eleP

  !----------------------------------------

	Ne = Mesh%Ne ; Ndim = Data%Ndim ; Nv = Sol%Nv ; lengdif=0 ; leng_recherche = 0
	vertcs_original = vertcs ; dist = 0.d0
	degre_rsct = 6 ; taille_stencil = 12 !nombre de noeuds pour la résolution 
	lengTab = 0 ; tour = 0 ; Test_ele = .FALSE.
	xx = mesh%Vertex(vertcs_original)%X(1)  ; yy =  mesh%Vertex(vertcs_original)%X(2)
  
  
	lambda = Icing_Lambda1 ! la reconstruction s'effectue dans la zone 1 

	ALLOCATE(TabRecherch(Nv)) ; TabRecherch = 0 ; leng_recherche = 0
	  
	!on enregistre tous les noeuds qui sont dans la zone 1 
	!car on reconstruit une valeur à gauche 
	DO ie= 1, Ne 
		DO i = 1, mesh%ele(ie)%nNodesPerElement 
			IF (mesh%ele(ie)%ref == IcingZone1_Ref .AND. mesh%ele(ie)%vertex(i)/=vertcs_original) THEN ! si on se trouve dans la bonne zone et pas sur le noeud ou il manque une valeur
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
	
	
	! Mainteant on va trier le tableau par distance par rapport a vertcs_original  
	! Le nouveau tableau sera TabTridist
	ALLOCATE(TabTridist(leng_Recherche)) ; TabTridist = 0 
	x= mesh%Vertex(vertcs_original)%X(1) ; y = mesh%Vertex(vertcs_original)%X(2)
	chang = 0 
	DO WHILE (chang /= leng_Recherche) !on doit faire autant de changement qu'il y a de valeur dans le tableau 
		mintab = 10000  ;  numb_min = 0
		DO k = 1, leng_Recherche
			IF(TabRecherch(k)/= 0) THEN !on met à 0 les valeurs que l'on a déja enlevé 
				x_k= mesh%Vertex(TabRecherch(k))%X(1) ; y_k = mesh%Vertex(TabRecherch(k))%X(2)
				IF ( mintab > sqrt((x-x_k)**2+(y-y_k)**2) ) THEN 
					mintab = sqrt((x-x_k)**2+(y-y_k)**2)
					numb_min = k ! position du noeud dansle tableau 
				END IF
			END IF 
		END DO  
		chang = chang + 1
		TabTridist(chang) = TabRecherch(numb_min)
		TabRecherch(numb_min) = 0 
	END DO 
	
	
	!on va enlever de Tabtridist les noeuds ou il manque aussi une valeur 
	Do j = 1, leng_Recherche
		IF (TabTridist(j)/=0) THEN 
			Nu_j = 	TabTridist(j)	   
			
			choix = 1 ! table precedente 
			CALL LookForComp(Mesh,InIt1,Nu_j,valcomp1,choix) 
			choix = 0 ! table courante 
			CALL LookForComp(Mesh,InIt2,Nu_j,valcomp2,choix) 
			IF(InIt1 .AND. InIt2) THEN 
				IF(valcomp1==Nu_j .AND. valcomp2/=Nu_j) THEN !alors on a un nvx cmpt de bord
					IF(mesh%vertex(Nu_j)%EmptyComp) THEN ! la valeur n'a pas etait reconstruite
						TabTridist(j) = 0 
					END IF 				
				END IF
			ELSE IF (.NOT. InIt1 .AND. InIt2) THEN !nouveaux noeuds à l'interface 
				IF(mesh%vertex(Nu_j)%EmptyComp) THEN ! la valeur n'a pas etait reconstruite
					TabTridist(j) = 0 
				END IF 
			END IF
		END IF 
	END DO 
	! maintenant on peut definir le stencil 
	

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
	ALLOCATE(A(taille_stencil,6)) ; ALLOCATE(B(taille_stencil)) 
	A = 0.d0 ; B = 0.d0 

	!reconstruction quadratique pour la température, reconstruction linéaire pour le flux 
	 
	DO i = 1, taille_stencil 
		Nu_i = Tab_stencil(i)
		x = mesh%Vertex(Nu_i)%X(1)  ; y =  mesh%Vertex(Nu_i)%X(2)
			
		!matrice du systeme lineaire reconstruction pour T 
		A(i,1)= 1        	 ; A(i,2)= x    			; A(i,3)= y
		A(i,4)= x*y    	     ; A(i,5)= x**2 			; A(i,6)= y**2

			! second membre pour la temperature 
		B(i)=Sol%pcg(Nu_i)
			
	END DO 
            
    ALLOCATE(Xa(6))
	Xa = 0 ! initial guess 
	
	CALL SolveSystemQR(A,b,Xa,taille_stencil,6)
	xx = mesh%Vertex(vertcs_original)%X(1)  ; yy =  mesh%Vertex(vertcs_original)%X(2)

	val1 = Xa(1) +Xa(2)*xx + Xa(3)*yy + Xa(4)*xx*yy  + Xa(5)*(xx**2)  + Xa(6)*(yy**2)
	val2 = -lambda*( Xa(2) + Xa(4)*yy  + 2*Xa(5)*xx)
	val3 = -lambda*( Xa(3) + Xa(4)*xx  + 2*Xa(6)*yy)
	
	ALLOCATE(B1(taille_stencil),B2(taille_stencil)) 
	B1 = 0.d0 ; B2 = 0.d0 

	!reconstruction quadratique pour la température, reconstruction linéaire pour le flux 
	!autre methode avec reconstruction indépendante sur Beta 
	DO i = 1, taille_stencil 
		Nu_i = Tab_stencil(i)
		x = mesh%Vertex(Tab_stencil(i))%X(1)  ; y =  mesh%Vertex(Tab_stencil(i))%X(2)
		B1(i)=Sol%bcg(Nu_i,1)  ;  B2(i)=Sol%bcg(Nu_i,2)

	END DO 
            
    ALLOCATE(Xa1(3),Xa2(3))
	Xa1 = 0 ; Xa2 = 0 ! initial guess 
	CALL SolveSystemQR(A(:,1:3),b1,Xa1,taille_stencil,3)
	!val2 = Xa1(1) + Xa1(2)*xx  + Xa1(3)*yy

	CALL SolveSystemQR(A(:,1:3),b2,Xa2,taille_stencil,3)
	!val3 = Xa2(1) + Xa2(2)*xx  + Xa2(3)*yy

	IF(icas == 001) THEN 
		val3 = 0.d0
	END IF 
	 
	! on a reconstruit la valeur, elle n'est plus manquante 
	IFoundYou = .FALSE. ; i = 0 
	DO WHILE (.NOT. IFoundYou .AND. i < Mesh%Nv ) 
		i = i + 1 
		IF(Mesh%vertex(i)%num == vertcs_original) THEN 
			Mesh%vertex(i)%EmptyComp = .FALSE. 
			IFOUNDYOU  = .TRUE. 
		END IF 
	END DO 	
	
	DEALLOCATE(Xa,Xa1,Xa2,A,B,b1,b2,Tab_stencil)

END SUBROUTINE PolyConstrQR
!========================================!
	
!=====================================================================================!
 SUBROUTINE JumpTest(edg, edgM, edgP, eleM, eleP, Data, xo,yo,SOl,vertcs_original,ILS,val1)
!=====================================================================================!

    IMPLICIT NONE
    type(Edge)              , INTENT(IN)    :: Edg
    type(Edge), DIMENSION(:), INTENT(IN)    :: edgM, edgP
    type(element)           , INTENT(IN)    :: eleM, eleP
    type(DataStructure)     , INTENT(IN)    :: Data
    real*8					, intent(out)   :: VAL1 
    real*8					, INTENT(IN)    :: xo,yo 
	type(SolStructure)      , INTENT(IN)  	:: Sol
	INTEGER					, INTENT(IN) 	:: vertcs_original, ILS

    !----------------------------------------------------------
    type(edge) :: edgIntM, edgIntP
    type(element) :: ELE_inter
    !----------------------------------------------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Np, Nm, Xq_M, Xq_P, Exa_p, Exa_m
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: XqM, XqP, X_M, X_P
    REAL*8, DIMENSION(:,:)  , ALLOCATABLE :: XeM, XeP
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: XxqP, XxqM
    REAL*8, DIMENSION(data%nDim)   :: d_qP, n_qP, d_qM, n_qM
    REAL*8, DIMENSION(data%nDim-1) :: tntM, tntP
    REAL*8, DIMENSION(data%nDim,data%nDim-1) :: t_qP, t_qM
    REAL*8, DIMENSION(data%nDim+1) :: exa1_q, exa2_q, exa1_q_prec, exa2_q_prec
    REAL*8, DIMENSION(data%nDim)   :: x_qP, x_qM
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: Exa_p_prec, Exa_M_prec
    !-------------------------------------------------------------
    REAL*8  :: Nrm, wq, jt, pD_M, pD_P, nnt
    REAL*8  :: a, ttM, ttP,res
    REAL*8  :: LS_qP, LS_qM, lambdaM, lambdaP,res1
    REAL*8  :: L, int,res2
    REAL*8  :: SF_iM, SF_iP
    !------------------------------------------------------------
    INTEGER :: i,iq, i1, j1, j, id, jd, ii
    INTEGER :: iM, iP
    INTEGER :: nDim, nNodesPerElement, nQuad, nEdgesPerElement, nVar
    !--------------------------------------------------------

    nDim = Data%nDim ; nNodesPerElement = eleM%nNodesPerElement
    nVar = nDim + 1 ; nEdgesPerElement = eleM%nEdgesPerElement ; jt =0
    
    
    ALLOCATE ( Np(nDim), Nm(nDim), Xq_M(nDim), Xq_P(nDim))
    ALLOCATE (  XxqP(nDim), XxqM(nDim) ) ; ALLOCATE ( Exa_p(nVar), Exa_M(nVar)  )
    ALLOCATE ( Exa_p_prec(nVar), Exa_M_prec(nVar)  )
    ALLOCATE ( XqM(nNodesPerElement,nDim), XqP(nNodesPerElement,nDim) )
    ALLOCATE ( X_M(nNodesPerElement,nDim), X_P(nNodesPerElement,nDim) )
    ALLOCATE ( XeM(nNodesPerElement,nDim), XeP(nNodesPerElement,nDim) )

    XeM = eleM%Coor ; XeP = eleP%Coor

    CALL IdentifyIFace(eleM, edg, iM, nDim) ; CALL IdentifyIFace(eleP, edg, iP, nDim)

    edgIntM = edgM(iM) ; edgIntP = edgP(iP)
    lambdaM=eleM%lambda ; lambdaP= eleP%lambda
    Np = -eleP%N(iP,:)
    Nrm = 0.d0
    DO id = 1, nDim
       Nrm = Nrm + Np(id)**2
    END DO
    Nrm = sqrt(Nrm) ; L=Nrm
    Np = Np/Nrm
    Nm = -Np

    CALL QuadratureOnEdge(xqM, eleM, edg, iM, nDim)
    CALL QuadratureOnEdge(xqP, eleP, edg, iP, nDim)

    nnt = 0.d0
    CALL GetNQuad(edg%eleType, nQuad)
    DO iq = 1, nQuad !boucle sur le nombre de points de quadrature
    
		IF (xo==XeM(iq,1) .AND. yo == XeM(iq,2))THEN 
		
		   xq_M = xqM(iq,:)  ;  xq_P = xqP(iq,:)
		   x_qP = 0.d0       ;  x_qM = 0.d0
		   
		   DO i = 1, nNodesPerElement
			  CALL SFEval(eleM%eleType, SF_iM, xq_M, i) ; CALL SFEval(eleP%eleType, SF_iP, xq_P, i)
			  DO id = 1, nDim
				 x_qP(id) = x_qP(id) + SF_iP*eleP%Coor(i,id)
				 x_qM(id) = x_qM(id) + SF_iM*eleM%Coor(i,id)
			  END DO
		   END DO
		   
		   CALL computeNodalDistance(x_qP, d_qP, LS_qP, Data, iLS)
		   CALL getNormalAndTangents(d_qP, n_qP, t_qP, nP, nnt, tntP, nDim)

		   d_qM = -d_qP

		   CALL ExactSolJumpt(eleP%L, exa1_q, x_qP + d_qP, Data%T-Data%dt, eleP, nDim)
		   CALL ExactSolJumpt(eleM%L, exa2_q, x_qM - d_qM, Data%T-Data%dt, eleM, nDim)
		   pD_P = exa1_q(nDim+1)!-(1.d0/lambdaP)*(exa1_q(1)*d_qP(1)+exa1_q(2)*d_qP(2))
		   pD_M = exa2_q(nDim+1)!-(1.d0/lambdaM)*(exa2_q(1)*d_qm(1)+exa2_q(2)*d_qm(2))
		   
		   x_qp(1)=xo ; x_qp(2)=yo
		   CALL ExactSolJumpt(eleP%L, exa1_q, x_qP, Data%T-Data%dt, eleP, nDim)
		   CALL ExactSolJumpt(eleM%L, exa2_q, x_qP, Data%T-Data%dt, eleM, nDim)
		   pD_P = pD_P+(1.d0/lambdaP)*(exa1_q(1)*d_qP(1)+exa1_q(2)*d_qP(2))
		   pD_M = pD_M+(1.d0/lambdaM)*(exa2_q(1)*d_qm(1)+exa2_q(2)*d_qm(2))
		   
		   
		   jt = pD_P - pD_M
		   
		END IF 
	end do 
	val1 = jt + Sol%Pcg(vertcs_original)
	
	! a mettre dans le code ou l'on veut utiliser la routine 
		
	!ib = 0 ; TEST_edg = .FALSE. 
	!DO WHILE (ib <mesh%nSurB .AND. .NOT. TEST_edg) !mesh%nSurB est le nombre d'arete sur la surrogate
     !     ib = ib+1 
      !    iLS= mesh%EmbBoundary(ib)%EmbLS
       !   edg = mesh%embBoundary(ib) 
        ! IF(edg%vertex(1)==vertcs_original .OR. edg%vertex(2)==vertcs_original) THEN
		!	TEST_edg = .TRUE.
        ! END IF 
    !END DO 
	
	 !IF (mesh%ele(edg%tri(2))%ref == IcingZone2_ref) THEN
		!eleP = mesh%ele(edg%tri(1))
		!eleM = mesh%ele(edg%tri(2))
	 !ELSEReconstructionVariableBx(mesh,Sol,vertcs_original,val1)
		!eleP = mesh%ele(edg%tri(2))
		!eleM = mesh%ele(edg%tri(1))
	!END IF
	
  END SUBROUTINE JumpTest
  !======================================================================!
  
  
 !======================================================================!
 SUBROUTINE ReconstructionVariableT(mesh,Sol,vertcs_original,val1,val2,val3)
 !======================================================================!

  
  IMPLICIT NONE
  
    type(MeshType)     		, INTENT(IN)    :: mesh
    REAL*8					, INTENT(OUT)   :: val1,val2,val3
	type(SolStructure)      , INTENT(IN)  	:: Sol
	INTEGER					, INTENT(IN) 	:: vertcs_original

    !----------------------------------------------------------
    REAL*8, DIMENSION(:,:)   , ALLOCATABLE :: A, AAT , Am1 
    REAL*8, DIMENSION(:)     , ALLOCATABLE :: B, AtB, Xa
    INTEGER, DIMENSION(:)    , ALLOCATABLE :: TabRef1, TabDist, Stencil 
	LOGICAL :: Test , InIt1, InIt2 
    !-------------------------------------------------------------
    REAL*8  :: dist , mindist, xx , yy , x_k, y_k , xi, yi 
    REAL*8  :: LAMBDA
 
    !------------------------------------------------------------
    INTEGER :: ie,j, Nu_j , i, k , post , choix , len_stencil , rempli 
    INTEGER :: Len_TabRef1, len_tabDist, valcomp1, valcomp2 

    !--------------------------------------------------------


	lambda = Icing_Lambda1
	  
	ALLOCATE(TabRef1(Mesh%Nv)) ; TabRef1 = 0 ; Len_TabRef1 = 0 
	
	Do ie = 1, Mesh%Ne 
		IF(mesh%ele(ie)%ref ==1) THEN 
			Do j = 1,mesh%ele(ie)%nNodesPerElement
				Nu_j = mesh%ele(ie)%Vertex(j) 
				IF(len_TabRef1 /= 0) THEN 
					TEST = .FALSE. 
					DO k =1, len_TabRef1
						IF(TabRef1(k)==Nu_j) THEN
							TEST = .TRUE.
						END IF 
					END DO 
					IF(.NOT. TEST .AND. Nu_j /= vertcs_original) THEN 	
						Len_TabRef1 = Len_TabRef1 + 1 
						TabRef1(Len_TabRef1)=Nu_j
					END IF 
				ELSE
					Len_TabRef1 = Len_TabRef1 + 1 
					TabRef1(Len_TabRef1)=Nu_j
				END IF 
			END DO 
		END IF 
	END DO 

	! On va maintenant trier TabRef1 par rapport à sa distance avec vertcs_original 
	ALLOCATE(TabDist(Len_TabRef1)) ; TabDist = 0 ; len_tabDist = 0
	xx=  Mesh%Vertex(vertcs_original)%X(1) ; yy = Mesh%Vertex(vertcs_original)%X(2)

	DO j = 1 , Len_TabRef1
	mindist = 100000
	  DO k = 1, Len_TabRef1
		IF( TabRef1(k) /= 0) THEN 
			x_k=  Mesh%Vertex(TabRef1(k))%X(1) ; y_k = Mesh%Vertex(TabRef1(k))%X(2)
			dist = sqrt((xx-x_k)**2+(yy-y_k)**2)
			IF(dist<mindist) THEN 
				mindist = dist 
				post = k 
			END IF 
		END IF 
	  END DO 
	len_tabDist = len_tabDist + 1 
	TabDist(len_tabDist) = TabRef1(post)
	TabRef1(post) = 0
  END DO 
  
  !on fait attention que les noeuds choisi non pas aussi une valeur manquante sur l'interface 
	Do j = 1, len_tabDist
		IF (TabDist(j)/=0) THEN 
			Nu_j = 	TabDist(j)	   
			choix = 1 ! table precedente 
			CALL LookForComp(Mesh,InIt1,Nu_j,valcomp1,choix) !renvoie la valeur du complementaire
			choix = 0 ! table courante 
			CALL LookForComp(Mesh,InIt2,Nu_j,valcomp2,choix) !renvoie la valeur du complementaire
			IF (.NOT. InIt1 .AND. InIt2) THEN ! i.e. que c'est un nouveau noeud à l'interface
				!IF(mesh%vertex(Nu_j)%EmptyComp) THEN ! la valeur n'a pas etait reconstruite
					TabDist(j) = 0 
				!END IF 
			END IF
			IF (InIt2) THEN ! i.e. que c'est un nouveau noeud à l'interface
				TabDist(j) = 0 
			END IF
		END IF 
	 END DO 
	  
	  ! maintenant on choisi le nombre de noeuds que l'on veut prendre 
	  len_stencil = 12
	  ALLOCATE(STENCIL(len_stencil)) ; STENCIL = 0 ; k=0
	  rempli = 0
	  
	  DO WHILE (rempli < len_stencil .AND. k < Len_TabRef1) 
	  k=k+1
		  IF(TabDist(k)/=0) THEN 
			rempli = rempli + 1
			STENCIL(rempli) = TabDist(k)
		  ENd IF 
	  END DO 
	  
	  IF(rempli < 6) THEN 
		PRINT*,"Error PolyReconstruction"
		PRINT*,"taille du stencil recherché non disponible"
		STOP
	  END IF 
	  len_stencil= rempli
	  
	  !creation du systeme linéaire 
	  ALLOCATE(A(3*len_stencil, 6)) ; ALLOCATE(B(3*len_stencil))
	  A = 0.d0 ; B =0.d0
	  
	  DO i = 1, len_stencil 
	  
		xi=  Mesh%Vertex(STENCIL(i))%X(1) ; yi = Mesh%Vertex(STENCIL(i))%X(2)
		
		A(i,1) = 1     ; A(i,2) = xi    ; A(i,3) = yi
		A(i,4) = xi*yi ; A(i,5) = xi**2 ; A(i,6) = yi**2
		
		IF(lambda > 5 ) THEN 

			A(i+len_stencil,1) = 0  ; A(i+len_stencil,2) = 1    ; A(i+len_stencil,3) = 0
			A(i+len_stencil,4) = yi ; A(i+len_stencil,5) = 2*xi ; A(i+len_stencil,6) = 0
			
			
			A(i+2*len_stencil,1) = 0   ; A(i+2*len_stencil,2) =0    ; A(i+2*len_stencil,3) = 1
			A(i+2*len_stencil,4) = xi  ; A(i+2*len_stencil,5) = 0   ; A(i+2*len_stencil,6) = 2*yi
			
			
			B(i) = Sol%Pcg(STENCIL(i))
			B(i+len_stencil) = -1.d0/lambda*Sol%bcg(STENCIL(i),1)
			b(i+2*len_stencil) = -1.d0/lambda*Sol%bcg(STENCIL(i),2)
		
		ELSE 
		
			A(i+len_stencil,1) = 0  ; A(i+len_stencil,2) = -lambda    ; A(i+len_stencil,3) = 0
			A(i+len_stencil,4) = -lambda*yi ; A(i+len_stencil,5) = -lambda*2*xi ; A(i+len_stencil,6) = 0
			
			
			A(i+2*len_stencil,1) = 0   ; A(i+2*len_stencil,2) =0    ; A(i+2*len_stencil,3) = -lambda
			A(i+2*len_stencil,4) = -lambda*xi  ; A(i+2*len_stencil,5) = 0   ; A(i+2*len_stencil,6) = -lambda*2*yi
			
			
			B(i) = Sol%Pcg(STENCIL(i))
			B(i+len_stencil) = Sol%bcg(STENCIL(i),1)
			B(i+2*len_stencil) = Sol%bcg(STENCIL(i),2)
		
		END IF 
		
	  END DO 
	  
	  !Calcul de la multiplication par la transposé 
	  ALLOCATE(AAt(6, 6)) ; ALLOCATE(AtB(6))
	  AAt = 0.d0 ; atB = 0.d0
	
	 Do i = 1 , 6
		Do j = 1, 6  
			Do k =1 , 3*len_stencil
				AAt(i,j) = AAt(i,j) + A(k,i)*A(k,j) 
			END DO
		END DO
	 END DO
	  
	Do i = 1 , 6
		Do k =1,3*len_stencil
			AtB(i) = AtB(i) + A(k,i)*B(k) 
		END DO
	END DO
	
	!calcul de l'inverse
	ALLOCATE(Am1(6,6)) ; Am1 = 0.d0
    CALL InverseMatrix_maxPivot(AAt,Am1)
    
    !calcul des coeff du polynome reconstruit 
    ALLOCATE(XA(6)) ; Xa = 0.d0 
    Do i = 1, 6
		Do k =1, 6
			Xa(i) = Xa(i) + Am1(i,k)*AtB(k) 
		END DO
	END DO
	
	xx=  Mesh%Vertex(vertcs_original)%X(1) ; yy = Mesh%Vertex(vertcs_original)%X(2)
	val1 = Xa(1) + Xa(2)*xx + Xa(3)*yy + Xa(4)*xx*yy + Xa(5)*xx*xx + Xa(6)*yy*yy
	val2 = -lambda*(Xa(2) + Xa(4)*yy + Xa(5)*xx*2)
	val3 = -lambda*(Xa(3) + Xa(4)*xx + Xa(6)*yy*2)


 END SUBROUTINE ReconstructionVariableT
 !======================================================================!

  
 !======================================================================!
 SUBROUTINE ReconstructionVariableB(mesh,Sol,vertcs_original,valBx,valBy)
 !======================================================================!

  
  IMPLICIT NONE
  
    type(MeshType)     		, INTENT(IN)    :: mesh
    REAL*8					, INTENT(OUT)   :: valBx,valBy
	type(SolStructure)      , INTENT(IN)  	:: Sol
	INTEGER					, INTENT(IN) 	:: vertcs_original

    !----------------------------------------------------------
    REAL*8, DIMENSION(:,:)   , ALLOCATABLE :: A, AAT , Am1 
    REAL*8, DIMENSION(:)     , ALLOCATABLE :: B, AtB, Xa
    INTEGER, DIMENSION(:)    , ALLOCATABLE :: TabRef1, TabDist, Stencil 
	LOGICAL :: Test , InIt1, InIt2 
    !-------------------------------------------------------------
    REAL*8  :: dist , mindist, xx , yy , x_k, y_k , xi, yi 
    REAL*8  :: LAMBDA
 
    !------------------------------------------------------------
    INTEGER :: ie,j, Nu_j , i, k , post , choix , len_stencil , rempli 
    INTEGER :: Len_TabRef1, len_tabDist, valcomp1, valcomp2 

    !--------------------------------------------------------


	lambda = Icing_Lambda1
	ALLOCATE(TabRef1(Mesh%Nv)) ; TabRef1 = 0 ; Len_TabRef1 = 0 
	
	Do ie = 1, Mesh%Ne 
		IF(mesh%ele(ie)%ref ==1) THEN 
			Do j = 1,mesh%ele(ie)%nNodesPerElement
				Nu_j = mesh%ele(ie)%Vertex(j) 
				IF(len_TabRef1 /= 0) THEN 
					TEST = .FALSE. 
					DO k =1, len_TabRef1
						IF(TabRef1(k)==Nu_j) THEN
							TEST = .TRUE.
						END IF 
					END DO 
					IF(.NOT. TEST .AND. Nu_j /= vertcs_original) THEN 	
						Len_TabRef1 = Len_TabRef1 + 1 
						TabRef1(Len_TabRef1)=Nu_j
					END IF 
				ELSE
					Len_TabRef1 = Len_TabRef1 + 1 
					TabRef1(Len_TabRef1)=Nu_j
				END IF 
			END DO 
		END IF 
	END DO 

	! On va maintenant trier TabRef1 par rapport à sa distance avec vertcs_original 
	ALLOCATE(TabDist(Len_TabRef1)) ; TabDist = 0 ; len_tabDist = 0
	xx=  Mesh%Vertex(vertcs_original)%X(1) ; yy = Mesh%Vertex(vertcs_original)%X(2)

	DO j = 1 , Len_TabRef1
	mindist = 100000
	  DO k = 1, Len_TabRef1
		IF( TabRef1(k) /= 0) THEN 
			x_k=  Mesh%Vertex(TabRef1(k))%X(1) ; y_k = Mesh%Vertex(TabRef1(k))%X(2)
			dist = sqrt((xx-x_k)**2+(yy-y_k)**2)
			IF(dist<mindist) THEN 
				mindist = dist 
				post = k 
			END IF 
		END IF 
	  END DO 
	len_tabDist = len_tabDist + 1 
	TabDist(len_tabDist) = TabRef1(post)
	TabRef1(post) = 0
  END DO 
  
  !on fait attention que les noeuds choisi non pas aussi une valeur manquante sur l'interface 
	Do j = 1, len_tabDist
		IF (TabDist(j)/=0) THEN 
			Nu_j = 	TabDist(j)	   
			choix = 1 ! table precedente 
			CALL LookForComp(Mesh,InIt1,Nu_j,valcomp1,choix) !renvoie la valeur du complementaire
			choix = 0 ! table courante 
			CALL LookForComp(Mesh,InIt2,Nu_j,valcomp2,choix) !renvoie la valeur du complementaire
			IF (.NOT. InIt1 .AND. InIt2) THEN ! i.e. que c'est un nouveau noeud à l'interface
					TabDist(j) = 0 
			END IF
		END IF 
	 END DO 
	  
	  ! maintenant on choisi le nombre de noeuds que l'on veut prendre 
	  len_stencil = 24
	  ALLOCATE(STENCIL(len_stencil)) ; STENCIL = 0 ; k=0
	  rempli = 0
	  
	  DO WHILE (rempli < len_stencil) 
	  k=k+1
		  IF(TabDist(k)/=0) THEN 
			rempli = rempli + 1
			STENCIL(rempli) = TabDist(k)
		  ENd IF 
	  END DO 
	  
	  !creation du systeme linéaire 
	  ALLOCATE(A(2*len_stencil,5)) ; ALLOCATE(B(2*len_stencil))
	  A = 0.d0 ; B =0.d0
	  
	  DO i = 1, len_stencil 
	  
		xi=  Mesh%Vertex(STENCIL(i))%X(1) ; yi = Mesh%Vertex(STENCIL(i))%X(2)
		
		A(i,1) = 1  ; A(i,2) = 0    ; A(i,3) = yi
		A(i,4) = 2*xi  ; A(i,5) = 0    
		
		
		A(i+len_stencil,1) = 0  ; A(i+len_stencil,2) = 1   ; A(i+len_stencil,3) = xi
		A(i+len_stencil,4) = 0  ; A(i+len_stencil,5) = 2*yi    

		B(i) = Sol%Bcg(STENCIL(i),1)
		B(i+len_stencil) = Sol%Bcg(STENCIL(i),2)

	  END DO 
	  
	  !Calcul de la multiplication par la transposé 
	  ALLOCATE(AAt(5, 5)) ; ALLOCATE(AtB(5))
	  AAt = 0.d0 ; atB = 0.d0
	
	 Do i = 1 , 5
		Do j = 1, 5  
			Do k =1 , 2*len_stencil
				AAt(i,j) = AAt(i,j) + A(k,i)*A(k,j) 
			END DO
		END DO
	 END DO
	  
	Do i = 1 , 5
		Do k =1,2*len_stencil
			AtB(i) = AtB(i) + A(k,i)*B(k) 
		END DO
	END DO
	
	!calcul de l'inverse
	ALLOCATE(Am1(5,5)) ; Am1 = 0.d0
    CALL InverseMatrix_maxPivot(AAt,Am1)
    
    !calcul des coeff du polynome reconstruit 
    ALLOCATE(XA(5)) ; Xa = 0.d0 
    Do i = 1, 5
		Do k =1,5
			Xa(i) = Xa(i) + Am1(i,k)*AtB(k) 
		END DO
	END DO
	
	xx=  Mesh%Vertex(vertcs_original)%X(1) ; yy = Mesh%Vertex(vertcs_original)%X(2)
	valBx = Xa(1) + 2*xx*Xa(4) + Xa(3)*yy
	valBy = Xa(2) + 2*yy*Xa(5) + Xa(3)*xx
	

 END SUBROUTINE ReconstructionVariableB
 !======================================================================!
 
END MODULE PolyReconstruction
