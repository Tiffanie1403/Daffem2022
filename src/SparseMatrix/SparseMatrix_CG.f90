MODULE SparseMatrix_CG

  USE Types
  USE PublicVar
  USE SparseMatrix_DG

  IMPLICIT NONE

CONTAINS

  !===================================================!
  SUBROUTINE AssembleRHSCG(F, Floc, ele, Nv, Ne, nVar)

  !===================================================!
  !*************************************************************************
  !   Adds the local matrix to the global matrix discretization of the second member
  !   Continuous Galerkin case
  !
  !  Parameters:
  !
  !    Input,Output, REAL*8, DIMENSION(:) F, Structure for the global second member
  !    Input, REAL*8, DIMENSION(:,:) Floc, Local matrix for the second member
  !    Input, type(element) ele, element considerated
  !    Input, INTEGER  Nv, Number of nodes (vertices)
  !    Input, INTEGER  Ne, Number of elements
  !    Input, INTEGER  nVar, Number of variables
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: F
    REAL*8, DIMENSION(:), INTENT(IN)    :: Floc
    type(element)       , INTENT(IN)    :: ele
    INTEGER             , INTENT(IN)    :: Nv, Ne, nVar
    !---------------------------------------------------
    INTEGER :: nNodesPerElement
    INTEGER :: i, iG, iL, NU_i, ivar
    !---------------------------------

    nNodesPerElement = ele%nNodesPerElement

    DO i = 1, nNodesPerElement
       NU_i = ele%SurVertex(i)
       DO ivar = 1, nVar
          iG = NU_i + (ivar-1)*Nv
          iL = i + (ivar-1)*nNodesPerElement
          F(iG) = F(iG) + Floc(iL)
       END DO
    END DO

  END SUBROUTINE AssembleRHSCG
  !===================================================!

  !========================================================!
  SUBROUTINE AssembleSparseMatrixCG(Mat, Aloc, ele, Nv, Ne)
    !========================================================!

    !*************************************************************************
    !   Adds the local matrix to the global matrix discretization
    !   Continuous Galerkin case
    !
    !  Parameters:
    !
    !    Input,Output, type(SparseMatrix) Mat, Structure for the global matrix
    !    Input, REAL*8, DIMENSION(:,:) Aloc, Local matrix
    !    Input, type(element) ele, element considerated
    !    Input, INTEGER  Nv, Number of nodes (vertices)
    !    Input, INTEGER  Ne, Number of elements
    !
    !*************************************************************************


    IMPLICIT NONE
    type(SparseMatrix)    , INTENT(INOUT) :: Mat
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: Aloc
    type(element)         , INTENT(IN)    :: ele
    INTEGER               , INTENT(IN)    :: Nv, Ne
    !-------------------------------------------------
    INTEGER :: i, j, NU_i, NU_j, Pos_i, Pos_j, k, Pos_k
    INTEGER :: Nei, Nej, NvLoc, Nmax
    INTEGER :: l, m, iG, jG, iL, jL, ivar, jvar
    INTEGER :: nNodesPerElement, nVar, Ns
    !---------------------------------------------------
    LOGICAL :: found
    !------------------

    nMax = mat%nMax ; nNodesPerElement = ele%nNodesPerElement ; nVar = Mat%nVar

    DO i = 1, nNodesPerElement
       Pos_i = ele%SurVertex(i) ! numerotation globale du noeud
       DO j = 1, nNodesPerElement
          Pos_j = ele%SurVertex(j) ! numerotation globale du noeud
          IF ( Pos_i == Pos_j ) THEN ! contribution sur le noeud ayant le numero Pos_i
             DO ivar = 1, nVar
                DO jvar = 1, nVar
                   iG = Pos_i + (ivar - 1)*Nv
                   jG = 1 + (jvar - 1)*NMax !!je comprends pas cette ligne pourquoi 1 et pas Pos_j
                   iL = i + (ivar - 1)*nNodesPerElement
                   jL = j + (jvar - 1)*nNodesPerElement
                   Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                END DO
             END DO
          ELSE
             NvLoc = Mat%Ind(Pos_i,1) ! nombre de voisin du noeud numero Pos_i
             DO k = 1, NvLoc
                Pos_k = Mat%Ind(Pos_i,1+k) ! on prend la numerotation du voisin dans la liste
                IF ( Pos_j == Pos_k ) THEN ! si le voisin choisi dans la liste et le noeud courant de l'élément
                   DO ivar = 1, nVar
                      DO jvar = 1, nVar
                         iG = Pos_i + (ivar-1)*Nv
                         jG = 1 + k + (jvar-1)*nMax ! k est la position du noeud voisin dans le tableau de voisin Mat%Ind
                         iL = i + (ivar-1)*nNodesPerElement
                         jL = j + (jvar-1)*nNodesPerElement
                         Mat%Coeff(iG,jG) = Mat%Coeff(iG,jG) + Aloc(iL,jL)
                      END DO
                   END DO
                   EXIT
                END IF
             END DO
          END IF
       END DO
    END DO
    
  END SUBROUTINE AssembleSparseMatrixCG
  !========================================================!

  !==============================!
  SUBROUTINE RechercheComplementaire(Mesh,NU_j,Comp)
  !==============================!

  !*************************************************************************
  !   Seeks for the number of the complementary node of the node Nu_j
  !
  !  Parameters:
  !
  !    Input,type(MeshType) Mesh, Mesh structure
  !    Input, INTEGER NU_j, Current node
  !    Input,Output, INTEGER Comp, complementary node
  !
  !*************************************************************************


    IMPLICIT NONE

    type(MeshType)    , INTENT(IN)    :: Mesh
    INTEGER, INTENT(IN)               :: NU_j
    INTEGER, INTENT(INOUT)            :: Comp
    !*********************************************
    INTEGER                           :: i
    !*********************************************


    Comp=0
    DO i=1,mesh%NarEmb+1+mesh%PtEmb
       IF(NU_j == mesh%Tab_NoeudComp(i,1)) THEN
          Comp=mesh%Tab_NoeudComp(i,2)
       END IF
    END DO

  END SUBROUTINE RechercheComplementaire
  !==============================!

  !==============================!
  SUBROUTINE FillIndCG(Mat, Mesh)
  !==============================!

  !*************************************************************************
  !   Fills Mat%Ind for the Storage of the discretized matrix
  !   Continuous version
  !
  !  Parameters:
  !
  !    Input,Output, type(MeshType) Mesh, Mesh structure
  !    Input,Output, type(SparseMatrix) Mat, Structure for the matrices of the resolution
  !
  !*************************************************************************


    IMPLICIT NONE
    type(SparseMatrix), INTENT(INOUT)    :: Mat
    type(MeshType)    , INTENT(IN)       :: Mesh
    !-----------------------------------------
    type(element) :: ele, Adj
    !-------------------------
    LOGICAL :: found
    !----------------
    INTEGER :: i, j, NMax, ie, NU_i, NU_j
    INTEGER :: k, l, AdjNum, NU_k, m, NU_m
    INTEGER :: nNodesPerElement, Comp, indice
    INTEGER :: nDDL,mt
    !--------------------------------------

    Mat%Ind = 0 ; NMax = Mat%NMax
    PRINT*, " -- Sparse matrix storage preparation -- "
    DO ie = 1, Mesh%Ne ! boucle sur le nombre d 'element

       ele = Mesh%Ele(ie) !on regarde l'element courant

       IF ( ele%Solve ) THEN

          nNodesPerElement = ele%nNodesPerElement

          DO i = 1, nNodesPerElement
             NU_i = mesh%mappingReal2Sur(ele%Vertex(i)) ! on se place sur un point de l'element courant
             DO j = 1, nNodesPerElement !typiquement ici on a du 3
                NU_j = mesh%mappingReal2Sur(ele%Vertex(j)) !ddl du triangle regarde
                IF ( i /= j ) THEN !si le noeud regardé et different du noeud courant de l'element

                   IF ( Mat%Ind(NU_i,1) == 0 ) THEN !si le nombre de point voisin est pour l'instant nulle
                      Mat%Ind(NU_i,1) = Mat%Ind(NU_i,1) + 1
                      l = Mat%Ind(NU_i,1) + 1 ! variable pour se déplacer sur les colonnes et ajouté la référence des voisins 
                      Mat%Ind(NU_i,l) = NU_j ! on rajoute le voisin comme deuxieme ddl voisin associe au point couran t
                      IF(mesh%NarEmb /= 0) THEN !alors on est dans le cas Icing
                         !on regarde si ce voisin n'a pas un point complémentaire
                         CALL RechercheComplementaire(mesh,NU_j,Comp) !Comp renvoie 0 si pas de complémentaire sinon on renvoie la valeur du complémentaire de NU_j
                         IF (Comp /= 0 .AND. Comp/=Nu_j) THEN !i.e si on se trouve pas sur un noeud de bord d'interface
                            Mat%Ind(NU_i,1) = Mat%Ind(NU_i,1) + 1 ! MAJ de la valeur du nombre de voisin
                            l = Mat%Ind(NU_i,1) + 1
                            Mat%Ind(NU_i,l) = Comp ! ajout de la numerotation du noeud dans le tableau
                         END IF
                      END IF

                   ELSE !cas ou le nombre de voisin n'est pas nulle alors il faut rechercher si un point est deja présent avant de l'ajouter

                      found = .FALSE. !si la valeur n'est pas nulle on va rechercher si le ddl regarder n'est pas deja dans la liste
                      DO k = 1, Mat%Ind(NU_i,1)
                         l = k + 1
                         IF ( Mat%Ind(NU_i, l) == NU_j ) THEN
                            found = .TRUE.
                            EXIT
                         END IF
                      END DO

                      IF ( found .eqv. .FALSE. ) THEN ! le voisin n'est pas encore présent 
                         Mat%Ind(NU_i,1) = Mat%Ind(NU_i,1) + 1 ! +1 sur le nombre de voisin
                         l = Mat%Ind(NU_i,1) + 1
                         Mat%Ind(NU_i,l) = NU_j ! ajout de la numerotation du noeud dans la liste
                         !on regarde si le point ajouter a un complementaire dans le tableau
                         IF(mesh%NarEmb /= 0) THEN
                            CALL RechercheComplementaire(Mesh,NU_j,Comp) !Comp renvoie 0 si pas de complémentaire sinon on renvoie la valeur du complémentaire de NU_j
                            IF (Comp /= 0 .AND. Nu_j /= Comp) THEN
                               Mat%Ind(NU_i,1) = Mat%Ind(NU_i,1) + 1
                               l = Mat%Ind(NU_i,1) + 1
                               Mat%Ind(NU_i,l) = Comp
                            END IF
                         END IF
                      END IF

                   END IF
                END IF
             END DO ! fin sur l'element courant

             !On ajoute aussi des valeurs voisines plus éloigné que les voisins par aretes directe
             !AJOUT   DE VALEURS  voisines
             DO j = 1, nNodesPerElement
                IF ( ele%Adj(j) /= -1 ) THEN ! si on est pas sur le bord

                   AdjNum = ele%Adj(j)/nNodesPerElement
                   Adj = Mesh%Ele(AdjNum)
                   IF ( adj%solve ) THEN

                      DO k = 1, nNodesPerElement
                         NU_k = mesh%mappingReal2Sur(Adj%Vertex(k))
                         IF ( NU_k /= NU_i ) THEN
                            found = .FALSE.
                            DO l = 1, Mat%Ind(NU_i,1)
                               m = l + 1
                               NU_m = Mat%Ind(NU_i,m)
                               IF ( NU_k == NU_m ) THEN
                                  found = .TRUE.
                                  EXIT
                               END IF
                            END DO

                            IF ( found .eqv. .FALSE. ) THEN

                               Mat%Ind(NU_i,1) = Mat%Ind(NU_i,1) + 1
                               m = Mat%Ind(NU_i,1) + 1
                               Mat%Ind(NU_i,m) = NU_k
                               IF( mesh%NarEmb /= 0) THEN
                                  CALL RechercheComplementaire(Mesh,NU_k,Comp) !Comp renvoie 0 si pas de complémentaire sinon on renvoie la valeur du complémentaire de NU_j
                                  IF (Comp /= 0 .AND. Nu_k /= Comp ) THEN
                                     Mat%Ind(NU_i,1) = Mat%Ind(NU_i,1) + 1
                                     m = Mat%Ind(NU_i,1) + 1
                                     Mat%Ind(NU_i,m) = Comp
                                  END IF
                               END IF

                            END IF

                         END IF
                      END DO
                   END IF
                END IF
             END DO
          END DO
       END IF
    END DO ! fin sur la boucle d'element

    !si le noeud a un complementaire on l'ajout aussi au tableau en position initiale debut
    IF (mesh%NarEmb /= 0 .AND. Mesh%NB_VertInterf /= 0) THEN !utile seulement si on considere qu'il y a des neouds compelmentaires
       DO ie=1, Mesh%NarEmb + mesh%PtEmb + 1
		IF(Mesh%Tab_NoeudComp(ie,1) /= Mesh%Tab_NoeudComp(ie,2) ) THEN !on ajoute seuelement si on est pas sur un noeud du bord 
			  Mat%Ind(Mesh%Tab_NoeudComp(ie,1),1)=Mat%Ind(Mesh%Tab_NoeudComp(ie,1),1)+1 !ajout de un dans le nombre des voisins

			  !on va ajouter la valeur au debut don on va decaler toutes les valeurs deja existante du tableau
			  DO mt=1,Mat%Ind(Mesh%Tab_NoeudComp(ie,1),1)
				 Mat%Ind(Mesh%Tab_NoeudComp(ie,1),Mat%Ind(Mesh%Tab_NoeudComp(ie,1),1)+2-mt)= Mat%Ind(Mesh%Tab_NoeudComp(ie,1) &
					  ,Mat%Ind(Mesh%Tab_NoeudComp(ie,1),1)-mt +1)
			  END DO
			  Mat%Ind(Mesh%Tab_NoeudComp(ie,2),:)=Mat%Ind(Mesh%Tab_NoeudComp(ie,1),:) ! copie des deux tableaux
			  Mat%Ind(Mesh%Tab_NoeudComp(ie,1),2)=Mesh%Tab_NoeudComp(ie,2)
			  Mat%Ind(Mesh%Tab_NoeudComp(ie,2),2)=Mesh%Tab_NoeudComp(ie,1)
		END IF 
       END DO

    END IF


    Mat%NZ = 0
    
    ! L: if not icing, there were a +1 that screwed the loop
    IF ( mesh%naremb /= 0 ) THEN
       nDDL = mesh%NvActive + Mesh%NB_VertInterf
    ELSE
       nDDL = mesh%nvActive
    END IF

    DO i = 1, nDDL
       Mat%NZ = Mat%NZ + Mat%Ind(i,1) + 1
    END DO
    Mat%NZ = Mat%nVar*Mat%nVar*Mat%NZ
    PRINT*, " -- Sparse matrix storage done -- "
    PRINT*, " "

  END SUBROUTINE FillIndCG
  !==============================!

  !=================================================!
  SUBROUTINE SparseMatrixVectorProduct_CG(Mat, U, B)
    !=================================================!

    !*************************************************************************
    !   Product of two vectors
    !
    !  Parameters:
    !
    !    Input, type(SparseMatrix) Mat, Structure for the matrices of the resolution
    !    Input,Output, REAL*8, DIMENSION(:) U, vector
    !    Input,Output, REAL*8, DIMENSION(:) B, vector
    !
    !*************************************************************************



    IMPLICIT NONE
    type(SparseMatrix)  , INTENT(IN)    :: Mat
    REAL*8, DIMENSION(:), INTENT(IN)    :: U
    REAL*8, DIMENSION(:), INTENT(INOUT) :: B
    !------------------------------------------
    INTEGER :: N, i, i1, ivar, j, NU_j,nVar, nLoc
    !---------------------------------------------

    nVar = Mat%nVar ; N = SIZE(U)/nVar
    B = 0.d0
    DO ivar = 1, nVar
       DO i = 1, N
          i1 = i + (ivar-1)*N
          B(i1) = B(i1) + Mat%Coeff(i,1)*U(i1)
          nLoc = Mat%Ind(i,1)
          DO j = 1, nLoc
             NU_j = Mat%ind(i,nLoc)
             B(i) = B(i) + Mat%Coeff(i,j)*U(NU_j)
          END DO
       END DO
    END DO

  END SUBROUTINE SparseMatrixVectorProduct_CG
  !=================================================!


  !=================================================!
  SUBROUTINE Affichage_Tab(Data, Sol, Mat, Mesh)
  !=================================================!

  !*************************************************************************
  !   Prints the Mat%Ind Tab
  !
  !  Parameters:
  !
  !    Input, type(DataStructure) Data, Data structure from the dom.data file
  !    Input, type(SolStructure) Sol, Strucuture of the solution
  !    Input, type(SparseMatrix) Mat, Structure of the matrices of the resolution
  !    Input, type(MeshType) Mesh, Mesh structure
  !
  !*************************************************************************

    IMPLICIT NONE
    type(DataStructure), INTENT(IN) :: Data
    type(SolStructure) , INTENT(IN) :: Sol
    type(SparseMatrix) , INTENT(IN) :: Mat
    type(MeshType)     , INTENT(IN) :: Mesh

    !------------------------------------------
    INTEGER :: i,j,Nmax,nNodesPerElement,NeActive,NvActive,nDim,Nv,Ns
    !---------------------------------------------

    Nmax = 10 !changer la valeur selon la taille ici Nmax vaut 70 mais pour ne pas
    !avoir une tonne de 0 a l'affichage d'un exemple tres simple seulement 10 ici

    IF (Data%Icing .eqv. .TRUE.) THEN
       Ns=Mesh%NB_VertInterf
    ELSE
       Ns=0
    END IF

    nDim = Data%nDim;  Nv = Sol%Nv+Ns
    nNodesPerElement = nDim + 1
    NeActive = Sol%NeActive ;NvActive = Sol%NvActive +Ns

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG" )
       DO i=1,NvActive
          WRITE(*,*) "T", i
          DO j=1,NMax
             WRITE(*,*) Mat%Ind(i,j)
          END DO
       END DO
    CASE ( "CG-Primal" )
       DO i=1,NvActive
          WRITE(*,*) "T", i
          DO j=1,NMax
             WRITE(*,*) Mat%Ind(i,j)
          END DO
       END DO
    CASE ( "DG" )
       DO i=1,nNodesPerElement*NeActive
          WRITE(*,*) "T", i
          DO j=1,NMax
             WRITE(*,*) Mat%Ind(i,j)
          END DO
       END DO
    CASE ( "DG-Primal" )
       DO i=1,nNodesPerElement*NeActive
          WRITE(*,*) "T", i
          DO j=1,NMax
             WRITE(*,*) Mat%Ind(i,j)
          END DO
       END DO
    CASE ( "MSEG", "MSEG-Stab" )
       DO i=1,Nv
          WRITE(*,*) "T", i
          DO j=1,NMax
             WRITE(*,*) Mat%Ind(i,j)
          END DO
       END DO
    CASE DEFAULT
       CALL PrintError("Error in type CASE")
    END SELECT

    PRINT*, '  --  End of TABLE  --  '

  END SUBROUTINE Affichage_Tab
  !=================================================!


  !=================================================!
  SUBROUTINE AssembleSparseMatrixInterface(Mat, Amm, Amp, Apm, App, EleM, EleP, Nv, Ne)
    !=================================================!

    !*************************************************************************
    !  Assembles of the different local matrices for the interface terms
    !
    !  Parameters:
    !
    !    Input,Output, type(SparseMatrix) Mat, Structure of the matrices of the resolution
    !    Input,REAL*8, DIMENSION(:,:) Amm, Amp, Apm, App, local matrices at the interface
    !    Input, type(MeshType) EleM, EleP, elements at the interface
    !    Input, INTEGER Nv, vertices number without the complementary nodes
    !    Input, INTEGER  Ne, elements number
    !
    !*************************************************************************


    IMPLICIT NONE
    type(SparseMatrix)    , INTENT(INOUT) :: Mat
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: Amm, Amp, Apm, App
    type(element)         , INTENT(IN)    :: EleM, EleP
    INTEGER               , INTENT(IN)    :: Nv, Ne
    !-------------------------------------------------

    CALL AssembleSparseMatrixDG2(Mat, Amm, EleM, EleM, Nv, Ne)
    CALL AssembleSparseMatrixDG2(Mat, Amp, EleM, EleP, Nv, Ne)
    CALL AssembleSparseMatrixDG2(Mat, Apm, EleP, EleM, Nv, Ne)
    CALL AssembleSparseMatrixDG2(Mat, App, EleP, EleP, Nv, Ne)

  END SUBROUTINE AssembleSparseMatrixInterface
  !=================================================!

END MODULE SparseMatrix_CG
