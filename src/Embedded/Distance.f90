MODULE Distance

  USE Types
  USE PublicVar
  USE CircleDistance
  USE BoxDistance
  USE SphereDistance

  IMPLICIT NONE

CONTAINS

  !==========================================!
  SUBROUTINE computeDistance(PhysicalInterface, Data, Mesh, Sol)
  !==========================================!

    IMPLICIT NONE
    type(VectorInterface), INTENT(IN)  :: PhysicalInterface 
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    !---------------------------------------
    REAL*8, DIMENSION(data%nDim) :: x_center
    !----------------------------------------
    REAL*8 :: r0, x0, y0
    !------------------------
    INTEGER :: i
    !------------

    DO i = 1, Data%N_LS !nombre de condition embedded
       x_center = data%LS_center(i,:)
       SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(i,1))) )
       CASE ( "Circle","circle" )
          r0 = Data%LS_L(i,1) !valeur d'un rayon
          SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(i,2))) )
          CASE ( "In", "in" )
             CALL InnerCircle(Sol, Mesh, x_center(:), r0, i)
          CASE ( "Out", "out" )
             CALL OuterCircle(Sol, Mesh, x_center(:), r0, i) ! a regarder pour le cas test 2
          CASE DEFAULT
             PRINT*, ' ** ERROR in the embedded geometry...: select "in" or "out"'
             STOP
          END SELECT
       CASE ( "box", "Box" )
          SELECT CASE ( Data%nDim)
          CASE ( 2 )
             x0 = Data%LS_L(i,1) ; y0 = Data%LS_L(i,2) !valeurs d'un point de repère
             SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(i,2))) )
             CASE("in","In")
                CALL InnerBox(PhysicalInterface, Sol, Mesh,Data, x_center, x0, y0, i)
             CASE("out","Out")
                CALL OuterBox(Sol, Mesh, x_center, x0, y0, i)
             CASE DEFAULT
                PRINT*, ' ** ERROR in the embedded geometry...: select "in" or "out"'
                STOP
             END SELECT
          CASE ( 3 )
             PRINT* , ' ** ERROR in embedded geometry, box not implemented in 3d'
             STOP
          END SELECT
       CASE ( "Sphere","sphere" )
          SELECT CASE ( data%nDim )
          CASE ( 3 )
             r0 = Data%LS_L(i,1)
             SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(i,2))) )
             CASE ( "In", "in" )
                CALL InnerSphere(Sol, Mesh, r0, i)
             CASE ( "Out", "out" )
                CALL OuterSphere(Sol, Mesh, r0, i)
             CASE DEFAULT
                PRINT*, ' ** ERROR in the embedded geometry...: select "in" or "out"'
                STOP
             END SELECT
          CASE DEFAULT
             PRINT*, ' ** ERROR in computeDistace...: Do not use Sphere and 2D'
             STOP
          END SELECT
       CASE DEFAULT
          PRINT*, ' ** ERROR computeDistance, unknown embedded geometry...'
          PRINT*, '    Available :'
          PRINT*, '      + circle'
          PRINT*, '      + box'
          STOP
       END SELECT
    END DO

  END SUBROUTINE computeDistance
  !==========================================!

  !===================================================!
  SUBROUTINE ComputeNodalDistance(X, d, LS, Data, iLS)
  !===================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: X
    REAL*8, DIMENSION(:), INTENT(OUT) :: d
    REAL*8              , INTENT(OUT) :: LS
    type(DataStructure) , INTENT(IN)  :: Data
    INTEGER             , INTENT(IN)  :: iLS
    !------------------------------------------
    REAL*8, DIMENSION(data%nDim) :: x_center
    !----------------------------------------
    REAL*8 :: r0, x0, y0, xc, yc
    !----------------------------

    x_center = data%LS_center(iLS,:)

    SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,1))) )
    CASE("Circle","circle")
       r0 = Data%LS_L(iLS,1)
       xc = x_center(1) ; yc = x_center(2)
       SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,2))) )
       CASE("in","In")
          CALL NodalInnerCircleDist(X, xc, yc, r0, d, LS)
       CASE("out","Out")
          CALL NodalOuterCircleDist(X, xc, yc, r0, d, LS)
       END SELECT
    CASE("Box","box")
       SELECT CASE ( Data%nDim )
       CASE ( 2 )
          x0 = Data%LS_L(iLS,1) ; y0 = Data%LS_L(iLS,2)
          xc = x_center(1) ; yc = x_center(2)
          SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,2))) )
          CASE("in","In")
			CALL NodalInnerBoxDist(X, xc, yc, x0, y0, d, LS)
          CASE("out","Out")
             CALL NodalOuterBoxDist(X, xc, yc, x0, y0, d, LS)
          END SELECT
       CASE ( 3 )
          PRINT*, "ERROR in computeNodalDistance, Box not implemented in 3D..."
          STOP
       END SELECT
    CASE ( "Sphere","sphere" )
       SELECT CASE ( data%nDim )
       CASE ( 3 )
          r0 = Data%LS_L(iLS,1)
          SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,2))) )
          CASE ( "In", "in" )
             CALL NodalInnerSphere(X, r0, d, LS)
          CASE ( "Out", "out" )
             CALL NodalOuterSphere(X, r0, d, LS)
          CASE DEFAULT
             PRINT*, ' ** ERROR in the embedded geometry...: select "in" or "out"'
             STOP
          END SELECT
       CASE DEFAULT
          PRINT*, ' ** ERROR in computeDistance...: Do not use Sphere and 2D'
          STOP
       END SELECT

    END SELECT

  END SUBROUTINE ComputeNodalDistance
  !======================================================!
  
    !===================================================!
  SUBROUTINE ComputeNodalDistancePhysicalinter(PhysicalInterface, X, d, LS, Data, iLS)
  !===================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:) , INTENT(IN)  :: X
    REAL*8, DIMENSION(:) , INTENT(OUT) :: d
    REAL*8               , INTENT(OUT) :: LS
    type(DataStructure)  , INTENT(IN)  :: Data
    INTEGER              , INTENT(IN)  :: iLS
    type(VectorInterface), INTENT(IN)  :: PhysicalInterface 
    !-------------------------------------------------------
    REAL*8, DIMENSION(data%nDim) :: x_center
    !-------------------------------------------------------
    REAL*8 :: r0, x0, y0, xc, yc
    !-------------------------------------------------------

    x_center = data%LS_center(iLS,:)

    SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,1))) )
    CASE("Circle","circle")
       r0 = Data%LS_L(iLS,1)
       xc = x_center(1) ; yc = x_center(2)
       SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,2))) )
       CASE("in","In")
          CALL NodalInnerCircleDist(X, xc, yc, r0, d, LS)
       CASE("out","Out")
          CALL NodalOuterCircleDist(X, xc, yc, r0, d, LS)
       END SELECT
    CASE("Box","box")
       SELECT CASE ( Data%nDim )
       CASE ( 2 )
          x0 = Data%LS_L(iLS,1) ; y0 = Data%LS_L(iLS,2)
          xc = x_center(1) ; yc = x_center(2)
          SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,2))) )
          CASE("in","In")
			IF (.NOT. Data%DefInterfaceNode) THEN
				CALL NodalInnerBoxDist(X, xc, yc, x0, y0, d, LS)
			ELSE
				CALL OrthoProjectionDistance(PhysicalInterface, X, d, LS)
			END IF 
          CASE("out","Out")
             CALL NodalOuterBoxDist(X, xc, yc, x0, y0, d, LS)
          END SELECT
       CASE ( 3 )
          PRINT*, "ERROR in computeNodalDistance, Box not implemented in 3D..."
          STOP
       END SELECT
    CASE ( "Sphere","sphere" )
       SELECT CASE ( data%nDim )
       CASE ( 3 )
          r0 = Data%LS_L(iLS,1)
          SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(iLS,2))) )
          CASE ( "In", "in" )
             CALL NodalInnerSphere(X, r0, d, LS)
          CASE ( "Out", "out" )
             CALL NodalOuterSphere(X, r0, d, LS)
          CASE DEFAULT
             PRINT*, ' ** ERROR in the embedded geometry...: select "in" or "out"'
             STOP
          END SELECT
       CASE DEFAULT
          PRINT*, ' ** ERROR in computeDistance...: Do not use Sphere and 2D'
          STOP
       END SELECT

    END SELECT

  END SUBROUTINE ComputeNodalDistancePhysicalinter
  !======================================================!
  

  !======================================!
  SUBROUTINE TagElement2Solve(Data, Mesh)
  !======================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(INOUT) :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    !-------------------------------------------
    type(element) :: ele, adj, ele1, ele2
    !-------------------------------------------
    INTEGER :: ie, i, j, NU_i, cntV, cntE, cntP
    INTEGER :: nNodesPerElement, nDim, Pos_i
    !-------------------------------------------
    CHARACTER(len=99) :: Tmp
    !-------------------------------------------

    nDim = data%nDim


    IF (Data%tour == 0 ) THEN
      ALLOCATE ( Data%TreatedNode(Mesh%Nv) )
      ALLOCATE ( Data%TreatedPos((nDim+1)*mesh%Ne) )
      ALLOCATE ( Data%TmpEmbEdg(mesh%Ned) )
    END IF
    Data%TreatedNode = .FALSE.  ;  Mesh%Vertex(:)%Active = .FALSE.
    Mesh%NvActive = 0 ; Mesh%NeActive = 0 ; mesh%NSurB = 0


    DO ie = 1, Mesh%Ne
       ele = mesh%ele(ie) ; nNodesPerElement  = ele%nNodesPerElement
       Mesh%Ele(ie)%Solve = .TRUE.
       Mesh%Ele(ie)%Id = .FALSE.
       ! Defintion of active elements
       DO i = 1, Data%N_LS
          IF ( Mesh%ele(ie)%Status(i) /= -1 ) THEN
             Mesh%Ele(ie)%Solve = .FALSE.
          END IF
       END DO
       ! Definition of active nodes
       DO i = 1, mesh%ele(ie)%nNodesPerElement
          IF ( mesh%ele(ie)%solve ) THEN
             IF ( Data%TreatedNode(mesh%ele(ie)%vertex(i)) .eqv. .FALSE. ) THEN
                Mesh%NvActive = Mesh%NvActive + 1
                Data%TreatedNode(mesh%ele(ie)%vertex(i)) = .TRUE.
             END IF
             IF ( i == 1 ) mesh%NeActive = mesh%NeActive + 1
          END IF
       END DO
    END DO
    IF ( data%immersed .OR. data%icing ) THEN
       mesh%NeActive = mesh%Ne
       mesh%NvActive = mesh%Nv
       mesh%vertex(:)%active = .TRUE.
    END IF
    IF ( data%embedded ) THEN
       ! Definition of the surrogate boundary
       DO i = 1, Mesh%Ned
          IF ( mesh%edg(i)%Tri(1) /= -1 .AND. mesh%edg(i)%Tri(2) /= -1 ) THEN
             ele1 = mesh%Ele(mesh%Edg(i)%Tri(1))
             ele2 = mesh%ele(mesh%edg(i)%Tri(2))
             IF ( ele1%solve .AND. (.NOT. ele2%solve) ) THEN
                DO j = 1, data%N_LS
                   IF ( ele2%status(j) == 0 ) THEN
                      mesh%NSurB = mesh%NSurB + 1
                      Data%tmpEmbEdg(mesh%NSurB) = mesh%edg(i)
                      Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j)
                      Data%tmpEmbEdg(mesh%NSurB)%num = i
                      Data%tmpEmbEdg(mesh%NSurB)%EmbLS = j
                   END IF
                END DO
             ELSE IF ( (.NOT. ele1%solve) .AND. ele2%solve ) THEN
                DO j = 1, data%N_LS
                   IF ( ele1%status(j) == 0 ) THEN
                      mesh%NSurB = mesh%NSurB + 1
                      Data%tmpEmbEdg(mesh%NSurB) = mesh%edg(i)
                      Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j)
                      Data%tmpEmbEdg(mesh%NSurB)%num = i
                      Data%tmpEmbEdg(mesh%NSurB)%EmbLS = j
                   END IF
                END DO
             END IF
          END IF
       END DO

       IF(data%tour == 0 ) THEN !donc si on est a l'étape d'initialisation
         ALLOCATE(mesh%EmbBoundary(2*mesh%NSurB))
       END IF
       mesh%EmbBoundary(1:mesh%NSurB) = Data%TmpEmbEdg(1:mesh%NSurB)

       PRINT*, " "
       PRINT*, " For the embedded simulation :: "
       IF ( .NOT. data%icing ) THEN
          WRITE(Tmp,*) mesh%NeActive
        !  PRINT*, " Found "//TRIM(ADJUSTL(Tmp))//" active elements"
          WRITE(Tmp,*) mesh%NvActive
        !  PRINT*, " Found "//TRIM(ADJUSTL(Tmp))//" active nodes"
       END IF
       WRITE(Tmp,*) mesh%NSurB
       PRINT*, " Found "//TRIM(ADJUSTL(Tmp))//" surrogate boundaries"
       PRINT*, " "

       ALLOCATE ( mesh%mappingSur2RealPos((nDim+1)*mesh%NeActive) )
       ALLOCATE ( Mesh%MappingSur2Real(Mesh%NvActive) )
       ALLOCATE ( Mesh%ActiveElement(Mesh%NeActive) )
       ALLOCATE ( Mesh%MappingReal2Sur(Mesh%Nv) )
       ALLOCATE ( mesh%mappingReal2SurPos(mesh%Ne*(nDim+1)) )

       Data%TreatedNode = .FALSE.  ;  Data%TreatedPos = .FALSE.
       cntE = 0 ; cntV = 0 ; cntP = 0
       ! Define new numbering on surrogate domain
       DO ie = 1, Mesh%Ne

          ele = mesh%ele(ie) ; nNodesPerElement = ele%nNodesPerElement
          IF ( ele%Solve ) THEN
             cntE = cntE + 1
             mesh%activeElement(cntE) = ie
             DO i = 1, nNodesPerElement
                NU_i = ele%Vertex(i)
                Pos_i = ele%Pos(i)
                ! If node is not already seen, increment, get new global number
                ! give the mapping and define active node
                IF ( Data%TreatedNode(NU_i) .eqv. .FALSE. ) THEN
                   Data%TreatedNode(NU_i) = .TRUE.
                   cntV = cntV + 1
                   Mesh%MappingSur2Real(cntV) = NU_i
                   Mesh%MappingReal2Sur(NU_i) = cntV
                   mesh%ele(ie)%surVertex(i) = cntV
                   mesh%vertex(NU_i)%active = .TRUE.
                ELSE
                   mesh%ele(ie)%surVertex(i) = mesh%mappingReal2Sur(NU_i)
                   mesh%vertex(NU_i)%active = .TRUE.
                END IF
                IF ( Data%TreatedPos(Pos_i) .eqv. .FALSE. ) THEN
                   Data%TreatedPos(Pos_i) = .TRUE.
                   cntP = cntP + 1
                   mesh%mappingSur2RealPos(cntP) = Pos_i
                   mesh%mappingReal2SurPos(Pos_i) = cntP
                   mesh%ele(ie)%surPos(i) = cntP
                END IF
             END DO
          END IF
       END DO

    END IF

    IF ( data%icing ) mesh%ele(:)%solve = .TRUE.

  END SUBROUTINE TagElement2Solve
  !======================================!

!  !==============================================!
!  SUBROUTINE TagElement2SolveImmersed(Data, Mesh)
!  !==============================================!

 !   IMPLICIT NONE
  !  type(DataStructure), INTENT(INOUT) :: Data
   ! type(MeshType)     , INTENT(INOUT) :: Mesh
    !------------------------------------------
 !   type(element) :: ele, adj, ele1, ele2
  !  !----------------------------------------------------
   ! INTEGER :: ie, i, j, NU_i, cntV, cntE, cntP, K1, k2, ll, taille_actuel, taille_WthBord
   ! INTEGER :: nNodesPerElement, nDim, Pos_i, id, TypeCond, iBC , prev_position
    !-------------------------------------------
   ! CHARACTER(len=99) :: Tmp
   ! LOGICAL           :: deja_present, IFound , TEST_BORD , deja_present_bord, testvertcs,Test 
   ! INTEGER, DIMENSION(3) :: Test_v, v 
    !-------------------------

   ! nDim = data%nDim 
   ! taille_WthBord = 0 ! nombre de noeud complementaire qui ne sont pas sur le bord dans le tableau
   ! mesh%edg(:)%InterfaceBound = .FALSE. ! savoir si une arete se trouve à l'interface et sur le bord
   ! mesh%Nb_EmbBord = 0 ; Mesh%NB_VertInterf = 0
   ! mesh%vertex(:)%OnSurrogate = .FALSE. !remise à 0 pour recommencer l'identification de la surrogate 
    
   ! IF(Data%tour == 0 ) ALLOCATE ( Data%TmpEmbEdg(mesh%Ned) )
    
   ! Mesh%Vertex(:)%Active = .FALSE.
   ! Mesh%Vertex(:)%OnOutBond = .FALSE. ! pour savoir si un noeud est aussi sur le bord gauche 
   ! mesh%NSurB = 0 ! on recalcul le nombre d'arete a l'interface en embbeded 

	! on commence par modifier le tag des éléments 
   ! DO ie = 1, Mesh%Ne
     !  ele = mesh%ele(ie) ; nNodesPerElement  = ele%nNodesPerElement
    !   Mesh%Ele(ie)%Solve = .TRUE. ; Mesh%Ele(ie)%Id = .FALSE. ; Mesh%Ele(ie)%Cut = .FALSE.        
       ! Defintion of active elements
      ! DO i = 1, Data%N_LS !Data%N_Ls nb de condition de bord embedded
	!	  IF( Mesh%ele(ie)%Status(i) /= -1 .AND.  Mesh%ele(ie)%Status(i) /= 2) THEN 
	!		Mesh%Ele(ie)%Cut = .TRUE. !element coupé par l'interface 
	!	  END IF 
     !     IF ( Mesh%ele(ie)%Status(i) /= -1 ) THEN !possede des points sur la veritable interface
      !       Mesh%Ele(ie)%Solve = .FALSE.
       !   END IF
        ! IF ( mesh%ele(ie)%status(i) == -1 ) THEN !lorsque tous les noeuds de l'elements sont à gauche de l'interface avec possiblement un noeud exactement sur l'interface
         !    mesh%ele(ie)%ref = IcingZone1_ref !changement des ref et du lambda associe bh
          !   mesh%ele(ie)%lambda  = Icing_Lambda1
           !  DO id = 1, nDim
            !    mesh%ele(ie)%L(id,id) = Icing_Lambda1
             !   mesh%ele(ie)%Lm1(id,id)= 1.d0/Icing_Lambda1
!             END DO
 !         ELSE
  !          mesh%ele(ie)%ref = IcingZone2_ref
   !         mesh%ele(ie)%lambda = Icing_Lambda2
    !         DO id = 1, nDim
     !           mesh%ele(ie)%L(id,id) = Icing_Lambda2
      !          mesh%ele(ie)%Lm1(id,id)= 1.d0/Icing_Lambda2
       !      END DO
        !  END IF
       !END DO
    !END DO
    
        
!	IF (Data%Resistance) THEN ! this elements are not cut by the interface 
!	!retag the element 
!	!to get the solid airfoil 
!	    DO ie = 1, Mesh%Ne
!			ele = mesh%ele(ie) ; nNodesPerElement  = ele%nNodesPerElement
!			IF (ele%ref == IcingZone1_ref) THEN 
!			    ! maybe we are in the airfoil solid 
!				CALL CheckSolidAirfoil(ele,Test,Data) 
!				IF(TEST) THEN 
!					mesh%ele(ie)%ref = IcingZone3_ref !changement des ref et du lambda associe bh
!					mesh%ele(ie)%lambda  = Icing_Lambda3
!					DO id = 1, nDim
!						mesh%ele(ie)%L(id,id) = Icing_Lambda3
!						mesh%ele(ie)%Lm1(id,id)= 1.d0/Icing_Lambda3
!					END DO
!				END IF 
!			END IF 
!		END DO 
!	END IF 
 !   !on a finit le retag des elements maintenant il faut refaire le tag des aretes
!
 !   ! Definition of the surrogate boundary
  !  DO i = 1, Mesh%Ned !Mesh%Ned est le nombre d'arretes du maillage
   !     IF ( mesh%edg(i)%Tri(1) /= -1 .AND. mesh%edg(i)%Tri(2) /= -1 ) THEN ! si on a pas une arete de bord externe
!            ele1 = mesh%Ele(mesh%Edg(i)%Tri(1)) ; ele2 = mesh%ele(mesh%Edg(i)%Tri(2))
 !           IF ( ele1%solve .AND. (.NOT. ele2%solve) ) THEN ! si les elements sont dans deux zones différentes on est sur la surrogate
	!			mesh%vertex(mesh%edg(i)%Vertex(1))%OnSurrogate = .TRUE.
	!			mesh%vertex(mesh%edg(i)%Vertex(2))%OnSurrogate = .TRUE.
	!			DO j = 1, data%N_LS
	!				IF ( ele2%status(j) == 0 ) THEN
	!				   mesh%NSurB = mesh%NSurB + 1
	!				   Data%tmpEmbEdg(mesh%NSurB) = mesh%edg(i)
	!				   Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j)
	!				   Data%tmpEmbEdg(mesh%NSurB)%num = i
	!				   Data%tmpEmbEdg(mesh%NSurB)%EmbLS = j
	!				END IF
	!			END DO
	!		ELSE IF ( (.NOT. ele1%solve) .AND. ele2%solve ) THEN
	!			mesh%vertex(mesh%edg(i)%Vertex(1))%OnSurrogate = .TRUE.
	!			mesh%vertex(mesh%edg(i)%Vertex(2))%OnSurrogate = .TRUE.
	!			DO j = 1, data%N_LS
	!				IF ( ele1%status(j) == 0 ) THEN
	!				   mesh%NSurB = mesh%NSurB + 1 !calcul du nombre d'arte sur la nouvelle interface
	!				   Data%tmpEmbEdg(mesh%NSurB) = mesh%edg(i) !enregistrement des aretes de la nouvelles interface, on obtient une nouvelle numerotation pour l'arrete
	!				   Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j) ! type de condition sur cette arete
	!				   Data%tmpEmbEdg(mesh%NSurB)%num = i ! numerotation globale de l'arete de l'interface dans le maillage
	!				   Data%tmpEmbEdg(mesh%NSurB)%EmbLS = j ! numerotation de la condition embedded
	!				END IF
	!			END DO          
	!		END IF
     !     !cas ou on à une patie de l'interface qui fait partie du bord du domaine
      !  ELSE IF ( mesh%edg(i)%Tri(1) == -1 .OR. mesh%edg(i)%Tri(2) == -1 ) THEN !si on est sur une arete de bord 
		!	PRINT*,"Need to check the tag function"
			!STOP	
!			IF (mesh%edg(i)%Tri(1) == -1) THEN 
!				ele1 = mesh%Ele(mesh%Edg(i)%Tri(2)) ! on se place sur le suele element associé à l'arete 
!			ELSE 
!				ele1 = mesh%Ele(mesh%Edg(i)%Tri(1))
!			END IF
!			IF (Mesh%ele(ele1%Num)%Cut) THEN !on regarde si cette element de bord est aussi croisé par l'interface 
!				! on est dans la situation ele%Status = 0 
!				!on va chercher sur qu'elle bord du domaine on se trouve 
!				j = 0 ; IFOUND = .FALSE.
!				DO WHILE(  j < Mesh%Nb .AND. .NOT. IFound) !
!					j=j+1
!					TEst_v(1) = Mesh%Boundary(j)%vertex(1) ; v(1) = mesh%edg(i)%Vertex(1)
!					TEst_v(2) = Mesh%Boundary(j)%vertex(2) ; v(2) = mesh%edg(i)%Vertex(2)
!					CALL OrderTable(TEst_v) ; CALL OrderTable(v)
!					!l'arete peut ne pas etre definis du meme sens au niveaud de l'enregistrement
!					! des noeuds 
!					IF(Test_v(1) == v(1) .AND. Test_v(2) == v(2)) THEN 
!						IFound = .TRUE. 
!						TypeCond = Mesh%Boundary(j)%Ref
!					END IF
!				END DO 
!				IF(TypeCond == 1 ) THEN ! alors on se trouve sur une arete du bord gauche ou l'interface peut coincider 
!
!					mesh%vertex(mesh%edg(i)%Vertex(1))%OnSurrogate = .TRUE.
!					mesh%vertex(mesh%edg(i)%Vertex(2))%OnSurrogate = .TRUE.
!					
!					mesh%edg(i)%InterfaceBound = .TRUE. ! on definit que l'arete de la surroagte est aussi sur une bord du domaine
!					
!					! on definit si un noeud qui est sur l'interface est aussi sur le bord  
!					mesh%vertex(mesh%edg(i)%Vertex(1))%OnOutBond = .TRUE.
!					mesh%vertex(mesh%edg(i)%Vertex(2))%OnOutBond = .TRUE.
!					mesh%Nb_EmbBord = mesh%Nb_EmbBord + 1 ! variable sur le nombre d'arete embedded de bord 
!					
!					DO j = 1, data%N_LS
!						IF ( ele1%status(j) == 0 ) THEN
!						   mesh%NSurB = mesh%NSurB + 1 !nombre d'aretes sur la surrogate 
!						   Data%tmpEmbEdg(mesh%NSurB) = mesh%edg(i)
!						   Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j)
!						   Data%tmpEmbEdg(mesh%NSurB)%num = i
!						   Data%tmpEmbEdg(mesh%NSurB)%EmbLS = j
!						END IF
!					END DO
!				END IF
!			END IF 
!		END IF 
 !   END DO
!
 !   IF(Data%tour == 0) THEN !donc si on est a l'étape d'initialisation
  !    Data%MaxRowEmbedded = 2*(mesh%NsurB+1+mesh%PtEmb)
   !   ALLOCATE(mesh%EmbBoundary(2*mesh%NSurB)) ! permet de faire l'allocation une seule fois
    !  ! on definit le tableau des noeuds complementaires apres la definition de la surrogate
    !  ALLOCATE (mesh%Tab_NoeudComp(Data%MaxRowEmbedded,2)) ! initialisé qu'une seule fois
    !END IF
!    mesh%Tab_NoeudComp = 0
!
 !   taille_actuel=0 ! taille du tableau de complementaire
  !  ! à voir si on ne peut pas juste utiliser le tableau mesh%EmbBoundary
   ! mesh%EmbBoundary(1:mesh%NSurB) = Data%TmpEmbEdg(1:mesh%NSurB) ! tableau d'arete sur la surrogate

 !   mesh%ele(:)%solve = .TRUE. ; mesh%vertex(:)%active = .TRUE.
  !  taille_WthBord = 0 ! nb de noeud du tableau n'étant pas sur le bord seulement ( i.e si deux aretes associé au noeuds alors elles sont toutes les deux sur le bord
	
!    DO i=1,mesh%NsurB
!      j=mesh%EmbBoundary(i)%num ! numero de l'arete !
	  ! on se place sur les deux noeuds de l'arete 
!      k1=mesh%edg(j)%Vertex(1) !numerotation global du noeud
!      k2=mesh%edg(j)%Vertex(2) !numerotation globale du noeud
!      deja_present = .FALSE. ;  deja_present_bord = .FALSE. !initialisation
!      ll=1
!      !on va regarder separemment les deux noeuds de l'arete
!      DO WHILE (ll<=taille_actuel .AND. .NOT. deja_present) ! on commence par le noeud k1 
!      ! on regarde si le noeud est deja dans le tableau de complementaire 
!        IF(k1==mesh%Tab_NoeudComp(ll,1)) THEN
!            deja_present = .TRUE.
!        END IF
    !    IF(k1==mesh%Tab_NoeudComp(ll,1) .AND. mesh%Tab_NoeudComp(ll,1) == mesh%Tab_NoeudComp(ll,2)) THEN
     !       ! ca veut dire que le noeud est present en tant que valeur de BORD 
      !      deja_present_bord = .TRUE. 
       !     prev_position = ll ! position du noeud dans le tableau 
        !END IF
 !           ll=ll+1
  !    END DO
   !   IF(deja_present .eqv. .FALSE. )THEN
    !    ! on l'ajoute au tableau
     !   taille_actuel=taille_actuel+1
      !  mesh%Tab_NoeudComp(taille_actuel,1)= k1
!
 !       IF(mesh%edg(j)%InterfaceBound) THEN ! on verifie comment ajouter le noeud dans le tableau 
	!		testvertcs = .FALSE. 
     !       CALL VertexBInterverif(Mesh, k1, testvertcs) ! ca veut dire que le noeud est rataché à deux aretes du bord 
      !      IF(testvertcs) THEN 
		!		mesh%Tab_NoeudComp(taille_actuel,2) = k1
		!	ELSE 
		!	    taille_WthBord = taille_WthBord+1
		!		mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
		!	END IF
        !ELSE ! si il ne s'agit pas d'uen arete de bord 
         !!   taille_WthBord = taille_WthBord+1
			!mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
	!	END IF 
	!	mesh%Vertex(k1)%Complementaire = mesh%Tab_NoeudComp(taille_actuel,2)
	 ! END IF
      !on fait la meme chose pour le deuxieme noeud de l'arete 
  !    deja_present = .FALSE. ;  deja_present_bord = .FALSE. !initialisation
   !   ll=1
    !  DO WHILE (ll<=taille_actuel .AND. .NOT. deja_present)
     !   IF(k2==mesh%Tab_NoeudComp(ll,1)) THEN
      !      deja_present = .TRUE.
       ! END IF
    !    IF(k2==mesh%Tab_NoeudComp(ll,1) .AND. mesh%Tab_NoeudComp(ll,1) == mesh%Tab_NoeudComp(ll,2)) THEN
     !       ! ca veut dire que le noeud est present en tant que valeur de BORD 
      !      deja_present_bord = .TRUE. 
       !     prev_position = ll 
       ! END IF
       ! ll=ll+1
      !END DO
      !IF(deja_present .eqv. .FALSE. )THEN
       ! on l'ajoute au tableau
       ! taille_actuel=taille_actuel+1
        !mesh%Tab_NoeudComp(taille_actuel,1)= k2
        !IF(mesh%edg(j)%InterfaceBound) THEN
		!	testvertcs = .FALSE. 
         !   CALL VertexBInterverif(Mesh, k2, testvertcs)
          !  IF(testvertcs) THEN 
			!	mesh%Tab_NoeudComp(taille_actuel,2)= k2
			!ELSE
			!	taille_WthBord = taille_WthBord+1
			!	mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
			!END IF  
		!ELSE
		 !  taille_WthBord = taille_WthBord+1
		  ! mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
!		END IF 
!		mesh%Vertex(k2)%Complementaire = mesh%Tab_NoeudComp(taille_actuel,2)
 !     END IF
  !  END DO
    
  !  !avec le tableau de noeud complementaire on definit le nombre de noeud à l'interface possedant un complementaire 
   ! Mesh%NB_VertInterf = 0 
    !DO i = 1, taille_actuel
	!	IF(mesh%Tab_NoeudComp(i,1)/= mesh%Tab_NoeudComp(i,2)) THEN 
	!!		Mesh%NB_VertInterf = Mesh%NB_VertInterf + 1  
		!END IF 
    !END DO 

!	PRINT*,"Affichage du tableau de complementaire "
!	DO I=1,taille_actuel
!		PRINT*,mesh%Tab_NoeudComp(i,1),mesh%Tab_NoeudComp(i,2)
!	END DO

 !   PRINT*, " "
  !  PRINT*, " For the immersed simulation :: "
   ! WRITE(Tmp,*) mesh%NSurB
   ! PRINT*, " Found "//TRIM(ADJUSTL(Tmp))//" surrogate boundaries"
   ! PRINT*, " "

  !END SUBROUTINE TagElement2SolveImmersed
  !==============================================!
  
  
  !==============================================!
  SUBROUTINE TagElement2SolveImmersed(Data, Mesh)
  !==============================================!

    IMPLICIT NONE
    
    type(DataStructure), INTENT(INOUT) :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    !-------------------------------------------------------------------------------------------------
    type(element) 		  :: ele, adj, ele1, ele2
    !-------------------------------------------------------------------------------------------------
    INTEGER       		  :: ie, i, j, NU_i, cntV, cntE, cntP, K1, k2, ll, taille_actuel, taille_WthBord
    INTEGER       		  :: nNodesPerElement, nDim, Pos_i, id, TypeCond, iBC , prev_position
    INTEGER				  :: Nu_1, Nu_2, Nu_3, k
    !-------------------------------------------------------------------------------------------------
    CHARACTER(len=99)     :: Tmp
    LOGICAL               :: deja_present, IFound , TEST_BORD , deja_present_bord, testvertcs,Test 
    LOGICAL				  :: TESTNode
    INTEGER, DIMENSION(3) :: Test_v, v 
    !-------------------------------------------------------------------------------------------------


    nDim									 = Data%nDim  
    taille_WthBord                           = 0         ! nombre de noeud complementaire qui ne sont pas sur le bord dans le tableau
    Mesh%edg(:)%InterfaceBound               = .FALSE.   ! savoir si une arete se trouve à l'interface et sur le bord
    Mesh%Nb_EmbBord = 0 ; Mesh%NB_VertInterf = 0
    Mesh%vertex(:)%OnSurrogate               = .FALSE.   ! remise à 0 pour recommencer l'identification de la surrogate 
    
    IF(Data%tour == 0 ) ALLOCATE ( Data%TmpEmbEdg(mesh%Ned) )
    
    Mesh%Vertex(:)%Active    = .FALSE.
    Mesh%Vertex(:)%OnOutBond = .FALSE. ! pour savoir si un noeud est aussi sur le bord gauche 
    Mesh%NSurB               = 0       ! on recalcul le nombre d'arete a l'interface en embbeded 

	! on commence par modifier le tag des éléments 
    DO ie = 1, Mesh%Ne
    
		ele = Mesh%ele(ie) ; nNodesPerElement  = ele%nNodesPerElement
		Mesh%Ele(ie)%Solve = .TRUE. ; Mesh%Ele(ie)%Id = .FALSE. ; Mesh%Ele(ie)%Cut = .FALSE.      
		  
		! Defintion of active elements
        DO i = 1, Data%N_LS !Data%N_LS nb de condition de bord embedded
			IF(Mesh%ele(ie)%Status(i) /= -1 .AND.  Mesh%ele(ie)%Status(i) /= 2) THEN 
				Mesh%Ele(ie)%Cut = .TRUE.          ! element coupé par l'interface 
			END IF 
			IF (Mesh%ele(ie)%Status(i) /= -1) THEN ! possede des points sur la veritable interface
				Mesh%Ele(ie)%Solve = .FALSE.
			END IF
			IF (Mesh%ele(ie)%status(i) == -1) THEN     ! lorsque tous les noeuds de l'elements sont à gauche de l'interface avec possiblement un noeud exactement sur l'interface
			
				IF (Data%Resistance) THEN ! this elements are not cut by the interface but have been tag IcingZone1_ref 
					!retag the element 
					CALL CheckSolidAirfoil(ele, Test, Data) ! the airfoil area is not define by the mesh 
					IF(TEST .OR. (.NOT. TEST .AND. Data%T==0) .OR. (.NOT. TEST .AND. mesh%ele(ie)%ref == IcingZone3_ref) ) THEN ! an element from area 3 stay 3 the one which is modified is area 2 
						mesh%ele(ie)%ref     = IcingZone3_ref !changement des ref et du lambda associe bh
						mesh%ele(ie)%lambda  = Icing_Lambda3
						DO id = 1, nDim
							mesh%ele(ie)%L(id,id)   = Icing_Lambda3
							mesh%ele(ie)%Lm1(id,id) = 1.d0/Icing_Lambda3
						END DO
					ELSE
						Mesh%ele(ie)%ref     = IcingZone1_ref  ! changement des ref et du lambda associe bh
						Mesh%ele(ie)%lambda  = Icing_Lambda1
						DO id = 1, nDim
							Mesh%ele(ie)%L(id,id)   = Icing_Lambda1
							Mesh%ele(ie)%Lm1(id,id) = 1.d0/Icing_Lambda1
						END DO
					END IF 
					
				ELSE
					Mesh%ele(ie)%ref     = IcingZone1_ref  ! changement des ref et du lambda associe bh
					Mesh%ele(ie)%lambda  = Icing_Lambda1
					DO id = 1, nDim
						Mesh%ele(ie)%L(id,id)   = Icing_Lambda1
						Mesh%ele(ie)%Lm1(id,id) = 1.d0/Icing_Lambda1
					END DO
				END IF 
			ELSE
				Mesh%ele(ie)%ref    = IcingZone2_ref
				Mesh%ele(ie)%lambda = Icing_Lambda2
				DO id = 1, nDim
					Mesh%ele(ie)%L(id,id)   = Icing_Lambda2
					Mesh%ele(ie)%Lm1(id,id) = 1.d0/Icing_Lambda2
				END DO
            END IF
        END DO
    END DO
	
    ! Definition of the surrogate boundary
    DO i = 1, Mesh%Ned !Mesh%Ned est le nombre d'arretes du maillage
        IF ( mesh%edg(i)%Tri(1)/= -1 .AND. mesh%edg(i)%Tri(2)/= -1 ) THEN ! si on a pas une arete de bord externe
            ele1 = mesh%Ele(mesh%Edg(i)%Tri(1)) ; ele2 = mesh%ele(mesh%Edg(i)%Tri(2))
            IF ( ele1%solve .AND. (.NOT. ele2%solve) ) THEN ! si les elements sont dans deux zones différentes on est sur la surrogate
				mesh%vertex(mesh%edg(i)%Vertex(1))%OnSurrogate = .TRUE.
				mesh%vertex(mesh%edg(i)%Vertex(2))%OnSurrogate = .TRUE.
				DO j = 1, data%N_LS
					IF ( ele2%status(j) == 0 ) THEN
					   mesh%NSurB                         = mesh%NSurB + 1
					   Data%tmpEmbEdg(mesh%NSurB)         = mesh%edg(i)
					   Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j)
					   Data%tmpEmbEdg(mesh%NSurB)%num     = i
					   Data%tmpEmbEdg(mesh%NSurB)%EmbLS   = j
					END IF
				END DO
			ELSE IF ( (.NOT. ele1%solve) .AND. ele2%solve ) THEN
				mesh%vertex(mesh%edg(i)%Vertex(1))%OnSurrogate = .TRUE.
				mesh%vertex(mesh%edg(i)%Vertex(2))%OnSurrogate = .TRUE.
				DO j = 1, data%N_LS
					IF ( ele1%status(j) == 0 ) THEN
					   mesh%NSurB                         = mesh%NSurB + 1 !calcul du nombre d'arte sur la nouvelle interface
					   Data%tmpEmbEdg(mesh%NSurB)         = mesh%edg(i) !enregistrement des aretes de la nouvelles interface, on obtient une nouvelle numerotation pour l'arrete
					   Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j) ! type de condition sur cette arete
					   Data%tmpEmbEdg(mesh%NSurB)%num     = i ! numerotation globale de l'arete de l'interface dans le maillage
					   Data%tmpEmbEdg(mesh%NSurB)%EmbLS   = j ! numerotation de la condition embedded
					END IF
				END DO          
			END IF
          !cas ou on à une patie de l'interface qui fait partie du bord du domaine
        ELSE IF ( mesh%edg(i)%Tri(1) == -1 .OR. mesh%edg(i)%Tri(2) == -1 ) THEN !si on est sur une arete de bord 
			!PRINT*,"Need to check the tag function"
			!STOP	
			IF (mesh%edg(i)%Tri(1) == -1) THEN 
				ele1 = mesh%Ele(mesh%Edg(i)%Tri(2)) ! on se place sur le seul element associé à l'arete 
			ELSE 
				ele1 = mesh%Ele(mesh%Edg(i)%Tri(1))
			END IF
			IF (Mesh%ele(ele1%Num)%Cut) THEN !on regarde si cette element de bord est aussi croisé par l'interface 
				! on est dans la situation ele%Status = 0 
				!on va chercher sur qu'elle bord du domaine on se trouve 
				j = 0 ; IFOUND = .FALSE.
				DO WHILE(  j < Mesh%Nb .AND. .NOT. IFound) !
					j=j+1
					TEst_v(1) = Mesh%Boundary(j)%vertex(1) ; v(1) = mesh%edg(i)%Vertex(1)
					TEst_v(2) = Mesh%Boundary(j)%vertex(2) ; v(2) = mesh%edg(i)%Vertex(2)
					CALL OrderTable(TEst_v) ; CALL OrderTable(v)
					!l'arete peut ne pas etre definis du meme sens au niveau de l'enregistrement
					! des noeuds 
					IF(Test_v(1) == v(1) .AND. Test_v(2) == v(2)) THEN 
						IFound = .TRUE. 
						TypeCond = Mesh%Boundary(j)%Ref
					END IF
				END DO 
 
				IF(TypeCond == 1 ) THEN ! alors on se trouve sur une arete du bord gauche ou l'interface peut coincider !!pb here 

					mesh%vertex(mesh%edg(i)%Vertex(1))%OnSurrogate = .TRUE.
					mesh%vertex(mesh%edg(i)%Vertex(2))%OnSurrogate = .TRUE.
					mesh%edg(i)%InterfaceBound                     = .TRUE. ! on definit que l'arete de la surroagte est aussi sur une bord du domaine
					
					! on definit si un noeud qui est sur l'interface est aussi sur le bord  
					mesh%vertex(mesh%edg(i)%Vertex(1))%OnOutBond = .TRUE.
					mesh%vertex(mesh%edg(i)%Vertex(2))%OnOutBond = .TRUE.
					mesh%Nb_EmbBord                              = mesh%Nb_EmbBord + 1 ! variable sur le nombre d'arete embedded de bord 
					
					DO j = 1, data%N_LS
						IF ( ele1%status(j) == 0 ) THEN
						   mesh%NSurB                         = mesh%NSurB + 1 !nombre d'aretes sur la surrogate 
						   Data%tmpEmbEdg(mesh%NSurB)         = mesh%edg(i)
						   Data%tmpEmbEdg(mesh%NSurB)%BC_Type = data%embeddedBCType(j)
						   Data%tmpEmbEdg(mesh%NSurB)%num     = i
						   Data%tmpEmbEdg(mesh%NSurB)%EmbLS   = j
						END IF
					END DO
				END IF
			END IF 
		END IF 
    END DO

    IF(Data%tour == 0) THEN !donc si on est a l'étape d'initialisation
		Data%MaxRowEmbedded = 2*(mesh%NsurB+1+mesh%PtEmb)
		ALLOCATE(mesh%EmbBoundary(2*mesh%NSurB)) ! permet de faire l'allocation une seule fois
		! on definit le tableau des noeuds complementaires apres la definition de la surrogate
		ALLOCATE (mesh%Tab_NoeudComp(Data%MaxRowEmbedded,2)) ! initialisé qu'une seule fois
    END IF
    
    mesh%Tab_NoeudComp = 0
    taille_actuel      = 0 ! taille du tableau de complementaire
    ! à voir si on ne peut pas juste utiliser le tableau mesh%EmbBoundary
    Mesh%EmbBoundary(1:mesh%NSurB) = Data%TmpEmbEdg(1:mesh%NSurB) ! tableau d'arete sur la surrogate

    Mesh%ele(:)%solve     = .TRUE. 
    Mesh%vertex(:)%active = .TRUE.
    taille_WthBord        = 0 ! nb de noeud du tableau n'étant pas sur le bord seulement ( i.e si deux aretes associé au noeuds alors elles sont toutes les deux sur le bord
	
    DO i = 1, mesh%NsurB
	
		j  = mesh%EmbBoundary(i)%num ! numero de l'arete 
		! on se place sur les deux noeuds de l'arete 
		k1 = mesh%edg(j)%Vertex(1) ! numerotation globale du noeud
		k2 = mesh%edg(j)%Vertex(2) ! numerotation globale du noeud
		deja_present = .FALSE. ;  deja_present_bord = .FALSE. !initialisation
		ll = 1
		!	on va regarder separemment les deux noeuds de l'arete
		DO WHILE (ll<=taille_actuel .AND. .NOT. deja_present) ! on commence par le noeud k1 
			! on regarde si le noeud est deja dans le tableau de complementaire 
			IF(k1==mesh%Tab_NoeudComp(ll,1)) THEN
				deja_present = .TRUE.
			END IF
			IF(k1==mesh%Tab_NoeudComp(ll,1) .AND. mesh%Tab_NoeudComp(ll,1) == mesh%Tab_NoeudComp(ll,2)) THEN
				! ca veut dire que le noeud est present en tant que valeur de BORD 
				deja_present_bord = .TRUE. 
				prev_position = ll ! position du noeud dans le tableau 
			END IF
            ll=ll+1
		END DO
		IF(deja_present .eqv. .FALSE. )THEN
			! on l'ajoute au tableau
			taille_actuel=taille_actuel+1
			mesh%Tab_NoeudComp(taille_actuel,1)= k1
			IF(mesh%edg(j)%InterfaceBound) THEN ! on verifie comment ajouter le noeud dans le tableau 
				testvertcs = .FALSE. 
				CALL VertexBInterverif(Mesh, k1, testvertcs) ! ca veut dire que le noeud est rataché à deux aretes du bord 
				IF(testvertcs) THEN 
					mesh%Tab_NoeudComp(taille_actuel,2) = k1
				ELSE 
					taille_WthBord = taille_WthBord+1
					mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
				END IF
			ELSE ! si il ne s'agit pas d'uen arete de bord 
				taille_WthBord = taille_WthBord+1
				mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
			END IF 
			mesh%Vertex(k1)%Complementaire = mesh%Tab_NoeudComp(taille_actuel,2)
		END IF
		!on fait la meme chose pour le deuxieme noeud de l'arete 
		deja_present = .FALSE. ;  deja_present_bord = .FALSE. !initialisation
		ll=1
		DO WHILE (ll<=taille_actuel .AND. .NOT. deja_present)
			IF(k2==mesh%Tab_NoeudComp(ll,1)) THEN
				deja_present = .TRUE.
			END IF
			IF(k2==mesh%Tab_NoeudComp(ll,1) .AND. mesh%Tab_NoeudComp(ll,1) == mesh%Tab_NoeudComp(ll,2)) THEN
				! ca veut dire que le noeud est present en tant que valeur de BORD 
				deja_present_bord = .TRUE. 
				prev_position = ll 
			END IF
			ll=ll+1
		END DO
		IF(deja_present .eqv. .FALSE. )THEN
			! on l'ajoute au tableau
			taille_actuel=taille_actuel+1
			mesh%Tab_NoeudComp(taille_actuel,1)= k2
			IF(mesh%edg(j)%InterfaceBound) THEN
				testvertcs = .FALSE. 
				CALL VertexBInterverif(Mesh, k2, testvertcs)
				IF(testvertcs) THEN 
					mesh%Tab_NoeudComp(taille_actuel,2)= k2
				ELSE
					taille_WthBord = taille_WthBord+1
					mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
				END IF  
			ELSE
				taille_WthBord = taille_WthBord+1
				mesh%Tab_NoeudComp(taille_actuel,2)= taille_WthBord + mesh%Nv
			END IF 
			mesh%Vertex(k2)%Complementaire = mesh%Tab_NoeudComp(taille_actuel,2)
		END IF
    END DO
    
	!avec le tableau de noeud complementaire on definit le nombre de noeud à l'interface possedant un complementaire 
	
	! we are going to check if a node really needs a complemntary value, it means for the problems with three areas 
	DO i =1, taille_actuel 
		Nu_i = mesh%Tab_NoeudComp(i,1)
		TestNode = .FALSE. 
		!Then we check if an element who has the node is an element ref 1 
		Do ie = 1, Mesh%Ne 
			ele = Mesh%ele(ie) 
			Nu_1 = ele%Vertex(1) ; Nu_2 = ele%Vertex(2) ; Nu_3 = ele%Vertex(3)
			IF ( (Nu_i == NU_1 .OR. Nu_i == NU_2 .OR. Nu_i == NU_3) .AND. ele%ref == IcingZone1_ref) THEN 
				TestNode = .TRUE.
			END IF   
		END DO 
		
		IF(.NOT. TESTNode) THEN 
		!it means we do not need a complementary value 
			 mesh%Tab_NoeudComp(i,2) =  mesh%Tab_NoeudComp(i,1)
			 !need to modify the following of the table 
			 DO k = i+1, taille_actuel 
				IF(mesh%Tab_NoeudComp(k,1)/=mesh%Tab_NoeudComp(k,2)) THEN 
					mesh%Tab_NoeudComp(k,2) = mesh%Tab_NoeudComp(k,2) -1 
				END IF 
			 END DO 
		END IF 
	END DO 
	
	!Modification of Mesh%vertex%complemntaire defined above 
	Do i = 1, Mesh%Nv 
		Nu_i = Mesh%vertex(i)%num
		DO k = 1, taille_actuel
			IF(nu_i == mesh%Tab_NoeudComp(k,1)) THEN 
				Mesh%vertex(i)%Complementaire = mesh%Tab_NoeudComp(k,2)
			END IF 
		END DO 
	END DO 
	
	Mesh%NB_VertInterf = 0 
	DO i = 1, taille_actuel
		IF(mesh%Tab_NoeudComp(i,1)/= mesh%Tab_NoeudComp(i,2)) THEN 
			Mesh%NB_VertInterf = Mesh%NB_VertInterf + 1  
		END IF 
	END DO 
	
	!PRINT*,"Affichage du tableau de complementaire "
	!DO I=1,taille_actuel
	!	PRINT*,Mesh%Tab_NoeudComp(i,1), Mesh%Tab_NoeudComp(i,2)
	!END DO

    PRINT*, " "
    PRINT*, " For the immersed simulation :: "
    WRITE(Tmp,*) mesh%NSurB
    PRINT*, " Found "//TRIM(ADJUSTL(Tmp))//" surrogate boundaries"
    PRINT*, " "

  END SUBROUTINE TagElement2SolveImmersed
  !==============================================!
  
  !==============================================!
  SUBROUTINE TagElementConforme(Data, Mesh)
  !==============================================!

    IMPLICIT NONE   
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    !------------------------------------------
    type(element) :: ele, adj, ele1, ele2
    !-------------------------------------------------
    INTEGER :: ie, i, j, NU_i, cntV, cntE, cntP, K1, k2
    INTEGER :: nNodesPerElement, nDim, Pos_i, id, ib
    !-------------------------------------------
    CHARACTER(len=99) :: Tmp
    LOGICAL           :: deja_present
    !-------------------------

    nDim = data%nDim

    DO ie = 1, Mesh%Ne
       IF(mesh%ele(ie)%ref == IcingZone2_ref ) THEN
           mesh%ele(ie)%ref = IcingZone1_ref !changement des ref et du lambda associe bh
           DO id = 1, nDim
              mesh%ele(ie)%L(id,id) = Icing_Lambda1
              mesh%ele(ie)%Lm1(id,id)= 1.d0/Icing_Lambda1
           END DO
        ELSE
          mesh%ele(ie)%ref = IcingZone2_ref
           DO id = 1, nDim
              mesh%ele(ie)%L(id,id) = Icing_Lambda2
              mesh%ele(ie)%Lm1(id,id)= 1.d0/Icing_Lambda2
           END DO
        END IF
    END DO

    DO ib = 1, mesh%nb
      IF (mesh%boundary(ib)%ref == IcingZone1_ref) THEN
        mesh%boundary(ib)%ref = IcingZone2_ref
      ELSE IF (mesh%boundary(ib)%ref == IcingZone2_ref) THEN
        mesh%boundary(ib)%ref = IcingZone1_ref
      END IF
    END DO

    !on a retag les elments maintenant on modifie celle des aretes

    PRINT*, " End : Tag Element for a Conform simulation  :: "

  END SUBROUTINE TagElementConforme
  
  !==============================================!

!==============================================!
  SUBROUTINE VertexBInterverif(Mesh, vertcs, test)
  !==============================================!

    IMPLICIT NONE   
    type(MeshType)     , INTENT(IN)    :: Mesh
    INTEGER 		   , INTENT(IN)    :: vertcs 
    LOGICAL 		   , INTENT(INOUT) :: test 
    !-------------------------------------------------
    INTEGER :: i,j, k1,k2
    !------------------------------------------
    LOGICAL           :: T1 , T2, TROUVER1 , TROUVER2
    !-------------------------

	T1 = .FALSE. ; T2 = .FALSE. ; test = .FALSE. ; TROUVER1 = .FALSE. ; TROUVER2 = .FALSE.
   DO i=1,mesh%NsurB
   j=mesh%EmbBoundary(i)%num
	  ! on se place sur les deux noeuds de l'arete 
      k1=mesh%edg(j)%Vertex(1) !numerotation global du noeud
      k2=mesh%edg(j)%Vertex(2) !numerotation globale du noeud
      
      IF(k1== vertcs) THEN 
      IF(.NOT. TROUVER1) THEN 
		TROUVER1 = .TRUE.
	  ELSE 
		TROUVER2 = .TRUE. 
	  END IF 
		IF( mesh%edg(j)%InterfaceBound) THEN ! si l'arete est une arete de bord 
			IF (.NOT. T1) THEN 
				T1 = .TRUE.
			ELSE IF (T1 .AND. .NOT. T2 ) THEN 
				T2 = .TRUE.
			END IF 
		END IF  
      ELSE IF(K2==vertcs) THEN 
       IF(.NOT. TROUVER1) THEN 
		TROUVER1 = .TRUE.
	  ELSE 
		TROUVER2 = .TRUE. 
	  END IF 
      IF( mesh%edg(j)%InterfaceBound) THEN ! si l'arete est une arete de bord 
			IF (.NOT. T1 ) THEN 
				T1 = .TRUE.
			ELSE IF (T1 .AND. .NOT. T2) THEN 
				T2 = .TRUE.
			END IF 
		END IF  
      
      END IF 
   END DO 

	IF (T1  .AND. T2 ) THEN
		test = .true. 
	else
		TEST = .FALSE. 
	END IF 
	
	IF(TROUVER1 .AND. .NOT. TROUVER2) THEN 
		TEst = .TRUE. 
	END IF 
	
	! si l'arete n'est associé qu'à un triangle on ne defini pas de complementaire 
  END SUBROUTINE VertexBInterverif
  !==============================================!
  
  !==============================================!
  SUBROUTINE CheckSolidAirfoil(ele,Test,Data)
  !==============================================!
  
  IMPLICIT NONE 
  type(Element), INTENT(IN)   		 :: ele 
  LOGICAL      , INTENT(OUT)		 :: Test 
  type(DataStructure), INTENT(IN)    :: Data
  !-------------------------------------------------
  REAL*8 							 :: x , y 
  INTEGER							 :: i, nNodesPerElement
  !------------------------------------------------
  
	nNodesPerElement = ele%nNodesPerElement
	Test = .TRUE. 
	! I need to check that the three coordinates of the element are in the airfoil solid 
	DO i = 1, nNodesPerElement
	  x = ele%Coor(i,1) ; y = ele%Coor(i,2)
	  IF (x>=  Data%Ha .AND. TEST ) THEN 
		TEST = .FALSE. 
	  END IF 
	END DO 
  
  !==============================================!
  END SUBROUTINE CheckSolidAirfoil
  
  
END MODULE Distance
