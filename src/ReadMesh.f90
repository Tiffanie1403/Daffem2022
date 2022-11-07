MODULE ReadMesh_MOD

  USE Types
  USE PublicVar
  USE Tools2D
  USE Tools3D
  USE Hash

  IMPLICIT NONE

CONTAINS


  !==============================!
  SUBROUTINE ReadMesh(Data, mesh)
  !==============================!

  !*************************************************************************
  !   Reads the dom.mesh file, definition of the mesh structure "Mesh"
  !
  !  Parameters:
  !
  !    Input, type(DataStructure) Data, Structure for the date from the dom.data file
  !    Output, type(MeshType) Mesh, Mesh structure from the dom.mesh file
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(IN)  :: Data
    type(MeshType)     , INTENT(OUT) :: Mesh
    !---------------------------------------
    CHARACTER(len=99) :: KW
    !-----------------------
    type(HashTable2D) :: Hash2D
    type(HashTable3D) :: Hash3D
    !-----------------------------
    type(Edge), DIMENSION(:), ALLOCATABLE :: BoundTmp
    !--------------------------------------------------
    LOGICAL, DIMENSION(:), ALLOCATABLE :: boundTmpIsBound
    !---------------------------------------------------------
    INTEGER, DIMENSION(3) :: Tmp
    !----------------------------
    LOGICAL :: found,deja_present
    !-----------------
    REAL*8 :: x1, x2, y1, y2
    !-------------------------
    INTEGER :: V1, V2, k
    INTEGER :: Nv, Ne, Nb, it, NbTmp
    INTEGER :: i, j, nDim,ll, taille_actuel,kk
    !--------------------------------
    taille_actuel=0
    mesh%nDim = data%nDim ; nDim = mesh%nDim

    !on commence par recuperer une des references de bords concernant l'interface
    !il peut y avoir plusieurs conditions sur l'interface
    mesh%NarEmb_ref=0
    mesh%NarEmb = 0
    IF(Data%Icing .AND. .NOT. Data%Embedded ) THEN ! cas conforme
      i=1
      DO WHILE (i<=Data%N_BC .AND. mesh%NarEmb_ref==0 )
        IF( TRIM(ADJUSTL(Data%BC(i)))=="NeumannJump") THEN
          mesh%NarEmb_ref=Data%BC_ref(i)
        END IF
        i=i+1
      END DO
    ELSE IF (Data%Icing .AND. Data%Embedded ) THEN ! cas embedded
      mesh%NarEmb_ref=Data%N_BC + 1
    END IF

    OPEN ( unit = 12, file = TRIM(ADJUSTL(Data%FileName))//'.mesh' )
    PRINT*, '  --  Reading the mesh  -- '
    DO
       READ(12,*) KW
       SELECT CASE ( TRIM(ADJUSTL(KW)) )

          ! Reading of the vertices
       CASE ( "Vertices" )
          READ(12,*) Nv
          mesh%Nv = Nv
          ALLOCATE(mesh%Vertex(Nv))
          DO i = 1, Nv
             ALLOCATE(mesh%Vertex(i)%X(nDim))
             ALLOCATE(mesh%Vertex(i)%Status(Data%N_LS))
             ALLOCATE(mesh%Vertex(i)%Ref(nDim))
             READ(12,*) mesh%Vertex(i)%X
             mesh%Vertex(i)%Boundary = 0
             mesh%Vertex(i)%Ref = 0
             mesh%Vertex(i)%Num = i
             mesh%VerteX(i)%Bound = .FALSE.
             mesh%Vertex(i)%Status = -1
             mesh%Vertex(i)%OnInterface = .FALSE.
             mesh%Vertex(i)%OnSurrogate = .FALSE.
             mesh%Vertex(i)%Complementaire = 0
             mesh%Vertex(i)%LogPtEmb = .FALSE.
             mesh%PtEmb = 0
          END DO
          ! Reading of the triangles
       CASE ( "Triangles" )
          SELECT CASE ( nDim )
          CASE ( 2 )
             READ(12,*) Ne
             mesh%Ne = Ne
             ALLOCATE ( mesh%Ele(Ne) )
             ! allocation/initialization of hash table
             Hash2D%SizeH = 2*mesh%Nv ; Hash2D%SizeD = 3*mesh%Nv
             ALLOCATE ( Hash2D%Item(Hash2D%SizeH + Hash2D%SizeD) )
             ALLOCATE ( Hash2D%EdgNum(Hash2D%SizeH + Hash2D%SizeD) )
             DO i = 1, Hash2D%sizeH + Hash2D%SizeD
                Hash2D%Item(i)%Nxt = 0
                Hash2D%Item(i)%Tri = -1
             END DO
             Hash2D%Compt = 1 ; mesh%Ned = 0
             ALLOCATE ( Hash2D%Edg(3*mesh%Nv) )
             ! Loop over possible number of edges to initialize
             ! both edge in hashtable and element related stuff
             DO i = 1, 3*mesh%Nv
                ALLOCATE ( Hash2D%edg(i)%vertex(2) )
                ALLOCATE ( Hash2D%edg(i)%N(2) )
                Hash2D%Edg(i)%Tri = -1
                IF ( i <= mesh%Ne ) THEN
                   ! Allocation of
                   mesh%Ele(i)%nNodesPerElement = 3
                   mesh%ele(i)%nEdgesPerElement = 3
                   mesh%Ele(i)%eleType = "TRI2D"
                   ALLOCATE ( mesh%Ele(i)%N(3,2) )
                   ALLOCATE ( mesh%Ele(i)%Coor(3,2) )
                   ALLOCATE ( mesh%Ele(i)%L(2,2), mesh%Ele(i)%Lm1(2,2) )
                   ALLOCATE ( mesh%Ele(i)%Vertex(3) )
                   ALLOCATE ( mesh%ele(i)%surVertex(3) )
                   ALLOCATE ( mesh%ele(i)%SurVertex_Sauvg(3) )
                   ALLOCATE ( mesh%ele(i)%surPos(3) )
                   ALLOCATE ( mesh%Ele(i)%Adj(3) )
                   ALLOCATE ( mesh%Ele(i)%Pos(3) )
                   ALLOCATE ( mesh%Ele(i)%Edg(3) )
                   ALLOCATE ( mesh%Ele(i)%Status(Data%N_LS) )
                   ALLOCATE ( mesh%Ele(i)%edgRef(3) )
                   ALLOCATE ( mesh%ele(i)%edgNum(3,2) )
                   mesh%Ele(i)%Edg = -1
                   mesh%Ele(i)%Adj = -1
                END IF
             END DO

             DO i = 1, Ne
                mesh%ele(i)%Num = i
                mesh%ele(i)%solve = .TRUE.
                mesh%ele(i)%Id = .FALSE.
                READ(12,*) mesh%ele(i)%Vertex, mesh%ele(i)%Ref
                IF(Data%Icing) THEN
                  IF(mesh%ele(i)%ref==IcingZone1_ref) THEN
                    mesh%ele(i)%lambda=Icing_Lambda1
                  ELSE IF (mesh%ele(i)%ref==IcingZone2_ref) THEN
                    mesh%ele(i)%lambda=Icing_Lambda2
                  END IF
                ELSE
                  mesh%ele(i)%lambda=0
                END IF
                DO j = 1, 3
                   mesh%ele(i)%Pos(j) = (i-1)*3 + j
                END DO
                DO j = 1, Data%N_LS
                   mesh%ele(i)%status(j) = -1
                END DO
                CALL FillHash2D(mesh, i, Hash2D)
                tmp = mesh%ele(i)%edg
                mesh%ele(i)%edg(1) = Tmp(3)
                mesh%ele(i)%edg(2) = Tmp(2)
                mesh%ele(i)%edg(3) = Tmp(1)
                DO j = 1, 3
                   mesh%vertex(mesh%ele(i)%vertex(j))%ele1 = i
                   mesh%ele(i)%coor(j,:) = mesh%vertex(mesh%ele(i)%vertex(j))%X
                END DO
                CALL Calc2DNormal(mesh%ele(i))
                CALL Calc2DArea(mesh%ele(i))

                DO j = 1, data%N_dom
                   IF ( mesh%ele(i)%ref == data%dom(j) ) THEN
                      mesh%ele(i)%Lm1(1,:) = Data%lambdam1(j,1:2)
                      mesh%ele(i)%Lm1(2,:) = Data%lambdam1(j,3:4)
                      mesh%ele(i)%L(1,:)   = Data%lambda(j,1:2)
                      mesh%ele(i)%L(2,:)   = Data%lambda(j,3:4)
                   END IF
                END DO
              !  mesh%ele(i)%lambda=mesh%ele(i)%L(1,1)
             END DO

             ! Filling edge info from hashtable
             ALLOCATE ( mesh%edg(mesh%Ned) )
             DO i = 1, mesh%ned
                ALLOCATE ( mesh%edg(i)%Vertex(2) )
                ALLOCATE ( mesh%edg(i)%N(2) )
                mesh%edg(i) = hash2D%edg(i)
                V1 = mesh%Edg(i)%Vertex(1) ; V2 = mesh%Edg(i)%Vertex(2)
                x1 = mesh%Vertex(V1)%X(1) ; y1 = mesh%Vertex(V1)%X(2)
                x2 = mesh%Vertex(V2)%X(1) ; y2 = mesh%Vertex(V2)%X(2)
                mesh%edg(i)%N(1) = y2 - y1
                mesh%edg(i)%N(2) = x1 - x2
                x1 = SQRT((y2-y1)**2 + (x1-x2)**2)
                mesh%edg(i)%N = mesh%edg(i)%N/x1
                mesh%edg(i)%treated = .FALSE.
                mesh%edg(i)%eleType = "EDG2D"
                mesh%edg(i)%nNodesPerEdge = 2
             END DO
             
          CASE ( 3 )
             ! Triangles = Boundary elements
             READ(12,*) NbTmp
             ALLOCATE ( boundTmp(NbTmp) )
             ALLOCATE ( boundTmpIsBound(NbTmp) )
             boundTmpIsBound = .FALSE.
             Nb = 0
             !mesh%Nb = Nb
             !ALLOCATE ( mesh%Boundary(mesh%Nb) )
             DO i = 1, NbTmp
                ALLOCATE ( boundTmp(i)%vertex(3) )
                ALLOCATE ( boundTmp(i)%N(3) )
                READ(12,*) boundTmp(i)%Vertex, boundTmp(i)%Ref
                found = .FALSE.
                DO j = 1, Data%N_BC
                   IF ( boundTmp(i)%Ref == data%BC_ref(j) ) THEN
                      boundTmp(i)%BC_Type = Data%BC(j)
                      found = .TRUE.
                      EXIT
                   END IF
                END DO
                IF ( found ) THEN
                   Nb = Nb + 1
                   boundTmpIsBound(i) = .TRUE.
                   mesh%Vertex(boundTmp(i)%Vertex)%Boundary = 1
                   DO j = 1, Ndim
                      IF ( mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(1) == 0 ) THEN
                         mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(1) = &
                              boundTmp(i)%Ref
                      ELSE IF ( mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(2) == 0 ) THEN
                         mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(2) = &
                              boundTmp(i)%Ref
                      ELSE
                         mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(3) = &
                              boundTmp(i)%Ref
                      END IF
                   END DO
                END IF
             END DO

             mesh%Nb = Nb
             ALLOCATE ( mesh%boundary(Nb) )
             it = 0
             DO i = 1, NbTmp
                IF ( boundTmpIsBound(i) ) THEN
                   it = it + 1
                   ALLOCATE ( mesh%boundary(it)%vertex(3) )
                   ALLOCATE ( mesh%boundary(it)%N(3) )
                   mesh%boundary(it) = boundTmp(i)
                END IF
             END DO

             DEALLOCATE ( boundTmpIsBound , boundTmp )

          END SELECT

       CASE ( "Tetrahedra" )

          READ(12,*) Ne
          mesh%Ne = Ne
          ALLOCATE ( mesh%ele(Ne) )
          ! allocation/initialization of hash table
          Hash3D%SizeH = 5*Ne ; Hash3D%SizeD = 5*Ne
          ALLOCATE ( Hash3D%Item(Hash3D%SizeH + Hash3D%SizeD) )
          ALLOCATE ( Hash3D%EdgNum(Hash3D%SizeH + Hash3D%SizeD) )
          DO i = 1, Hash3D%sizeH + Hash3D%SizeD
             Hash3D%Item(i)%Nxt = 0
             Hash3D%Item(i)%Tri = -1
          END DO
          Hash3D%Compt = 1 ; mesh%Ned = 0
          ALLOCATE ( Hash3D%Edg(10*mesh%Ne) )
          ! Loop over possible number of edges to initialize
          ! both edge in hashtable and element related stuff
          DO i = 1, 10*mesh%Ne
             ALLOCATE ( Hash3D%edg(i)%vertex(3) )
             ALLOCATE ( Hash3D%edg(i)%N(3) )
             Hash3D%Edg(i)%Tri = -1
             IF ( i <= mesh%Ne ) THEN
                ! Allocation of
                mesh%Ele(i)%nNodesPerElement = 4
                mesh%ele(i)%nEdgesPerElement = 6
                mesh%Ele(i)%eleType = "TET3D"
                ALLOCATE ( mesh%Ele(i)%N(4,3) )
                ALLOCATE ( mesh%Ele(i)%Coor(4,3) )
                ALLOCATE ( mesh%Ele(i)%L(3,3), mesh%Ele(i)%Lm1(3,3) )
                ALLOCATE ( mesh%Ele(i)%Vertex(4) )
                ALLOCATE ( mesh%ele(i)%surVertex(4) )
                ALLOCATE ( mesh%ele(i)%surPos(4) )
                ALLOCATE ( mesh%Ele(i)%Adj(4) )
                ALLOCATE ( mesh%Ele(i)%Pos(4) )
                ALLOCATE ( mesh%Ele(i)%Edg(4) )
                ALLOCATE ( mesh%Ele(i)%Status(Data%N_LS) )
                ALLOCATE ( mesh%Ele(i)%EdgRef(4) )
                ALLOCATE ( mesh%ele(i)%edgNum(4,3) )
                mesh%Ele(i)%Edg = -1
                mesh%Ele(i)%Adj = -1
                mesh%Ele(i)%EdgRef = 0
             END IF
          END DO

          DO i = 1, Ne
             mesh%ele(i)%Num = i
             mesh%ele(i)%solve = .TRUE.
             mesh%ele(i)%Id = .FALSE.
             READ(12,*) mesh%ele(i)%Vertex, mesh%ele(i)%Ref
             DO j = 1, 3
                mesh%ele(i)%Pos(j) = (i-1)*3 + j
             END DO
             DO j = 1, Data%N_LS
                mesh%ele(i)%status(j) = -1
             END DO
             CALL FillHash3D(mesh, i, Hash3D)
             DO j = 1, 4
                mesh%Vertex(mesh%Ele(i)%Vertex(j))%Ele1 = i
                mesh%Ele(i)%Coor(j,:) = mesh%Vertex(mesh%Ele(i)%Vertex(j))%X
                mesh%Ele(i)%Pos(j) = (i-1)*4 + j
             END DO
             CALL Calc3DNormal(mesh%Ele(i))
             CALL Calc3DArea(mesh%Ele(i))
             DO j = 1, data%n_dom
                IF ( mesh%ele(i)%ref == data%dom(j) ) THEN
                   mesh%ele(i)%Lm1(1,:) = data%lambdam1(j,1:3)
                   mesh%ele(i)%Lm1(2,:) = data%lambdam1(j,4:6)
                   mesh%ele(i)%Lm1(3,:) = data%lambdam1(j,7:9)
                   mesh%ele(i)%L(1,:)   = data%lambda(j,1:3)
                   mesh%ele(i)%L(2,:)   = data%lambda(j,4:6)
                   mesh%ele(i)%L(3,:)   = data%lambda(j,7:9)
                END IF
             END DO
          END DO
          ALLOCATE ( mesh%edg(mesh%Ned) )
          DO i = 1, mesh%Ned
             ALLOCATE ( mesh%edg(i)%vertex(3) )
             ALLOCATE ( mesh%edg(i)%N(3) )
             mesh%edg(i) = hash3D%edg(i)
             CALL CalcNormalFace3D(mesh%edg(i), mesh)
             mesh%edg(i)%treated = .FALSE.
             mesh%edg(i)%eleType = "TRI3D"
             mesh%edg(i)%nNodesPerEdge = 3
          END DO

       CASE ( "Edges" )
          IF ( nDim == 2 ) THEN
             READ(12,*) NbTmp !nombre d'edges
             !mesh%Nb = Nb
             Nb = 0
             ALLOCATE ( BoundTmp(NbTmp) )
             ALLOCATE ( boundTmpIsBound(NbTmp) )
             boundTmpIsBound = .FALSE.
             !ALLOCATE ( mesh%boundary(Nb) )
             DO i = 1, NbTmp
                ALLOCATE ( BoundTmp(i)%Vertex(2) )
                ALLOCATE ( BoundTmp(i)%N(2) )
                READ(12,*) BoundTmp(i)%Vertex, BoundTmp(i)%Ref
                !ALLOCATE ( mesh%boundary(i)%vertex(2) )
                !ALLOCATE ( mesh%boundary(i)%N(2) )
                !READ(12,*) mesh%Boundary(i)%Vertex, mesh%Boundary(i)%Ref
                found = .FALSE.
                !mesh%Vertex(mesh%boundary(i)%vertex(1))%bound = .TRUE.
                !mesh%Vertex(mesh%boundary(i)%vertex(2))%bound = .TRUE.



                DO j = 1, Data%N_BC
                   IF ( boundTmp(i)%Ref == Data%BC_ref(j) ) THEN
                   !IF ( mesh%Boundary(i)%Ref == Data%BC_ref(j) ) THEN
                      boundTmp(i)%BC_Type = Data%BC(j)
                      !mesh%Boundary(i)%BC_Type = Data%BC(j)
                      found = .TRUE.
                   END IF

                   IF ( (found .eqv. .TRUE.) .AND. boundTmp(i)%Ref == mesh%NarEmb_ref ) THEN
                     mesh%NarEmb=mesh%NarEmb +1
                     IF(mesh%Vertex(boundTmp(i)%Vertex(1))%LogPtEmb .eqv. .TRUE. .AND. .NOT. Data%Embedded ) THEN
                       mesh%PtEmb=mesh%PtEmb+1
                     END IF
                     IF(mesh%Vertex(boundTmp(i)%Vertex(2))%LogPtEmb .eqv. .TRUE. .AND. .NOT. Data%Embedded  ) THEN
                       mesh%PtEmb=mesh%PtEmb+1
                     END IF
                       mesh%Vertex(boundTmp(i)%Vertex(1))%LogPtEmb= .TRUE.
                       mesh%Vertex(boundTmp(i)%Vertex(2))%LogPtEmb= .TRUE.
                   END IF

                   IF(found .eqv. .TRUE.) THEN
                    EXIT
                   END IF
                END DO


                IF ( found .eqv. .TRUE. ) THEN
                   Nb = Nb + 1 !variable pour les aretes avec conditions
                   boundTmpIsBound(i) = .TRUE.
                   mesh%Vertex(BoundTmp(i)%Vertex)%Boundary = 1 !definit que l'on a une arete du bord
                   DO j = 1, Ndim
                      IF ( mesh%Vertex(BoundTmp(i)%Vertex(j))%Ref(1) == 0 ) THEN
                         mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(1) = BoundTmp(i)%Ref
                      ELSE
                         mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(2) = BoundTmp(i)%Ref
                      END IF
                   END DO
                ELSE
                END IF


                !meme si on est dans le cas embedded on peut calculer les elements à la veritable interface
                !il suffit d'avoir les edges de l'interface dans le fichier mesh avec comme ref la derniere definie
                IF (Data%Icing .AND. Data%Embedded) THEN
                  found = .FALSE.
                  IF (boundTmp(i)%Ref == Data%N_BC+ 1) THEN !on considere une seule condition a l'interface
                       found = .TRUE.
                  END IF

                  IF ( (found .eqv. .TRUE.) .AND. boundTmp(i)%Ref == mesh%NarEmb_ref ) THEN
                    mesh%NarEmb=mesh%NarEmb +1
                    IF(mesh%Vertex(boundTmp(i)%Vertex(1))%LogPtEmb .eqv. .TRUE. .AND. .NOT. Data%Embedded ) THEN
                      mesh%PtEmb=mesh%PtEmb+1
                    END IF
                    IF(mesh%Vertex(boundTmp(i)%Vertex(2))%LogPtEmb .eqv. .TRUE. .AND. .NOT. Data%Embedded  ) THEN
                      mesh%PtEmb=mesh%PtEmb+1
                    END IF
                      mesh%Vertex(boundTmp(i)%Vertex(1))%LogPtEmb= .TRUE.
                      mesh%Vertex(boundTmp(i)%Vertex(2))%LogPtEmb= .TRUE.
                  END IF

                  IF ( found .eqv. .TRUE. ) THEN
                     Nb = Nb + 1 !variable pour les aretes avec conditions
                     boundTmpIsBound(i) = .TRUE.
                     mesh%Vertex(BoundTmp(i)%Vertex)%Boundary = 1 !definit que l'on a une arete du bord
                     DO j = 1, Ndim
                        IF ( mesh%Vertex(BoundTmp(i)%Vertex(j))%Ref(1) == 0 ) THEN
                           mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(1) = BoundTmp(i)%Ref
                        ELSE
                           mesh%Vertex(boundTmp(i)%Vertex(j))%Ref(2) = BoundTmp(i)%Ref
                        END IF
                     END DO
                  ELSE
                  END IF

                END IF
             END DO !fin de la boucle sur le nombre d'arete


             IF(Data%Icing .AND. .NOT. Data%Embedded) THEN
               ALLOCATE ( mesh%Tab_NoeudComp(mesh%NarEmb+1+mesh%PtEmb,2))
             END IF

             !permet en fonction de la géométrie de savoir comment on lie le nombre d'arete
             ! à l'interface avec le nombre de noeuds à l'interface
             IF(mesh%PtEmb==mesh%NarEmb) THEN
              mesh%PtEmb=-1
             ELSE
               mesh%PtEmb=0
             END IF

             mesh%Nb = Nb
             ALLOCATE ( mesh%boundary(Nb) )
             it = 0

             DO i = 1, NbTmp
                IF ( boundTmpIsBound(i) ) THEN
                   it = it + 1
                   ALLOCATE ( mesh%boundary(it)%vertex(2) )
                   ALLOCATE ( mesh%boundary(it)%N(2) )
                   mesh%boundary(it) = boundTmp(i)

                   IF (Data%Icing .AND. .NOT. Data%Embedded) THEN
                      !si on se trouve sur une arete du bord
                      IF(mesh%boundary(it)%Ref == mesh%NarEmb_ref) THEN
                         deja_present = .FALSE. !initialisation
                         ll=1
                         !on va regarder separemment les deux noeuds de l'arete
                         DO WHILE (ll<=taille_actuel)
                            IF(mesh%boundary(it)%vertex(1)==mesh%Tab_NoeudComp(ll,1)) THEN
                               deja_present = .TRUE.
                            END IF
                            ll=ll+1
                         END DO
                         IF(deja_present .eqv. .FALSE. )THEN
                            ! on l'ajoute au tableau
                            taille_actuel=taille_actuel+1
                            mesh%Tab_NoeudComp(taille_actuel,1)= mesh%boundary(it)%vertex(1)
                            mesh%Tab_NoeudComp(taille_actuel,2)= taille_actuel + mesh%Nv
                            mesh%Vertex(mesh%boundary(it)%vertex(1))%OnInterface = .TRUE.
                            mesh%Vertex(mesh%boundary(it)%vertex(1))%Complementaire = mesh%Tab_NoeudComp(taille_actuel,2)
                         END IF
                         deja_present = .FALSE.
                         ll=1
                         DO WHILE (ll<=taille_actuel)
                            IF(mesh%boundary(it)%vertex(2)==mesh%Tab_NoeudComp(ll,1)) THEN
                               deja_present = .TRUE.
                            END IF
                            ll=ll+1
                         END DO
                         IF(deja_present .eqv. .FALSE. )THEN
                            ! on l'ajoute au tableau
                            taille_actuel=taille_actuel+1
                            mesh%Tab_NoeudComp(taille_actuel,1)= mesh%boundary(it)%vertex(2)
                            mesh%Tab_NoeudComp(taille_actuel,2)= taille_actuel + mesh%Nv
                            mesh%Vertex(mesh%boundary(it)%vertex(2))%OnInterface = .TRUE.
                            mesh%Vertex(mesh%boundary(it)%vertex(2))%Complementaire = mesh%Tab_NoeudComp(taille_actuel,2)
                         END IF
                    END IF

                END IF !FIN CAS ICING

              END IF
          END DO
          DEALLOCATE ( boundTmpIsBound , boundTmp )
       END IF ! end case dim 2

       CASE ( "End" )

          EXIT

       END SELECT
    END DO

    PRINT*, "      + nNodes      = ", mesh%Nv
    PRINT*, "      + nElements   = ", mesh%Ne
    PRINT*, "      + nEdges      = ", mesh%Ned
    PRINT*, "      + nBoundaries = ", mesh%Nb
    PRINT*, "      + nBoundaries true interface = ", mesh%NarEmb
    PRINT*, "      + mtpl nodes true interface = ", mesh%PtEmb

    Mesh%NarEmb_prec= Mesh%NarEmb
	!tableau des noeuds complémentaire dans le cas conforme 
  !  PRINT*, "  --  mesh read --  "
   ! PRINT*, "  --  Affichage du tableau de noeuds complémentaires  --  "
   ! i=1
   ! DO WHILE (i<=taille_actuel)
    !  WRITE(*,*) 'T =', mesh%Tab_NoeudComp(i,1), mesh%Tab_NoeudComp(i,2)
    !  i=i+1
    !END DO

    CLOSE ( 12 )

  END SUBROUTINE ReadMesh
  !==============================!

!definition of the boundary conditions
  !=====================!
  SUBROUTINE SetBC(Mesh)
  !=====================!

  !*************************************************************************
  !  Definition of the boundary conditions in the mesh structure
  !
  !  Parameters:
  !
  !    Input,Output, type(meshType) mesh, Mesh Structure
  !
  !*************************************************************************


    IMPLICIT NONE
    type(meshType), INTENT(INOUT) :: Mesh
    !----------------------------------
    type(element) :: ele
    type(edge)    :: Edg, EdgBC
    !--------------------------
    INTEGER, DIMENSION(:), ALLOCATABLE :: N, Nb
    !--------------------------------------------
    INTEGER :: ibc, Nbc, ied, Ned, Nele, i, ie, Ne
    INTEGER :: Nmin, Nmid, Nmax, ifound, it, j
    INTEGER :: nNodesPerEdge, nNodesPerElement
    !------------------------------------------------
    LOGICAL :: found, isEqual
    !--------------------------


    Nbc = mesh%Nb ; Ned = mesh%Ned ; Ne = mesh%Ne
    ALLOCATE ( N(mesh%nDim), Nb(mesh%nDim) )


    it = 0
    DO ied = 1, Ned ! nombre d'aretes
       Edg = mesh%Edg(ied) ! on se place sur une arete du maillage
       nNodesPerEdge = Edg%nNodesPerEdge
       IF ( Edg%Tri(2) == -1 ) THEN ! si on se trouve sur une arete de bord
          DO i = 1, nNodesPerEdge
             !it = it + 1
             N(i) = Edg%Vertex(i) ! on enregistre les deux noeuds associé à l'arete
          END DO
          CALL OrderTable(N) ! on a un vecteur contenant la reference des deux noeuds de l'aretes

          found = .FALSE.
          DO ibc = 1, Nbc ! nombre de conditions embedded
             EdgBC = mesh%Boundary(iBC)
             DO i = 1, nNodesPerEdge
                Nb(i) = edgBC%vertex(i)
             END DO
             CALL OrderTable(Nb)
             isEqual = .TRUE.
             DO i = 1, nNodesPerEdge
                IF ( Nb(i) /= N(i) ) THEN
                   isEqual = .FALSE.
                   EXIT
                END IF
             END DO
             IF ( isEqual ) THEN
                found = .TRUE.
                ifound = ibc
                EXIT
             END IF
          END DO
          IF ( found .eqv. .TRUE. ) THEN ! si l'arete definit une  condition de bord
             Nele = edg%tri(1) ! car associe a un seul elemen etant au bord
             ele = mesh%ele(Nele) ! on se place sur l'element associé 
             nNodesPerElement = ele%nNodesPerElement
             DO i = 1, nNodesPerElement
                edgBC = mesh%edg(mesh%ele(Nele)%edg(i))
                nNodesPerEdge = edgBC%nNodesPerEdge
                DO j = 1, nNodesPerEdge
                   Nb(j) = edgBC%vertex(j)
                END DO
                CALL OrderTable(Nb)
                isEqual = .TRUE.
                DO j = 1, nNodesPerEdge
                   IF ( N(j) /= Nb(j) ) THEN
                      isEqual = .FALSE.
                      !EXIT
                   END IF
                END DO
                IF ( isEqual ) THEN
                   mesh%Ele(Nele)%EdgRef(i) = mesh%Boundary(ifound)%Ref
                   mesh%ele(nele)%edgNum(i,:) = mesh%Boundary(ifound)%vertex
                   mesh%boundary(ifound)%tri = edgBC%tri
                   it = it+1
                   EXIT
                END IF
             END DO
             IF ( .NOT. isEqual ) THEN
                PRINT*, "ERROR in SetBC, ref not found..."
                STOP
             END IF
          ELSE
             PRINT*, N
             PRINT*, " ERROR in SetBC, no match for BC..."
             STOP
          END IF
       ELSE
          ! For internal edges such as NeumannJump
          DO i = 1, nNodesPerEdge
             N(i) = Edg%Vertex(i)
          END DO
          CALL OrderTable(N)

          found = .FALSE.
          DO ibc = 1, Nbc
             EdgBC = mesh%Boundary(iBC)
             DO i = 1, nNodesPerEdge
                Nb(i) = edgBC%vertex(i)
             END DO
             CALL OrderTable(Nb)
             isEqual = .TRUE.
             DO i = 1, mesh%nDim
                IF ( Nb(i) /= N(i) ) THEN
                   isEqual = .FALSE.
                   EXIT
                END IF
             END DO
             IF ( isEqual ) THEN
                ifound = iBC
                EXIT
             END IF
          END DO
          IF ( isEqual .eqv. .TRUE. ) THEN
             mesh%edg(ied)%BC_Type = mesh%boundary(ifound)%BC_Type
             mesh%edg(ied)%Ref = mesh%boundary(ifound)%Ref
             Mesh%boundary(iFound) = mesh%edg(ied)
             it = it+1
             !EXIT
          END IF
       END IF
    END DO

    DEALLOCATE ( N, Nb )

  END SUBROUTINE SetBC
  !========================!

  !==================================!
  SUBROUTINE GetFace(ele, edg, iFace)
  !==================================!

  !*************************************************************************
  !  Gives the corresponding number of the edg in the current ele
  !
  !  Parameters:
  !
  !    Input type(Element) ele, Current element
  !    Input, type(Edge) edg, Current edg
  !    Output, INTEGER  iFace, edg number in element ele
  !
  !*************************************************************************




    IMPLICIT NONE
    type(Element), INTENT(IN)  :: ele
    type(Edge)   , INTENT(IN)  :: edg
    INTEGER      , INTENT(OUT) :: iFace
    !-----------------------------------
    INTEGER, DIMENSION(4,3) :: PermutEdg
    !------------------------------------
    INTEGER :: i, nDim, id
    INTEGER, DIMENSION(SIZE(edg%vertex)) :: Nele, Nedg
    !---------------------------------------------------
    LOGICAL :: found
    !----------------

    nDim = SIZE(edg%vertex)
    PermutEdg(1,1) = 2 ; PermutEdg(1,2) = 3 ; PermutEdg(1,3) = 4
    PermutEdg(2,1) = 3 ; PermutEdg(2,2) = 1 ; PermutEdg(2,3) = 4;
    PermutEdg(3,1) = 1 ; PermutEdg(3,2) = 2 ; PermutEdg(3,3) = 4;
    PermutEdg(4,1) = 1 ; PermutEdg(4,2) = 2 ; PermutEdg(4,3) = 3;

    DO i = 1, ele%nNodesPerElement
       found = .TRUE.
       DO id = 1, nDim
          Nele(id) = ele%vertex(permutEdg(i,id))
          Nedg(id) = edg%vertex(id)
       END DO
       CALL OrderTable(Nele)
       CALL OrderTable(Nedg)

       found = .TRUE.
       DO id = 1, nDim
          IF ( Nele(id) /= Nedg(id) ) THEN
             found = .FALSE.
             EXIT
          END IF
       END DO

       IF ( found ) THEN
          iFace = i
          EXIT
       END IF
    END DO

  END SUBROUTINE GetFace
  !==================================!  

  !==================================!
  SUBROUTINE VertcsPos(Mesh, vertcs, Test, typeV)
  !==================================!

  !*************************************************************************
  !  returns if vertcs is a point on the boundary 
  ! if its the case the table typeV co,ntains the ref of the type conditon ( the vertex is contained by two edges)
  !
  !  Parameters:
  !
  !    Input type(MeshType) Mesh , Mesh structure 
  !    Input, INTEGER vertcs, numerotation of the current vertcs 
  !    Output, LOGICAL test, says if vertcs is on the boundary or not 
  !
  !*************************************************************************


    IMPLICIT NONE
    type(MeshType)    		, INTENT(IN)   		:: Mesh
    INTEGER            		, INTENT(IN) 		:: vertcs
	INTEGER, DIMENSION(2,2)	, INTENT(INOUT)   	:: typeV ! on va regarder le type des deux aretes associé 
	LOGICAL					, INTENT(inOUT) 	    :: Test 
    !---------------------------------------------------
    INTEGER :: i,j, k1, k2 
	LOGICAL :: TEST1, TEST2, IFOUND ! si il s'agit d'un noeud de bord il est au plus affilié à deux aretes 
    !---------------------------------------------------

	!la numerotation 5 coincide avec la numerotation de la surrogate du maillage de départ 
	
	Test = .FALSE. ; typeV = 0 ! there is no reference 0 
	TEST1 = .FALSE. ; TEST2 = .FALSE. 
	Do i=1, Mesh%Nb ! nombre d'arete sur le bord , comprend aussi les aretes sur la surrogate 
		k1 = Mesh%Boundary(i)%Vertex(1) ; k2 = Mesh%Boundary(i)%Vertex(2)
		IF(k1==vertcs .OR. k2==vertcs) THEN 
			IF (mesh%vertex(k1)%OnSurrogate .AND. mesh%vertex(k2)%OnSurrogate) THEN !cela veut dire qu'on est pas sur une arete de la surrogate 
				!on ne fait rien car on est sur une arte de la surrogate 	
			ELSE 
		
				IF(.NOT. test1) THEN 
					Test1 = .TRUE. ! le noeud est bien sur le bord 
					! on doit identifier sur qu'elle bord on se trouve 
					typeV(1,1) = Mesh%Boundary(i)%Ref
					typeV(1,2) = i !numero de l'arete de bord dans Mesh%nb 
				ELSE IF (.NOT. test2) THEN 
					TEST2= .TRUE. 
					typeV(2,1) = Mesh%Boundary(i)%Ref
					typeV(2,2) = i
				ELSE IF((typeV(2,1)==5 .AND. typeV(1,1)/=5) .OR. (typeV(2,1)/=5 .AND. typeV(1,1)==5)) THEN !il faut faire attention car les aretes de la surrogate d'origine sont toujours préesentes 
				
					IF( Mesh%Boundary(i)%Ref/= 5 ) THEN !on a bien un noeud de bord on remplace le type 5 par la valeur du nouveau noeud 
					
						IF(typeV(1,1)==5) THEN 
							typeV(1,1) = Mesh%Boundary(i)%Ref
							typeV(1,2) = i !numero de l'arete de bord dans Mesh%nb 
						ELSE 
							typeV(2,1) = Mesh%Boundary(i)%Ref
							typeV(2,2) = i
						END IF 
					ELSE 
						PRINT*,"Error in VertcsPos pos 1"
					END IF 
				
				ELSE IF (typeV(1,1)/=5 .AND. typeV(2,1)/=5 .AND. Mesh%Boundary(i)%Ref == 5) THEN 
				
					! on ne change rien 
				ELSE 
				
					PRINT*,"Error in VertcsPos pos 2"

				END IF 
			END IF 
		END IF 
	END DO
	
	IF( TEST1 .AND. Test2 .AND. (typeV(1,1)==5 .AND. typeV(2,1)==5) ) THEN 
		test = .FALSE. ! on a un ancien noeud de la surrogate du maillage
		
	ELSE IF(TEST1 .AND. Test2) THEN ! on a bien un noeud de bord
		test=.TRUE. 
	END IF 

  END SUBROUTINE  VertcsPos
  !==================================!

  


END MODULE ReadMesh_MOD
