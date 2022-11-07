MODULE WriteSol

  USE Types
  USE ExactSolutiont
  USE ExactFunctionJumpt

  IMPLICIT NONE

CONTAINS

  !==============================================!
  SUBROUTINE WriteSolution(Mesh, Sol, Data, char)
  !==============================================!

  !*************************************************************************
  !   Writes the solution on a *.tec file
  !
  !  Parameters:
  !
  !    Input, type(MeshType)  Mesh, mesh structure from the dom.mesh file
  !    Input, type(SolStructure) Sol, Solution structure
  !    Input, type(DataStructure) Data, data structure from the dom.data file
  !    Input, CHARACTER(len=*)   char, name of the structure used to identified the files (here "dom")
  !
  !*************************************************************************



    IMPLICIT NONE
    type(MeshType)     , INTENT(IN) :: Mesh
    type(SolStructure) , INTENT(IN) :: Sol
    type(DataStructure), INTENT(IN) :: Data
    CHARACTER(len=*)   , INTENT(IN) :: char
    !---------------------------------------
    type(element) :: ele
    !--------------------
    CHARACTER(len=99) :: Ne_Char, Nv_Char, ind
    CHARACTER(len=99) :: SpSc, WriteDx, WriteDy, WriteDz, WriteSt, WriteLS
    !------------------------------------------------------------------------
    REAL*8, DIMENSION(Data%nDim+1) :: Exact
    !---------------------------------------
    REAL*8 :: NrmB, NrmBe, conv 
    !---------------------
    INTEGER :: Nv, Ne, i, NU_i, j, NU_j, Pos_i, Pos_j
    INTEGER :: nNodesPerElement, nDim, id, ie, surPos_j
    !-------------------------------------------------
    LOGICAL :: eleSolve
    !--------------------
    
    PRINT*, " "
    PRINT*, " -- Writing the solution"

    SpSc = Data%Space_Scheme

    nDim = Data%nDim

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG", "CG-Primal" )
       IF ( .NOT. data%Icing ) THEN
          OPEN(unit = 12, file = TRIM(ADJUSTL(Data%FileName))//'.'&
               //TRIM(ADJUSTL(SpSc))//'.'//char//'.tec')

          WRITE(12,*) 'TITLE = DarcyFlow'

          IF ( .NOT. Data%Embedded .AND. .Not. Data%Immersed) THEN

             WRITE(Nv_char,*) Mesh%Nv
             WRITE(Ne_char,*) Mesh%Ne

             SELECT CASE ( nDim )
             CASE ( 2 )
                WRITE(12,*) 'VARIABLES = "x","y","Coor_x","Coor_y","NodeID","Bx","By","NrmB","p"&
                     ,"Bxe","Bye","NrmBe","pe","eBx","eBy","eNrmB","eP"'
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TRIANGLE'
             CASE ( 3 )
                WRITE(12,*) 'VARIABLES = "x","y","z","Coor_x","Coor_y","Coor_z","NodeID","Bx","By","Bz","NrmB","p","Bxe","Bye",&
                     "Bze","NrmBe","pe","eBx","eBy","eBz","eNrmB","eP"'
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TETRAHEDRON'
             END SELECT
             DO i = 1, Mesh%Nv
                IF ( .NOT. data%icing ) THEN
                   CALL ExactSolt(mesh%ele(1)%L, Exact,Mesh%Vertex(i)%X,Data%T, Mesh%Ele(Mesh%Vertex(i)%Ele1), nDim)
                ELSE
                   CALL ExactSolJumpt(mesh%ele(1)%L, Exact,Mesh%Vertex(i)%X,Data%T, Mesh%Ele(Mesh%Vertex(i)%Ele1), nDim)
                END IF
                NrmB = 0.d0 ; NrmBe = 0.d0
                DO id = 1, nDim
                   NrmB = NrmB + Sol%Bcg(i,id)**2
                   NrmBe = NrmBe + Exact(id)**2
                END DO
                NrmB = sqrt(NrmB)
                NrmBe = sqrt(NrmBe)
                WRITE(12,*) Mesh%vertex(i)%X, Mesh%vertex(i)%X, i,&
                     Sol%Bcg(i,:), NrmB, Sol%pcg(i), &
                     Exact(1:nDim), NrmBe, Exact(nDim+1), &
                     abs(Sol%Bcg(i,:)-Exact(1:nDim)), abs(NrmB-NrmBe), abs(Sol%pcg(i)-Exact(nDim+1))
             END DO
             DO i = 1, Mesh%Ne
                WRITE(12,*) Mesh%Ele(i)%Vertex
             END DO

          ELSE ! Embedded
             WRITE(Nv_char,*) (Data%nDim+1)*Mesh%Ne
             WRITE(Ne_char,*) Mesh%Ne
             SELECT CASE ( Data%nDim )
             CASE ( 2 )
                WriteDx = ',"dist_x1"'
                WriteDy = ',"dist_y1"'
                WriteSt = ',"EleStatus1"'
                WriteLS = ',"LS1"'
                IF ( Data%N_LS > 1 ) THEN
                   DO i = 2, Data%N_LS
                      WRITE(ind,*) i
                      WriteDx = TRIM(ADJUSTL(WriteDx))//'"dist_x'//TRIM(ADJUSTL(ind))//'"'
                      WriteDy = TRIM(ADJUSTL(WriteDy))//'"dist_y'//TRIM(ADJUSTL(ind))//'"'
                      WriteSt = TRIM(ADJUSTL(WriteSt))//',"EleStatus'//TRIM(ADJUSTL(ind))//'"'
                      WriteLS = TRIM(ADJUSTL(WriteLS))//',"LS'//TRIM(ADJUSTL(ind))//'"'
                   END DO
                END IF
                WRITE(12,*) 'VARIABLES = "x", "y","Coor_x","Coor_y","NodeID","Bx","By","NrmB","p"&
                     ,"Bxe","Bye","NrmBe","pe","eBx","eBy","eNrmB","eP"'//&
                     TRIM(ADJUSTL(WriteDx))//TRIM(ADJUSTL(WriteDy))//&
                     TRIM(ADJUSTL(WriteLS))//TRIM(ADJUSTL(WriteSt))
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TRIANGLE'
             CASE ( 3 )
                WriteDx = ',"dist_x1"'
                WriteDy = ',"dist_y1"'
                WriteDz = ',"dist_z1"'
                WriteSt = ',"EleStatus1"'
                WriteLS = ',"LS1"'
                IF ( Data%N_LS > 1 ) THEN
                   DO i = 2, Data%N_LS
                      WRITE(ind,*) i
                      WriteDx = TRIM(ADJUSTL(WriteDx))//',"dist_x'//TRIM(ADJUSTL(ind))//'"'
                      WriteDy = TRIM(ADJUSTL(WriteDy))//',"dist_y'//TRIM(ADJUSTL(ind))//'"'
                      WriteDz = TRIM(ADJUSTL(WriteDz))//',"dist_z'//TRIM(ADJUSTL(ind))//'"'
                      WriteSt = TRIM(ADJUSTL(WriteSt))//',"EleStatus'//TRIM(ADJUSTL(ind))//'"'
                      WriteLS = TRIM(ADJUSTL(WriteLS))//',"LS'//TRIM(ADJUSTL(ind))//'"'
                   END DO
                END IF
                WRITE(12,*) 'VARIABLES = "x", "y","z","Coor_x","Coor_y","Coor_z", "NodeID","Bx","By","Bz",&
                     "NrmB","p","Bxe","Bye","Bze","NrmBe","pe","eBx","eBy","eBz","eNrmB","eP"'//&
                     TRIM(ADJUSTL(WriteDx))//TRIM(ADJUSTL(WriteDy))//TRIM(ADJUSTL(WriteDz))//&
                     TRIM(ADJUSTL(WriteLS))//TRIM(ADJUSTL(WriteSt))
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TETRAHEDRON'
             END SELECT

             DO ie = 1, Mesh%Ne
                ele = mesh%ele(ie)
                eleSolve = .TRUE.
                DO i = 1, Data%N_LS
                   IF ( ele%Status(i) /= -1 ) eleSolve = .FALSE.
                END DO

                nNodesPerElement = ele%nNodesPerElement
                DO i = 1, nNodesPerElement
                   NU_i = ele%Vertex(i) ; Pos_i = ele%Pos(i)
                   Exact = 0.d0
                   IF ( EleSolve .OR. data%icing) THEN
                      IF ( .NOT. data%icing ) THEN
                         CALL ExactSolt(ele%L, Exact, Mesh%Vertex(NU_i)%X,Data%T, ele, nDim)
                      ELSE
                         CALL ExactSolJumpt(ele%L, Exact, Mesh%Vertex(NU_i)%X,Data%T, ele, nDim)
                      END IF
                   END IF

                   NrmB = 0.d0 ; NrmBe = 0.d0
                   DO id = 1, nDim
                      NrmB = NrmB + Sol%Bcg(NU_i,id)**2
                      NrmBe = NrmBe + Exact(id)**2
                   END DO
                   NrmB = sqrt(NrmB)
                   NrmBe = sqrt(NrmBe)
                   SELECT CASE ( Data%nDim )

                   CASE ( 2 )
                      WRITE(12,*) ele%Coor(i,:), ele%Coor(i,:), NU_i,&
                           Sol%Bcg(NU_i,:), NrmB, Sol%pcg(NU_i), &
                           Exact(1:nDim), NrmBe, Exact(nDim+1), &
                           abs(Sol%Bcg(NU_i,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                           abs(Sol%pcg(NU_i)-Exact(nDim+1)), &
                           Sol%Dist(:,NU_i,1), Sol%Dist(:,NU_i,2), &
                           Sol%LS(:,NU_i), ele%status(:)!sol%eleStatus(:,Pos_i)
                   CASE ( 3 )
                      WRITE(12,*) ele%Coor(i,:), ele%Coor(i,:), NU_i,&
                           Sol%Bcg(NU_i,:), NrmB, Sol%pcg(NU_i), &
                           Exact(1:nDim), NrmBe, Exact(nDim+1), &
                           abs(Sol%Bcg(NU_i,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                           abs(Sol%pcg(NU_i)-Exact(nDim+1)), &
                           Sol%Dist(:,NU_i,1), Sol%Dist(:,NU_i,2), sol%dist(:,NU_i,3),&
                           Sol%LS(:,NU_i), ele%status(:)
                   END SELECT
                END DO
             END DO

             DO i = 1, Mesh%Ne
                WRITE(12,*) Mesh%Ele(i)%Pos
             END DO
          END IF

       ELSE ! data%icing, we can have discontinuous solution -> as in DG

          OPEN(unit = 12, file = TRIM(ADJUSTL(Data%FileName))//'.'&
               //TRIM(ADJUSTL(SpSc))//'.'//char//'.tec')

          IF ( .NOT. data%embedded .AND. .NOT. data%immersed ) THEN

             WRITE(Nv_char,*) (nDim+1)*Mesh%Ne
             WRITE(Ne_char,*) Mesh%Ne

             WRITE(12,*) 'TITLE = Stefan'
             SELECT CASE ( nDim )
             CASE ( 2 )
                WRITE(12,*) 'VARIABLES = "x","y","Coor_x","Coor_y","NodeID","NodeDGID","Bx","By","NrmB","T","Bxe","Bye",&
                     "NrmBe","Te","eBx","eBy","eNrmB","eT"'
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TRIANGLE'
             CASE ( 3 )
                WRITE(12,*) 'VARIABLES = "x","y","z","Coor_x","Coor_y","Coor_z","NodeID","NodeDGID",&
                     "Bx","By","Bz","NrmB","p","Bxe","Bye",&
                     "Bze","NrmBe","pe","eBx","eBy","eBz","eNrmB","eP"'
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TETRAHEDRON'
             END SELECT
             DO i = 1, Mesh%Ne
                ele = Mesh%ele(i)
                nNodesPerElement = ele%nNodesPerElement
                DO j = 1, nNodesPerElement
                   Pos_j = ele%Pos(j)
                   NU_j = ele%Vertex(j)
                   surPos_j = ele%surVertex(j)
                   CALL ExactSolJumpt(ele%L, Exact, Mesh%Vertex(NU_j)%X, Data%T, ele, nDim)
                   NrmB = 0.d0
                   DO id = 1, nDim
                      NrmB = NrmB + Sol%Bcg(surPos_j,id)**2
                      NrmBe = NrmBe + Exact(id)**2
                   END DO
                   NrmB = sqrt(NrmB)
                   NrmBe = sqrt(NrmBe)
                   WRITE(12,*) Mesh%Vertex(NU_j)%X, Mesh%Vertex(NU_j)%X, NU_j, Pos_j,&
                        Sol%Bcg(surPos_j,:), NrmB, Sol%pcg(surPos_j), &
                        Exact(1:nDim), NrmBe, Exact(nDim+1), &
                        abs(Sol%Bcg(surPos_j,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                        abs(Sol%pcg(surPos_j)-Exact(nDim+1))
                END DO
             END DO
             DO i = 1, Mesh%Ne
                WRITE(12,*) Mesh%Ele(i)%Pos
             END DO

          ELSE

             WRITE(Nv_char,*) (Data%nDim+1)*Mesh%Ne
             WRITE(Ne_char,*) Mesh%Ne
             SELECT CASE ( Data%nDim )
             CASE ( 2 )
                WriteDx = ',"dist_x1"'
                WriteDy = ',"dist_y1"'
                WriteSt = ',"EleStatus1"'
                WriteLS = ',"LS1"'
                IF ( Data%N_LS > 1 ) THEN
                   DO i = 2, Data%N_LS
                      WRITE(ind,*) i
                      WriteDx = TRIM(ADJUSTL(WriteDx))//'"dist_x'//TRIM(ADJUSTL(ind))//'"'
                      WriteDy = TRIM(ADJUSTL(WriteDy))//'"dist_y'//TRIM(ADJUSTL(ind))//'"'
                      WriteSt = TRIM(ADJUSTL(WriteSt))//',"EleStatus'//TRIM(ADJUSTL(ind))//'"'
                      WriteLS = TRIM(ADJUSTL(WriteLS))//',"LS'//TRIM(ADJUSTL(ind))//'"'
                   END DO
                END IF
                WRITE(12,*) 'VARIABLES = "x", "y","Coor_x","Coor_y","NodeID","Bx","By","NrmB","T"&
                     ,"Bxe","Bye","NrmBe","Te","eBx","eBy","eNrmB","eT"'//&
                     TRIM(ADJUSTL(WriteDx))//TRIM(ADJUSTL(WriteDy))//&
                     TRIM(ADJUSTL(WriteLS))//TRIM(ADJUSTL(WriteSt))
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TRIANGLE'
             CASE ( 3 )
                WriteDx = ',"dist_x1"'
                WriteDy = ',"dist_y1"'
                WriteDz = ',"dist_z1"'
                WriteSt = ',"EleStatus1"'
                WriteLS = ',"LS1"'
                IF ( Data%N_LS > 1 ) THEN
                   DO i = 2, Data%N_LS
                      WRITE(ind,*) i
                      WriteDx = TRIM(ADJUSTL(WriteDx))//',"dist_x'//TRIM(ADJUSTL(ind))//'"'
                      WriteDy = TRIM(ADJUSTL(WriteDy))//',"dist_y'//TRIM(ADJUSTL(ind))//'"'
                      WriteDz = TRIM(ADJUSTL(WriteDz))//',"dist_z'//TRIM(ADJUSTL(ind))//'"'
                      WriteSt = TRIM(ADJUSTL(WriteSt))//',"EleStatus'//TRIM(ADJUSTL(ind))//'"'
                      WriteLS = TRIM(ADJUSTL(WriteLS))//',"LS'//TRIM(ADJUSTL(ind))//'"'
                   END DO
                END IF
                WRITE(12,*) 'VARIABLES = "x", "y","z","Coor_x","Coor_y","Coor_z", "NodeID","Bx","By","Bz",&
                     "NrmB","p","Bxe","Bye","Bze","NrmBe","pe","eBx","eBy","eBz","eNrmB","eP"'//&
                     TRIM(ADJUSTL(WriteDx))//TRIM(ADJUSTL(WriteDy))//TRIM(ADJUSTL(WriteDz))//&
                     TRIM(ADJUSTL(WriteLS))//TRIM(ADJUSTL(WriteSt))
                WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                     TRIM(ADJUSTL(Ne_char))//', ET=TETRAHEDRON'
             END SELECT

             DO ie = 1, Mesh%Ne
                ele = mesh%ele(ie)
                eleSolve = .TRUE.
                !DO i = 1, Data%N_LS
                !   IF ( ele%Status(i) /= -1 ) eleSolve = .FALSE.
                !END DO
				
				
                nNodesPerElement = ele%nNodesPerElement
                DO i = 1, nNodesPerElement
                   NU_i = ele%surVertex(i) ; Pos_i = ele%Pos(i)
                   Exact = 0.d0
                   IF ( eleSolve ) CALL ExactSolJumpt(ele%L, Exact, Mesh%Vertex(ele%Vertex(i))%X,Data%T, ele, nDim)
                   NrmB = 0.d0 ; NrmBe = 0.d0
                   DO id = 1, nDim
                      NrmB = NrmB + Sol%Bcg(NU_i,id)**2
                      NrmBe = NrmBe + Exact(id)**2
                   END DO
                   NrmB = sqrt(NrmB)
                   NrmBe = sqrt(NrmBe)
                   
                   SELECT CASE ( Data%nDim )
                    

                   CASE ( 2 )
                      WRITE(12,*) ele%Coor(i,:), ele%Coor(i,:), NU_i,&
                           Sol%Bcg(NU_i,:), NrmB, Sol%pcg(NU_i), &
                           Exact(1:nDim), NrmBe, Exact(nDim+1), &
                           abs(Sol%Bcg(NU_i,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                           abs(Sol%pcg(NU_i)-Exact(nDim+1)), &
                           Sol%Dist(:,NU_i,1), Sol%Dist(:,NU_i,2), &
                           Sol%LS(:,NU_i), ele%status(:)!sol%eleStatus(:,Pos_i)
                   CASE ( 3 )
                      WRITE(12,*) ele%Coor(i,:), ele%Coor(i,:), NU_i,&
                           Sol%Bcg(NU_i,:), NrmB, Sol%pcg(NU_i), &
                           Exact(1:nDim), NrmBe, Exact(nDim+1), &
                           abs(Sol%Bcg(NU_i,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                           abs(Sol%pcg(NU_i)-Exact(nDim+1)), &
                           Sol%Dist(:,NU_i,1), Sol%Dist(:,NU_i,2), &
                           sol%dist(:,NU_i,3),&
                           Sol%LS(:,Pos_i), sol%eleStatus(:,Pos_i)
                   END SELECT
                END DO
             END DO

             DO i = 1, Mesh%Ne
                WRITE(12,*) Mesh%Ele(i)%Pos
             END DO


          END IF

       END IF

       CLOSE(12)

    CASE ( "DG", "DG-Primal", "MSEG", "MSEG-Stab")

       OPEN(unit = 12, file = TRIM(ADJUSTL(Data%FileName))//'.'&
            //TRIM(ADJUSTL(SpSc))//'.'//char//'.tec')

       IF ( .NOT. data%embedded .AND. .NOT. data%immersed ) THEN

          WRITE(Nv_char,*) (nDim+1)*Mesh%Ne
          WRITE(Ne_char,*) Mesh%Ne

          WRITE(12,*) 'TITLE = DarcyFlow'
          SELECT CASE ( nDim )
          CASE ( 2 )
             WRITE(12,*) 'VARIABLES = "x","y","Coor_x","Coor_y","NodeID","NodeDGID","Bx","By","NrmB","p","Bxe","Bye",&
                  "NrmBe","pe","eBx","eBy","eNrmB","eP"'
             WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                  TRIM(ADJUSTL(Ne_char))//', ET=TRIANGLE'
          CASE ( 3 )
             WRITE(12,*) 'VARIABLES = "x","y","z","Coor_x","Coor_y","Coor_z","NodeID","NodeDGID",&
                  "Bx","By","Bz","NrmB","p","Bxe","Bye",&
                  "Bze","NrmBe","pe","eBx","eBy","eBz","eNrmB","eP"'
             WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                  TRIM(ADJUSTL(Ne_char))//', ET=TETRAHEDRON'
          END SELECT
          DO i = 1, Mesh%Ne
             ele = Mesh%ele(i)
             nNodesPerElement = ele%nNodesPerElement
             DO j = 1, nNodesPerElement
                Pos_j = ele%Pos(j)
                NU_j = ele%Vertex(j)
                CALL ExactSolt(ele%L, Exact, Mesh%Vertex(NU_j)%X,Data%T, Mesh%Ele(Mesh%Vertex(NU_j)%Ele1), nDim)
                NrmB = 0.d0
                DO id = 1, nDim
                   NrmB = NrmB + Sol%Bdg(Pos_j,id)**2
                   NrmBe = NrmBe + Exact(id)**2
                END DO
                NrmB = sqrt(NrmB)
                NrmBe = sqrt(NrmBe)
                WRITE(12,*) Mesh%Vertex(NU_j)%X, Mesh%Vertex(NU_j)%X, NU_j, Pos_j,&
                     Sol%Bdg(Pos_j,:), NrmB, Sol%pdg(Pos_j), &
                     Exact(1:nDim), NrmBe, Exact(nDim+1), &
                     abs(Sol%Bdg(Pos_j,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                     abs(Sol%pdg(Pos_j)-Exact(nDim+1))
             END DO
          END DO
          DO i = 1, Mesh%Ne
             WRITE(12,*) Mesh%Ele(i)%Pos
          END DO

       ELSE

          WRITE(Nv_char,*) (Data%nDim+1)*Mesh%Ne
          WRITE(Ne_char,*) Mesh%Ne
          SELECT CASE ( Data%nDim )
          CASE ( 2 )
             WriteDx = ',"dist_x1"'
             WriteDy = ',"dist_y1"'
             WriteSt = ',"EleStatus1"'
             WriteLS = ',"LS1"'
             IF ( Data%N_LS > 1 ) THEN
                DO i = 2, Data%N_LS
                   WRITE(ind,*) i
                   WriteDx = TRIM(ADJUSTL(WriteDx))//'"dist_x'//TRIM(ADJUSTL(ind))//'"'
                   WriteDy = TRIM(ADJUSTL(WriteDy))//'"dist_y'//TRIM(ADJUSTL(ind))//'"'
                   WriteSt = TRIM(ADJUSTL(WriteSt))//',"EleStatus'//TRIM(ADJUSTL(ind))//'"'
                   WriteLS = TRIM(ADJUSTL(WriteLS))//',"LS'//TRIM(ADJUSTL(ind))//'"'
                END DO
             END IF
             WRITE(12,*) 'VARIABLES = "x", "y","Coor_x","Coor_y","NodeID","Bx","By","NrmB","p"&
                  ,"Bxe","Bye","NrmBe","pe","eBx","eBy","eNrmB","eP"'//&
                  TRIM(ADJUSTL(WriteDx))//TRIM(ADJUSTL(WriteDy))//&
                  TRIM(ADJUSTL(WriteLS))//TRIM(ADJUSTL(WriteSt))
             WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                  TRIM(ADJUSTL(Ne_char))//', ET=TRIANGLE'
          CASE ( 3 )
             WriteDx = ',"dist_x1"'
             WriteDy = ',"dist_y1"'
             WriteDz = ',"dist_z1"'
             WriteSt = ',"EleStatus1"'
             WriteLS = ',"LS1"'
             IF ( Data%N_LS > 1 ) THEN
                DO i = 2, Data%N_LS
                   WRITE(ind,*) i
                   WriteDx = TRIM(ADJUSTL(WriteDx))//',"dist_x'//TRIM(ADJUSTL(ind))//'"'
                   WriteDy = TRIM(ADJUSTL(WriteDy))//',"dist_y'//TRIM(ADJUSTL(ind))//'"'
                   WriteDz = TRIM(ADJUSTL(WriteDz))//',"dist_z'//TRIM(ADJUSTL(ind))//'"'
                   WriteSt = TRIM(ADJUSTL(WriteSt))//',"EleStatus'//TRIM(ADJUSTL(ind))//'"'
                   WriteLS = TRIM(ADJUSTL(WriteLS))//',"LS'//TRIM(ADJUSTL(ind))//'"'
                END DO
             END IF
             WRITE(12,*) 'VARIABLES = "x", "y","z","Coor_x","Coor_y","Coor_z", "NodeID","Bx","By","Bz",&
                  "NrmB","p","Bxe","Bye","Bze","NrmBe","pe","eBx","eBy","eBz","eNrmB","eP"'//&
                  TRIM(ADJUSTL(WriteDx))//TRIM(ADJUSTL(WriteDy))//TRIM(ADJUSTL(WriteDz))//&
                  TRIM(ADJUSTL(WriteLS))//TRIM(ADJUSTL(WriteSt))
             WRITE(12,*) 'ZONE T="P1", F=FEPOINT, N='//TRIM(ADJUSTL(Nv_char))//', E='//&
                  TRIM(ADJUSTL(Ne_char))//', ET=TETRAHEDRON'
          END SELECT

          DO ie = 1, Mesh%Ne
             ele = mesh%ele(ie)
             eleSolve = .TRUE.
             DO i = 1, Data%N_LS
                IF ( ele%Status(i) /= -1 ) eleSolve = .FALSE.
             END DO

             nNodesPerElement = ele%nNodesPerElement
             DO i = 1, nNodesPerElement
                NU_i = ele%Vertex(i) ; Pos_i = ele%Pos(i)
                Exact = 0.d0
                IF ( eleSolve ) CALL ExactSolt(ele%L, Exact,Mesh%Vertex(NU_i)%X,Data%T, ele, nDim)
                NrmB = 0.d0 ; NrmBe = 0.d0
                DO id = 1, nDim
                   NrmB = NrmB + Sol%Bdg(Pos_i,id)**2
                   NrmBe = NrmBe + Exact(id)**2
                END DO
                NrmB = sqrt(NrmB)
                NrmBe = sqrt(NrmBe)
                SELECT CASE ( Data%nDim )

                CASE ( 2 )
                   WRITE(12,*) ele%Coor(i,:), ele%Coor(i,:), NU_i,&
                        Sol%Bdg(Pos_i,:), NrmB, Sol%pdg(Pos_i), &
                        Exact(1:nDim), NrmBe, Exact(nDim+1), &
                        abs(Sol%Bdg(Pos_i,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                        abs(Sol%pdg(Pos_i)-Exact(nDim+1)), &
                        Sol%Dist(:,Pos_i,1), Sol%Dist(:,Pos_i,2), &
                        Sol%LS(:,Pos_i), ele%status(:)!sol%eleStatus(:,Pos_i)
                CASE ( 3 )
                   WRITE(12,*) ele%Coor(i,:), ele%Coor(i,:), NU_i,&
                        Sol%Bdg(Pos_i,:), NrmB, Sol%pdg(Pos_i), &
                        Exact(1:nDim), NrmBe, Exact(nDim+1), &
                        abs(Sol%Bdg(Pos_i,:)-Exact(1:nDim)), abs(NrmB-NrmBe), &
                        abs(Sol%pdg(Pos_i)-Exact(nDim+1)), &
                        Sol%Dist(:,Pos_i,1), Sol%Dist(:,Pos_i,2), sol%dist(:,Pos_i,3),&
                        Sol%LS(:,Pos_i), sol%eleStatus(:,Pos_i)
                END SELECT
             END DO
          END DO

          DO i = 1, Mesh%Ne
             WRITE(12,*) Mesh%Ele(i)%Pos
          END DO
       END IF

       CLOSE(12)

    CASE DEFAULT

       CALL PrintError("WriteSol")

    END SELECT

    PRINT*, " -- Solution wrote "
    PRINT*, " "

  END SUBROUTINE WriteSolution
  !=============================!


END MODULE WriteSol
