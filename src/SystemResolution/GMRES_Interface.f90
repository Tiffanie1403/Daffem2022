MODULE GMRES_Interface

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !=====================================================!
  SUBROUTINE GMRES_Init(GMRES_Mat, Sol, NDof, Mat, Data)
  !=====================================================!

  !*************************************************************************
  !   Creation and allocation of variable to handle the resolution by GMRES
  !
  !  Parameters:
  !
  !    Input,Output, type(GMRESMat) GMRES_Mat, Structure of the matrix for the resolution by GMRES
  !    Input,Output, type(SolStructure) Sol, Solution strucutre
  !    Input, type(SparseMatrix) Mat, Structure for the resolution matrice
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !    Input, INTEGER NDof, the number of degrees of freedom
  !
  !*************************************************************************



    IMPLICIT NONE
    type(GMRESMat)     , INTENT(INOUT) :: GMRES_Mat
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(IN)    :: Mat
    type(DataStructure), INTENT(IN)    :: Data
    INTEGER            , INTENT(IN)    :: NDof
    !------------------------------------------------
    CHARACTER(len=99) :: ST
    !-----------------------
    INTEGER :: nDim, nVar
    !---------------------

    nDim = Data%nDim ; nVar = Mat%nVar
    ST = Data%SolverType

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( 'CG', 'DG', 'MSEG', 'MSEG-Stab')
       GMRES_Mat%N = nVar*NDof !nombre de variable fois le nombre de noeud
    CASE ( 'CG-Primal', 'DG-Primal' )
       GMRES_Mat%N = NDof
    CASE DEFAULT
       CALL PrintError("GMRES_Init")
    END SELECT

    GMRES_Mat%NZ = Mat%NZ

    IF ( TRIM(ADJUSTL(ST)) == "GMRES" ) THEN
       ALLOCATE( GMRES_Mat%IA( GMRES_Mat%NZ ) )
    ELSE
       ALLOCATE( GMRES_Mat%IA( GMRES_Mat%N+1 ) )
    END IF
    ALLOCATE( GMRES_Mat%JA( GMRES_Mat%NZ ) )
    ALLOCATE( GMRES_Mat%A( GMRES_Mat%NZ ) )
    ALLOCATE( GMRES_Mat%RHS( GMRES_Mat%N ) )

  END SUBROUTINE GMRES_Init
  !===========================================================!

  !==============================================================!
  SUBROUTINE GMRES_FillMat(GMRES_Mat, Mesh, Sol, NDof, Ns, Mat, Data)
  !==============================================================!

  !*************************************************************************
  !  Sparse storage for the resolution by GMRES
  !
  !  Parameters:
  !
  !    Input,Output, type(GMRESMat) GMRES_Mat, Structure of the matrix for the resolution by GMRES
  !    Input,Output, type(SolStructure) Sol, Solution strucutre
  !    Input,type(MeshType)  Mesh, Mesh structure from dom.mesh
  !    Input,type(SparseMatrix) Mat, Structure for the resolution matrice
  !    Input,type(DataStructure) Data, Data structure from dom.data
  !    Input,INTEGER NDof, the number of degrees of freedom
  !    Input,INTEGER Ns, the number of complementary nodes
  !
  !*************************************************************************



    IMPLICIT NONE
    type(GMRESMat)     , INTENT(INOUT) :: GMRES_Mat
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SolStructure) , INTENT(IN)    :: Sol
    type(SparseMatrix) , INTENT(IN)    :: Mat
    type(DataStructure), INTENT(IN)    :: Data
    INTEGER            , INTENT(IN)    :: NDof, Ns
    !--------------------------------------------------
    type(element) :: ele
    !---------------------
    INTEGER :: Nmax, i, j, k, l, m, cpt, NnL, ie, ied
    INTEGER :: iG, jG, iM, jM, OSt, Nn, Pos_i, Pos_k
    INTEGER :: iv, NedL, nnLoc, ntLoc, nTri, NU_i, NU_j, OffSet
    INTEGER :: nvLoc, OSv, NU_k, jd, kd, NnTmp, NtriTmp, NtriLoc, ivar, jvar
    INTEGER :: OS, NVar, Nv, Ne, Ned, nNodesPerElement, NeActive
    !------------------------------------------------------------------------

    Nv = Mesh%NvActive + Ns ; Ne = Mesh%Ne ; Ned = Mesh%Ned ; nVar = Mat%nVar
    NeActive = mesh%NeActive
    cpt = 0
    NMax = Mat%NMax

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG" , "CG-Primal", "MSEG", "MSEG-Stab" )
       DO i = 1, NDof
          DO l = 1, nVar
             DO m = 1, nVar
                cpt = cpt + 1
                !iG = i + (l-1)*Nv
                iG = i + (l-1)*(Mesh%NvActive + Data%MaxRowEmbedded)
                jG = 1 + (m-1)*Nmax ; iM = i + (l-1)*Nv ; jM = i + (m-1)*Nv
                GMRES_Mat%IA(cpt) = iM ; GMRES_Mat%JA(cpt) = jM
                GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
             END DO
          END DO
          NnL = Mat%Ind(i,1)
          DO j = 1, NnL
             NU_j = Mat%Ind(i,j+1)
             DO l = 1, nVar
                DO m = 1, nVar
                   cpt = cpt+1
                  ! iG = i + (l-1)*Nv
                   iG = i + (l-1)*(Mesh%NvActive + Data%MaxRowEmbedded)
                   jG = 1 + j + (m-1)*NMax
                   iM = i + (l-1)*Nv ; jM = Mat%Ind(i,j+1) + (m-1)*Nv
                   GMRES_Mat%IA(cpt) = iM ; GMRES_Mat%JA(cpt) = jM
                   GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
                END DO
             END DO
          END DO
       END DO

    CASE ( "DG" , "DG-Primal" )

       Nmax = Mat%Nmax
       DO ie = 1, Ne
          ele = Mesh%Ele(ie)
          IF ( ele%solve ) THEN
             nNodesPerElement = ele%nNodesPerElement
             DO i = 1, nNodesPerElement
                Pos_i = ele%surPos(i)
                DO l = 1, nVar
                   DO m = 1, nVar
                      cpt = cpt + 1
                      iG = Pos_i + (l-1)*NDof ; jG = 1 + (m-1)*Nmax
                      iM = Pos_i + (l-1)*NDof ; jM = Pos_i + (m-1)*NDof
                      GMRES_Mat%IA(cpt) = iM ; GMRES_Mat%JA(cpt) = jm
                      GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
                   END DO
                END DO
                NvLoc = Mat%Ind(Pos_i,1)
                DO k = 1, NvLoc
                   Pos_k = Mat%Ind(Pos_i,1+k)
                   DO l = 1, nVar
                      DO m = 1, nVar
                         cpt = cpt + 1
                         iG = Pos_i + (l-1)*NDof ; jG = 1 + k + (m-1)*Nmax
                         iM = Pos_i + (l-1)*NDof ; jM = Pos_k + (m-1)*NDof
                         GMRES_Mat%IA(cpt) = iM ; GMRES_Mat%JA(cpt) = jM
                         GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
                      END DO
                   END DO
                END DO
             END DO
          END IF
       END DO

    CASE DEFAULT
       CALL PrintError("GMRES_FillMat")
    END SELECT

  END SUBROUTINE GMRES_FillMat
  !==============================================================!

  !=================================================================!
  SUBROUTINE GMRESILU_FillMat(GMRES_Mat, Mesh, Sol, NDof, Mat, Data)
  !=================================================================!

  !*************************************************************************
  !  Sparse storage for the resolution by GMRES with LU decomposition
  !
  !  Parameters:
  !
  !    Input,Output, type(GMRESMat) GMRES_Mat, Structure of the matrix for the resolution by GMRES
  !    Input,Output, type(SolStructure) Sol, Solution strucutre
  !    Input,type(MeshType)  Mesh, Mesh structure from dom.mesh
  !    Input,type(SparseMatrix) Mat, Structure for the resolution matrice
  !    Input,type(DataStructure) Data, Data structure from dom.data
  !    Input,INTEGER NDof, the number of degrees of freedom
  !
  !*************************************************************************



    IMPLICIT NONE
    type(GMRESMat)     , INTENT(INOUT) :: GMRES_Mat
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SolStructure) , INTENT(IN)    :: Sol
    type(SparseMatrix) , INTENT(IN)    :: Mat
    type(DataStructure), INTENT(IN)    :: Data
    INTEGER            , INTENT(IN)    :: NDof
    !--------------------------------------------------
    type(element) :: ele
    !---------------------
    INTEGER :: Nmax, i, j, k, l, m, cpt, NnL, ie, ied
    INTEGER :: iG, jG, iM, jM, OSt, Nn, Pos_i, Pos_k, MaxRowEmbedded
    INTEGER :: iv, NedL, nnLoc, ntLoc, nTri, NU_i, NU_j, OffSet,iG_bis
    INTEGER :: nvLoc, OSv, NU_k, jd, kd, NnTmp, NtriTmp, NtriLoc, ivar, jvar
    INTEGER :: OS, NVar, Nv, Ne, Ned, nNodesPerElement, NeActive, Ns, Nvv
    !------------------------------------------------------------------------

    IF (Data%Icing .eqv. .TRUE.) then
      Ns=Mesh%NB_VertInterf
    ELSE
      Ns=0
    END IF
    
    Nv = Mesh%NvActive+Ns !contains the complementary nodes 
    Ne = Mesh%Ne ; Ned = Mesh%Ned ; nVar = Mat%nVar
    Nvv = Mesh%NvActive ; MaxRowEmbedded = Data%MaxRowEmbedded
    NeActive = mesh%NeActive
    cpt = 0
    NMax = Mat%NMax

    !GMRES_Mat%RHS = Sol%RHS ! il faut changer la copie ici
    GMRES_Mat%RHS(1:Nv) = Sol%RHS(1:Nv)
    IF(Data%Space_Scheme=="CG") THEN !only the mixed form has an explicit definition of B 
		GMRES_Mat%RHS(Nv+1:2*Nv) = Sol%RHS(Nvv + MaxRowEmbedded + 1 : Nvv + MaxRowEmbedded + Nv )
		GMRES_Mat%RHS(2*Nv+1:3*Nv) = Sol%RHS(2*(Nvv + MaxRowEmbedded)+1:2*(Nvv + MaxRowEmbedded) + Nv)
	END IF 
	
    GMRES_Mat%IA = 0
    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG" , "CG-Primal", "MSEG", "MSEG-Stab")
        GMRES_Mat%IA(1) = 1
        DO l = 1, nVar
            DO i = 1, NDof
				DO m = 1, nVar
					cpt = cpt + 1
				  !  iG = i + (l-1)*Nv
					iG = i + (l-1)*(Mesh%NvActive + Data%MaxRowEmbedded)
					iG_bis =  i + (l-1)*Nv
					jG = 1 + (m-1)*Nmax
					iM = i + (l-1)*Nv ; jM = i + (m-1)*Nv
					!GMRES_Mat%IA(iG + 1) = GMRES_Mat%IA(iG+1) + 1
					GMRES_Mat%IA(iG_bis + 1) = GMRES_Mat%IA(iG_bis+1) + 1
					GMRES_Mat%JA(cpt) = jM
					GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
					!WRITE(*,*) "GMRES_Mat%IA ",  GMRES_Mat%IA(iG + 1) ," GMRES_Mat%JA ", &
					! GMRES_Mat%JA(cpt), "GMRES_Mat%A ",GMRES_Mat%A(cpt)
				END DO
				NnL = Mat%Ind(i,1)
				DO j = 1, NnL
					NU_j = Mat%Ind(i,j+1)
					DO m = 1, nVar
					   cpt = cpt+1
					   !iG = i + (l-1)*Nv
					   iG_bis = i + (l-1)*Nv
					   iG = i + (l-1)*(Mesh%NvActive + Data%MaxRowEmbedded)
					   jG = 1 + j + (m-1)*NMax
					   iM = i + (l-1)*Nv ; jM = Mat%Ind(i,j+1) + (m-1)*Nv
					   GMRES_Mat%IA(iG_bis + 1) = GMRES_Mat%IA(iG_bis+1) + 1
					   GMRES_Mat%JA(cpt) = jM ; GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
					  ! WRITE(*,*) "GMRES_Mat%IA ",  GMRES_Mat%IA(iG + 1) ," GMRES_Mat%JA ",&
					  !   GMRES_Mat%JA(cpt), "GMRES_Mat%A ",GMRES_Mat%A(cpt)
					END DO
				END DO
			END DO
		END DO
    CASE ( "DG", "DG-Primal" )

       Nmax = Mat%Nmax ; NVar = Mat%NVar
       GMRES_Mat%IA(1) = 1
       DO l = 1, NVar
          DO ie = 1, Ne
             ele = Mesh%Ele(ie)
             IF ( ele%solve ) THEN
                nNodesPerElement = ele%nNodesPerElement
                DO i = 1, nNodesPerElement
                   Pos_i = ele%surPos(i)
                   DO m = 1, nVar
                      cpt = cpt + 1
                      iG = Pos_i + (l-1)*NDof ; jG = 1 + (m-1)*Nmax
                      iM = Pos_i + (l-1)*NDof ; jM = Pos_i + (m-1)*NDof
                      GMRES_Mat%IA(iG+1) = GMRES_Mat%IA(iG+1) + 1
                      GMRES_Mat%JA(cpt) = jM
                      GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
                   END DO
                   NvLoc = Mat%Ind(Pos_i,1)
                   DO k = 1, NvLoc
                      Pos_k = Mat%Ind(Pos_i,1+k)
                      DO m = 1, NVar
                         cpt = cpt + 1
                         iG = Pos_i + (l-1)*NDof ; jG = 1 + k + (m-1)*Nmax
                         iM = Pos_i + (l-1)*NDof ; jM = Pos_k + (m-1)*NDof
                         GMRES_Mat%IA(iG+1) = GMRES_Mat%IA(iG+1) + 1
                         GMRES_Mat%JA(cpt) = jM
                         GMRES_Mat%A(cpt) = Mat%Coeff(iG,jG)
                      END DO
                   END DO
                END DO
             END IF
          END DO
       END DO
    CASE DEFAULT
       CALL PrintError("GMRES_FillMat")
    END SELECT
    DO i = 1, SIZE(GMRES_Mat%IA)-1
       GMRES_Mat%IA(i+1) = GMRES_Mat%IA(i+1) + GMRES_Mat%IA(i)
    END DO

  END SUBROUTINE GMRESILU_FillMat
  !==============================================================!


END MODULE GMRES_Interface
