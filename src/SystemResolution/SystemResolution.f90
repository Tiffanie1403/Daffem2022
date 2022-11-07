MODULE SystemResolution

  USE Types
  USE MUMPS_Interface
  USE GMRES_Interface
  USE GMRES

  IMPLICIT NONE

CONTAINS

  !======================================================!
  SUBROUTINE SolveSystem(Data, Sol, Mesh, Mat, SolveType)
  !======================================================!

  !*************************************************************************
  !  Solves the system
  !
  !  Parameters:
  !
  !    Input,Output, type (SparseMatrix) Mat, Structure for the global matrix
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !    Input, CHARACTER(len=*) SolveType      IF(Mesh%NarEmb==Mesh%Nb_EmbBord ) THEN 
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    CHARACTER(len=*)   , INTENT(IN)    :: SolveType
    !-----------------------------------------------
#ifdef MUMPS
    type(DMUMPS_STRUC) :: Id
#endif
    type(GMRESMat)    :: GMRES_Mat
    !------------------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: U, Test
    !---------------------------------------------
    INTEGER :: Nv, Ne, NDof, i, Ned, nVar,Ns, k, Nvv
    !-------------------------------------

    IF (Data%Icing .eqv. .TRUE.) then
      Ns=Mesh%NB_VertInterf
    ELSE
      Ns=0
    END IF
	
    Nv = Mesh%NvActive + Ns ; Nvv = Mesh%NvActive
    !Nv = Sol%Nv + Data%MaxRowEmbedded
    Ne = Mesh%NeActive ; Ned = Mesh%Ned ; nVar = Mat%nVar
    !allocation du vecteur solution U
    CALL ResolutionInitialization( NDof, U, Nv, Ne, Ned, Data, nVar)


    SELECT CASE ( TRIM(ADJUSTL(Data%SolverType)) )

    CASE ( "MUMPS" )
#ifdef MUMPS
        PRINT*, '    ****    MUMPS Resolution    ****'

        CALL CPU_TIME(t1_matrixTransfer)
        CALL MUMPS_Init(Id, Mesh, Sol, NDof, Mat, Data, SolveType)
        CALL CPU_TIME(t2_matrixTransfer)
        CALL CPU_TIME(t1_solve)
        CALL DMUMPS(id)

        PRINT*, '    ****    MUMPS Over **** '

        Sol%RHS = id%RHS

        CALL MUMPS_Finalize(Id)
        CALL CPU_TIME(t2_solve)
#else
        PRINT*, " "
        PRINT*, "    !!!!    !!!!    !!!!    !!!!    !!!!    !!!!    !!!!    !!!!    "
        PRINT*, "SOLVER ASKED IS MUMPS BUT COMPILATION WITHOUT INCLUDING LIBRARY !!!!"
        PRINT*, "    !!!!    !!!!    !!!!    !!!!    !!!!    !!!!    !!!!    !!!!    "
        PRINT*, " "
        STOP
#endif

    CASE ( "GMRES" )

        PRINT*, '    ****    GMRES Resolution    ****'
        CALL CPU_TIME(t1_matrixTransfer)
        CALL GMRES_Init( GMRES_Mat, Sol, NDof, Mat, Data )
        CALL GMRES_FillMat(GMRES_Mat, Mesh, Sol, NDof, Ns, Mat, Data)
        U = Sol%RHS
        CALL CPU_TIME(t2_matrixTransfer)

        CALL CPU_TIME(t1_solve)
        CALL mgmres_st(GMRES_Mat%N, GMRES_Mat%NZ, GMRES_Mat%IA, GMRES_Mat%JA, &
            GMRES_Mat%A, U, GMRES_Mat%RHS, Data%IterMax, Data%NKrilov, &
            Data%TolAbs, Data%TolRel) 
        Sol%RHS = U
        PRINT*, '    ****    GMRES Over **** '

    CASE ( "GMRES-ILU" )
       !cas qui nous interesse
        PRINT*, '    ****    GMRES-ILU Resolution    ****'
        CALL GMRES_Init( GMRES_Mat, Sol, NDof, Mat, Data )
        CALL GMRESILU_FillMat(GMRES_Mat, Mesh, Sol, NDof, Mat, Data)

        U = Sol%RHS 

        CALL pmgmres_ilu_cr(GMRES_Mat%N, GMRES_Mat%NZ, GMRES_Mat%IA, GMRES_Mat%JA, &
            GMRES_Mat%A, U, GMRES_Mat%RHS, Data%IterMax, Data%NKrilov, &
            Data%TolAbs, Data%TolRel)
        ! la solution du système linéaire se trouve dans le vecteur U
        ! Sol%RHS = U
        !on laisse des espace de 0 dans le vecteurs 
        Sol%RHS = 0.d0 ; Sol%RHS(1:Nv) = U(1:Nv)
        IF(Data%Space_Scheme=="CG") THEN !only in mixed form the flux is explicitly defined 
		   Sol%RHS(Nvv + Data%MaxRowEmbedded + 1 : Nvv + Data%MaxRowEmbedded + Nv ) = U(Nv+1:2*Nv)
		   Sol%RHS(2*(Nvv + Data%MaxRowEmbedded)+1:2*(Nvv + Data%MaxRowEmbedded) + Nv) = U(2*Nv+1:3*Nv)
		END IF 

        CALL CPU_TIME(t2_solve)
        PRINT*, '    ****    GMRES-ILU Over **** '

    END SELECT

  END SUBROUTINE SolveSystem
  !================================!

  !=====================================================================!
  SUBROUTINE ResolutionInitialization( NDof, U, Nv, Ne, Ned, Data, nVar)
  !=====================================================================!

  !*************************************************************************
  !  Defines the number of degrees of freedom depending of the type of resolution
  !  'CG', 'MSDG', 'MSEG', 'MSEG-Stab','CG-Primal', 'MSDG-Primal','DG', 'DG-Primal'
  !
  !  Parameters:
  !
  !    Output,INTEGER NDof, the number of degrees of freedom
  !    Output, REAL*8, DIMENSION(:) U, vector solution
  !    Input, INTEGER Nv, numver of vertices
  !    Input, INTEGER  Ne, number of elements
  !    Input, INTEGER  Ned, number of edges
  !    Input, INTEGER  nVar, number of variables
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************


    IMPLICIT NONE
    INTEGER                          , INTENT(OUT) :: NDof
    REAL*8, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: U
    INTEGER                          , INTENT(IN)  :: Nv, Ne, Ned, nVar
    type(DataStructure)              , INTENT(IN)  :: Data
    !-------------------------------------------------------------------
    INTEGER :: nNodesPerElement, nDim
    !----------------------------------

    nDim = Data%nDim ; nNodesPerElement = nDim + 1

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( 'CG', 'MSDG', 'MSEG', 'MSEG-Stab' )
       NDof = Nv !nombre de noeud de la triangulation sans ajout des noeuds complementaires à l'interface 
       ALLOCATE(U(nVar*NDof))
    CASE ( 'CG-Primal', 'MSDG-Primal' )
       NDof = Nv
       ALLOCATE( U(NDof) )
    CASE ( 'DG', 'DG-Primal' )
       NDof = Ne*nNodesPerElement
       ALLOCATE( U(nVar*NDof) )
    CASE DEFAULT
       CALL PrintError("ResolutionInitialization")
    END SELECT

  END SUBROUTINE ResolutionInitialization
  !================================================================!

END MODULE SystemResolution
