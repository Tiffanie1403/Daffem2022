MODULE EmbeddedVariableInitialization

  USE Types
  USE PublicVar

CONTAINS

  !================================================!
  SUBROUTINE AllocEmbeddedVariables(Mesh,Data,Sol,Sol_prec,Mat)
  !================================================!

  !*************************************************************************
  !  allocates variables in the solution structure for embedded computation
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Solution structure for current time
  !    Input,Output, type(SolStructure) Sol_prec, Solution structure for previous time
  !    Input, type(SparseMatrix) Mat, structure for the matrix of the resolution
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !    Input, type(MeshType) Mesh, Mesh structure from dom.mesh
  !
  !*************************************************************************

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(SolStructure) , INTENT(INOUT) :: Sol, Sol_prec
    type(SparseMatrix) , INTENT(IN)    :: Mat
    type(MeshType)     , INTENT(IN)    :: Mesh
    !----------------------------------------
    INTEGER :: Ne, nDim, nNodesPerElement
    !------------------------------------------

     Ne = sol%Ne ; nDim = Data%nDim ; nNodesPerElement = nDim + 1

    PRINT*, ' '
    PRINT*, '  -- Allocation of necessary variables for embedded computations --  '
    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( 'CG', 'CG-Primal', 'DG', 'DG-Primal' )
       ALLOCATE ( Sol%dist(Data%N_LS, Sol%Nv, nDim) ) ! Data%N_LS est le nombre de condition embedded
       ALLOCATE ( Sol_prec%dist(Data%N_LS, Sol_prec%Nv, nDim) ) ! Data%N_LS est le nombre de condition embedded
       ALLOCATE ( Sol%LS( Data%N_LS, Sol%Nv) ) ; ALLOCATE ( Sol_prec%LS( Data%N_LS, Sol_prec%Nv) )
       ALLOCATE ( Sol%NodeStatus( Data%N_LS,nNodesPerElement*Ne ) )
       ALLOCATE ( Sol%EleStatus( Data%N_LS,nNodesPerElement*Ne ) )
       ALLOCATE ( Sol_prec%NodeStatus( Data%N_LS,nNodesPerElement*Ne ) )
       ALLOCATE ( Sol_prec%EleStatus( Data%N_LS,nNodesPerElement*Ne ) )
    CASE DEFAULT
       CALL PrintError("AllocEmbeddedVariables")
    END SELECT

    PRINT*, '  -- Embedded variables allocated  --  '

  END SUBROUTINE AllocEmbeddedVariables
  !=================================================!

  !================================================!
  SUBROUTINE AllocEmbeddedVariables_Tinit(Data,Sol,Mesh)
  !================================================!

    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(MeshType)     , INTENT(IN)    :: Mesh

    !----------------------------------------
    INTEGER :: Ne, nDim, nNodesPerElement,Nv
    !------------------------------------------

     Ne = Mesh%Ne ; nDim = Data%nDim ; nNodesPerElement = nDim + 1
     Nv = Mesh%Nv 

    PRINT*, ' '
    PRINT*, '  -- Allocation of necessary variables for embedded computations at Tinit --  '
    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( 'CG', 'CG-Primal', 'DG', 'DG-Primal' )
       ALLOCATE ( Sol%dist(Data%N_LS,Nv, nDim))
       ALLOCATE ( Sol%LS( Data%N_LS, Nv) )
       ALLOCATE ( Sol%NodeStatus( Data%N_LS,nNodesPerElement*Ne ) )
       ALLOCATE ( Sol%EleStatus( Data%N_LS,nNodesPerElement*Ne ))
       ALLOCATE ( Sol%N(Nv,nDim))

    CASE DEFAULT
       CALL PrintError("AllocEmbeddedVariables_Tinit")
    END SELECT

    PRINT*, '  -- Embedded variables allocated  --  '

  END SUBROUTINE AllocEmbeddedVariables_Tinit
  !=================================================!


END MODULE EmbeddedVariableInitialization
