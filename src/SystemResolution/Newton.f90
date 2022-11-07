MODULE Newton

  USE Types
  USE SparseMatrix_CG
  USE SparseMatrix_DG

  IMPLICIT NONE
  REAL*8, PUBLIC :: NewtonResidual, normalizedNewtonResidual, NewtonIterationCounter
  REAL*8, PRIVATE :: NewtonErr0
CONTAINS

  !===============================================!
  SUBROUTINE ComputeResidual(data, sol, mesh, mat)
  !===============================================!

  !*************************************************************************
  !  Calculates the residual of the newton iterations 
  !
  !  Parameters:
  !
  !    Input, type (SparseMatrix) Mat, Structure for the global matrix
  !    Input,Output, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************


    IMPLICIT NONE
    type(dataStructure), INTENT(IN)    :: Data
    type(solStructure) , INTENT(INOUT) :: Sol
    type(meshType)     , INTENT(IN)    :: Mesh
    type(SparseMatrix) , INTENT(IN)    :: Mat
    !--------------------------------------------
    REAL*8 :: aii, aij, ui, uj
    REAL*8 :: residual, res2
    !--------------------------
    INTEGER :: i, j, nAdj, NU_j
    INTEGER :: NDof, Nv, Ne, nNodesPerElement, nDim,Ns
    !-----------------------------------------------
    CHARACTER(len=99) :: resC, nresC
    !---------------------------------


    IF (Data%Icing .eqv. .TRUE.) then
      Ns=Mesh%NarEmb + 1 + mesh%PtEmb
    ELSE
      Ns=0
    END IF

    Nv = mesh%NvActive+Ns ; Ne = mesh%NeActive
    nDim = data%nDim ; nNodesPerElement = nDim + 1
    sol%rhsN = 0.d0

    SELECT CASE ( TRIM(ADJUSTL(data%space_Scheme)) )
    CASE ( 'CG-Primal','CG' )
       nDof = Nv
    CASE ( 'DG-Primal' )
       nDof = Ne*nNodesPerElement
    END SELECT

    residual = 0.d0
    DO i = 1, nDof
       residual = residual + sol%rhsTmp(i)**2
    END DO

    residual = sqrt(residual/nDof)

    NewtonResidual = residual
    IF ( NewtonIterationCounter == 1 ) NewtonErr0 = residual
    normalizedNewtonResidual = NewtonResidual/NewtonErr0

    WRITE(resC,*) residual
    WRITE(nresC,*) normalizedNewtonResidual

    PRINT*, " "
    PRINT*, " ** residual of Newton Iteration = "//TRIM(ADJUSTL(resC))//"  "//TRIM(ADJUSTL(nresC))
    PRINT*, " "

  END SUBROUTINE ComputeResidual
  !===============================================!

END MODULE Newton
