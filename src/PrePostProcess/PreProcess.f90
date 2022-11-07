MODULE PreProcess

  USE Types
  USE GradientReconstruction

  IMPLICIT NONE

CONTAINS

  !===========================================!
  SUBROUTINE PostProcess(Data,Mesh,Sol,Mat)
  !===========================================!

  !*************************************************************************
  !   Calls the appropriate routine to calculate the flux in the case
  !   of a problem in primal form
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
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SparseMatrix) , INTENT(IN)    :: Mat

    CALL GreenGaussReconstruction_Primal(Sol, Mesh, data)

  END SUBROUTINE PostProcess
  !===========================================!

END MODULE PreProcess
