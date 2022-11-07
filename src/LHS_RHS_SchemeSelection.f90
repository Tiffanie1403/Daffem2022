
MODULE LHS_RHS_SchemeSelection

  USE Types
  USE PublicVar
  USE Data_mod
  USE Compute_CG
  USE Compute_CG_Explicit
  USE Compute_CG_Primal
  USE Compute_DG
  USE Compute_DG_Primal
  USE Compute_MSEG

  IMPLICIT NONE

CONTAINS

  !=======================================================!
  SUBROUTINE RHS_LHS_SchemeSelection(Data, Mesh, Sol, Sol_prec, Sol_2prec, Mat, PhysicalInterface)
  !========================================================!

  !*************************************************************************
  !   Sends to the correct routine for the computation of the matrix
  !   The choice depends on Data%Space_Scheme
  !  Parameters:
  !
  !    Input, type(DataStructure) Data, Structure for the date from the dom.data file
  !    Output, type(MeshType) Mesh, Mesh structure
  !    Input, type(SolStructure) Sol, Sol structure for the current time
  !    Input, type(SolStructure) Sol_prec, Sol structure for the previous time
  !    Input, type(SparseMatrix) Mat, Structure for the global matrix
  !
  !*************************************************************************

    IMPLICIT NONE
    type(DataStructure), INTENT(INOUT) :: Data
    type(MeshType)     , INTENT(INOUT) :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(SolStructure) , INTENT(IN)    :: Sol_prec, Sol_2prec 
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    type(VectorInterface), INTENT(IN)  :: PhysicalInterface 
    !-----------------------------------------

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( 'CG' )
		CALL Compute_LHS_RHS_CG(Data, Mesh, Sol, Sol_prec, Sol_2prec, Mat, PhysicalInterface) !modif en temps seulement sur le cas CG pour le moment
    CASE ( 'CG-Primal' )
		IF (Data%Space_Type == "EXPLICIT") THEN 
			CALL Compute_LHS_RHS_CG_Explicit(Data, Mesh, Sol, Sol_prec, Mat) !primal form but explicit formulation
		ELSE 
		   CALL Compute_LHS_RHS_CG_Primal(Data, Mesh, Sol, Mat) !stay intact for implicit formulation 
		END IF 
    CASE ( 'DG' )
       CALL Compute_LHS_RHS_DG(Data, Mesh, Sol, Mat)
    CASE ( 'DG-Primal' )
       CALL Compute_LHS_RHS_DG_Primal(Data, Mesh, Sol, Mat)
    CASE ( 'MSEG', 'MSEG-Stab' )
       CALL Compute_LHS_RHS_MSEG(data, mesh, sol, mat)
    CASE DEFAULT
       CALL PrintError("LHS_RHS_SchemeSelection")
    END SELECT

  END SUBROUTINE RHS_LHS_SchemeSelection
  !==========================================================!

END MODULE LHS_RHS_SchemeSelection
