MODULE MUMPS_Interface

  USE Types
  USE PublicVar

  IMPLICIT NONE

#ifdef MUMPS
  include 'dmumps_struc.h'
#endif

CONTAINS

#ifdef MUMPS
  !========================================================!
  SUBROUTINE MUMPS_Init(Id, Mesh, Sol, NDof, Mat, Data, ST)
  !========================================================!

  !*************************************************************************
  !  Creation of the Mumps structure for the resolution by MUMPS
  !
  !  Parameters:
  !
  !    Input,Output, type(DMUMPS_STRUC) Id, Mumps structure for the resolution
  !    Input,type(MeshType)  Mesh, Mesh structure from dom.mesh
  !    Input,type(SolStructure)  Sol, Solution structure
  !    Input,INTEGER NDof, the number of degrees of freedom
  !    Input,Output,type(SparseMatrix) Mat, Structure for the resolution matrice
  !    Input,type(DataStructure) Data, Data structure from dom.data
  !    Input,CHARACTER(len=*) ST, ???
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DMUMPS_STRUC) , INTENT(INOUT) :: Id
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SolStructure) , INTENT(IN)    :: Sol
    INTEGER            , INTENT(IN)    :: NDof
    type(SparseMatrix) , INTENT(INOUT) :: Mat
    type(DataStructure), INTENT(IN)    :: Data
    CHARACTER(len=*)   , INTENT(IN)    :: ST
    !------------------------------------------------
    type(element) :: ele
    !--------------------
    INTEGER :: Nmax, i, j, cpt, ied, iv
    INTEGER :: Nv, Ne, Num
    INTEGER :: k, l, NU_l, m, NU_k
    INTEGER :: Ntri, Nn, Ntloc, NnLoc
    INTEGER :: OSv, OSt
    INTEGER :: iG, jG, iM, jM, NU_j
    INTEGER :: ie, Pos_i, Pos_k, NvLoc, NnL, NedL, Ned
    !----------------------------------------------------

    CALL MUMPS_Parameters(id, Data)

    PRINT*, '       Fill MUMPS matrix...'

    Nv = Mesh%NvActive ; Ne = Mesh%NeActive ; Ned = Mesh%Ned
    CALL MUMPS_Alloc(Id, Mat, NDof, Nv, Ne, Ned, Data)
    CALL MUMPS_ResolutionFillMat(id, NDof, Nv, Ne, Ned, Mat, Mesh, Data)
    IF (.NOT. data%NewtonIterations )  DEALLOCATE ( Mat%Coeff )
    id%RHS = Sol%RHS

    PRINT*, '      MUMPS Matrix filled'


  END SUBROUTINE MUMPS_Init
  !===========================================!

  !====================================!
  SUBROUTINE MUMPS_Parameters(id, Data)
  !====================================!

  !*************************************************************************
  !  Parameters for the MUMPS Structure for the resolution
  !
  !  Parameters:
  !
  !    Input,Output, type(DMUMPS_STRUC) Id, Mumps structure for the resolution
  !    Input,type(DataStructure) Data, Data structure from dom.data
  !
  !*************************************************************************

    IMPLICIT NONE
    type(DMUMPS_STRUC) , INTENT(INOUT) :: id
    type(DataStructure), INTENT(IN)    :: Data
    !------------------------------------------

    Id%SYM = 0
    Id%PAR = 1
    Id%JOB = -1
    id%ICNTL(1:3) = 3 ! prevents from printing...
    CALL DMUMPS(id)

    id%ICNTL(1) = 6
    id%ICNTL(2) = 0
    id%ICNTL(3) = 0
!!$
    Id%ICNTL(7) = 7
    id%ICNTL(14) = Data%MUMPS_Option
    id%ICNTL(11) = 2
    id%JOB = 6

  END SUBROUTINE MUMPS_Parameters
  !==============================!

  !=======================================================!
  SUBROUTINE MUMPS_Alloc(Id, Mat, NDof, Nv, Ne, Ned, Data)
  !=======================================================!

  !*************************************************************************
  !  Allocation of the MUMPS Structure for the resolution
  !
  !  Parameters:
  !
  !    Input,Output, type(DMUMPS_STRUC) Id, Mumps structure for the resolution
  !    Input,type(SparseMatrix) Mat, Structure for the resolution matrice
  !    Input,type(DataStructure) Data, Data structure from dom.data
  !    Input,INTEGER NDof, the number of degrees of freedom
  !    Input,INTEGER Nv, the number of vertices
  !    Input,INTEGER Ne, the number of elements
  !    Input,INTEGER Ned, the number of edges
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DMUMPS_STRUC) , INTENT(INOUT) :: id
    type(SparseMatrix) , INTENT(IN)    :: Mat
    type(DataStructure), INTENT(IN)    :: Data
    INTEGER            , INTENT(IN)    :: NDof, Nv, Ne, Ned
    !--------------------------------------------------------
    INTEGER :: nVar
    !-----------------

    nVar = Mat%nVar

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( 'CG', 'DG', 'MSEG', 'MSEG-Stab' )
       Id%N = nVar*NDof
    CASE ( 'CG-Primal', 'DG-Primal' )
       Id%N = NDof
    CASE DEFAULT
       CALL PrintError("MUMPS_Alloc")
    END SELECT

    Id%NZ = Mat%NZ

    ALLOCATE( Id%IRN( Id%NZ ) )
    ALLOCATE( Id%JCN( Id%NZ ) )
    ALLOCATE( Id%A( Id%NZ ) )
    ALLOCATE( Id%RHS( Id%N ) )

  END SUBROUTINE MUMPS_Alloc
  !==================================================!

  !==========================================================================!
  SUBROUTINE MUMPS_ResolutionFillMat(id, NDof, Nv, Ne, Ned, Mat, Mesh, Data)
  !==========================================================================!

  !*************************************************************************
  !  Sparse storage for the resolution by MUMPS
  !
  !  Parameters:
  !
  !    Input,type(MeshType)  Mesh, Mesh structure from dom.mesh
  !    Input,type(SparseMatrix) Mat, Structure for the resolution matrice
  !    Input,type(DataStructure) Data, Data structure from dom.data
  !    Input,INTEGER NDof, the number of degrees of freedom
  !    Input,INTEGER Nv, the number of vertices
  !    Input,INTEGER Ne, the number of elements
  !    Input,INTEGER Ned, the number of edges
  !    Input,Output, type(DMUMPS_STRUC) id, Mumps structure for the resolution
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DMUMPS_STRUC) , INTENT(INOUT) :: id
    type(SparseMatrix) , INTENT(IN)    :: Mat
    INTEGER            , INTENT(IN)    :: NDof, Nv, Ne, Ned
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(DataStructure), INTENT(IN)    :: Data
    !-------------------------------------------------------
    type(element) :: ele
    !-----------------------
    INTEGER :: Nmax, i, j, k, l, m, cpt, NnL, ie, ied
    INTEGER :: iG, jG, iM, jM, OSt, Nn, Pos_i, Pos_k
    INTEGER :: iv, NedL, nnLoc, ntLoc, nTri, NU_i, NU_j, OffSet
    INTEGER :: nvLoc, OSv, NU_k, jd, kd, NnTmp, NtriTmp, NtriLoc, ivar, jvar
    INTEGER :: OS, NVar, nNodesPerElement
    !------------------------------------------------------------------------

    cpt = 0
    NMax = Mat%NMax ; nVar = Mat%Nvar

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG", "CG-Primal", "MSEG", "MSEG-Stab" )
       DO i = 1, NDof
          DO l = 1, nVar
             DO m = 1, nVar
                cpt = cpt + 1
                iG = i + (l-1)*Nv ; jG = 1 + (m-1)*Nmax
                iM = i + (l-1)*Nv ; jM = i + (m-1)*Nv
                Id%IRN(cpt) = iM ; Id%JCN(cpt) = jM
                Id%A(cpt) = Mat%Coeff(iG,jG)
             END DO
          END DO
          NnL = Mat%Ind(i,1)
          DO j = 1, NnL
             NU_j = Mat%Ind(i,j+1)
             DO l = 1, nVar
                DO m = 1, nVar
                   cpt = cpt+1
                   iG = i + (l-1)*Nv ; jG = 1 + j + (m-1)*NMax
                   iM = i + (l-1)*Nv ; jM = Mat%Ind(i,j+1) + (m-1)*Nv
                   Id%IRN(cpt) = iM ; Id%JCN(cpt) = jM
                   Id%A(cpt) = Mat%Coeff(iG,jG)
                END DO
             END DO
          END DO
       END DO

    CASE ( "DG-Primal", "DG" )

       Nmax = Mat%Nmax
       DO ie = 1, Ne
          ele = Mesh%Ele(ie)
          nNodesPerElement = ele%nNodesPerElement
          DO i = 1, nNodesPerElement
             Pos_i = ele%Pos(i)
             DO l = 1, nVar
                DO m = 1, nVar
                   cpt = cpt + 1
                   iG = Pos_i + (l-1)*NDof ; jG = 1 + (m-1)*Nmax
                   iM = Pos_i + (l-1)*NDof ; jM = Pos_i + (m-1)*NDof
                   Id%IRN(cpt) = iM ; Id%JCN(cpt) = jm
                   Id%A(cpt) = Mat%Coeff(iG,jG)
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
                      Id%IRN(cpt) = iM ; Id%JCN(cpt) = jM
                      Id%A(cpt) = Mat%Coeff(iG,jG)
                   END DO
                END DO
             END DO
          END DO
       END DO

    CASE DEFAULT
       CALL PrintError("MUMPS_ResolutionFillMat_3D")
    END SELECT

  END SUBROUTINE MUMPS_ResolutionFillMat
  !=============================================!

  !============================!
  SUBROUTINE MUMPS_FINALIZE(id)
  !============================!

  !*************************************************************************
  !  Resolution by MUMPS
  !
  !  Parameters:
  !
  !    Input,Output, type(DMUMPS_STRUC) Id, Mumps structure for the resolution
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DMUMPS_STRUC), INTENT(INOUT) :: id
    !---------------------------------------

    Id%JOB = -2
    CALL DMUMPS(id)

  END SUBROUTINE MUMPS_FINALIZE
  !============================!

#endif

END MODULE MUMPS_Interface
