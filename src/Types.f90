  MODULE Types

  !*************************************************************************
  ! Defines the different variable types
  !*************************************************************************

  IMPLICIT NONE

  !*****************************************!
  !********        Mesh Type        ********!

  !** Vertex **!
  type Vertex
     REAL*8 , DIMENSION(:), ALLOCATABLE :: X
     INTEGER, DIMENSION(:), ALLOCATABLE :: Ref
     INTEGER, DIMENSION(:), ALLOCATABLE :: Status
     INTEGER							:: Boundary
     INTEGER							:: Ele1
     INTEGER							:: Num
     LOGICAL							:: Bound
     LOGICAL							:: Active,LogPtEmb
     LOGICAL							:: OnOutBond
     LOGICAL							:: OnInterface
     LOGICAL							:: OnSurrogate
     INTEGER							:: Complementaire
     LOGICAL							:: EmptyComp

  end type Vertex

  !** Boundary elements -> Edges or triangles **!
  type Edge
     INTEGER, DIMENSION(:), ALLOCATABLE :: Vertex
     REAL*8,  DIMENSION(:), ALLOCATABLE :: N
     INTEGER, DIMENSION(2)              :: Tri
     INTEGER                            :: Ref
     INTEGER                            :: BC_Pos
     LOGICAL                            :: Treated,InterfaceBound
     CHARACTER(len=20)                  :: BC_Type
     CHARACTER(len=5)                   :: eleType
     INTEGER                            :: nNodesPerEdge
     INTEGER                            :: num
     INTEGER                            :: EmbLS
  end type Edge

  !** Element -> Triangle or Tetrahedrons **!
    type Element
		 REAL*8 , DIMENSION(:,:), ALLOCATABLE :: N
		 REAL*8 , DIMENSION(:,:), ALLOCATABLE :: Coor
		 REAL*8 , DIMENSION(:,:), ALLOCATABLE :: L, Lm1
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: Vertex
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: SurVertex
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: SurVertex_Sauvg
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: Adj
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: Pos
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: surPos
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: Edg
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: Status
		 INTEGER, DIMENSION(:)  , ALLOCATABLE :: EdgRef
		 INTEGER, DIMENSION(:,:), ALLOCATABLE :: EdgNum
		 REAL*8                               :: A
		 REAL*8                               :: lambda
		 INTEGER                              :: Num
		 INTEGER                              :: Ref
		 INTEGER                              :: RefEdg
		 LOGICAL                              :: Solve, Id,Cut
		 CHARACTER(len=5)                     :: eleType
		 INTEGER                              :: nNodesPerElement, nEdgesPerElement !nombre de ddl et d'aretes par elements
  end type Element

  !** Mesh **!
  type MeshType
     INTEGER      , DIMENSION(:),   ALLOCATABLE :: MappingSur2Real
     INTEGER      , DIMENSION(:),   ALLOCATABLE :: MappingReal2Sur
     INTEGER      , DIMENSION(:),   ALLOCATABLE :: mappingSur2RealPos
     INTEGER      , DIMENSION(:),   ALLOCATABLE :: mappingReal2SurPos
     INTEGER      , DIMENSION(:),   ALLOCATABLE :: ActiveElement
     type(Vertex) , DIMENSION(:),   ALLOCATABLE :: Vertex
     type(Element), DIMENSION(:),   ALLOCATABLE :: Ele
     type(Edge)   , DIMENSION(:),   ALLOCATABLE :: Boundary, embBoundary
     type(Edge)   , DIMENSION(:),   ALLOCATABLE :: Edg
     INTEGER      , DIMENSION(:,:), ALLOCATABLE :: Tab_NoeudComp, Tab_NoeudComp_prec, Tab_NoeudComp_prec_inter
     INTEGER                                    :: Nv, Ne, Nb, Ndim, Ned, NvActive, NeActive
     INTEGER                                    :: NSurB, NarEmb,NarEmb_ref,PtEmb,NarEmb_prec,NarEmb_prec_inter
	 INTEGER                                    :: Nb_EmbBord
	 INTEGER 									:: NB_VertInterf
  end type MeshType

  !** Item -> For filling Hashtable for Mesh construction **!
  type Item2D
     integer               :: Tri    ! 3*Num_tri + numero arete local
     integer, dimension(2) :: Vertex ! Numéro global des noeuds
     integer               :: Nxt    ! pour dépacement
  end type Item2D
  
  type VectorInterface
     REAL*8, dimension(:,:), ALLOCATABLE :: Coor ! Numéro global des noeuds
     REAL*8, dimension(:,:), ALLOCATABLE :: Coor_prec ! Numéro global des noeuds
     REAL*8, dimension(:,:), ALLOCATABLE :: Coor_2prec ! Numéro global des noeuds
  end type VectorInterface

  type Item3D
     integer               :: Tri
     integer, dimension(3) :: Vertex
     integer               :: Nxt
  end type Item3D

  !** Hashtable **!
  type HashTable2D
     type(Item2D), DIMENSION(:), ALLOCATABLE :: Item
     INTEGER                                 :: SizeH ! Taille de la HashTable
     INTEGER                                 :: SizeD ! Taille du dépacement
     INTEGER                                 :: Compt ! Compteur de Nxt
     INTEGER     , DIMENSION(:), ALLOCATABLE :: EdgNum
     type(Edge)  , DIMENSION(:), ALLOCATABLE :: Edg
  end type HashTable2D

  type HashTable3D
     type(Item3D), DIMENSION(:), ALLOCATABLE :: Item
     INTEGER                                 :: SizeH ! Taille de la HashTable
     INTEGER                                 :: SizeD ! Taille du dépacement
     INTEGER                                 :: Compt ! Compteur de Nxt
     INTEGER     , DIMENSION(:), ALLOCATABLE :: EdgNum
     type(Edge)  , DIMENSION(:), ALLOCATABLE :: Edg
  end type HashTable3D

  !********    All For Mesh Above   ********!
  !*****************************************!

  !*****************************************!
  !********        Data Type        ********!

  type DataStructure
     INTEGER           								:: Ndim
     INTEGER           								:: icas
     CHARACTER(len=99)  							:: FileName, Inversion
     CHARACTER(len=99)  							:: Space_Scheme, Space_Type 
     REAL*8              							:: MUMPS_Option
     INTEGER             							:: MUMPS_Verbose
     REAL*8               							:: alphaP, alphaB, alphaIP1, alphaIP2
     REAL*8               							:: alphaS, C2 
     REAL*8, DIMENSION(:), ALLOCATABLE 				:: alphaD, betaN, betaT
     REAL*8   , DIMENSION(:,:), ALLOCATABLE			:: lambda
     REAL*8   , DIMENSION(:,:), ALLOCATABLE 		:: lambdam1
     CHARACTER(len=20), DIMENSION(:), ALLOCATABLE 	:: BC
     INTEGER, DIMENSION(:), ALLOCATABLE 			:: BC_ref
     INTEGER, DIMENSION(:)  , ALLOCATABLE 			:: Dom
     REAL*8   , DIMENSION(:)  , ALLOCATABLE 		:: BCVal, BCEmbVal
     INTEGER 										:: N_BC, N_dom, SetEquation
     LOGICAL 										:: Embedded, immersed, StabTransfer, PressureEnrichment, testFunctionEnrichment
     LOGICAL 										:: Output, ComputeError, EnrichStab  
     INTEGER 										:: N_LS, MaxRowEmbedded, tour
     CHARACTER(len=20), DIMENSION(:,:), ALLOCATABLE :: LSDef
     CHARACTER(len=20), DIMENSION(:)  , ALLOCATABLE :: EmbeddedBCType
     REAL*8 , DIMENSION(:,:), ALLOCATABLE 			:: LS_L, LS_center,Tab_Prev2Jump,CoordResistance
     REAL*8 , DIMENSION(:)  , ALLOCATABLE 			:: Emb_alphaD, Emb_BetaN, Emb_BetaT
     type(Edge), DIMENSION(:), ALLOCATABLE 			:: TmpEmbEdg
     LOGICAL, DIMENSION(:), ALLOCATABLE 			:: TreatedNode, TreatedPos
     LOGICAL 										:: PermAnisotropy
     CHARACTER(len = 20) 							:: SolverType, TScheme, InterfaceScheme
     INTEGER 										:: NKrilov, IterMax, KPrint,TourMax,nNewtonIter, IterJump
     REAL*8  										:: ErrorT, ErrorBx, ErrorBy,PrevJumpSurrogate,Jump
     REAL*8											:: PrevJumpSurrogate_x, PrevJumpSurrogate_y
     REAL*8  										:: TolAbs, TolRel ,NewtonTolerance, SvgSpeed2Prec
     LOGICAL 										:: NewtonIterations, DtChoiceCFL, SpeedCst, Icing, NewtonProcess
     REAL*8  										:: Tmax , h, CFL, Dt, T, Speed, Icing_Prev_Velocity, Tinit 
     REAL*8											:: ResNewt, ResNewtD, InitResNewt, Speed_prec, MeanT
     LOGICAL										:: DefInterfaceNode
     LOGICAL										:: ResisTance
     REAL*8 										:: Hl, Ha, Hi, Lh, L
     REAL*8											:: InitialTime, FinalTime 
     
  end type DataStructure

  !********    All For Data Above   ********!
  !*****************************************!

  !******************************************!
  !******** Sparse Matrix Storing    ********!
  type SparseMatrix
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: Ind
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: IndTri
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: IndEdg
     REAL*8 , DIMENSION(:,:), ALLOCATABLE :: Coeff
     REAL*8 , DIMENSION(:,:), ALLOCATABLE :: CoeffTri
     REAL*8 , DIMENSION(:,:), ALLOCATABLE :: CoeffEdg
     INTEGER                              :: NMax
     INTEGER                              :: Nn   ! Directly connected nodes
     INTEGER                              :: NTri ! Triangle the node belongs to (connected or not)
     INTEGER                              :: Ned  ! Edges with whom the node interact
     INTEGER                              :: NZ, NVar
  end type SparseMatrix
  !******************************************!

  !***********************************************!
  !********    GMRES Library interface    ********!
  type GMRESMat
     REAL*8 , DIMENSION(:), ALLOCATABLE :: A, RHS
     INTEGER, DIMENSION(:), ALLOCATABLE :: IA
     INTEGER, DIMENSION(:), ALLOCATABLE :: JA
     INTEGER                            :: NZ, N
  END type GMRESMat
  !***********************************************!

	!**********************************!
	!********    Solution Type ********!
	type SolStructure
		REAL*8   , DIMENSION(:,:), ALLOCATABLE 	    :: Bcg, Beg, BcgRecons, Bdg, N,Bcg_inter,Bcg_Newton, BcgN,Bcg_prec
		REAL*8   , DIMENSION(:)  , ALLOCATABLE 	    :: pcg, peg, pcgRecons, pdg,Pcg_inter, PcgN,Pcg_prec
		REAL*8   , DIMENSION(:,:), ALLOCATABLE 	    :: DBcg, DBeg, DBdg
		REAL*8   , DIMENSION(:)  , ALLOCATABLE 	    :: Dpcg, Dpeg, Dpdg
		REAL*8   , DIMENSION(:)  , ALLOCATABLE 	    :: Brt, pcgNp1,DualArea
		REAL*8   , DIMENSION(:,:), ALLOCATABLE 	    :: grad_pcg
		REAL*8   , DIMENSION(:,:,:), ALLOCATABLE 	:: Dist
		REAL*8   , DIMENSION(:,:)  , ALLOCATABLE 	:: LS
		INTEGER  , DIMENSION(:,:)  , ALLOCATABLE 	:: NodeStatus, EleStatus
		REAL*8   , DIMENSION(:)  , ALLOCATABLE 		:: RHS, rhsN, rhs0, rhsTmp
		INTEGER                              		:: Nv, Ne, Ned, Nb, NZ, nVar
		INTEGER                              		:: NvActive, NeActive
		INTEGER, DIMENSION(:,:), ALLOCATABLE 		:: Indices
	end type SolStructure
	!******** All for Solution above ********!
	!****************************************!

END MODULE Types
