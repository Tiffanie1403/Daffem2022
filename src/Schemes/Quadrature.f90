MODULE Quadrature

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !==================================!
  SUBROUTINE GetNQuad(eleType, nQuad)
  !==================================!

    IMPLICIT NONE
    CHARACTER(len=5), INTENT(IN)  :: eleType
    INTEGER         , INTENT(OUT) :: nQuad
    !-----------------------------------------

    SELECT CASE ( eleType )
    CASE ( "TET3D" )
       nQuad = 4
    CASE ( "TRI3D", "TRI2D" )
       nQuad = 3
    CASE ( "EDG2D" )
       nQuad = 3
    CASE ( "DEFAULT" )
       PRINT*, "ERROR in GetNQuad, unkwown element"
       STOP
    END SELECT

  END SUBROUTINE GetNQuad
  !==================================!

  !===================================!
  SUBROUTINE CoorQuad(eleType, Xq, iq)
  !===================================!

    IMPLICIT NONE
    CHARACTER(len=5)    , INTENT(IN)    :: eleType
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Xq
    INTEGER             , INTENT(IN)    :: iq
    !-----------------------------------------------

    SELECT CASE ( eleType )
    CASE ( "TET3D" )
       CALL TetCoorQuad(Xq, iq)
    CASE ( "TRI2D" )
       CALL TriCoorQuad(Xq, iq)
    CASE ( "DEFAULT" )
       PRINT*, "ERROR in CoorQuad, unkwown element"
       STOP
    END SELECT

  END SUBROUTINE CoorQuad
  !==========================!

  !===========================================!
  SUBROUTINE FaceCoorQuad(EdgType, Xq, iq, ii)
  !===========================================!

    IMPLICIT NONE
    CHARACTER(len=5)    , INTENT(IN)    :: EdgType
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Xq
    INTEGER             , INTENT(IN)    :: iq, ii
    !-----------------------------------------------

    SELECT CASE (EdgType)
    CASE ( "TRI3D" )
       CALL TetFaceCoorQuad(Xq, iq, ii)
    CASE ( "EDG2D" )
       CALL TriFaceCoorQuad(Xq, iq, ii)
    CASE DEFAULT
       PRINT*, "ERROR in FaceCoorQuad, unknown case"
       STOP
    END SELECT


  END SUBROUTINE FaceCoorQuad
  !===========================================!

  !=====================================!
  SUBROUTINE TetFaceCoorQuad(Xq, iq, ii)
  !=====================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Xq
    INTEGER             , INTENT(IN)    :: iq, ii
    !----------------------------------------------

    SELECT CASE ( ii )
    CASE ( 1 )
       SELECT CASE ( iq )
       CASE ( 1 )
          Xq(1) = 2.d0/3.d0
          Xq(2) = 1.d0/6.d0
          Xq(3) = 1.d0/6.d0
       CASE ( 2 )
          Xq(1) = 1.d0/6.d0
          Xq(2) = 2.d0/3.d0
          Xq(3) = 1.d0/6.d0
       CASE ( 3 )
          Xq(1) = 1.d0/6.d0
          Xq(2) = 1.d0/6.d0
          Xq(3) = 2.d0/3.d0
       END SELECT
    CASE ( 2 )
       SELECT CASE ( iq )
       CASE ( 1 )
          Xq(1) = 0.d0
          Xq(2) = 1.d0/6.d0
          Xq(3) = 1.d0/6.d0
       CASE ( 2 )
          Xq(1) = 0.d0
          Xq(2) = 2.d0/3.d0
          Xq(3) = 1.d0/6.d0
       CASE ( 3 )
          Xq(1) = 0.d0
          Xq(2) = 1.d0/6.d0
          Xq(3) = 2.d0/3.d0
       END SELECT
    CASE ( 3 )
       SELECT CASE ( iq )
       CASE ( 1 )
          Xq(1) = 1.d0/6.d0
          Xq(2) = 0.d0
          Xq(3) = 1.d0/6.d0
       CASE ( 2 )
          Xq(1) = 2.d0/3.d0
          Xq(2) = 0.d0
          Xq(3) = 1.d0/6.d0
       CASE ( 3 )
          Xq(1) = 1.d0/6.d0
          Xq(2) = 0.d0
          Xq(3) = 2.d0/3.d0
       END SELECT
    CASE ( 4 )
       SELECT CASE ( iq )
       CASE ( 1 )
          Xq(1) = 1.d0/6.d0
          Xq(2) = 1.d0/6.d0
          Xq(3) = 0.d0
       CASE ( 2 )
          Xq(1) = 2.d0/3.d0
          Xq(2) = 1.d0/6.d0
          Xq(3) = 0.d0
       CASE ( 3 )
          Xq(1) = 1.d0/6.d0
          Xq(2) = 2.d0/3.d0
          Xq(3) = 0.d0
       END SELECT
    END SELECT

  END SUBROUTINE TetFaceCoorQuad
  !=====================================!

  !=====================================!
  SUBROUTINE TriFaceCoorQuad(Xq, iq, ii)
  !=====================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Xq
    INTEGER             , INTENT(IN)    :: iq, ii
    !----------------------------------------------

    SELECT CASE ( ii )
    CASE ( 1 )
       SELECT CASE ( iq )
       CASE ( 1 )
          Xq(1) = 1.d0
          Xq(2) = 0.d0
       CASE ( 2 )
          Xq(1) = 0.5d0
          Xq(2) = 0.5d0
       CASE ( 3 )
          Xq(1) = 0.d0
          Xq(2) = 1.d0
       END SELECT
    CASE ( 2 )
       SELECT CASE ( iq )
       CASE ( 1 )
          Xq(1) = 0.d0
          Xq(2) = 1.d0
       CASE ( 2 )
          Xq(1) = 0.d0
          Xq(2) = 0.5d0
       CASE ( 3 )
          Xq(1) = 0.d0
          Xq(2) = 0.d0
       END SELECT
    CASE ( 3 )
       SELECT CASE ( iq )
       CASE ( 1 )
          Xq(1) = 0.d0
          Xq(2) = 0.d0
       CASE ( 2 )
          Xq(1) = 0.5d0
          Xq(2) = 0.d0
       CASE ( 3 )
          Xq(1) = 1.d0
          Xq(2) = 0.d0
       END SELECT
    END SELECT

  END SUBROUTINE TriFaceCoorQuad
  !=====================================!

  !=============================!
  SUBROUTINE TetCoorQuad(Xq, iq)
  !=============================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Xq
    INTEGER             , INTENT(IN)    :: iq
    !------------------------------------------
    REAL*8 :: c1, c2
    !-----------------

    c1 = 0.5854101966249685 ; c2 =0.1381966011250105
    SELECT CASE (iq)
    CASE (1)
       Xq(1) = c1
       Xq(2) = c2
       Xq(3) = c2
    CASE (2)
       Xq(1) = c2
       Xq(2) = c2
       Xq(3) = c2
    CASE(3)
       Xq(1) = c2
       Xq(2) = c2
       Xq(3) = c1
    CASE ( 4 )
       Xq(1) = c2
       Xq(2) = c1
       Xq(3) = c2
    END SELECT

  END SUBROUTINE TetCoorQuad
  !==========================!

  !=============================!
  SUBROUTINE TriCoorQuad(Xq, iq)
  !=============================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: Xq
    INTEGER             , INTENT(IN)    :: iq
    !------------------------------------------

    SELECT CASE (iq)
    CASE (1)
       Xq(1) = 1.d0/6.d0 ; Xq(2) = 1.d0/6.d0
    CASE (2)
       Xq(1) = 2.d0/3.d0 ; Xq(2) = 1.d0/6.d0
    CASE(3)
       Xq(1) = 1.d0/6.d0 ; Xq(2) = 2.d0/3.d0
    END SELECT

  END SUBROUTINE TriCoorQuad
  !==========================!

  !=====================================!
  SUBROUTINE WeightQuad(eleType, wq, iq)
  !=====================================!

    IMPLICIT NONE
    CHARACTER(len=5), INTENT(IN)  :: eleType
    REAL*8          , INTENT(OUT) :: wq
    INTEGER         , INTENT(IN)  :: iq
    !------------------------------------------

    SELECT CASE ( eleType )
    CASE ( "TET3D" )
       wq = 1.d0/4.d0
    CASE ( "TRI2D" )
       wq = 1.d0/3.d0
    CASE ( "EDG2D" )
       SELECT CASE ( iq )
       CASE ( 1 )
          wq = 1.d0/6.d0
       CASE ( 2 )
          wq = 2.d0/3.d0
       CASE ( 3 )
          wq = 1.d0/6.d0
       END SELECT
       !wq = 0.5d0*wq
    CASE ( "TRI3D" )
       wq = 1.d0/3.d0
    CASE ( "DEFAULT" )
       PRINT*, "ERROR in WeightQuad, unknown element"
       STOP
    END SELECT

  END SUBROUTINE WeightQuad
  !=====================================!

  !====================================!
  SUBROUTINE SFEval(eleType, SF, Xq, i)
  !====================================!

    IMPLICIT NONE
    CHARACTER(len=5)    , INTENT(IN)  :: eleType
    REAL*8              , INTENT(OUT) :: SF
    REAL*8, DIMENSION(:), INTENT(IN)  :: Xq
    INTEGER             , INTENT(IN)  :: i
    !--------------------------------------------
    REAL*8 :: x, y, z
    !------------------

    SELECT CASE ( eleType )
    CASE ( "TET3D" )
       x = Xq(1) ; y = Xq(2) ; z = Xq(3)
       SELECT CASE ( i )
       CASE ( 1 )
          SF = 1.d0 - x - y - z
       CASE ( 2 )
          SF = x
       CASE ( 3 )
          SF = y
       CASE ( 4 )
          SF = z
       CASE ( 5 )
          SF = 1.d0
       END SELECT
    CASE ( "TRI2D" )
       x = Xq(1) ; y = Xq(2)
       SELECT CASE ( i )
       CASE ( 1 )
          SF = 1.d0 - x - y
       CASE ( 2 )
          SF = x
       CASE ( 3 )
          SF = y
       CASE ( 4 )
          SF = 1.d0
       END SELECT
    END SELECT

  END SUBROUTINE SFEval
  !====================================!

  !==================================!
  SUBROUTINE GSFEval(GSF, xq, i, ele)
  !==================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: GSF
    REAL*8, DIMENSION(:), INTENT(IN)    :: xq
    INTEGER             , INTENT(IN)    :: i
    type(element)       , INTENT(IN)    :: ele
    !-------------------------------------------

    SELECT CASE ( ele%eleType )
    CASE ( "TET3D" )
       SELECT CASE ( i )
       CASE ( 1, 2, 3, 4 )
          GSF = ele%N(i,:)/(3.d0*ele%A)
       CASE ( 5 )
          GSF = 0.d0
       END SELECT
    CASE ( "TRI2D" )
       SELECT CASE ( i )
       CASE ( 1, 2, 3 )
          GSF = ele%N(i,:)/(2.d0*ele%A)
       CASE ( 4 )
          GSF = 0.d0
       END SELECT
    CASE DEFAULT
       PRINT*, "ERROR in GSFEval, unknown element"
       STOP
    END SELECT

  END SUBROUTINE GSFEval
  !==================================!

  !===========================================!
  SUBROUTINE IdentifyIFace(ele, edg, iOK, nDim)
  !===========================================!

    IMPLICIT NONE
    type(element), INTENT(IN)  :: ele
    type(Edge)   , INTENT(IN)  :: Edg
    INTEGER      , INTENT(OUT) :: iOK
    INTEGER      , INTENT(IN)  :: nDim
    !------------------------------------
    INTEGER, DIMENSION(4,3)            :: Permut
    INTEGER, DIMENSION(:), ALLOCATABLE :: N1, N2
    !------------------------------------------------
    INTEGER :: i, j, nNodesPerElement
    INTEGER :: Ne1, Ne2, Nt1, Nt2
    !----------------------------------
    LOGICAL :: found
    !-----------------

    Permut(1,1) = 2 ; Permut(1,2) = 3 ; Permut(1,3) = 4
    Permut(2,1) = 3 ; Permut(2,2) = 1 ; Permut(2,3) = 4
    Permut(3,1) = 1 ; Permut(3,2) = 2 ; Permut(3,3) = 4
    Permut(4,1) = 1 ; Permut(4,2) = 2 ; Permut(4,3) = 3

    ALLOCATE ( N1(nDim), N2(nDim) )

    nNodesPerElement = ele%nNodesPerElement
    N1 = Edg%Vertex
    CALL OrderTable(N1)

    DO i = 1, nNodesPerElement
       N2 = ele%Vertex(Permut(i,1:nDim))
       CALL OrderTable(N2)
       found = .TRUE.
       DO j = 1, SIZE(N2)
          IF ( N1(j) /= N2(j) ) THEN
             found = .FALSE.
             EXIT
          END IF
       END DO
       IF ( found .eqv. .TRUE. ) THEN
          iOK = i
          EXIT
       END IF
    END DO

    IF ( found .eqv. .FALSE. ) THEN
       PRINT*, "Error in IdentifyI_BC, did not found local edge number"
       STOP
    END IF

    DEALLOCATE ( N1, N2 )

  END SUBROUTINE IdentifyIFace
  !=====================================!

  !=========================================================!
  SUBROUTINE QuadratureOnEdge(xq, ele, edg, iL, nDim)
  !=========================================================!

    IMPLICIT NONE
    type(element)         , INTENT(IN)  :: ele
    type(Edge)            , INTENT(IN)  :: edg
    REAL*8, DIMENSION(:,:), INTENT(OUT) :: xq
    INTEGER               , INTENT(IN)  :: iL, nDim
    !--------------------------------------------------
    REAL*8, DIMENSION(SIZE(xq,2)) :: xTmp
    !--------------------------------------------------
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: confQuad, conf
    INTEGER, DIMENSION(:)  , ALLOCATABLE :: N1, N2
    !-------------------------------------------------------
    INTEGER, DIMENSION(4,3) :: Permut
    !---------------------------------
    INTEGER :: i, j, nNodesPerEdge, id, iConf
    INTEGER :: nQuad
    !----------------------------------------
    LOGICAL :: foundConf
    !--------------------

    nNodesPerEdge = edg%nNodesPerEdge

    IF ( nDim == 2 ) THEN
       Permut(1,1) = 2 ; Permut(1,2) = 3 ; Permut(1,3) = 4
       Permut(2,1) = 3 ; Permut(2,2) = 1 ; Permut(2,3) = 4
       Permut(3,1) = 1 ; Permut(3,2) = 2 ; Permut(3,3) = 4
       Permut(4,1) = 1 ; Permut(4,2) = 2 ; Permut(4,3) = 3
       ALLOCATE ( conf(2,2) )
       ALLOCATE ( confQuad(2,3) )
       Conf(1,1) = 1 ; Conf(1,2) = 2
       Conf(2,1) = 2 ; Conf(2,2) = 1
       confQuad(1,1) = 1 ; confQuad(1,2) = 2 ; confQuad(1,3) = 3
       confQuad(2,1) = 3 ; confQuad(2,2) = 2 ; confQuad(2,3) = 1
       nQuad = 3
    ELSE
       Permut(1,1) = 2 ; Permut(1,2) = 3 ; Permut(1,3) = 4
       Permut(2,1) = 1 ; Permut(2,2) = 3 ; Permut(2,3) = 4
       Permut(3,1) = 1 ; Permut(3,2) = 2 ; Permut(3,3) = 4
       Permut(4,1) = 1 ; Permut(4,2) = 2 ; Permut(4,3) = 3
       ALLOCATE ( conf(6,3) )
       ALLOCATE ( confQuad(6,3) )
       Conf(1,1) = 1 ; Conf(1,2) = 2 ; Conf(1,3) = 3
       Conf(2,1) = 1 ; Conf(2,2) = 3 ; Conf(2,3) = 2
       Conf(3,1) = 2 ; Conf(3,2) = 1 ; Conf(3,3) = 3
       Conf(4,1) = 2 ; Conf(4,2) = 3 ; Conf(4,3) = 1
       Conf(5,1) = 3 ; Conf(5,2) = 1 ; Conf(5,3) = 2
       Conf(6,1) = 3 ; Conf(6,2) = 2 ; Conf(6,3) = 1
       confQuad = conf
       nQuad = 3
    END IF

    ALLOCATE ( N1(nDim), N2(nDim) )
    
    !on verifie que l'arete recherché fait bient partie de l'élément envoyé 
    N1 = edg%Vertex !reference for the nodes of the edge 
    N2 = ele%Vertex(Permut(iL,1:nDim))    ! numerotation d'ordre des deux noeuds de l'arete 3-2, 3-1,1-2
    DO id = 1, SIZE(confQuad,1)       
       foundConf = .TRUE.
       DO i = 1, nNodesPerEdge
          IF ( N1(conf(id,i)) /= N2(i) ) THEN
             foundConf = .FALSE.
             EXIT
          END IF
       END DO
       IF ( foundConf ) THEN
          iConf = id
          EXIT
       END IF
    END DO
    IF ( foundConf .eqv. .FALSE. ) THEN
       PRINT*, "Error in QuadratureOnEdge, no conf found"
       PRINT*, N1
       PRINT*, N2
       STOP
    END IF

	
    DO i = 1, nQuad
       CALL FaceCoorQuad(edg%eleType, xTmp, i, iL)
       Xq(confQuad(iConf,i),:) = xTmp
    END DO

    DEALLOCATE ( N1, N2, conf, confQuad )

	!si on est sur l'arete 1 il resort (1,0), (0.5,0.5), (0,1) 
	! si on est sur l'arête 2 il resort (0,1), (0,0.5), (0,0)
	!si on est sur l'arete 3 il resort (0,0) (0.5,0), (1,0)
  END SUBROUTINE QuadratureOnEdge
  !============================================!


END MODULE Quadrature
