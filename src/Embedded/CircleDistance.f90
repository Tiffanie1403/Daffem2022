MODULE CircleDistance

  USE Types

  IMPLICIT NONE

CONTAINS

  !===================================================!
  SUBROUTINE InnerCircle(Sol, Mesh, x_center, r0, iLS)
  !===================================================!

    IMPLICIT NONE
    type(SolStructure)  , INTENT(INOUT) :: Sol
    type(MeshType)      , INTENT(INOUT) :: Mesh
    REAL*8, DIMENSION(:), INTENT(IN)    :: x_center
    REAL*8              , INTENT(IN)    :: r0
    INTEGER             , INTENT(IN)    :: iLS
    !---------------------------------------------------
    type(element) :: ele
    !--------------------
    REAL*8, DIMENSION(3) :: LocDist
    !-------------------------------
    INTEGER, DIMENSION(3) :: NU, Pos
    !--------------------------------
    REAL*8 :: x, y, r, nx, ny, Nrm, dist, x0, y0
    REAL*8 :: min_dist
    !---------------------------------------------
    INTEGER :: ie, Ne, i, cpt, Nv
    !------------------------------

    Nv = Mesh%Nv ; Ne = Mesh%Ne
    x0 = x_center(1) ; y0 = x_center(2)

    DO ie = 1, Ne

       ele = Mesh%Ele(ie)
       NU = ele%Vertex ; Pos = ele%Pos
       DO i = 1, 3
          CALL NodalInnerCircleDist(ele%Coor(i,:), x0, y0, r0, Sol%Dist(iLS,NU(i),:), Sol%LS(iLS,NU(i)))
       END DO

       min_dist = 10000.
       DO i = 1, 3
          dist = Sol%LS(iLS,NU(i))
          LocDist(i) = Sol%LS(iLS,NU(i))
          min_dist = MIN(min_dist,dist)
          IF ( dist > 0.d0 ) THEN
             Sol%NodeStatus(iLS,ele%Pos(i)) = 1
             Mesh%Vertex(ele%Vertex(i))%Status(iLS) = 1
          ELSE
             Sol%NodeStatus(iLS,ele%Pos(i)) = 0
             Mesh%Vertex(ele%Vertex(i))%Status(iLS) = 0
          END IF
       END DO

       Cpt = SUM(Sol%NodeStatus(iLS,ele%Pos))
       IF ( cpt == 1 ) THEN
          Mesh%Ele(ie)%Status(iLS) = 0
       ELSE IF ( cpt == 2 ) THEN
          Sol%NodeStatus(iLS, ele%Pos) = 1
          Mesh%Ele(ie)%Status(iLS) = 1
       ELSE IF ( cpt == 0 ) THEN
          Mesh%Ele(ie)%Status(iLS) = -1
       ELSE
          Mesh%Ele(ie)%Status(iLS) = 2
       END IF

       Sol%EleStatus(iLS,ele%Pos) = Mesh%Ele(ie)%Status(iLS)

    END DO

  END SUBROUTINE InnerCircle
  !====================================!

  !===================================================!
  SUBROUTINE OuterCircle(Sol, Mesh, x_center, r0, iLS)
  !===================================================!

    IMPLICIT NONE
    type(SolStructure)  , INTENT(INOUT) :: Sol
    type(MeshType)      , INTENT(INOUT) :: Mesh
    REAL*8, DIMENSION(:), INTENT(IN)    :: x_center
    REAL*8              , INTENT(IN)    :: r0
    INTEGER             , INTENT(IN)    :: iLS
    !--------------------------------------------
    type(element) :: ele
    !-------------------------------
    INTEGER, DIMENSION(3) :: NU, Pos
    !--------------------------------
    REAL*8 :: x, y, r, nx, ny, Nrm, dist
    REAL*8 :: min_dist, x0, y0
    !------------------------------------
    INTEGER :: ie, i, cpt, Ne
    !------------------------------

    Ne = Mesh%Ne
    x0 = x_center(1) ; y0 = x_center(2)
    DO ie = 1, Ne

       ele = Mesh%Ele(ie)
       NU = ele%Vertex ; Pos = ele%Pos
       DO i = 1, 3
          CALL NodalOuterCircleDist(ele%Coor(i,:), x0, y0, r0, Sol%Dist(iLS,NU(i),:), Sol%LS(iLS,NU(i)))
       END DO

       min_dist = 10000.
       DO i = 1, 3
          dist = Sol%LS(iLS,NU(i))
          min_dist = MIN(min_dist,dist)
          IF ( dist > 0. ) THEN
             Sol%NodeStatus(iLS,ele%Pos(i)) = 1
             Mesh%Vertex(ele%Vertex(i))%Status(iLS) = 1
          ELSE
             Sol%NodeStatus(iLS,ele%Pos(i)) = 0
             Mesh%Vertex(ele%Vertex(i))%Status(iLS) = 0
          END IF
       END DO

       cpt = SUM(Sol%NodeStatus(iLS,ele%Pos))
       IF ( cpt == 2 ) THEN
          Mesh%Ele(ie)%Status(iLS) = 0
       ELSE IF ( cpt == 1 ) THEN
          Sol%NodeStatus(iLS,ele%Pos) = 1
          Mesh%Ele(ie)%Status(iLS) = 1
       ELSE IF ( cpt == 3 ) THEN
          Mesh%Ele(ie)%Status(iLS) = -1
       ELSE
          Mesh%Ele(ie)%Status(iLS) = 2
       END IF

       Sol%EleStatus(iLS,ele%Pos) = Mesh%Ele(ie)%Status(iLS)

    END DO

  END SUBROUTINE OuterCircle
  !====================================!


  !=====================================================!
  SUBROUTINE NodalInnerCircleDist(Xi, x0, y0, r0, d, LS)
  !=====================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: Xi
    REAL*8, DIMENSION(2), INTENT(OUT) :: d
    REAL*8              , INTENT(IN)  :: x0, y0, r0
    REAL*8              , INTENT(OUT) :: LS
    !----------------------------------------
    REAL*8 :: r, dist, nx, ny, Nrm, x, y
    !------------------------------------

    x = Xi(1) ; y = Xi(2)
    r = sqrt((x-x0)**2 + (y-y0)**2)

    dist = r - r0
    nx = -(x-x0)/r ; ny = -(y-y0)/r
    nrm = sqrt(nx*nx + ny*ny)
    nx = nx/nrm ; ny = ny/nrm

    d(1) = dist*nx ; d(2) = dist*ny ; LS = dist

  END SUBROUTINE NodalInnerCircleDist
  !================================================!

  !=====================================================!
  SUBROUTINE NodalOuterCircleDist(Xi, x0, y0, r0, d, LS)
  !=====================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: Xi
    REAL*8, DIMENSION(2), INTENT(OUT) :: d
    REAL*8              , INTENT(IN)  :: x0, y0, r0
    REAL*8              , INTENT(OUT) :: LS
    !----------------------------------------
    REAL*8 :: r, dist, nx, ny, Nrm, x, y
    !------------------------------------

    x = Xi(1) ; y = Xi(2)
    r = sqrt((x-x0)**2 + (y-y0)**2)

    dist = r - r0
    nx = -x/r ; ny = -y/r
    nrm = sqrt(nx*nx + ny*ny)
    nx = nx/nrm ; ny = ny/nrm

    d(1) = dist*nx
    d(2) = dist*ny
    LS = dist

  END SUBROUTINE NodalOuterCircleDist
  !============================================!

  !=====================================================!
  SUBROUTINE NodalInnerCircleDist2(Xi, x0, y0, r0, d, LS,n)
  !=====================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: Xi
    REAL*8, DIMENSION(2), INTENT(OUT) :: d
    REAL*8              , INTENT(IN)  :: x0, y0, r0
    REAL*8              , INTENT(OUT) :: LS
    REAL*8, DIMENSION(2), INTENT(OUT) :: n
    !----------------------------------------
    REAL*8 :: r, dist, nx, ny, Nrm, x, y
    REAL*8 :: x1,y1,x2,y2,d1,d2
    !------------------------------------

    x = Xi(1) ; y = Xi(2)
    r = sqrt((x-x0)**2 + (y-y0)**2)
    dist = r - r0
    x1=(x*r0)/r ; y1=(y*r0)/r
    x2=-(x*r0)/r ; y2=-(y*r0)/r
    d1 = sqrt((x1-x)**2 + (y1-y)**2)
    d2 = sqrt((x2-x)**2 + (y2-y)**2)
    IF (d1 < d2) THEN
      d(1)=x1-x ; d(2)=y1-y
    ELSE
      d(1)=x2-x ; d(2)=y2-y
    END IF
    LS=dist
    n(1)=d(1)/sqrt(d(1)**2+d(2)**2) ; n(2)=d(2)/sqrt(d(1)**2+d(2)**2)

  END SUBROUTINE NodalInnerCircleDist2
  !============================================!

END MODULE CircleDistance
