MODULE SphereDistance

  USE Types

  IMPLICIT NONE

CONTAINS

  !=========================================!
  SUBROUTINE InnerSphere(Sol, Mesh, r0, iLS)
  !=========================================!

  !*************************************************************************
  !  Handle the calculation of the composantes for the normals and distance vectors
  !  Case of the surrogate inside the true domain
  !  3d case
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Solution structure
  !    Input,Output, REAL*8, DIMENSION(3)  Mesh, Mesh structure from dom.mesh
  !    Input, REAL*8  r0, radius of the cercle from the true interface
  !    Input, INTEGER iLS, ????
  !
  !*************************************************************************


    IMPLICIT NONE
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(MeshType)     , INTENT(INOUT) :: Mesh
    REAL*8             , INTENT(IN)    :: r0
    INTEGER            , INTENT(IN)    :: iLS
    !---------------------------------------------------
    type(element) :: ele
    !--------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: LocDist
    !---------------------------------------------
    INTEGER, DIMENSION(:), ALLOCATABLE :: NU, Pos
    !----------------------------------------------
    REAL*8 :: x, y, r, nx, ny, Nrm, dist
    REAL*8 :: min_dist
    !------------------------------------
    INTEGER :: ie, Ne, i, cpt, Nv, nNodesPerElement
    !------------------------------------------------

    Nv = Mesh%Nv ; Ne = Mesh%Ne

    DO ie = 1, Ne

       ele = Mesh%Ele(ie)
       nNodesPerElement = ele%nNodesPerElement
       ALLOCATE ( LocDist(nNodesPerElement) )
       ALLOCATE ( NU(nNodesPerElement), Pos(nNodesPerElement) )
       NU = ele%Vertex ; Pos = ele%Pos
       DO i = 1, nNodesPerElement
          CALL NodalInnerSphere(ele%Coor(i,:), r0, Sol%Dist(iLS,NU(i),:), Sol%LS(iLS,NU(i)))
          !IF ( ele%Vertex(i) == 1 ) STOP
       END DO

       min_dist = 10000.
       DO i = 1, 4
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

       cpt = SUM(Sol%NodeStatus(iLS,ele%Pos))
!!$       IF ( cpt == 1 ) THEN
!!$          Mesh%Ele(ie)%Status(iLS) = 0
!!$       ELSE IF ( cpt == 2 ) THEN
!!$          Sol%NodeStatus(iLS, ele%Pos) = 1
!!$          Mesh%Ele(ie)%Status(iLS) = 1
!!$       ELSE IF ( cpt == 0 ) THEN
!!$          Mesh%Ele(ie)%Status(iLS) = -1
!!$       ELSE
!!$          Mesh%Ele(ie)%Status(iLS) = 2
!!$       END IF
       IF ( cpt == 0 ) THEN
          mesh%ele(ie)%status(iLS) = -1
       ELSE IF ( cpt == 3 ) THEN
          sol%nodeStatus(iLS,ele%Pos) = 1
          mesh%ele(ie)%status = 1
       ELSE IF ( cpt == 2 .OR. cpt == 1 ) THEN
          Mesh%Ele(ie)%Status(iLS) = 0
       ELSE
          Mesh%Ele(ie)%Status(iLS) = 2
       END IF
       Sol%EleStatus(iLS,ele%Pos) = Mesh%Ele(ie)%Status(iLS)

       DEALLOCATE ( LocDist, NU, Pos )

    END DO

  END SUBROUTINE InnerSphere
  !===============================!

  !=========================================!
  SUBROUTINE OuterSphere(Sol, Mesh, r0, iLS)
  !=========================================!

  !*************************************************************************
  !  Handle the calculation of the composantes for the normals and distance vectors
  !  Case of the surrogate outside the true domain
  !  3d case
  !
  !  Parameters:
  !
  !    Input,Output, type(SolStructure) Sol, Solution structure
  !    Input,Output, REAL*8, DIMENSION(3)  Mesh, Mesh structure from dom.mesh
  !    Input, REAL*8  r0, radius of the cercle from the true interface
  !    Input, INTEGER iLS, ????
  !
  !*************************************************************************

    IMPLICIT NONE
    type(SolStructure) , INTENT(INOUT) :: Sol
    type(MeshType)     , INTENT(INOUT) :: Mesh
    REAL*8             , INTENT(IN)    :: r0
    INTEGER            , INTENT(IN)    :: iLS
    !--------------------------------------------
    type(element) :: ele
    !--------------------
    REAL*8, DIMENSION(4) :: LocDist
    !-------------------------------
    INTEGER, DIMENSION(4) :: NU, Pos
    !--------------------------------
    REAL*8 :: x, y, r, nx, ny, Nrm, dist
    REAL*8 :: min_dist
    !------------------------------------
    INTEGER :: ie, i, cpt, Nv, Ne
    !------------------------------

    Ne = Mesh%Ne ; Nv = Mesh%Nv

    DO ie = 1, Ne

       ele = Mesh%Ele(ie)
       NU = ele%Vertex ; Pos = ele%Pos
       DO i = 1, 4
          CALL NodalOuterSphere(ele%Coor(i,:), r0, Sol%Dist(iLS,NU(i),:), Sol%LS(iLS,NU(i)))
       END DO

       min_dist = 10000.
       DO i = 1, 4
          dist = Sol%LS(iLS,NU(i))
          LocDist(i) = Sol%LS(iLS,NU(i))
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
              IF ( cpt == 0 ) THEN
          mesh%ele(ie)%status(iLS) = -1
       ELSE IF ( cpt == 3 ) THEN
          sol%nodeStatus(iLS,ele%Pos) = 1
          mesh%ele(ie)%status = 1
       ELSE IF ( cpt == 2 .OR. cpt == 1 ) THEN
          Mesh%Ele(ie)%Status(iLS) = 0
       ELSE
          Mesh%Ele(ie)%Status(iLS) = 2
       END IF
       Sol%EleStatus(iLS,ele%Pos) = Mesh%Ele(ie)%Status(iLS)

    END DO

  END SUBROUTINE OuterSphere
  !================================!


  !===============================================!
  SUBROUTINE NodalInnerSphere(Xi, r0, d, LS)
  !===============================================!

  !*************************************************************************
  !  Defines the distance vector in the case of a surrogate chosen
  !  inside the cercle difined by the true interface
  !  3d case
  !
  !  Parameters:
  !
  !    Input,REAL*8, DIMENSION(:) Xi, Current point
  !    Output, REAL*8, DIMENSION(3) d, distance vector
  !    Input, REAL*8  r0, radius of the cercle from the true interface
  !    Input, REAL*8  LS, distance value
  !
  !*************************************************************************

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: Xi
    REAL*8, DIMENSION(3), INTENT(OUT) :: d
    REAL*8              , INTENT(IN)  :: r0
    REAL*8              , INTENT(OUT) :: LS
    !----------------------------------------
    REAL*8 :: r, dist, nx, ny, nz, Nrm, x, y, z
    !-------------------------------------------

    x = Xi(1) ; y = Xi(2) ; z = Xi(3)
    r = sqrt(x*x + y*y + z*z)
    dist = r - r0
    nx = -x/r ; ny = -y/r ; nz = -z/r
    nrm = sqrt(nx*nx + ny*ny + nz*nz)
    nx = nx/nrm ; ny = ny/nrm ; nz = nz/nrm
    d(1) = dist*nx ; d(2) = dist*ny ; d(3) = dist*nz ; LS = dist

  END SUBROUTINE NodalInnerSphere
  !================================================!

  !=============================================!
  SUBROUTINE NodalOuterSphere(Xi, r0, d, LS)
  !=============================================!

  !*************************************************************************
  !  Defines the distance vector in the case of a surrogate chosen
  !  outside the cercle difined by the true interface
  !  3d case
  !
  !  Parameters:
  !
  !    Input,REAL*8, DIMENSION(:) Xi, Current point
  !    Output, REAL*8, DIMENSION(3) d, distance vector
  !    Input, REAL*8  r0, radius of the cercle from the true interface
  !    Input, REAL*8  LS, distance value (with an opposite sign)
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: Xi
    REAL*8, DIMENSION(3), INTENT(OUT) :: d
    REAL*8              , INTENT(IN)  :: r0
    REAL*8              , INTENT(OUT) :: LS
    !----------------------------------------
    REAL*8 :: r, dist, nx, ny, nz, Nrm, x, y, z
    !-------------------------------------------

    x = Xi(1) ; y = Xi(2) ; z = Xi(3)
    r = sqrt(x*x + y*y + z*z)

    dist = r - r0
    nx = -x/r ; ny = -y/r ; nz = -z/r
    nrm = sqrt(nx*nx + ny*ny + nz*nz)
    nx = nx/nrm ; ny = ny/nrm ; nz = nz/nrm

    d(1) = -dist*nx ; d(2) = -dist*ny ; d(3) = -dist*nz ; LS = -dist

  END SUBROUTINE NodalOuterSphere
  !============================================!

END MODULE SphereDistance
