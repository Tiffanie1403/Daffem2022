MODULE BoxDistance

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !====================================================!
  SUBROUTINE InnerBox(PhysicalInterface, Sol, Mesh, Data, x_center, x0, y0, iLS)
  !====================================================!

    IMPLICIT NONE
    
    type(VectorInterface), INTENT(IN)    :: PhysicalInterface 
    type(SolStructure)   , INTENT(INOUT) :: Sol
    type(MeshType)       , INTENT(INOUT) :: Mesh
    type(DataStructure)  , INTENT(IN)    :: Data
    REAL*8, DIMENSION(:) , INTENT(IN)    :: x_center
    REAL*8               , INTENT(IN)    :: x0, y0
    INTEGER              , INTENT(IN)    :: iLS
    !---------------------------------------------------
    type(element)		  :: ele
    !---------------------------------------------------
    REAL*8, DIMENSION(3)  :: LocDist
    !---------------------------------------------------
    INTEGER, DIMENSION(3) :: NU, Pos
    !---------------------------------------------------
    REAL*8 				  :: x, y, r, nx, ny, Nrm, dist
    REAL*8 				  :: min_dist, x_c, y_c
    !---------------------------------------------------
    INTEGER 			  :: ie, Ne, i, cpt, Nv
    !---------------------------------------------------

    Nv = Mesh%Nv ; Ne = Mesh%Ne
    x_c = x_center(1) ; y_c = x_center (2) ! typiquement (0,0)
    DO ie = 1, Ne

       ele = Mesh%Ele(ie)
       NU = ele%Vertex !vecteur conteneant les noeud de l'element  
       Pos = ele%Pos
       
       DO i = 1, 3
        IF(Data%DefInterfaceNode) THEN 
		  CALL OrthoProjectionDistance(PhysicalInterface, ele%Coor(i,:), Sol%Dist(iLS,NU(i),:), Sol%LS(iLS,NU(i)))
        ELSE 
          CALL NodalInnerBoxDist(ele%Coor(i,:), x_c, y_c, x0, y0, Sol%Dist(iLS,NU(i),:), Sol%LS(iLS,NU(i)))
		  !remplit pour les trois noeuds de l'elements Sol%dist qui ets le vecteur distance non normalisé, et Sol%Ls qui est la distance signé
		END IF 
       END DO

       min_dist = 10000.
       DO i = 1, 3
          dist = Sol%LS(iLS,NU(i)) !distance signé 
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

  END SUBROUTINE InnerBox
  !====================================!

  !=============================================!
  SUBROUTINE OuterBox(Sol, Mesh, x_center, x0, y0, iLS)
  !=============================================!

    IMPLICIT NONE
    type(SolStructure)  , INTENT(INOUT) :: Sol
    type(MeshType)      , INTENT(INOUT) :: Mesh
    REAL*8, DIMENSION(:), INTENT(IN)    :: x_center
    REAL*8              , INTENT(IN)    :: x0, y0
    INTEGER             , INTENT(IN)    :: iLS
    !--------------------------------------------
    type(element) :: ele
    !--------------------
    REAL*8, DIMENSION(3) :: LocDist
    !-------------------------------
    INTEGER, DIMENSION(3) :: NU, Pos
    !--------------------------------
    REAL*8 :: x, y, r, nx, ny, Nrm, dist
    REAL*8 :: min_dist, xc, yc
    !------------------------------------
    INTEGER :: ie, i, cpt, Ne
    !------------------------------

    Ne = Mesh%Ne
    xc = x_center(1) ; yc = x_center(2)
    DO ie = 1, Ne

       ele = Mesh%Ele(ie)
       NU = ele%Vertex ; Pos = ele%Pos !Nu numerotation locale
       DO i = 1, 3
          CALL NodalOuterBoxDist(ele%Coor(i,:), xc, yc, x0, y0, Sol%Dist(iLS,NU(i),1:2), Sol%LS(iLS,NU(i)))
          !PRINT*, iLS, NU(i), ele%coor(i,:), sol%dist(1,NU(i),2)
       END DO

       min_dist = 10000.
       DO i = 1, 3
          dist = Sol%LS(iLS,NU(i)) !on recupere la distance calcule par NodalOuterBoxDist
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

    !PRINT*, sol%dist(1,:,2)

  END SUBROUTINE OuterBox
  !====================================!

  !======================================================!
  SUBROUTINE NodalInnerBoxDist(Xi, xc, yc, x0, y0, d, LS)
  !======================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: Xi
    REAL*8, DIMENSION(2), INTENT(OUT) :: d
    REAL*8              , INTENT(IN)  :: x0, y0, xc, yc
    REAL*8              , INTENT(OUT) :: LS
    !-------------------------------------------
    REAL*8 :: r, dist, nx, ny, Nrm, x, y
    REAL*8 :: x1, y1
    !------------------------------------

    CALL NodalBoxDist(Xi, xc, yc, x0, y0, d, LS)

  END SUBROUTINE NodalInnerBoxDist
  !================================================!

  !==============================================!
  SUBROUTINE NodalOuterBoxDist(Xi, xc, yc, x0, y0, d, LS)
  !==============================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)     :: Xi
    REAL*8              , INTENT(IN)     :: xc, yc, x0, y0
    REAL*8, DIMENSION(2), INTENT(OUT)    :: d
    REAL*8              , INTENT(OUT)    :: LS
    !-----------------------NodalOuter----------------------
    REAL*8 :: r, dist, nx, ny, Nrm, x, y, x1, y1
    !--------------------------------------------

    CALL NodalBoxDist(Xi, xc, yc, x0, y0, d, LS)

  END SUBROUTINE NodalOuterBoxDist
  !================================================!

  !==============================================!
  SUBROUTINE NodalBoxDist(Xi, xc, yc, x0, y0, d, LS)
  !==============================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)     :: Xi
    REAL*8              , INTENT(IN)     :: xc, yc, x0, y0
    REAL*8, DIMENSION(2), INTENT(OUT)    :: d
    REAL*8              , INTENT(OUT)    :: LS

    !-----------------------NodalOuter----------------------
    REAL*8 :: r, dist, nx, ny, Nrm, x, y, x1, y1, x2, y2, r0
    REAL*8 :: xc2, yc2, xc3, yc3
    !----------------------------------------------------
    
    x = Xi(1) ; y = Xi(2)
    
    !avant modif 
    x1 = xc + x0/2 ; y1 = yc + y0/2
    x2 = xc - x0/2 ; y2 = yc - y0/2

    IF ( x > x1 ) THEN

       IF ( ( y < y1 ) .AND. ( y > y2 ) ) THEN
          dist = x - x1
          nx = -1.d0 ; ny = 0.d0

       ELSE IF ( y >= y1 ) THEN
          r = sqrt((x-x1)**2 + (y-y1)**2)
          dist = r
          nx = -(x - x1)/r ; ny = -(y-y1)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       ELSE

          r = sqrt((x-x1)**2 + (y-y2)**2)
          dist = r
          nx = -(x-x1)/r ; ny = -(y-y2)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       END IF

    ELSE IF ( x < x2 ) THEN

       IF ( ( y < y1 ) .AND. ( y > y2 ) ) THEN
          dist = -(x - x2)
          nx = 1.d0 ; ny = 0.d0

       ELSE IF ( y >= y1 ) THEN
          r = sqrt((x-x2)**2 + (y-y1)**2)
          dist = r
          nx = - (x-x2)/r ; ny = -(y - y1)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       ELSE
          r = sqrt((x-x2)**2 + (y-y2)**2)
          dist = r
          nx = - (x-x2)/r ; ny = -(y-y2)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       END IF

    ELSE

       IF ( y >= y1 ) THEN

          dist = y - y1
          nx = 0.d0 ; ny = -1.d0

       ELSE IF ( y <= y2 ) THEN

          dist = -(y - y2)
          nx = 0.d0 ; ny = 1.d0

       ELSE

          dist = MIN(x1-x, x-x2)
          dist = MIN(dist, y1-y)
          dist = MIN(dist, y-y2)
          IF ( abs(dist - (x1-x)) < 1.d-10 ) THEN
             nx = 1.d0 ; ny = 0.d0
          ELSE IF ( abs(dist - (x-x2)) < 1.d-10 ) THEN
             nx = -1.d0 ; ny = 0.d0
          ELSE IF ( abs(dist - (y1-y)) < 1.d-10 ) THEN
             nx = 0.d0 ; ny = 1.d0
          ELSE
             nx = 0.d0 ;
             ny = -1.d0
          END IF
          dist = -dist
       END IF
    END IF

    d(1) = abs(dist)*nx ; d(2) = abs(dist)*ny ; LS = dist
    
  END SUBROUTINE NodalBoxDist
  !================================================!

  !==============================================!
  SUBROUTINE NodalBoxDist2(Xi, xc, yc, x0, y0, d, LS)
  !==============================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)     :: Xi
    REAL*8              , INTENT(IN)     :: xc, yc, x0, y0
    REAL*8, DIMENSION(2), INTENT(OUT)    :: d
    REAL*8              , INTENT(OUT)    :: LS

    !-----------------------NodalOuter----------------------
    REAL*8 :: r, dist, nx, ny, Nrm, x, y, x1, y1, x2, y2, r0
    REAL*8 :: xc2, yc2, xc3, yc3
    !----------------------------------------------------

    x = Xi(1) ; y = Xi(2)
    x1 = xc + x0 ; y1 = yc + y0/2.d0
    x2 = xc - x0 ; y2 = yc - y0/2.d0


    IF ( x > x1 ) THEN

       IF ( ( y < y1 ) .AND. ( y > y2 ) ) THEN
          dist = x - x1
          nx = -1.d0 ; ny = 0.d0

       ELSE IF ( y >= y1 ) THEN
          r = sqrt((x-x1)**2 + (y-y1)**2)
          dist = r
          nx = -(x - x1)/r ; ny = -(y-y1)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       ELSE

          r = sqrt((x-x1)**2 + (y-y2)**2)
          dist = r
          nx = -(x-x1)/r ; ny = -(y-y2)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       END IF

    ELSE IF ( x < x2 ) THEN

       IF ( ( y < y1 ) .AND. ( y > y2 ) ) THEN
          dist = -(x - x2)
          nx = 1.d0 ; ny = 0.d0

       ELSE IF ( y >= y1 ) THEN
          r = sqrt((x-x2)**2 + (y-y1)**2)
          dist = r
          nx = - (x-x2)/r ; ny = -(y - y1)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       ELSE
          r = sqrt((x-x2)**2 + (y-y2)**2)
          dist = r
          nx = - (x-x2)/r ; ny = -(y-y2)/r
          nrm = sqrt(nx*nx + ny*ny)
          nx = nx/nrm ; ny = ny/nrm

       END IF

    ELSE

       IF ( y >= y1 ) THEN

          dist = y - y1
          nx = 0.d0 ; ny = -1.d0

       ELSE IF ( y <= y2 ) THEN

          dist = -(y - y2)
          nx = 0.d0 ; ny = 1.d0

       ELSE

          dist = MIN(x1-x, x-x2)
          dist = MIN(dist, y1-y)
          dist = MIN(dist, y-y2)
          IF ( abs(dist - (x1-x)) < 1.d-10 ) THEN
             nx = 1.d0 ; ny = 0.d0
          ELSE IF ( abs(dist - (x-x2)) < 1.d-10 ) THEN
             nx = -1.d0 ; ny = 0.d0
          ELSE IF ( abs(dist - (y1-y)) < 1.d-10 ) THEN
             nx = 0.d0 ; ny = 1.d0
          ELSE
             nx = 0.d0 ;
             ny = -1.d0
          END IF
          dist = -dist
       END IF
    END IF

    d(1) = abs(dist)*nx ; d(2) = abs(dist)*ny ; LS = dist
  END SUBROUTINE NodalBoxDist2
  !================================================!
  
  !==============================================!
  SUBROUTINE OrthoProjectionDistance(PhysicalInterface, Xi, d, LS)
  !==============================================!
	!function define only if the surrogate is define before the physical interface, on the left side 
    IMPLICIT NONE
    
    type(VectorInterface), INTENT(IN)  :: PhysicalInterface 
    REAL*8, DIMENSION(:), INTENT(IN)   :: Xi
    REAL*8, DIMENSION(2), INTENT(OUT)  :: d
    REAL*8              , INTENT(OUT)  :: LS

    !-----------------------NodalOuter----------------------
    REAL*8 								:: dist, nx, ny, Nrm, x, y, x1, y1, x2, y2, minDis, x_w, y_w 
    LOGICAL 							:: TEST
    INTEGER								:: i , NumDis
    REAL*8, DIMENSION(2)				:: u, v, w, e
    REAL*8								:: MinDist, xw, yw , cur_dist
    !----------------------------------------------------

    x = Xi(1) ; y = Xi(2)
    !what I need first is to define x1,y1
    !We do the orthogonal projection on every segments defining the Physical interface 
    MinDist = -1.d0 
    xw = 0 ; yw = 0 
    DO i = 1,Icing_NodePhysicalInterface-1
		x1 = PhysicalInterface%Coor(i,1) ; x2 = PhysicalInterface%Coor(i+1,1) 
		y1 = PhysicalInterface%Coor(i,2) ; y2 = PhysicalInterface%Coor(i+1,2) 
		!coordonnées associé à une edge de la vraie interface 
		
		u(1) = x  - x1 ; u(2) = y  - y1 ! projected vector 
		v(1) = x2 - x1 ; v(2) = y2 - y1 ! domain receiving the projection  

		!results of the projection 
		w(1)  = (v(1)*u(1)+v(2)*u(2))/(v(1)**2+v(2)**2)*v(1) 
		w(2)  = (v(1)*u(1)+v(2)*u(2))/(v(1)**2+v(2)**2)*v(2) 
		e(1)  = u(1) - w(1) ; e(2)  = u(2) - w(2)            !projected vector 
		
		x_w = x-e(1) ; y_w = y-e(2) ! coordinate of the the node from the orthogonale projection  
		
		!we now check if the node is on the edge or not  
		TEST = .FALSE. 
		CALL CheckPointOnEdge(TEST,x_w,y_w,x1,x2,y1,y2) 
		IF(TEST) THEN 
			cur_dist = sqrt((x-x_w)**2+(y-y_w)**2)
		
			IF(MinDist<0) THEN !first iteration 
				MinDist = cur_dist
				xw = x_w 
				yw = y_w
			ELSE IF (MinDist > cur_dist) THEN 
				MinDist = cur_dist
				xw = x_w 
				yw = y_w
			END IF 
		END IF 
    END DO 
    
    !PRINT*,"Value of the projection", x, y, xw, yw, MinDist
    !the projection has been done on every edges of the physical interface 
    !at this point I have the distance and the point on the physical interface 
    
    !I have now to define the sign distance and the normal vector 
    IF((xw-x)>0) THEN 
		dist =  -MinDist
		nx   = (xw-x)/abs(dist) ; ny = (yw-y)/abs(dist) 
    ELSE IF ((xw-x)<= 0) THEN	
		dist = MinDist
		nx   =  -(xw-x)/abs(dist) ; ny = -(yw-y)/abs(dist)
	ELSE 
		PRINT*,"The points are on the same abscisse, it shouldn't"
		STOP 
    END IF
	
    d(1) = abs(dist)*nx ; d(2) = abs(dist)*ny 
    LS = dist !signed distance 
    
    END SUBROUTINE OrthoProjectionDistance
    !================================================!

   !================================================!
	SUBROUTINE CheckPointOnEdge(TEST,x_w,y_w,x1,x2,y1,y2)
   !================================================!
	! we check if a point is between two points (x1,y1) et (x2,y2) 
	! we already do the assumption that (x_w,y_w) is on the line defined by the two points above 
	IMPLICIT NONE
    
    LOGICAL, INTENT(INOUT) :: TEST 
    REAL*8 , INTENT(IN)    :: x_w, y_w, x1, x2, y1, y2
    !----------------------------------------------------
    REAL*8                 :: K1, K2
    !----------------------------------------------------
    
    K1 = (x2-x1)*(x_w-x1) + (y2-y1)*(y_w-y1)
    K2 = (x2-x1)*(x2-x1)  + (y2-y1)*(y2-y1)

	IF(K1 < 0 .OR. K1 > K2 ) THEN 
		TEST = .FALSE. 
	ELSE 
		TEST = .TRUE.
	END IF    
	
	END SUBROUTINE CheckPointOnEdge
    !================================================!


END MODULE BoxDistance
