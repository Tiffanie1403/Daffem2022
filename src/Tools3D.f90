MODULE Tools3D

  USE Types

  IMPLICIT NONE

CONTAINS

  !===============================!
  SUBROUTINE CalcNormalFace3D(edg, Mesh)
  !===============================!

  !*************************************************************************
  !   Calculates the interior normal vectors of the faces of the element
  !   3d VERSION
  !  Parameters:
  !
  !    Input,Output, type(Element) Ele, Current element
  !
  !*************************************************************************


    IMPLICIT NONE
    type(MeshType), INTENT(IN)    :: Mesh
    type(edge)    , INTENT(INOUT) :: edg
    !----------------------------------
    REAL*8, DIMENSION(3) :: UU, VV
    !-------------------------------

    UU = Mesh%Vertex(edg%Vertex(2))%X - Mesh%Vertex(edg%Vertex(1))%X
    VV = Mesh%Vertex(edg%VerteX(3))%X - Mesh%Vertex(edg%Vertex(1))%X
    edg%N(1) = UU(2)*VV(3) - UU(3)*VV(2)
    edg%N(2) = UU(3)*VV(1) - UU(1)*VV(3)
    edg%N(3) = UU(1)*VV(2) - UU(2)*VV(1)

  END SUBROUTINE CalcNormalFace3D
  !===============================!

  !===========================!
  SUBROUTINE Calc3DNormal(Ele)
  !===========================!

  !*************************************************************************
  !   Calculates the interior normal vectors of the edges of the element
  !   3d VERSION
  !  Parameters:
  !
  !    Input,Output, type(Element) Ele, Current element
  !
  !*************************************************************************


    IMPLICIT NONE
    type(Element), INTENT(INOUT) :: Ele
    !-------------------------------------
    REAL*8, DIMENSION(3) :: UU, VV
    !----------------------------

    UU = Ele%Coor(3,:) - Ele%Coor(2,:)
    VV = Ele%Coor(4,:) - Ele%Coor(2,:)

    Ele%N(1,1) = UU(2)*VV(3) - UU(3)*VV(2)
    Ele%N(1,2) = UU(3)*VV(1) - UU(1)*VV(3)
    Ele%N(1,3) = UU(1)*VV(2) - UU(2)*VV(1)

    UU = Ele%Coor(4,:) - Ele%Coor(1,:)
    VV = Ele%Coor(3,:) - Ele%Coor(1,:)

    Ele%N(2,1) = UU(2)*VV(3) - UU(3)*VV(2)
    Ele%N(2,2) = UU(3)*VV(1) - UU(1)*VV(3)
    Ele%N(2,3) = UU(1)*VV(2) - UU(2)*VV(1)

    UU = Ele%Coor(1,:) - Ele%Coor(4,:)
    VV = Ele%Coor(2,:) - Ele%Coor(4,:)

    Ele%N(3,1) = UU(2)*VV(3) - UU(3)*VV(2)
    Ele%N(3,2) = UU(3)*VV(1) - UU(1)*VV(3)
    Ele%N(3,3) = UU(1)*VV(2) - UU(2)*VV(1)

    UU = Ele%Coor(1,:) - Ele%Coor(2,:)
    VV = Ele%Coor(3,:) - Ele%Coor(2,:)

    Ele%N(4,1) = UU(2)*VV(3) - UU(3)*VV(2)
    Ele%N(4,2) = UU(3)*VV(1) - UU(1)*VV(3)
    Ele%N(4,3) = UU(1)*VV(2) - UU(2)*VV(1)

    Ele%N = - 0.5*Ele%N

  END SUBROUTINE Calc3DNormal
  !===========================!

  !=========================!
  SUBROUTINE Calc3DArea(Ele)
  !=========================!

  !*************************************************************************
  !   Calculates the area of the element
  !   3d VERSION
  !  Parameters:
  !
  !    Input,Output, type(Element) Ele, Current element
  !
  !*************************************************************************


    IMPLICIT NONE
    type(Element), INTENT(INOUT) :: Ele
    !-----------------------------------
    REAL*8, DIMENSION(3,3) :: A1
    !--------------------------
    REAL*8 :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
    REAL*8 :: det
    !-------------------------------------------------------
    REAL*8 :: a, b, c, d, e, f
    REAL*8 :: D2, E2, F2
    REAL*8 :: P, Q, R, area
    !---------------------------

    x1 = ele%Coor(1,1) ; y1 = ele%Coor(1,2) ; z1 = ele%Coor(1,3)
    x2 = ele%Coor(2,1) ; y2 = ele%Coor(2,2) ; z2 = ele%Coor(2,3)
    x3 = ele%Coor(3,1) ; y3 = ele%Coor(3,2) ; z3 = ele%Coor(3,3)
    x4 = ele%Coor(4,1) ; y4 = ele%Coor(4,2) ; z4 = ele%Coor(4,3)

    A1(1,1) = x2 - x1 ; A1(1,2) = x3 - x1 ; A1(1,3) = x4 - x1
    A1(2,1) = y2 - y1 ; A1(2,2) = y3 - y1 ; A1(2,3) = y4 - y1
    A1(3,1) = z2 - z1 ; A1(3,2) = z3 - z1 ; A1(3,3) = z4 - z1

    det = A1(1,1)*A1(2,2)*A1(3,3) + A1(2,1)*A1(3,2)*A1(1,3) + A1(3,1)*A1(1,2)*A1(2,3) &
         - A1(3,1)*A1(2,2)*A1(1,3) - A1(2,1)*A1(1,2)*A1(3,3) - A1(1,1)*A1(3,2)*A1(2,3)

    ele%A = ABS(det)/6.


  END SUBROUTINE Calc3DArea
  !=========================!

  !=========================!
  SUBROUTINE CheckAreas(ele)
  !=========================!
  !*************************************************************************
  !   ??????
  !  Parameters:
  !
  !    Input,Output, type(Element) Ele, Current element
  !
  !*************************************************************************


    IMPLICIT NONE
    type(element), INTENT(IN) :: ele
    !----------------------------------
    REAL*8 :: x1, x2, x3, y1, y2, y3, z1, z2, z3
    REAL*8 :: Nrm, a, b, c, p, L
    !--------------------------------------------
    REAL*8, DIMENSION(3) :: N
    !---------------------------
    INTEGER, DIMENSION(4,3) :: Permut
    !----------------------------------
    INTEGER :: i, n1, n2, n3
    !----------------------------

    Permut(1,1) = 2 ; Permut(1,2) = 3 ; Permut(1,3) = 4
    Permut(2,1) = 1 ; Permut(2,2) = 3 ; Permut(2,3) = 4
    Permut(3,1) = 1 ; Permut(3,2) = 2 ; Permut(3,3) = 4
    Permut(4,1) = 1 ; Permut(4,2) = 2 ; Permut(4,3) = 3

    DO i = 1, 4
       n1 = Permut(i,1)
       n2 = Permut(i,2)
       n3 = Permut(i,3)
       x1 = ele%Coor(n1,1) ; y1 = ele%Coor(n1,2) ; z1 = ele%Coor(n1,3)
       x2 = ele%Coor(n2,1) ; y2 = ele%Coor(n2,2) ; z2 = ele%Coor(n2,3)
       x3 = ele%Coor(n3,1) ; y3 = ele%Coor(n3,2) ; z3 = ele%Coor(n3,3)

       a = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
       b = sqrt((x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2)
       c = sqrt((X2-x3)**2 + (y2-y3)**2 + (z2-z3)**2)
       p = 0.5d0*(a+b+c)
       L = sqrt(p*(p-a)*(p-b)*(p-c))

       N = ele%N(i,:)
       Nrm = sqrt(N(1)**2 + N(2)**2 + N(3)**2)
    END DO

  END SUBROUTINE CheckAreas
  !=========================!

  !==========================!
  SUBROUTINE CheckEdges(Mesh)
  !==========================!

  !*************************************************************************
  !   ??????
  !  Parameters:
  !
  !    Input,Output, type(Element) Ele, Current element
  !
  !*************************************************************************


    IMPLICIT NONE
    type(MeshType), INTENT(IN) :: Mesh
    !------------------------------------
    type(element) :: ele
    !----------------------
    INTEGER :: Ne, Ned, ie, AdjNum, it, i
    !-------------------------------------

    Ne = Mesh%Ne ; Ned = Mesh%Ned
    it = 0
    DO ie = 1, Ne
       ele = Mesh%Ele(ie)
       DO i = 1, 4
          AdjNum = ele%Adj(i)/4
          IF ( AdjNum > ie .OR. AdjNum == 0 ) it = it + 1
       END DO
    END DO

  END SUBROUTINE CheckEdges
  !=========================!

END MODULE Tools3D
