MODULE Tools2D

  USE Types

  IMPLICIT NONE

CONTAINS

  !===========================!
  SUBROUTINE Calc2DNormal(Ele)
  !===========================!

  !*************************************************************************
  !   Calculates the interior normal vectors of the element
  !   2d VERSION
  !  Parameters:
  !
  !    Input,Output, type(Element) Ele, Current element
  !
  !*************************************************************************


    IMPLICIT NONE
    type(Element), INTENT(INOUT) :: Ele
    !-------------------------------------

    Ele%N(1,1) = Ele%Coor(2,2) - Ele%Coor(3,2)
    Ele%N(1,2) = Ele%Coor(3,1) - Ele%Coor(2,1)

    Ele%N(2,1) = Ele%Coor(3,2) - Ele%Coor(1,2)
    Ele%N(2,2) = Ele%Coor(1,1) - Ele%Coor(3,1)

    Ele%N(3,1) = Ele%Coor(1,2) - Ele%Coor(2,2)
    Ele%N(3,2) = Ele%Coor(2,1) - Ele%Coor(1,1)

  END SUBROUTINE Calc2DNormal
  !===========================!

  !=========================!
  SUBROUTINE Calc2DArea(Ele)
  !=========================!

  !*************************************************************************
  !   Calculates the area of an element
  !
  !  Parameters:
  !
  !    Input,Output, type(Element) Ele, Current element
  !
  !*************************************************************************


    IMPLICIT NONE
    type(Element), INTENT(INOUT) :: Ele
    !-----------------------------------
    REAL*8 :: a, b, c, p
    !---------------------

    a = SQRT( (Ele%Coor(1,1) - Ele%Coor(2,1))**2 + (Ele%Coor(1,2) - Ele%Coor(2,2))**2 )
    b = SQRT( (Ele%Coor(1,1) - Ele%Coor(3,1))**2 + (Ele%Coor(1,2) - Ele%Coor(3,2))**2 )
    c = SQRT( (Ele%Coor(3,1) - Ele%Coor(2,1))**2 + (Ele%Coor(3,2) - Ele%Coor(2,2))**2 )
    p = 0.5d0*(a + b + c)
    Ele%A = SQRT(p*(p - a)*(p - b)*(p - c)) !attributes for the area on the element

  END SUBROUTINE Calc2DArea
  !=========================!


END MODULE Tools2D
