MODULE ExactFunctionJump

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !==============================================!
  SUBROUTINE ExactSolJump(L, Exa, Xin, ele, nDim)
  !==============================================!

  !*************************************************************************
  !   Gives the exact value of the solutions on the element ele
  !   at the point Xin, subroutine for non time dependent problem
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(:,:) L, permeability matrix associated to the current element
  !    Output, REAL*8, DIMENSION(:) exa, analytical solution at the point Xin
  !    Input, REAL*8, DIMENSION(:) Xin, Current point on the element ele
  !    Input, type(Element) ele, Current element
  !    Input, INTEGER  nDim, Dimension of the problem (2 or 3)
  !

  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN)  :: L
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: exa
    type(Element)         , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:)  , INTENT(IN)  :: xin
    INTEGER               , INTENT(IN)  :: nDim
    !--------------------------------------------
    REAL*8 :: x, y, z, lambda1, lambda2, R1,k1,k2
    REAL*8 :: r, A, l1, l2, R0, R2, q, T0, T2, lambda
    REAL*8 :: alpha_s , alpha_l 

    !---------------------------------------

    x = Xin(1) ; y = Xin(2)
    r = sqrt(x*x + y*y)
    IF ( nDim == 3 ) THEN
       z = Xin(3) ; r = sqrt(x*x + y*y + z*z)
    END IF


    T0 = Icing_T1 ; T2 = Icing_T2
	
    SELECT CASE ( icas )

    CASE ( 101 )
       L1 = Icing_L1 ; L2 = Icing_L2
       q = - (T2 - T0)/(L1/lambda1 + L2/lambda2)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          exa(1) = q
          exa(2) = 0
          exa(3) = T0 - q/lambda1*x
       ELSE
          exa(1) = q
          exa(2) = 0
          exa(3) = T0 - (L1/lambda1 + (x-L1)/lambda2)*q
       END IF
    CASE ( 102 )
      R0 = Icing_L1 ; R1 = Icing_L3 ; R2 = Icing_L2
       A = (T2-T0)/(log(R2)/lambda2 - log(R0)/lambda1 + &
       log(R1)*(1.d0/lambda1 - 1.d0/lambda2))
       exa(1) = (-A*x)/r**2
       exa(2) = (-A*y)/r**2
       IF ( ele%ref == IcingZone1_Ref ) THEN
          exa(3) = (A/lambda1)*log(r/R0) + T0 +5
       ELSE
          exa(3) = (A/lambda2)*log(r/R2) + T2
       END IF
     CASE ( 122 )
       lambda = ele%lambda
       IF ( ele%ref == IcingZone1_Ref ) THEN
          exa(3) = log(r)/log(3.d0) + 3
       ELSE
          exa(3) = log(3.d0*r)/log(3.d0)
       END IF
       Exa(1) = -lambda*x/(log(3.d0)*r**2)
       Exa(2) = -lambda*y/(log(3.d0)*r**2)
    CASE ( 103 )
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          exa(1) = -lambda1
          exa(2) = 0.d0
          exa(3) = x
       ELSE
          exa(1) = -10.d0*lambda2
          exa(2) = 0.d0
          exa(3) = 10.d0*x
       END IF
     CASE ( 104 )
        IF ( ele%Ref == IcingZone1_Ref ) THEN
           exa(1) = -2*Icing_Lambda1*x
           exa(2) = 0
           exa(3) = x**2
        ELSE
           exa(1) = -20*Icing_Lambda2*x
           exa(2) = 0
           exa(3) = 10*x**2 ! -9*Icing_L1**2
        END IF
      CASE ( 105 )
        Exa(1) = 2.d0*pi*sin(2.d0*pi*x)*cos(2.d0*PI*y)*sin(2.d0*pi*x)*sin(2.d0*pi*y) - &
                 2.d0*pi*cos(2.d0*pi*x)*cos(2.d0*PI*y)*cos(2.d0*pi*x)*sin(2.d0*pi*y)
        Exa(2) = 2.d0*pi*sin(2.d0*pi*y)*cos(2.d0*PI*x)*sin(2.d0*pi*x)*sin(2.d0*pi*y) - &
                 2.d0*pi*cos(2.d0*pi*x)*cos(2.d0*PI*y)*cos(2.d0*pi*y)*sin(2.d0*pi*x)
        Exa(3) = cos(2.d0*pi*x)*cos(2.d0*pi*y)*sin(2.d0*pi*x)*sin(2.d0*pi*y)
      CASE ( 106 )
         IF ( ele%Ref == IcingZone1_Ref ) THEN
            exa(1) = -3*Icing_Lambda1*x**2
            exa(2) = 0.d0
            exa(3) = x**3
         ELSE
            exa(1) = -10*3*Icing_Lambda2*x**2
            exa(2) = 0.d0
            exa(3) = 10*x**3-9*Icing_L1**3
         END IF
       CASE ( 107 )
          IF ( ele%ref == IcingZone1_Ref ) THEN
             exa(1) = (-lambda1*x)/r**2
             exa(2) = (-lambda1*y)/r**2
             exa(3) = log(r)
          ELSE
            exa(1) = (-lambda2*x)/r**2
            exa(2) = (-lambda2*y)/r**2
            exa(3) = log(r)
          END IF
    CASE DEFAULT
       PRINT*, " ** ERROR in ExactSolJump :: unknown test case"
       PRINT*, " "
       STOP
    END SELECT

  END SUBROUTINE ExactSolJump
  !==============================================!

  !==============================================!
  SUBROUTINE  Finverse(ele,x,y,n,nDim,L)
  !==============================================!

  !*************************************************************************
  !   Sends from the current element ele to the reference element
  !
  !  Parameters:
  !
  !    Input,Output, REAL*8, L, value on the referent element
  !    Input, REAL*8  x, value on the x axe
  !    Input, REAL*8  y, value on the y axe
  !    Input, INTEGER  n, wanted point number
  !    Input, type(Element) ele, Current element
  !    Input, INTEGER  nDim, Dimension of the problem (2 or 3)
  !

  !*************************************************************************


    IMPLICIT NONE
    REAL*8, INTENT(INOUT)      :: L
    REAL*8, INTENT(IN)         :: x
    REAL*8, INTENT(IN)         :: y
    INTEGER, INTENT(IN)        :: n
    type(Element), INTENT(IN)  :: ele
    INTEGER      , INTENT(IN)  :: nDim


    !--------------------------------------------
    REAL*8, DIMENSION(:)    , ALLOCATABLE :: X1,X2,X3
    !--------------------------------------------
    REAL*8  ::xx,yy, A
    !--------------------------------------------
    INTEGER :: nVar
    !--------------------------------------------

    nVar = nDim + 1
    ALLOCATE ( X1(nDim), X2(nDim), X3(nDim))
    X1=ele%Coor(1,:) ; X2=ele%Coor(2,:) ; X3=ele%Coor(3,:)

    A = (X2(1)-X1(1))*(X3(2)-X1(2))-(X2(2)-X1(2))*(X3(1)-X1(1))
    xx=1/A*(x-X1(1))*(X3(2)-X1(2))+(y-X1(2))*(X1(1)-X3(1))
    yy =1/A*(x-X1(1))*(X1(2)-X2(2))+(y-X1(2))*(X2(1)-X1(1))

    SELECT CASE ( n )
     CASE ( 1 )
       L=1-xx-yy
     CASE ( 2 )
       L=xx
     CASE ( 3 )
       L=yy
    CASE DEFAULT
       PRINT*, " ** ERROR in Finverse **"
       PRINT*, " "
       STOP
    END SELECT




  END SUBROUTINE Finverse
  !==============================================!


END MODULE ExactFunctionJump
