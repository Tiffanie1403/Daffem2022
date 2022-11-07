MODULE ExactSolution

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !==========================================!
  SUBROUTINE ExactSol(L, Exa, Xin, Ele, nDim)
  !==========================================!

  !*************************************************************************
  !  Analytical solutions of the problem
  !  at the point Xin, subroutine for non time dependent problem
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
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: L
    REAL*8, DIMENSION(:)  , INTENT(INOUT) :: Exa
    REAL*8, DIMENSION(:)  , INTENT(IN)    :: Xin
    type(element)         , INTENT(IN)    :: Ele
    INTEGER               , INTENT(IN)    :: nDim
    !--------------------------------------------
    REAL*8 :: x, y, z
    REAL*8 :: r, theta, c, s, Br, Bt, mu1, mu2, denom
    !-------------------------------------------------

    SELECT CASE ( nDim )
    CASE ( 2 )
       x = Xin(1) ; y = Xin(2)
    CASE ( 3 )
       x = Xin(1) ; y = Xin(2) ; z = Xin(3)
    END SELECT
    SELECT CASE ( icas )
    CASE ( 1 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = 0.d0
          Exa(2) = 0.d0
          Exa(3) = 1.d0
       CASE ( 3 )
          Exa(1) = 0.d0
          Exa(2) = 0.d0
          Exa(3) = 0.d0
          Exa(4) = 1.d0
       END SELECT
       WRITE(*,*) 12

    CASE ( 2 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = -L(1,1)
          Exa(2) = 0.d0
          Exa(3) = x
       CASE ( 3 )
          Exa(1) = -L(1,1)
          Exa(2) = 0.d0
          Exa(3) = 0.d0
          Exa(4) = x
       END SELECT
       WRITE(*,*) 22

    CASE ( 3 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = 0.d0
          Exa(2) = -L(2,2)
          Exa(3) = y
       CASE ( 3 )
          Exa(1) = 0.d0
          Exa(2) = -L(2,2)
          Exa(3) = 0.d0
          Exa(4) = y
       END SELECT
       WRITE(*,*) 32

    CASE ( 4 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = -L(1,1)
          Exa(2) = -L(2,2)
          Exa(3) = x + y
       CASE ( 3 )
          Exa(1) = 0.d0
          Exa(2) = 0.d0
          Exa(3) = -L(3,3)
          Exa(4) = z
       END SELECT
       WRITE(*,*) 42

    CASE ( 5 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = -L(1,1)
          Exa(2) = -L(1,1)
          Exa(3) = x + y
       CASE ( 3 )
          Exa(1) = -L(1,1)
          Exa(2) = -L(2,2)
          Exa(3) = -L(3,3)
          Exa(4) = x + y + z
       END SELECT
       WRITE(*,*) 52

    CASE ( 6 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = 2.d0*pi*sin(2.d0*pi*x)*cos(2.d0*PI*y)*sin(2.d0*pi*x)*sin(2.d0*pi*y) - &
                   2.d0*pi*cos(2.d0*pi*x)*cos(2.d0*PI*y)*cos(2.d0*pi*x)*sin(2.d0*pi*y)
          Exa(2) = 2.d0*pi*sin(2.d0*pi*y)*cos(2.d0*PI*x)*sin(2.d0*pi*x)*sin(2.d0*pi*y) - &
                   2.d0*pi*cos(2.d0*pi*x)*cos(2.d0*PI*y)*cos(2.d0*pi*y)*sin(2.d0*pi*x)
          Exa(3) = cos(2.d0*pi*x)*cos(2.d0*pi*y)*sin(2.d0*pi*x)*sin(2.d0*pi*y)
       CASE ( 3 )
          Exa(1) = 2*pi*sin(2*pi*x)**2*cos(2*pi*y)*sin(2*pi*y)*cos(2*pi*z)*sin(2*pi*z&
               & )-2*pi*cos(2*pi*x)**2*cos(2*pi*y)*sin(2*pi*y)*cos(2*pi*z)*sin(2&
               & *pi*z)
          Exa(2) = 2*pi*cos(2*pi*x)*sin(2*pi*x)*sin(2*pi*y)**2*cos(2*pi*z)*sin(2*pi*z&
               &)-2*pi*cos(2*pi*x)*sin(2*pi*x)*cos(2*pi*y)**2*cos(2*pi*z)*sin(2&
               &*pi*z)
          Exa(3) = 2*pi*cos(2*pi*x)*sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*y)*sin(2*pi*z)**&
               &2-2*pi*cos(2*pi*x)*sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*y)*cos(2*pi&
               &*z)**2
          Exa(4) = cos(2*pi*x)*sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*y)*&
               cos(2*pi*z)*sin(2*pi*z)
       END SELECT
       WRITE(*,*) 62

    CASE ( 7 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          r = sqrt(x*x + y*y)
          theta = atan2(y,x)
          Exa(3) = 1.d0 + (r + 1.d0/r)*cos(theta)
          Br = -(1.d0 - 1.d0/r**2)*cos(theta)
          Bt = (r + 1.d0/r)*sin(theta)/r
          c = cos(theta) ; s = sin(theta)
          Exa(1) = c*Br - s*Bt
          Exa(2) = s*Br + c*Bt
       CASE ( 3 )
          r = sqrt(x*x + y*y)
          theta = atan2(y,x)
          Exa(4) = 1.d0 + (r + 1.d0/r)*cos(theta)
          Br = -(1.d0 - 1.d0/r**2)*cos(theta)
          Bt = (r + 1.d0/r)*sin(theta)/r
          c = cos(theta) ; s = sin(theta)
          Exa(1) = c*Br - s*Bt
          Exa(2) = s*Br + c*Bt
          Exa(3) = 0.d0
       END SELECT
       WRITE(*,*) 72
    CASE ( 8 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = -2.5d-1*(cos(x)*sinh(y) + 2.d0*x/(y**2+x**2) - 2.d0*x)
          Exa(2) = -2.5d-1*(sin(x)*cosh(y) + 2.d0*y/(y**2+x**2) - 2.d0*y)
          Exa(3) = 2.5d-1*(log(y**2 + x**2)+ sin(x)*sinh(y) - y**2 - x**2 - &
               2.d0*log(3.d0) + 9.d0)
       CASE ( 3 )
          PRINT*, "ERROR in ExactSol, case not implemented in 3D"
          STOP
       END SELECT
       WRITE(*,*) 82

    CASE ( 9 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = 2.d0*lxy*pi*cos(2.d0*pi*x)*sin(2.d0*pi*x)*sin(2.d0*pi*y)**2+&
               (2.d0*lxx*pi*sin(2.d0*pi&
               &*x)**2-2.d0*lxx*pi*cos(2.d0*pi*x)**2)*cos(2.d0*pi*y)*sin(2.d0*pi*y)-&
               2.d0*lxy*p&
               &i*cos(2.d0*pi*x)*sin(2.d0*pi*x)*cos(2.d0*pi*y)**2
          Exa(2) = 2.d0*lyy*pi*cos(2.d0*pi*x)*sin(2.d0*pi*x)*sin(2.d0*pi*y)**2+&
               (2.d0*lyx*pi*sin(2.d0*pi&
               &*x)**2-2.d0*lyx*pi*cos(2.d0*pi*x)**2)*cos(2.d0*pi*y)*sin(2.d0*pi*y)-&
               2.d0*lyy*p&
               &i*cos(2.d0*pi*x)*sin(2.d0*pi*x)*cos(2.d0*pi*y)**2
          Exa(3) = cos(2.d0*pi*x)*sin(2.d0*pi*x)*cos(2.d0*pi*y)*sin(2.d0*pi*y)
       CASE ( 3 )
          PRINT*, "ERROR in ExactSol, case not implemented in 3D"
          STOP
       END SELECT
       WRITE(*,*) 92

    CASE ( 10 )
       mu1 = mu_1
       mu2 = mu_2
       denom = 0.5d0*(mu1+mu2) + 4.d0*mu1*mu2
       Exa(1) = -mu1*mu2/denom
       Exa(2) = 0.d0
       IF ( ele%ref == 4 ) THEN
          Exa(3) = 1.d0/denom*(mu2*(x+0.5d0) + 2.d0*mu1*mu2)
       ELSE IF ( ele%ref == 5 ) THEN
          Exa(3) = 1.d0/denom*(mu1*(x+0.5d0) + 2.d0*mu1*mu2 + 0.5d0*(mu2-mu1))
       ELSE
          PRINT*, " ** ERROR :: for test case 10, references for permeability regions are 4 and 5"
          STOP
       END IF
       WRITE(*,*) 102

    CASE ( 11 )
       mu1 = mu_1 ; mu2 = mu_2
       denom = 0.5d0*(mu1 + mu2) + 4.d0*mu1*mu2
       Exa(1) = -mu1*mu2/denom * (x + 0.5d0)
       Exa(2) = 0.d0
       IF ( ele%ref == 4 ) THEN
          Exa(3) = 1.d0/denom*(mu2/2.d0*(x + 0.5d0)**2 + 2.d0*mu1*mu2)
       ELSE IF ( ele%ref == 5 ) THEN
          Exa(3) = 1.d0/denom*(mu1/2.d0*(x + 0.5d0)**2 + &
               2.d0*mu1*mu2 + 1.d0/8.d0*(mu2-mu1))
       ELSE
          PRINT*, " ** ERROR :: for test case 11, references for permeability regions are 4 and 5"
          STOP
       END IF
       WRITE(*,*) 112

    CASE ( 12 )
       mu1 = mu_1 ; mu2 = mu_2
       denom = 0.5d0*(mu1 + mu2) + 4.d0*mu1*mu2
       Exa(1) = -mu1*mu2/denom * (x + 0.5d0)**2
       Exa(2) = 0.d0
       IF ( ele%ref == 4 ) THEN
          Exa(3) = 1.d0/denom*(mu2/3.d0*(x + 0.5d0)**3 + 2.d0*mu1*mu2)
       ELSE
          Exa(3) = 1.d0/denom*(mu1/3.d0*(x + 0.5d0)**3 + &
               2.d0*mu1*mu2 + 1.d0/24.d0*(mu2-mu1))
       END IF
       WRITE(*,*) 122

    CASE ( 0 )
       SELECT CASE ( nDim )
       CASE ( 2 )
          Exa(1) = ele%L(1,1)
          Exa(2) = 0.d0
          Exa(3) = 1.d0 - x
       CASE ( 3 )
          Exa(1) = ele%L(1,1)
          Exa(2) = 0.d0
          Exa(3) = 0.d0
          Exa(4) = 1.d0 - x
       END SELECT
       WRITE(*,*) 0

    CASE DEFAULT
       PRINT*, " ERROR in ExactSol :: unknown test case"
       STOP
    END SELECT

  END SUBROUTINE ExactSol
  !==================================!

END MODULE ExactSolution
