MODULE SourceTerm

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS


  !======================================================!
  SUBROUTINE Compute_SourceTerm(phi, Xq, T, ele, icas, nDim)
  !======================================================!

  !*************************************************************************
  !   Gives the value of the source term phi for the point Xq
  !   11* embedded resolution for time dependent cases
  !   10* embedded resolution for non-time dependent cases

  !  Parameters:
  !
  !    Input,Output, REAL*8, DIMENSION(:) phi, Value of the solution
  !    Input,Output, REAL*8, DIMENSION(:) Xq, Current point
  !    Input, type(element) ele, Current element
  !    Input, INTEGER   icas, Case number from the dom.data file
  !    Input, INTEGER   nDim, Dimension of the problem (2 or 3)
  !    Input, REAL*8 T, Current time
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(INOUT) :: phi
    REAL*8, DIMENSION(:), INTENT(IN)    :: Xq
    type(element)       , INTENT(IN)    :: ele
    INTEGER             , INTENT(IN)    :: icas, nDim
    REAL*8              , INTENT(IN)    :: T
    !-------------------------------------------------
    REAL*8 :: x, y, z, denom, mu1, mu2,T0,T2,R0,R1,R2,lambda1,lambda2
    REAL*8 :: A,r,k1,k2,q, L1
    !-----------------------------------

    phi = 0.d0
    lambda1 = Icing_Lambda1 ; lambda2 = Icing_Lambda2

    SELECT CASE ( nDim )
    CASE ( 2 )
       x = Xq(1) ; y = Xq(2)
    CASE ( 3 )
       x = Xq(1) ; y = Xq(2) ; z = Xq(3)
    END SELECT

    r=sqrt(x*x+y*y)

    SELECT CASE ( icas )
	 CASE (2,3,4,5,7,102,101,107,134,001)
       phi =0.d0
     CASE(100012)
		
	  IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = Icing_rho*icing_c1*6*t -2*lambda1
       ELSE
          phi = Icing_rho*icing_c2*6*t -2*lambda2
       END IF
     CASE(201)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 3*t**2
       ELSE
         phi = 3*t**2
       END IF
     CASE(202)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 3*t**2
       ELSE
         phi = 2
       END IF
     CASE(111)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 3*t**2+lambda1/((x+1)**2)
       ELSE
          phi = 3*t**2+lambda2/((x+1)**2)
       END IF
       CASE(1111)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 3*x*y-Icing_Lambda1*(12*x**2)
       ELSE
          phi = 3*t**2+lambda2/((x+1)**2)
       END IF
     CASE(1044)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = -Icing_Lambda1*20
       ELSE
         phi = -Icing_Lambda2*20
       END IF
     CASE(777)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = -2*Icing_Lambda1
       ELSE
          phi = -2*Icing_Lambda2
       END IF
     CASE(773)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = x**2+2.d0/3.d0*y*t-2*lambda1*t
       ELSE
         phi = x
       END IF
     CASE(112)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 1
       ELSE
         phi = x*y
       END IF
     CASE(132)
         phi = 1
     CASE(115)
           phi = 2*t
     CASE(113)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 0.5d0*(1.d0/(t+1.d0))-2*lambda1
       ELSE
         phi = 0.5d0*(1.d0/(t+1.d0))-2*lambda2
       END IF
     CASE(114)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi =1
       ELSE
         phi =1
       END IF
     CASE(137)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 3*x
       ELSE
         phi = x
       END IF
     CASE(138)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = x*y
       ELSE
         phi = x
       END IF
    CASE(106)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = -6*Icing_Lambda1*x
      ELSE
         phi = -6*Icing_Lambda2*x
      END IF
    CASE(166)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
        phi = -6*Icing_Lambda1*(x+y)
      ELSE
        phi = -6*Icing_Lambda2*(x+y)
      END IF
    CASE(104)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = -20*Icing_Lambda1
      ELSE
         phi = -2*Icing_Lambda2
      END IF
    CASE(108)
      k1=log(3.d0)/lambda2
      k2=1.d0/lambda2-1.d0/lambda1
      IF ( ele%Ref == IcingZone1_Ref ) THEN
        phi = (144*k1*k2)/(r**2*(k2*2*log(r)+2*k1)**3)
      ELSE
        phi = (144*k1*k2+144*k2**2*log(3.d0))/(r**2*(k2*2*log(r)+2*k1)**3)
      END IF
     CASE(109)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = -2*Icing_Lambda1
       ELSE
          phi = -18*Icing_Lambda2
       END IF
     CASE (20211,20212,20213,20214) !T polynome de degre 0 continue
          phi = 0
     CASE(20215)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = -Icing_Lambda1*(6+2)
       ELSE
          phi = -Icing_Lambda2*(6+2)
       END IF
     CASE(20216)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = -Icing_Lambda1*(6+2)
       ELSE
          phi = -Icing_Lambda2*(8)
       END IF
     CASE(20217)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = -Icing_Lambda1*(6+6*x+2+18*x)
       ELSE
          phi = -Icing_Lambda2*(6+6*x+2+18*x)
       END IF
     CASE(20218)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = -Icing_Lambda1*(6+6*x+2+18*x)
       ELSE
          phi = -Icing_Lambda2*(8+2*y**2+2*x**2+6*y)
       END IF
     CASE(20219)
          phi = 7
     CASE(202110)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 7
       ELSE
          phi = 3
       END IF
     CASE(202111)
          phi = 7+6*t
     CASE(202112)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 7+6*t
       ELSE
          phi = 3 + 22*t
       END IF
     CASE(202113)
          phi = 7+6*t + 33*t**2
     CASE(202114)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 7+6*t + 33*t**2
       ELSE
          phi = 3 + 22*t + 1.5d0*t**2
       END IF
     CASE(202115)
          phi = 5
     CASE(202116)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 5
       ELSE
          phi = 3
       END IF
     CASE(202117)
          phi = 5 + 46*t
     CASE(202118)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 5 + 46*t
       ELSE
          phi = 3 + 18*t
       END IF
     CASE(202119)
          phi = 5 + 46*t+ 3*t**2
     CASE(202120)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          phi = 5 + 46*t + 3*t**2
       ELSE
          phi = 3 + 18*t
       END IF
     CASE(202121)
          IF ( ele%Ref == IcingZone1_Ref ) THEN
            phi = 11 - Icing_Lambda1*42
          ELSE
            phi = 11 - Icing_Lambda2*42
          END IF
     CASE(202122)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 11 - Icing_Lambda1*42
       ELSE
          phi = 2 - Icing_Lambda2*2
       END IF
     CASE(202123)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 11 + 12*t - Icing_Lambda1*42
       ELSE
         phi = 11 + 12*t - Icing_Lambda2*42
       END IF
     CASE(202124)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 11 + 12*t - Icing_Lambda1*42
       ELSE
          phi = 2-16*t - Icing_Lambda2*2
       END IF
     CASE(202125)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 11 + 12*t-42*t**2 - Icing_Lambda1*42
       ELSE
         phi = 11 + 12*t-42*t**2 - Icing_Lambda2*42
       END IF
     CASE(202126)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 11 + 12*t-42*t**2 - Icing_Lambda1*42
       ELSE
          phi = 2-16*t-3*t**2 - Icing_Lambda2*2
       END IF
     CASE(202127)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 7-Icing_Lambda1*42*y
       ELSE
          phi =7 -Icing_Lambda2*42*y
       END IF
     CASE(202128)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 7-Icing_Lambda1*42*y
       ELSE
          phi = -10 -Icing_Lambda2*(6*x+10*y)
       END IF
     CASE(202129)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 7-Icing_Lambda1*42*y-6*T
       ELSE
          phi =-Icing_Lambda2*42*y + 5*t**4
       END IF
     CASE(202130)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 7-Icing_Lambda1*42*y-6*T
       ELSE
          phi = -10 -Icing_Lambda2*(6*x+4*y+6*y)
       END IF
     CASE(202131)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 7-Icing_Lambda1*42*y-6*T+ 36*t**2
       ELSE
          phi =7 -Icing_Lambda2*42*y-6*T+ 36*t**2
       END IF
     CASE(202132)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 7-Icing_Lambda1*42*y-6*T
       ELSE
          phi = -10 -Icing_Lambda2*(6*x+4*y+6*y) + 36*t**2
       END IF
     CASE(202133)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = -11 + 12*x -3*y + 8*x*y +18*t -Icing_Lambda1*(42-30)
       ELSE
         phi = -11 + 12*x -3*y + 8*x*y +18*t -Icing_Lambda2*(42-30)
       END IF
     CASE(202134)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = -11 + 12*x -3*y +18*t -Icing_Lambda1*(126*x-60*y)
       ELSE
          phi = 66*x*y*t
       END IF
	CASE(202135)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 5*t**4 - Icing_Lambda1*(48*x*y**2+42*y**5+16*x**3)
       ELSE
          phi = 5*t**4 - Icing_Lambda2*(48*x*y**2+42*y**5+16*x**3)
       END IF
       
    CASE(202136)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
         phi = 1.d0/(2.d0*(t+1)) + 3*t**2 - Icing_Lambda1*(9*y**2*exp(3*x*y)+9*x**2*exp(3*x*y))
       ELSE
          phi = 1.d0/sqrt(t) - Icing_Lambda2*(9*y**2*exp(3*x*y)+9*x**2*exp(3*x*y))
       END IF
       
    END SELECT


  END SUBROUTINE Compute_SourceTerm
  !========================================================!


    !======================================================!
    SUBROUTINE Compute_SourceTermGrad(Gphi, Xq, T, ele, icas, nDim)
    !======================================================!

    !*************************************************************************
    !   Gives the value of the gradient of thhe source term phi for the point Xq


    !  Parameters:
    !
    !    Input,Output, REAL*8, DIMENSION(:) Gphi, Value of the gradient solution
    !    Input,Output, REAL*8, DIMENSION(:) Xq, Current point
    !    Input, type(element) ele, Current element
    !    Input, INTEGER   icas, Case number from the dom.data file
    !    Input, INTEGER   nDim, Dimension of the problem (2 or 3)
    !    Input, REAL*8 T, Current time
    !
    !*************************************************************************


      IMPLICIT NONE
      REAL*8, DIMENSION(:), INTENT(INOUT) :: Gphi
      REAL*8, DIMENSION(:), INTENT(IN)    :: Xq
      type(element)       , INTENT(IN)    :: ele
      INTEGER             , INTENT(IN)    :: icas, nDim
      REAL*8              , INTENT(IN)    :: T
      !-------------------------------------------------
      REAL*8 :: x, y, z, denom, mu1, mu2,T0,T2,R0,R1,R2,lambda1,lambda2
      REAL*8 :: A,r,k1,k2,q
      !-----------------------------------

      Gphi = 0.d0
      lambda1 = Icing_Lambda1 ; lambda2 = Icing_Lambda2

      SELECT CASE ( nDim )
      CASE ( 2 )
         x = Xq(1) ; y = Xq(2)
      CASE ( 3 )
         x = Xq(1) ; y = Xq(2) ; z = Xq(3)
      END SELECT

      r=sqrt(x*x+y*y)

      SELECT CASE ( icas )
      CASE (1,2,3,4,5,7,102,101,107,134,1044,114,104,109)
         Gphi =0
       CASE(132)
         IF ( ele%Ref == IcingZone1_Ref ) THEN
            Gphi = 0
         ELSE
           Gphi  = 0
         END IF
       CASE(137)
         IF ( ele%Ref == IcingZone1_Ref ) THEN
            Gphi(1) = 3 ; Gphi(2)=0
         ELSE
           Gphi(1) = 1 ; Gphi(2)=0
         END IF
       CASE(111)
         IF ( ele%Ref == IcingZone1_Ref ) THEN
            Gphi(1) =-lambda1*((x+1)/((x+1)**4)) ; Gphi(2) = 0
         ELSE
           Gphi(1) = 0 ; Gphi(2)= 3
         END IF
       CASE(112)
         IF ( ele%Ref == IcingZone1_Ref ) THEN
           Gphi(1) = 0 ; Gphi(2)= 0
         ELSE
           Gphi(1) = y ; Gphi(2)= x
         END IF
       CASE(138)
         IF ( ele%Ref == IcingZone1_Ref ) THEN
           Gphi(1) = y ; Gphi(2)= x
         ELSE
           Gphi(1) = 1 ; Gphi(2)=0
         END IF
       CASE(113)
         IF ( ele%Ref == IcingZone1_Ref ) THEN
            Gphi(1) = 9*x**2-18*Icing_Lambda1*t ; Gphi(2)= 0
         ELSE
           Gphi(1) = 0 ; Gphi(2)= -24*Icing_Lambda2*y
         END IF
      CASE(106)
        IF ( ele%Ref == IcingZone1_Ref ) THEN
          Gphi(1) = -6*Icing_Lambda1 ; Gphi(2)= 0
        ELSE
          Gphi(1) = -10*6*Icing_Lambda2 ; Gphi(2)= 0
        END IF
       CASE DEFAULT
         Gphi = 0
      END SELECT


    END SUBROUTINE Compute_SourceTermGrad
    !========================================================!
END MODULE SourceTerm
