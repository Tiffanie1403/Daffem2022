MODULE ExactFunctionJumpt

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !==============================================!
  SUBROUTINE ExactSolJumpt(L, Exa, Xin, t, ele, nDim)
  !==============================================!

  !*************************************************************************
  !   Gives the exact value of the solutions on the element ele
  !   at the point Xin, subroutine for time dependent problem
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(:,:) L, permeability matrix associated to the current element
  !    Output, REAL*8, DIMENSION(:) exa, analytical solution at the point Xin
  !    Input, REAL*8, DIMENSION(:) Xin, Current point on the element ele
  !    Input, type(Element) ele, Current element
  !    Input, INTEGER  nDim, Dimension of the problem (2 or 3)
  !    Input, INTEGER  t, Current time value
  !

  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN)  :: L    ! lambda associe a la zone ou se trouve l'element
    REAL*8, DIMENSION(:)  , INTENT(OUT) :: exa  ! resultat solution exacte
    REAL*8, DIMENSION(:)  , INTENT(IN)  :: xin  !coordonnees des nodes sur l'element courant
    type(Element)         , INTENT(IN)  :: ele  !element etudie
    INTEGER               , INTENT(IN)  :: nDim !dimension du probleme
    REAL*8                , INTENT(IN)  :: t    ! variable temporelle

    !--------------------------------------------
    REAL*8 :: x, y, z, lambda1, lambda2, R1,k1,k2,res1
    REAL*8 :: r, A, l1, l2, R0, R2, q, T0, T2, lambda, chi
    REAL*8 :: alpha_l , alpha_s 

    !---------------------------------------

    x = Xin(1) ; y = Xin(2)
    r = sqrt(x*x + y*y)
    IF ( nDim == 3 ) THEN
       z = Xin(3) ; r = sqrt(x*x + y*y + z*z)
    END IF
	res1 = 0.d0

    T0 = Icing_T1 ; T2 = Icing_T2
    lambda1 = Icing_Lambda1 ; lambda2 = Icing_Lambda2
	
	chi = Icing_chi 
	alpha_l = lambda1/(Icing_rho*Icing_c1)
	alpha_s = lambda2/(Icing_rho*Icing_c2)
	
    SELECT CASE ( icas )

	 CASE ( 001 )
	 IF (x==Icing_L2 .OR. t==0) THEN !right boundary where we make the assumption of having a semi-infinite slab 
		
			exa(1) = 0
			exa(2) = 0
			Exa(3) = Icing_Tinit ! we only have ice at t=0 which means only area 2
		ElSE IF(x==0) THEN !left boundary 
		
			res1 = (Icing_Tw-Icing_Tm)/erf(chi)

			exa(1) = -lambda1*(-res1*1.d0/(sqrt(PI*alpha_l*t)))
			exa(2) = 0
			Exa(3) = Icing_Tw
		ELSE 
			! le problème est indépendant de la position en ordonnée 
			IF ( ele%Ref == IcingZone1_Ref ) THEN
				res1 = (Icing_Tw-Icing_Tm)/erf(chi)

				exa(1) = -lambda1*(-res1/(sqrt(PI*alpha_l*t))*exp(-x**2/(4.d0*alpha_l*t))) !Bx
				exa(2) = 0.d0 ! By
				exa(3) = Icing_Tw - res1*erf(x/(2.d0*sqrt(alpha_l*t))) ! T 
		    ELSE
				res1 = (Icing_Tm-Icing_Tinit)/erfc(chi*sqrt(alpha_l/alpha_s))
			  
				exa(1) = -lambda2*(-res1/(sqrt(PI*alpha_s*t))*exp(-x**2/(4.d0*alpha_s*t)))
				exa(2) = 0.d0
				exa(3) = Icing_Tinit + res1*erfc(x/(2.d0*sqrt(alpha_s*t)))
		    END IF
        END IF 
    CASE(100012) 
		IF ( ele%Ref == IcingZone1_Ref ) THEN

				exa(1) = -lambda1*(2*x+3*y) !Bx
				exa(2) = -lambda1*(3*x + 7) !By
				exa(3) = x**2+3*x*y+7*y+1+3*t**2 ! T 
		    ELSE
	
				exa(1) = -lambda2*(2*x+3*y) !Bx
				exa(2) = -lambda2*(3*x + 7) !By
				exa(3) = x**2+3*x*y+7*y+1+3*t**2 ! T 
		    END IF
    CASE ( 101 )
       L1 = 1 ; L2 = Icing_L2
       q = - (T2 - T0)/(L1/lambda1 + L2/lambda2)
       IF ( ele%Ref == IcingZone1_Ref ) THEN
          exa(1) = q
          exa(2) = 0
          exa(3) = (T0 - q/lambda1*x)
       ELSE
          exa(1) = q
          exa(2) = 0
          exa(3) = (T0 - (L1/lambda1 + (x-L1)/lambda2)*q)
       END IF
    CASE(111) ! cas de type (11_) dependant du temps
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         exa(1) = -lambda1/(x+1)
         exa(2) = -lambda1
         exa(3) = log(x+1)+t**3+y
      ELSE
          exa(1) = -lambda2/(x+1)
         exa(2) = -lambda2
         exa(3) = log(x+1)+t**3+y
      END IF
      CASE(1111) ! cas de type (11_) dependant du temps
      IF ( ele%Ref == IcingZone1_Ref ) THEN
          exa(1) = -lambda1*(4*x**3+3*y*t)
         exa(2) = -lambda1*(3*x*t)
         exa(3) = x**4+3*x*y*t+7
      ELSE
          exa(1) = -lambda2/(x+1)
         exa(2) = -lambda2
         exa(3) = log(x+1)+t**3+y
      END IF
    CASE(112)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         exa(1) = -lambda1
         exa(2) = 0
         exa(3) = x+t
      ELSE
        exa(1) = -lambda2*t*y
        exa(2) = -lambda2*t*x
        exa(3) = x*y*t
      END IF
    CASE(132)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         exa(1) = -lambda1
         exa(2) = 0
         exa(3) = x+t
      ELSE
        exa(1) = -lambda2
        exa(2) = 0
        exa(3) = x+t
      END IF
    CASE(137)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         exa(1) = -lambda1*3*t
         exa(2) = 0
         exa(3) = 3*x*t
      ELSE
        exa(1) = -lambda2*t
        exa(2) = 0
        exa(3) = x*t
      END IF
    CASE(138)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
        exa(1) = -lambda1*y*t
        exa(2) = -lambda1*x*t
        exa(3) = x*y*t
      ELSE
        exa(1) = -lambda2*t
        exa(2) = 0
        exa(3) = x*t
      END IF
    CASE(777)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
        exa(1) = -lambda1*2*x
        exa(2) = -lambda1*1.d0/3.d0
        exa(3) = x**2 + 1.d0/3.d0*y
      ELSE
        exa(1) = -lambda2*2*x
        exa(2) = -lambda2*1.d0/3.d0
        exa(3) = x**2 + 1.d0/3.d0*y
      END IF
    CASE(773)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
        exa(1) = -lambda1*2*x*t
        exa(2) = -lambda1*1.d0/3.d0*t**2
        exa(3) = x**2*t + 1.d0/3.d0*y*t**2
      ELSE
        exa(1) = -lambda2*t
        exa(2) = 0
        exa(3) = x*t
      END IF
    CASE(113)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
          exa(1) = 0
        exa(2) = -2*lambda1*y
        exa(3) = 0.5d0*log(t+1)+y**2
      ELSE
        exa(1) = 0
        exa(2) = -2*lambda2*y
        exa(3) = 0.5d0*log(t+1)+y**2
      END IF
    CASE(114)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         exa(1) = 0
         exa(2) = 0
         exa(3) = t
      ELSE
        exa(1) = 0
        exa(2) = 0
        exa(3) = t
      END IF
    CASE(115)
      IF ( ele%Ref == IcingZone1_Ref ) THEN
         exa(1) = 0
         exa(2) = 0
         exa(3) = t**2
      ELSE
        exa(1) = 0
        exa(2) = 0
        exa(3) = t**2
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
     CASE ( 104 )
        IF ( ele%Ref == IcingZone1_Ref ) THEN
          exa(1) = -20*Icing_Lambda1*x
          exa(2) = 0
          exa(3) = 10*x**2 ! -9*Icing_L1**2
        ELSE
          exa(1) = -2*Icing_Lambda2*x
          exa(2) = 0
          exa(3) = x**2
        END IF
      CASE ( 1044 )
         IF ( ele%Ref == IcingZone1_Ref ) THEN
           exa(1) = -Icing_Lambda1*(14*x+2*y)
           exa(2) = -Icing_Lambda1*(6*y+2*x)
           exa(3) = 7*x**2+3*y**2+2*x*y+7
         ELSE
           exa(1) = -Icing_Lambda2*(14*x+2*y)
           exa(2) = -Icing_Lambda2*(6*y+2*x)
           exa(3) = 7*x**2+3*y**2+2*x*y+7
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
            exa(1) = -3*Icing_Lambda2*x**2
            exa(2) = 0.d0
            exa(3) = x**3
         END IF
       CASE ( 166 )
          IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = -3*Icing_Lambda1*x**2
             exa(2) = -3*Icing_Lambda1*y**2
             exa(3) = x**3+y**3
          ELSE
            exa(1) = -3*Icing_Lambda2*x**2
            exa(2) = -3*Icing_Lambda2*y**2
            exa(3) = x**3+y**3
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
        CASE(201)
          IF ( ele%ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = t**3
          ELSE
            exa(1) = 0
            exa(2) = 0
            exa(3) = t**3
          END IF
        CASE(202)
          IF ( ele%ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = t**3
          ELSE
            exa(1) = 0
            exa(2) = 0
            exa(3) = 2*t+1
          END IF
        CASE ( 109 )
           IF ( ele%Ref == IcingZone1_Ref ) THEN
              exa(1) = -Icing_Lambda1*(2*x)
              exa(2) = -Icing_Lambda1*(7)
              exa(3) = x**2 +7*y + 3
           ELSE
             exa(1) = -Icing_Lambda2*(6*x+21*y)
             exa(2) = -Icing_Lambda2*(12*y+21*x)
             exa(3) = 3*x**2 +6*y**2 +21*x*y
           END IF !On commence par des polynomes sans dependances temporelle)
         CASE (20211) !T polynome de degre 0 continue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
              exa(1) = 0
              exa(2) = 0
              exa(3) = 7
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 7
           END IF
         CASE (20212) !T polynome de degre 0 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
              exa(1) = 0
              exa(2) = 0
              exa(3) = 7
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3
           END IF
         CASE (20213) !T polynome de degre 1 continue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
              exa(1) = -Icing_Lambda1*3
              exa(2) = -Icing_Lambda1*2
              exa(3) = 7 + 3*x + 2*y
           ELSE
             exa(1) = -Icing_Lambda2*3
             exa(2) = -Icing_Lambda2*2
             exa(3) = 7 + 3*x + 2*y
           END IF
         CASE (20214) !T polynome de degre 1 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
              exa(1) = -Icing_Lambda1*3
              exa(2) = -Icing_Lambda1*2
              exa(3) = 7 + 3*x + 2*y
           ELSE
             exa(1) = 0
             exa(2) = -Icing_Lambda2*6
             exa(3) = 3 + 6*y
           END IF
         CASE (20215) !T polynome de degre 2 continue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
              exa(1) = -Icing_Lambda1*(3+y+6*x)
              exa(2) = -Icing_Lambda1*(2+x+2*y)
              exa(3) = 7 + 3*x + 2*y + x*y+3*x**2+y**2
           ELSE
             exa(1) = -Icing_Lambda2*(3+y+6*x)
             exa(2) = -Icing_Lambda2*(2+x+2*y)
             exa(3) = 7 + 3*x + 2*y + x*y+3*x**2+y**2
           END IF
         CASE (20216) !T polynome de degre 2 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = -Icing_Lambda1*(3+y+6*x)
             exa(2) = -Icing_Lambda1*(2+x+2*y)
             exa(3) = 7 + 3*x + 2*y + x*y+3*x**2+y**2
           ELSE
             exa(1) = -Icing_Lambda2*(7*y+8*x)
             exa(2) = -Icing_Lambda2*(6+7*x)
             exa(3) = 3 + 6*y + 7*x*y + 4*x**2
           END IF
         CASE (20217) !T polynome de degre 3 continue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
              exa(1) = -Icing_Lambda1*(3+y+6*x+9*y**2+3*x**2)
              exa(2) = -Icing_Lambda1*(2+x+2*y+18*x*y)
              exa(3) = 7 + 3*x + 2*y + x*y+3*x**2+y**2+9*x*y**2+x**3
           ELSE
             exa(1) = -Icing_Lambda2*(3+y+6*x+9*y**2+3*x**2)
             exa(2) = -Icing_Lambda2*(2+x+2*y+18*x*y)
             exa(3) = 7 + 3*x + 2*y + x*y+3*x**2+y**2+9*x*y**2+x**3
           END IF
         CASE (20218) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = -Icing_Lambda1*(3+y+6*x+9*y**2+3*x**2)
             exa(2) = -Icing_Lambda1*(2+x+2*y+18*x*y)
             exa(3) = 7 + 3*x + 2*y + x*y+3*x**2+y**2+9*x*y**2+x**3
           ELSE
             exa(1) = -Icing_Lambda2*(7*y+8*x+2*x*y**2)
             exa(2) = -Icing_Lambda2*(6+7*x+2*y*x**2+3*y**2)
             exa(3) = 3 + 6*y + 7*x*y + 4*x**2+y**2*x**2+y**3
           END IF
           !debut des cas test avec dépendance en temps seulement
         CASE (20219) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = 7*t + 1
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 7*t + 1
           END IF
         CASE (202110) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = 7*t + 1
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t + 3
           END IF
         CASE (202111) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t**2 + 7*t + 1
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t**2 + 7*t + 1
           END IF
         CASE (202112) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t**2 + 7*t + 1
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t + 3 + 11*t**2
           END IF
         CASE (202113) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t**2 + 7*t + 1 + 11*t**3
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t**2 + 7*t + 1 + 11*t**3
           END IF
         CASE (202114) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t**2 + 7*t + 1 + 11*t**3
           ELSE
             exa(1) = 0
             exa(2) = 0
             exa(3) = 3*t + 3 + 11*t**2 + 0.5d0*t**3
           END IF
           !cas test avec x et t
         CASE (202115) !T polynome de degre 3 discontinue
           IF ( ele%Ref == IcingZone1_Ref ) THEN
             exa(1) = -Icing_Lambda1*2
             exa(2) = -Icing_Lambda1*(11)
             exa(3) = 3 + 2*x + 11 * y + 5*t
           ELSE
             exa(1) = -Icing_Lambda2*2
             exa(2) = -Icing_Lambda2*(11)
             exa(3) = 3 + 2*x + 11 * y + 5*t
           END IF
         CASE (202116) !T polynome de degre 3 discontinue
             IF ( ele%Ref == IcingZone1_Ref ) THEN
               exa(1) = -Icing_Lambda1*2
               exa(2) = -Icing_Lambda1*(11)
               exa(3) = 3 + 2*x + 11 * y + 5*t
             ELSE
               exa(1) = -Icing_Lambda2*7
               exa(2) = -Icing_Lambda2*7
               exa(3) = 4 + 7*x+7*y+3*t
             END IF
           CASE (202117) !T polynome de degre 3 discontinue
             IF ( ele%Ref == IcingZone1_Ref ) THEN
               exa(1) = -Icing_Lambda1*2
               exa(2) = -Icing_Lambda1*(11)
               exa(3) = 3 + 2*x + 11 * y + 5*t + 23*t**2
             ELSE
               exa(1) = -Icing_Lambda2*2
               exa(2) = -Icing_Lambda2*(11)
               exa(3) = 3 + 2*x + 11 * y + 5*t+23*t**2
             END IF
           CASE (202118) !T polynome de degre 3 discontinue
               IF ( ele%Ref == IcingZone1_Ref ) THEN
                 exa(1) = -Icing_Lambda1*2
                 exa(2) = -Icing_Lambda1*(11)
                 exa(3) = 3 + 2*x + 11 * y + 5*t + 23*t**2
               ELSE
                 exa(1) = -Icing_Lambda2*7
                 exa(2) = -Icing_Lambda2*7
                 exa(3) = 4 + 7*x+7*y+3*t + 9*t**2
               END IF
             CASE (202119) !T polynome de degre 3 discontinue
               IF ( ele%Ref == IcingZone1_Ref ) THEN
                 exa(1) = -Icing_Lambda1*2
                 exa(2) = -Icing_Lambda1*(11)
                 exa(3) = 3 + 2*x + 11 * y + 5*t + 23*t**2 + t**3
               ELSE
                 exa(1) = -Icing_Lambda2*2
                 exa(2) = -Icing_Lambda2*(11)
                 exa(3) = 3 + 2*x + 11 * y + 5*t+23*t**2 + t**3
               END IF
             CASE (202120) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*2
                   exa(2) = -Icing_Lambda1*(11)
                   exa(3) = 3 + 2*x + 11 * y + 5*t + 23*t**2 + t**3
                 ELSE
                   exa(1) = -Icing_Lambda2*7
                   exa(2) = -Icing_Lambda2*7
                   exa(3) = 4 + 7*x+7*y+3*t + 9*t**2
                 END IF
               CASE (202121) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(33*y)
                   exa(2) = -Icing_Lambda1*(33*x+42*y)
                   exa(3) =4 + 33*x*y + 21*y**2 + 11*t
                 ELSE
                   exa(1) = -Icing_Lambda2*(33*y)
                   exa(2) = -Icing_Lambda2*(33*x+42*y)
                   exa(3) =4 + 33*x*y + 21*y**2 + 11*t
                 END IF
               CASE (202122) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(33*y)
                   exa(2) = -Icing_Lambda1*(33*x+42*y)
                   exa(3) = 4 + 33*x*y + 21*y**2 + 11*t
                 ELSE
                   exa(1) = -Icing_Lambda2*(2*x)
                   exa(2) = 0
                   exa(3) = 77 + x**2 + 2*t
                 END IF
               CASE (202123) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(33*y)
                   exa(2) = -Icing_Lambda1*(33*x+42*y)
                   exa(3) = 4 + 33*x*y + 21*y**2 + 11*t + 6*t**2
                 ELSE
                   exa(1) = -Icing_Lambda2*(33*y)
                   exa(2) = -Icing_Lambda2*(33*x+42*y)
                   exa(3) = 4 + 33*x*y + 21*y**2 + 11*t + 6*t**2
                 END IF
               CASE (202124) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(33*y)
                   exa(2) = -Icing_Lambda1*(33*x+42*y)
                   exa(3) = 4 + 33*x*y + 21*y**2 + 11*t + 6*t**2
                 ELSE
                   exa(1) = -Icing_Lambda2*(2*x)
                   exa(2) = 0
                   exa(3) = 77 + x**2 + 2*t - 8*t**2
                 END IF
               CASE (202125) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(33*y)
                   exa(2) = -Icing_Lambda1*(33*x+42*y)
                   exa(3) =4 + 33*x*y + 21*y**2 + 11*t + 6*t**2 -14*t**3
                 ELSE
                   exa(1) = -Icing_Lambda2*(33*y)
                   exa(2) = -Icing_Lambda2*(33*x+42*y)
                   exa(3) =4 + 33*x*y + 21*y**2 + 11*t + 6*t**2 - 14*t**3
                 END IF
               CASE (202126) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(33*y)
                   exa(2) = -Icing_Lambda1*(33*x+42*y)
                   exa(3) =4 + 33*x*y + 21*y**2 + 11*t + 6*t**2 - 14*t**3
                 ELSE
                   exa(1) = -Icing_Lambda2*(2*x)
                   exa(2) = 0
                   exa(3) = 77 + x**2 + 2*t - 8*t**2 - t**3
                 END IF
               CASE (202127) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(7*y)
                   exa(2) = -Icing_Lambda1*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t
                 ELSE
                   exa(1) = -Icing_Lambda2*(7*y)
                   exa(2) = -Icing_Lambda2*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t
                 END IF
               CASE (202128) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(7*y)
                   exa(2) = -Icing_Lambda1*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t
                 ELSE
                   exa(1) = -Icing_Lambda2*(3*x**2+4*x*y)
                   exa(2) = -Icing_Lambda2*(3*y**2+2*x**2)
                   exa(3) = x**3 + 2*x**2*y+ y**3 - 10 * t
                 END IF
               CASE (202129) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(7*y)
                   exa(2) = -Icing_Lambda1*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t -3*t**2
                 ELSE
                   exa(1) = -Icing_Lambda2*(7*y)
                   exa(2) = -Icing_Lambda2*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + t**5
                 END IF
               CASE (202130) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(7*y)
                   exa(2) = -Icing_Lambda1*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t -3*t**2
                 ELSE
                   exa(1) = -Icing_Lambda2*(3*x**2+4*x*y)
                   exa(2) = -Icing_Lambda2*(3*y**2+2*x**2)
                   exa(3) = x**3 + 2*x**2*y+ y**3 - 10 * t
                 END IF
               CASE (202131) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(7*y)
                   exa(2) = -Icing_Lambda1*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t -3*t**2 + 12*t**3
                 ELSE
                   exa(1) = -Icing_Lambda2*(7*y)
                   exa(2) = -Icing_Lambda2*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t -3*t**2 + 12*t**3
                 END IF
               CASE (202132) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(7*y)
                   exa(2) = -Icing_Lambda1*(7*x+21*y**2)
                   exa(3) = 7 + 7*x*y + 7*y**3 + 7*t -3*t**2
                 ELSE
                   exa(1) = -Icing_Lambda2*(3*x**2+4*x*y)
                   exa(2) = -Icing_Lambda2*(3*y**2+2*x**2)
                   exa(3) = x**3 + 2*x**2*y+ y**3 - 10 * t+ 12*t**3
                 END IF
               CASE (202134) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(8 + 7*y + 12*t +63*x**2)
                   exa(2) = -Icing_Lambda1*(6+7*x-3*t -45*y**2)
                   exa(3) = 14 + 8*x + 6*y - 11*t + 7*x*y + 12*x*t -3*y*t + &
                   21*x**3 -15*y**3 + 9*t**2
                 ELSE
                   exa(1) = -Icing_Lambda2*(33*y*t**2)
                   exa(2) = -Icing_Lambda2*(33*x*t**2)
                   exa(3) = 33*x*y*t**2
                 END IF
                 CASE (202135) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(24*x**2*y**2)
                   exa(2) = -Icing_Lambda1*(16*x**3*y+7*y**6)
                   exa(3) = 14 + 8*x**3*y**2 + t**5+y**7
                 ELSE
                   exa(1) = -Icing_Lambda2*(24*x**2*y**2)
                   exa(2) = -Icing_Lambda2*(16*x**3*y+7*y**6)
                   exa(3) = 14 + 8*x**3*y**2 + t**5+y**7
                 END IF
                 CASE (202136) !T polynome de degre 3 discontinue
                 IF ( ele%Ref == IcingZone1_Ref ) THEN
                   exa(1) = -Icing_Lambda1*(3*y*exp(3*x*y))
                   exa(2) = -Icing_Lambda1*(3*x*exp(3*x*y))
                   exa(3) = exp(3*x*y) + log(t+1)/2.d0 + t**3
                 ELSE
                   exa(1) = -Icing_Lambda2*(3*y*exp(3*x*y))
                   exa(2) = -Icing_Lambda2*(3*x*exp(3*x*y))
                   exa(3) = exp(3*x*y) + sqrt(t)*2.d0
                 END IF
       
    CASE DEFAULT
       PRINT*, " ** ERROR in ExactSolJumpt :: unknown test case"
       PRINT*, " "
       STOP
    END SELECT

  END SUBROUTINE ExactSolJumpt
  !==============================================!

END MODULE ExactFunctionJumpt
