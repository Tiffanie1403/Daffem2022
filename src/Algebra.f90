MODULE Algebra
! module for matrix operations
  IMPLICIT NONE

CONTAINS

  !================================!
  SUBROUTINE InverseMatrix(A, Am1)
  !================================!

    !*************************************************************************
    !   Inversion of the matrix A return in the variable Am1
    !
    !  Parameters:
    !
    !    Input,Output, REAL*8, DIMENSION(:,:) A, matrix we want to inverse
    !    Output, REAL*8, DIMENSION(:,:) Am1, inverse of the matrix A
    !
    !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: A
    REAL*8, DIMENSION(:,:), INTENT(OUT)   :: Am1
    !-------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: A2, Id
    !----------------------------------------------
    REAL*8 :: aii, f
    !----------------------------
    INTEGER :: i, j, k, N
    !-------------------------

    N = SIZE(A,1)
    ALLOCATE(A2(N,N),Id(N,N))
    Id = 0.
    !creation d'une matrice identite
    DO i = 1, N
       Id(i,i) = 1.
    END DO
    A2 = A

    DO i = 1, N
       aii = A2(i,i) ! pivot
       DO j = i+1, N
          f = A2(j,i)/aii
          DO k = 1, N
             Id(j,k) = Id(j,k) - f*Id(i,k)
             A2(j,k) = A2(j,k) - f*A2(i,k)
          END DO
          Id(j,:) = Id(j,:)/A2(j,j)
          A2(j,:) = A2(j,:)/A2(j,j)
       END DO
    END DO
    Id(1,:) = Id(1,:)/A2(1,1)
    A2(1,:) = A2(1,:)/A2(1,1)

    DO i = N, 1, -1
       aii = A2(i,i)
       DO j = i-1, 1, -1
          f = A2(j,i)/aii

          DO k = N, 1, -1
             Id(j,k) = Id(j,k) - f*Id(i,k)
             A2(j,k) = A2(j,k) - f*A2(i,k)
          END DO
       END DO
    END DO

    Am1 = Id

    DEALLOCATE(A2, Id)

  END SUBROUTINE InverseMatrix
  !=====================================!


   !================================!
   SUBROUTINE GradCong(A,P,b,X,eps,leng)
   !================================!

     IMPLICIT NONE
     REAL*8, DIMENSION(:), INTENT(INOUT) :: X
     REAL*8, DIMENSION(:,:), INTENT(IN)  :: A, P
     REAL*8, DIMENSION(:), INTENT(IN)    :: b
     INTEGER,              INTENT(IN)    :: leng
     REAL*8,               INTENT(IN)    :: eps
     !-------------------------------------------
     REAL*8, DIMENSION(:), ALLOCATABLE :: G, H,qk, inter, interG ,Cgi,Cgip
     !----------------------------------------------
     REAL*8 :: rho, alpha, beta, inter1, inter2, ggamma, inter3 
     !----------------------------
     INTEGER :: k,i, iterMAx, j
     LOGICAL :: test
     !-------------------------

     ALLOCATE(G(leng), H(leng), inter(leng), interG(leng), Cgi(leng))
     ALLOCATE(CGip(leng))
     interG = 0.d0
     ! initialisation
     !Calcul de AX
     G = 0.d0
      Do i=1,leng
        inter1 = 0
        Do j =1,leng
          inter1 = inter1 + A(i,j)*X(j)
        END DO
        G(i) = inter1 - b(i)
      END DO

      H = 0
      DO i=1,leng
        DO j=1,leng
          H(i)=H(i)+ P(i,j)*G(j)
        END DO
      END DO
      H = - H

      TEST = .FALSE. ; k = 0 ; iterMAx = 100
      DO WHILE (k < leng .AND. .NOT. Test )

        k=k+1
        inter = 0.d0
        !Calcul de AH
        Do i= 1, leng
          DO j=1,leng
            inter(i)= inter(i)+A(i,j)*H(j)
          END DO
        END DO

        inter1 = 0.d0 ; inter2=0.d0
        DO i= 1,leng
          inter1 = inter1 + G(i)*H(i)
          inter2 = inter2 + H(i)*inter(i)
        END DO
        rho = - inter1/inter2


        X = X+rho*H

        interG = G ! on enregistre G on temps precedent
        G = G + rho*inter

        !Calcul de PGk+1 et PGk
        Cgip = 0.d0 ; Cgi = 0.d0
        DO i=1,leng
          DO j=1,leng
            Cgip(i)=Cgip(i) + P(i,j)*G(j)
            Cgi(i) = Cgi(i) + P(i,j)*inter(j)
          END DO
        END DO

        inter1 = 0.d0 ; inter2=0.d0 ; inter3 =0 
        DO i=1,leng
          inter1= inter1+ Cgip(i)*G(i)
          inter2= inter2+ Cgi(i)*inter(i)
          inter3 = inter3 + G(i)*G(i)
        END DO
        ggamma = inter1/inter2
        H=-Cgip+ggamma*H

        IF(inter1 < eps) THEN ! on a la tolerance verifié
    ! IF(sqrt(inter3) < eps) THEN
          TEST = .TRUE.
        END IF
      END DO


   END SUBROUTINE GradCong
   !=====================================!


! methode sans preconditionnement 
      !================================!
      SUBROUTINE GradCong2(A,b,X,eps,leng)
      !================================!

        IMPLICIT NONE
        REAL*8, DIMENSION(:), INTENT(INOUT) :: X ! initial guess at first 
        REAL*8, DIMENSION(:,:), INTENT(IN)  :: A
        REAL*8, DIMENSION(:), INTENT(IN)    :: b
        INTEGER,              INTENT(IN)    :: leng
        REAL*8,               INTENT(INOUT) :: eps
        !-------------------------------------------
        REAL*8, DIMENSION(:), ALLOCATABLE :: R0, P0, Qk,Tvef
        !----------------------------------------------
        REAL*8 :: res, rho, alpha, beta, test, init_rho, test_residu
        !----------------------------
        INTEGER :: iterMAx, k , i, j
        !-------------------------
        !initialisation
        ALLOCATE(R0(leng), P0(leng), Qk(leng), Tvef(leng)) ; res = 0 ; alpha = 0
        itermax = 100 ; beta = 0 ; test= 10000 ; R0 = 0.d0 ; P0 = 0.d0 
        eps = 0.000000001 ; test_residu = 10000
		!X = rand(0)*1
		
        DO i=1,leng
          res = 0
          DO j=1,leng
            res = res +A(i,j)*X(j)
          END DO
          R0(i) =b(i)-res
        END DO

        P0 = R0
        rho = 0
        Do i=1, leng
          rho = rho + R0(i)**2
        END DO
		init_rho = rho
		
        !Iterations
        k = 1
         Do WHILE ( k < iterMAx .AND. test_residu > eps )
        !Do WHILE ( k < iterMAx .AND. rho > eps )
      !  PRINT*,"voici rho", test_residu
          Qk =0
          DO i = 1 ,leng
            DO j = 1, leng
              Qk(i) = Qk(i) + A(i,j)*P0(j)
            END DO
          END DO
          res = 0
          DO i =1 , leng
            res = res + P0(i)*Qk(i)
          END DO
          alpha = rho / res
          X = X + alpha*P0
          R0 = R0 - alpha*Qk
          res = rho ! on garde rho au tour precedent
          rho = 0
          DO i = 1,leng
            rho = rho + R0(i)**2
          END DO
          beta = rho / res
          P0 = R0 + beta*P0

		  k = k+1 ! compteur pour le nombre d'itération 
		  test_residu = sqrt(rho) ! /sqrt(init_rho)
		  ! Calcul du résidu du système 
		  
          !calcul du residu du systeme lineaire
      !    Do i = 1 ,leng
       !     res = 0
        !    do j=1,leng
         !     res = res + A(i,j)*X(j)
          !  END DO
    !        Tvef(i)=b(i)-res
     !     END DO
      !    test = 0
       !   do i =1, leng
        !    test = test + Tvef(i)**2
         ! END DO
         ! test = sqrt(test)
         ! PRINT*,"Voici test", test , "au tour", k

        END DO

      END SUBROUTINE GradCong2
      !=====================================!

      !================================!
      SUBROUTINE GradCong2bis(A,Am1,b,X,eps,leng)
      !================================!

        IMPLICIT NONE
        REAL*8, DIMENSION(:), INTENT(INOUT) :: X
        REAL*8, DIMENSION(:,:), INTENT(IN)  :: A, Am1
        REAL*8, DIMENSION(:), INTENT(IN)    :: b
        INTEGER,              INTENT(IN)    :: leng
        REAL*8,               INTENT(IN)    :: eps
        !-------------------------------------------
        REAL*8, DIMENSION(:), ALLOCATABLE :: R0, P0, Qk,Tvef
        REAL*8, DIMENSION(leng,leng)      :: Acop
        REAL*8, DIMENSION(leng)           :: Bcop
        !----------------------------------------------
        REAL*8 :: res, rho, alpha, beta, test
        !----------------------------
        INTEGER :: iterMAx, k , i, j
        !-------------------------
        !initialisation
        Acop = 0.d0
        DO i = 1 , leng
          DO j=1 ,leng
            DO k=1, leng
              Acop(i,j)= Acop(i,j)+ Am1(i,k)*A(k,j)
            END DO
          END DO
        END DO

        Bcop = 0.d0
        Do i=1, leng
          Do j=1, leng
            Bcop(i) = Bcop(i) + Am1(i,j)*B(j)
          END DO
        END DO

        ALLOCATE(R0(leng), P0(leng), Qk(leng), Tvef(leng)) ; res = 0 ; alpha = 0
        itermax = 10000 ; beta = 0 ; test= 10000

        DO i=1, leng
          res = 0
          DO j=1,leng
            res = res +Acop(i,j)*X(j)
          END DO
          R0(i) =bcop(i)-res
        END DO

        P0 = R0
        rho = 0
        Do i=1, leng
          rho = rho + R0(i)**2
        END DO

        !Iterations
        k = 1
        Do WHILE ( k < iterMAx .AND. test > eps )
          Qk =0
          DO i = 1 ,leng
            DO j = 1, leng
              Qk(i) = Qk(i) + Acop(i,j)*P0(j)
            END DO
          END DO
          res = 0
          DO i =1 , leng
            res = res + P0(i)*Qk(i)
          END DO
          alpha = rho / res
          X = X + alpha*P0
          R0 = R0 - alpha*Qk
          res = rho ! on garde rho au tour precedent
          rho = 0
          DO i = 1,leng
            rho = rho + R0(i)**2
          END DO
          beta = rho / res
          P0 = R0 + beta*P0

          !calcul du residu du systeme lineaire
          Do i = 1 ,leng
            res = 0
            do j=1,leng
              res = res + Acop(i,j)*X(j)
            END DO
            Tvef(i)=bcop(i)-res
          END DO
          test = 0
          do i =1, leng
            test = test + Tvef(i)**2
          END DO
      !    PRINT*,"Voici test", test , "au tour", k
          k = k+1

        END DO


      END SUBROUTINE GradCong2bis
      !=====================================!

  !================================!
  SUBROUTINE CheckInversion(M, Mm1)
  !================================!

  !*************************************************************************
  !   Checks if the product of the matrix M and Mm1 gives Id
  !   Gives the possibility to verify if a matrix is the inverse of another one
  !  Parameters:
  !
  !    Input,Output, REAL*8, DIMENSION(:,:) M,  Matrix of the same size than Mm1
  !    Input,Output, REAL*8, DIMENSION(:,:) Mm1, Matrix of the same size than M
  !
  !*************************************************************************

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: M, Mm1
    !-----------------------------------------------
    REAL*8, DIMENSION(SIZE(M,1),SIZE(M,1)) :: Id
    !--------------------------------------
    INTEGER :: i, j, k
    !----------------

    !Id = MATMUL(M,Mm1)
    Id = 0.d0
    Do i = 1, size(M,1)
		Do j =1, size(M,1)
			DO k= 1, size(m,1)
				Id(i,j) = Id(i,j) + M(i,k)*Mm1(k,j)
			END DO  
		END DO 
    END DO 
    DO i = 1, SIZE(M,1)

       IF (ABS(Id(i,i) - 1 ) > 0.000001 ) THEN
          PRINT*, 'Error in matrix inversion'
          PRINT*, "On Diagonal :: ", i, i, Id(i,i)
          STOP
       END IF
       DO j = 1, SIZE(M,1)
          IF (i /= j) THEN
             IF ( ABS(Id(i,j)) > 0.000001 ) THEN
                PRINT*, 'Error in matrix inversion'
                PRINT*, "On extra diagonal ", i, j, Id(i,j)
                STOP
             END IF
          END IF
       END DO

    END DO

  END SUBROUTINE CheckInversion
  !================================!

  !==========================!
  SUBROUTINE Transpose(A, AT)
  !==========================!

  !*************************************************************************
  !   Gives the transpose matrix of a matrix a return in At
  !
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(:,:) A,  Matrix
  !    Output, REAL*8, DIMENSION(:,:) AT, Matrix for the transposition of A
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: A
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: AT
    !----------------------------------------
    INTEGER :: i
    !------------

    DO i = 1, SIZE(A,1)
       AT(:,i) = A(i,:)
    END DO

  END SUBROUTINE Transpose
  !==========================!

  !==================================!
  SUBROUTINE Inverse3x3Matrix(A, Am1)
  !==================================!

  !*************************************************************************
  !   Gives the inverse of 3x3 matrix
  !
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(3,3) A,  Matrix
  !    Output, REAL*8, DIMENSION(3,3) Am1, Matrix for the inverse of A
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(3,3), INTENT(IN)  :: A
    REAL*8, DIMENSION(3,3), INTENT(OUT) :: Am1
    !------------------------------------------
    REAL*8 :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    REAL*8 :: det
    !-------------------------------------------------------

    a11 = A(1,1) ; a12 = A(1,2) ; a13 = A(1,3)
    a21 = A(2,1) ; a22 = A(2,2) ; a23 = A(2,3)
    a31 = A(3,1) ; a32 = A(3,2) ; a33 = A(3,3)

    det = a11*a22*a33 + a21*a32*a13 + a31*a12*a23 - a31*a22*a13 - a21*a12*a33 - a11*a32*a23

    IF ( abs(det) < 0.00000000001 ) THEN
       PRINT*, " ** ERROR In Inverse3x3Matrix :: Singular matrix !!!!"
       STOP
    END IF

    Am1 = 0.d0
    !Am1(1,1) = a22*a33 - a32*a23
    !Am1(1,2) = - ( a21*a33 - a31*a23 )
    !Am1(1,3) = a21*a32 - a31*a22
  
    !Am1(2,1) = - ( a12*a33 - a32*a13 )
  
    !Am1(2,2) = a11*a33 - a31*a13
    !Am1(2,3) = - ( a11*a32 - a31*a12 )
    !Am1(3,1) = a12*a23 - a22*a13
    !Am1(3,2) = - ( a11*a23 - a21*a13 )
    !Am1(3,3) = a11*a22 - a21*a12
    
    !we take the transpose of the matrix of the different commatrix 
    Am1(1,1) = a22*a33 - a32*a23
    
    Am1(2,1) = - ( a21*a33 - a31*a23 )
    
    Am1(3,1) = a21*a32 - a31*a22
    
    Am1(1,2) = - ( a12*a33 - a32*a13 )
    
    Am1(2,2) = a11*a33 - a31*a13
    Am1(3,2) = - ( a11*a32 - a31*a12 )
    Am1(1,3) = a12*a23 - a22*a13
    Am1(2,3) = - ( a11*a23 - a21*a13 )
    Am1(3,3) = a11*a22 - a21*a12
    
    
    Am1 = Am1/det

  END SUBROUTINE Inverse3x3Matrix
  !==================================!

  !====================================!
  SUBROUTINE ComputeVectorNorm2(V, nrm)
  !====================================!

  !*************************************************************************
  !   Calculates the norm of the vector V
  !
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(:) V, norm of the vector V
  !    Output, REAL*8 nrm, norm of the vector V
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)  :: V
    REAL*8              , INTENT(OUT) :: nrm
    !----------------------------------------
    INTEGER :: i, N
    !----------------

    N = size(V)
    nrm = 0.d0
    DO i = 1, N
       nrm = nrm + V(i)**2
    END DO
    nrm = sqrt(nrm/N)

  END SUBROUTINE ComputeVectorNorm2
  !====================================!


 !================================!
  SUBROUTINE QrFactorization(A,Q,R,m,n)
  !================================!

    !*************************************************************************
  
    !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN) 	:: A   !matrix of the original system independ of the normal equation 
    REAL*8, DIMENSION(:,:), INTENT(OUT)	:: Q,R !matrix for the QR factorization
    INTEGER,				INTENT(IN)	:: m,n !size of the matrix A 
    !-------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: A2,A2_inter,Q_inter,Hk
    !----------------------------------------------
    INTEGER :: i, j,l, k
    REAL*8  :: norm
    REAL*8, DIMENSION(:), ALLOCATABLE	:: ak,vk,ek
    LOGICAL ::	TEST
    !-------------------------

	ALLOCATE(A2(m,n),A2_inter(m,n),Q_inter(m,m),ak(m),vk(m),ek(m),Hk(m,m))
	
	A2 = A ; A2_inter = 0.d0 ; Q = 0.d0
	ek = 0.d0 ; ek(1)=1.d0 !vecteur de la base canonique  
	
	!Initialization with the identity matrix 
	DO i = 1, m
		Q(i,i) = 1.d0
	END DO
	
	DO k =1,n 
		ak = 0.d0
		ak(k:m) = A2(k:m,k)
		
		norm= 0.d0
		DO j = k, m
			norm = norm + ak(j)**2
		END DO 
		norm = sqrt(norm)
		
		vk = 0.d0 
		vk(k:m) = ak(k:m) - ek(1:m-k+1)*norm
		
		TEST = .TRUE.
		DO j = k, m
			IF(ak(j)/=0) THEN
				TEST=.FALSE.
			END IF 
		END DO 
		
		IF (.NOT. TEST) THEN 
			norm= 0.d0
			DO j = k, m
				norm = norm + vk(j)**2
			END DO 
			
			!definition of the matrix Hk 
			Hk = 0.d0
			DO i=k,m
				DO j=k,m
					Hk(i,j) = -2.d0*vk(i)*vk(j)/norm
					IF(i==j) THEN 
						Hk(i,j) = Hk(i,j) + 1
					END IF 
				END DO 
			END DO 
			DO i=1,k-1
				Hk(i,i) = 1.d0
			END DO 
			
			
			!Calcul de Ak+1 
			A2_inter = A2 ; A2 =0.d0
			Do i = 1 ,m
				DO j =1,n 
					DO l = 1, m
						A2(i,j) = A2(i,j) + Hk(i,l)*A2_inter(l,j)
					END DO 
				END DO 
			END DO 
		ELSE 
			Hk = 0.d0
			DO j=1,m
				Hk(j,j)=1.d0
			END DO
		END IF
		!Calcul de Q 
		Q_inter = Q ; Q = 0.d0
		DO i = 1,m 
			Do j=1,m 
				DO l=1,m
					Q(i,j) =Q(i,j) + Q_inter(i,l)*Hk(j,l)
				END DO 
			END DO 
		END DO 
	END DO 
	R = A2 
	!PRINT*,"Affichage de la factorisation QR"
	!!PRINT*,"Affichage de Q"
	!DO i =1,m 
	!	PRINT*,Q(i,:)
	!END DO
	!PRINT*,"Affichage de R"
	!DO i =1,m 
	!	PRINT*,R(i,:)
	!END DO 

  END SUBROUTINE QrFactorization
  !=====================================!

!================================!
  SUBROUTINE SolveSystemQR(A,b,x,m,n)
  !================================!

    !*************************************************************************
  
    !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN) 	:: A   !matrix of the original system independ of the normal equation 
    REAL*8, DIMENSION(:), INTENT(IN) 	:: b   !second member of the system
    REAL*8, DIMENSION(:), INTENT(OUT)	:: x   !solution vector
    INTEGER,				INTENT(IN)	:: m,n !size of the matrix A 
    !-------------------------------------------
    REAL*8, DIMENSION(:,:),ALLOCATABLE			:: Q,R !matrix for the QR factorization
    REAL*8, DIMENSION(:), ALLOCATABLE			:: Beta
    !----------------------------------------------
    INTEGER :: i, j,l, k
    !-------------------------
    
    	
	ALLOCATE(Q(m,m),R(m,n))
	CALL QrFactorization(A,Q,R,m,n)
	
	ALLOCATE(Beta(m)) ; Beta = 0.d0 
	
	!Value of Q^t* b 
	DO k = 1, m
		Do j = 1,m
			Beta(k) = Beta(k) + Q(j,k)*B(j)
		END DO 
	END DO 
	
	X = 0.d0 
	DO i=1,n
		k = n-i+1
		IF(k==n) THEN 
			X(k) = Beta(k)/R(k,k)
		ELSE 
			X(k) = Beta(k)
			Do j=k+1,n
				X(k) = X(k) - R(k,j)*X(j)
			END DO 
			X(k) = X(k)/R(k,k)
		END IF 
	END DO

	!PRINT*,"Value de x",X
	
  END SUBROUTINE SolveSystemQR
  !=====================================!


END MODULE Algebra
