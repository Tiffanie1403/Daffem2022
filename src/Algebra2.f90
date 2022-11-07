MODULE Algebra2
  
  IMPLICIT NONE
  
CONTAINS
  
  !================================!
  SUBROUTINE InverseMatrix2(A, Am1)
  !================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: A
    REAL*8, DIMENSION(:,:), INTENT(OUT)   :: Am1
    !-------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: A2, Id
    !----------------------------------------------
    REAL*8 :: aii, f
    !----------------------------
    INTEGER :: i, j, k, N,m
    !-------------------------

    N = SIZE(A,1)
    ALLOCATE(A2(N,N),Id(N,N))
    Id = 0.
    DO i = 1, N
       Id(i,i) = 1.
    END DO
    A2 = A

    DO i = 1, N
       aii = A2(i,i)
       DO j = i+1, N
          f = A2(j,i)/aii
          DO k = 1, N
             Id(j,k) = Id(j,k) - f*Id(i,k)
             A2(j,k) = A2(j,k) - f*A2(i,k)          ! A2(i,:) ligne pivot    
          END DO
          IF ( ABS(A2(j,j)) < 1.d-10  ) THEN
             PRINT*, "ERROR in InverseMatrix, division by 0 !!!!"
             STOP
          END IF
          Id(j,:) = Id(j,:)/A2(j,j) ! pour definir un 1 sur la diagonale 
          A2(j,:) = A2(j,:)/A2(j,j)
       END DO
    END DO
    
    
    IF ( ABS(A2(1,1)) < 1.d-10  ) THEN
       PRINT*, "ERROR in InverseMatrix, division by 0 !!!!"
       STOP
    END IF
    Id(1,:) = Id(1,:)/A2(1,1)
    A2(1,:) = A2(1,:)/A2(1,1)

    DO i = N, 1, -1
       aii = A2(i,i)
       DO j = i-1, 1, -1
          IF ( ABS(aii) < 1.d-10 ) THEN
             PRINT*, "ERROR in InverseMatrix, division by 0 !!!!"
             STOP
          END IF
          f = A2(j,i)/aii

          DO k = N, 1, -1
             Id(j,k) = Id(j,k) - f*Id(i,k)
             A2(j,k) = A2(j,k) - f*A2(i,k)
          END DO
       END DO
    END DO

    Am1 = Id

    DEALLOCATE(A2, Id)

  END SUBROUTINE InverseMatrix2
  !=====================================!
  
   !================================!
  SUBROUTINE InverseMatrix_maxPivot(A, Am1)
  !================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN)    :: A
    REAL*8, DIMENSION(:,:), INTENT(OUT)   :: Am1
    !-------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: A2, Ident
    REAL*8, DIMENSION(:), ALLOCATABLE   :: inter 

    !----------------------------------------------
    REAL*8 :: aii, f, max_pivot 
    !----------------------------
    INTEGER :: i, j, k, N,m, ligne_max 
    !-------------------------

   N = SIZE(A,1) ! taille de la matrice 
    ! allocation 
    ALLOCATE(A2(N,N),Ident(N,N),inter(N))
    
    !creation de la matrice identité 
    Ident = 0.d0
    DO i = 1, N
       Ident(i,i) = 1.
    END DO
    
    !on copie la matrice que l'on veut inverser dans A2 
    A2 = A !permet de ne pas toucher à A 
    
    !vecteur de sauvegarde pour l'échange des lignes
    inter = 0.d0 
    
    DO i = 1, N
		!recherche du max sur la colonne i 
		max_pivot = A2(i,i) ; ligne_max=i
		IF(i/=N) THEN ! si on ne se trouve pas sur la derniére ligne 
			Do j=i+1,N 
				IF(abs(A2(j,i)) > max_pivot) THEN 
				!on enregistre la valeur pivot et la ligne associé 
					max_pivot = A2(j,i) ; ligne_max  = j
				END IF  
			END DO 
		END IF 
		
		
		IF(MAX_PIVOT ==0) THEN	 
		! ON PASSE AU i suivant 
		PRINT*,"PB max pivot nulle",MAX_PIVOT
		STOP
		!on ne fait rien 
		ELSE ! pivot non nul 
		
			IF(ligne_max/=i) THEN 
				! on doit échanger les lignes 
				inter= A2(i,:) ! on sauvegarde la ligne i 
				A2(i,:) = A2(ligne_max,:)
				A2(ligne_max,:) = inter 
				
				inter= Ident(i,:) ! on sauvegarde la ligne i 
				Ident(i,:) = Ident(ligne_max,:)
				Ident(ligne_max,:) = inter 
			END IF ! on a fini l'échange de ligne 
			
			!on ramene la valeur pivot à 1 sur la ligne pivot i 
			aii = A2(i,i) ! valeur du pivot 
			A2(i,:) = A2(i,:)/aii
			Ident(i,:) = Ident(i,:)/aii 
			
			!operation sur les lignes inférieur 
			If(i/=N) THEN 
				Do j =i+1,N 
					Ident(j,:) = Ident(j,:)-A2(j,i)*Ident(i,:) 
					A2(j,:) = A2(j,:)-A2(j,i)*A2(i,:) 
				END DO 
			END IF 

			!Operation sur les lignes supérieur 
			If(i > 1) THEN 
				DO k = 1,i-1
					Ident(k,:) = Ident(k,:)-A2(k,i)*Ident(i,:) 
					A2(k,:) = A2(k,:)-A2(k,i)*A2(i,:) 
				END DO 
			END IF 
		END IF 
	END DO

	Am1 = Ident ! inverse de A 
!	A2=0.d0
!	!on va verifier que le produit des deux fait bien l'identité 
!	Do i = 1, N 
!		Do j=1,n 
!			DO k=1,N
!				A2(i,j) = A2(i,j) + A(i,k)*Am1(k,j)
!			END DO 
!		END DO 
!	END DO 
	!PRINT*,"test matrice inverse"
	!DO i=1,n 
	!	PRINT*,A2(i,:)
	!END DO 
	
    DEALLOCATE(A2, Ident)

  END SUBROUTINE InverseMatrix_maxPivot
  !=====================================!

  !================================!
  SUBROUTINE CheckInversion2(M, Mm1)
  !================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(IN) :: M, Mm1
    !-----------------------------------------------
    REAL*8, DIMENSION(SIZE(M,1),SIZE(M,1)) :: Id
    !--------------------------------------
    INTEGER :: i, j
    !----------------

    Id = MATMUL(M,Mm1)
    
    DO i = 1, SIZE(M,1)

       IF (ABS(Id(i,i) - 1 ) > 0.000001 ) THEN
          PRINT*, 'Error in matrix inversion'
          PRINT*, "On Diagonal :: ", i, i, Id(i,i)
          DO j = 1, SIZE(M,1)
             PRINT*, Id(j,:)
          END DO
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
    
  END SUBROUTINE CheckInversion2
  !================================!
 
  
END MODULE Algebra2
