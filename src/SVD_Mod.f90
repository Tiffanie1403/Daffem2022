MODULE SVD_Mod
! module for SVD related functions
  IMPLICIT NONE

CONTAINS


  !==================================!
  SUBROUTINE Puissance(A,X,lambda,eps,n)
  !==================================!

  !*************************************************************************
  !   Power method used to find the highest eigen value of a matrix A 
  !   and its associated eigen vector 
  !		Remark : A is a matrix of the form M^t M homélie
  !*************************************************************************



    IMPLICIT NONE
    REAL*8    , DIMENSION(:,:), INTENT(IN)    :: A ! matrix where the eigen info are needed 
    REAL*8    , DIMENSION (:) , INTENT(INOUT) :: X ! X eigen vector , \Y vector such as || y||_2 = L the eigen value 
    REAL*8                    , INTENT(INOUT) :: eps ! tolerance pour la convergence 
    REAL*8                    , INTENT(OUT)   :: lambda ! tolerance pour la convergence 
    INTEGER					  , INTENT(IN)    :: n ! size of the system  
    !------------------------------------------------------
	INTEGER :: i , iterMax , cmpt, j
	REAL*8 :: res, normY
	REAL*8, DIMENSION (:) , ALLOCATABLE :: X_prec, Y  ! X eigen vector , \Y vector such as || y||_2 = L the eigen value 

	!------------------------------------------------------


	ALLOCATE(X_prec(n), Y(n)) ; Y = 0.d0  ;  X = 0.d0 
	
	!definition d'un vecteur de départ aléatoire 
	DO i = 1, n
		X(i) = rand(0)*1
	END DO 
	
	iterMAX = 10000 ! nombre d'iteration max 
	cmpt = 0 ! variable pour le comptage dans le tant que 
	res = 10000 ! valeur de la norme entre le vecteur x entre deux iterations
	eps= 0.00000000000000000001 ! precision voulu pour la convergence 
	
	DO WHILE ( cmpt < iterMAX .AND. res > eps) 
		cmpt = cmpt + 1 ! valeur de comptage de la boucle 
		!Y = AX produit matrice vecteur 
		X_prec = X ! valeu du vecteur au temps precedent 
		Y = 0.d0 
		DO i= 1, n 
			DO j=1,n 
				Y(i) = Y(i) + X(j)*A(i,j)
			END DO
		END DO 
		
		normY = 0.d0 
		DO i =1,n
			normY = normY + Y(i)**2
		END DO 
		normY = sqrt(normY) ! on prend la norme et pas son carré 
		
		X = Y/normY 
		
		res = 0
		Do i =1, n 
			res = res + (X(i)-X_prec(i))**2
		END DO
		res = sqrt(res) 
	END DO 

	lambda = normY
	!PRINT*,"Voici le reisidu de fin ", res , cmpt 

  END SUBROUTINE Puissance
  !==================================!
  
  !Permet de calculer les autres valeurs propres et vecteurs propres en utilisant puissance sur une suite de matrice 
  !==================================!
  SUBROUTINE Deflation(M,U,S,V,mn,nn)
  !==================================!

  !*************************************************************************
  !   Take in entries a Matrix to which the SVD is wanted
  !		E_val is the vector of eigen value in a decreasing order, V the orthonormal matrix built with the eigen vector from M^t M 
  !*************************************************************************


    IMPLICIT NONE
    REAL*8    , DIMENSION(:,:), INTENT(IN)    :: M 
    REAL*8    , DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: V , S , U 
    INTEGER					  , INTENT(IN)    :: nn, mn ! size of the system  which can be rectangular m * n 
    !------------------------------------------------------
	INTEGER :: i, j , k 
	REAL*8, DIMENSION (:,:) , ALLOCATABLE :: MtM, C , UU ! matrice pour le calcul de M^t * M 
	REAL*8, DIMENSION (:)   , ALLOCATABLE :: X ! matrice pour le calcul de M^t * M 
	REAL*8 :: eps , lambda , res 
	!------------------------------------------------------
	
	! definition des matrices de la decomposition SVD 
	ALLOCATE(V(nn,nn),U(mn,mn)) ; ALLOCATE(S(mn,nn))
	V = 0.d0 ; S = 0.d0 ; U =0.d0
	
	
	!Deux cas se presente en fonction de la taille de la matrice M 
	
	IF(nn>=mn) THEN ! on va reconstruire V par deflation 
		ALLOCATE(MtM(nn,nn)) ; MtM = 0.d0 
		
		ALLOCATE(C(nn,nn), UU(nn,nn))
		C = 0.d0 ; UU = 0.d0 
		
		! CAlcul de M^t * M 
		Do i = 1 , nn 
			Do j=1, nn 
				DO k =1, mn
					MtM(i,j) = MtM(i,j) + M(k,i)*M(k,j)
				END DO 
			END DO 
		END DO 
		C = MtM 
		
		!on commence par calculer la premiere valeur propre et le vecteur propre associé 
		!valeur de tolerance pour la fonction puissance 
		eps = 0.000000000001 ; ALLOCATE(X(nn)) 
		CALL Puissance(C,X,lambda,eps,nn)
		V(:,1) = X(:) ; S(1,1) = lambda
		! on commence par faire le calcul de M^t M 
		
		DO i = 2,nn

			
			!calcul de la norme de V(:,i-1)
			res = 0 
			DO j = 1 ,nn 
				res = res + V(j,i-1)**2
			END DO 
			
			UU = 0.d0
			Do j = 1 ,nn 
				DO k=1, nn 
					UU(j,k) = V(j,i-1)*V(k,i-1)
				END DO 
			END DO 
			
			C = C - (S(i-1,i-1)/res)*UU
			X = 0 ; lambda = 0 
			CALL Puissance(C,X,lambda,eps,nn)
			V(:,i) = X 
			If(i<=mn) THEN 
				S(i,i) = lambda
			END IF 
		END DO 
		
		dO i=1, mn 
			S(i,i) = sqrt(S(i,i))
		END DO 
	ELSE 
		!on va sur le calcul de U est non de V 
		ALLOCATE(MtM(mn,mn)) ; MtM = 0.d0 ! calcul de M M^ t 
		
		
		ALLOCATE(C(mn,mn), UU(mn,mn))
		C = 0.d0 ; UU = 0.d0 

		Do i = 1 ,mn 
			Do j=1 ,mn 
				DO k =1, nn
					MtM(i,j) = MtM(i,j) + M(i,k)*M(j,k)
				END DO 
			END DO 
		END DO 
		C = MtM 

		!valeur de tolerance pour la fonction puissance 
		eps = 0.00000000001 ; ALLOCATE(X(mn)) 
		CALL Puissance(C,X,lambda,eps,mn)
		U(:,1) = X ; S(1,1) = lambda 
		! on commence par faire le calcul de M^t M 
		
		DO i = 2,mn
			
			!calcul de la norme de U(:,i-1)
			res = 0 
			DO j = 1 ,mn 
				res = res + U(j,i-1)**2
			END DO 
			
			UU = 0.d0
			Do j = 1 ,mn 
				DO k=1, mn 
					UU(j,k) = U(j,i-1)*U(k,i-1)
				END DO 
			END DO 
			
			C = C - (S(i-1,i-1)/res)*UU
			X = 0 ; lambda = 0 
			CALL Puissance(C,X,lambda,eps,mn)
			U(:,i) = X 
			IF(i<= nn) THEN 
				S(i,i) = lambda
			END IF 
		END DO
		
		dO i=1, nn 
			S(i,i) = sqrt(S(i,i))
		END DO 
	END IF 

	
	
	

  END SUBROUTINE Deflation
  !==================================!

  
  !==================================!
  SUBROUTINE SVD_facto(A,U,S,V,m,n)
  !==================================!
  
  
    !*************************************************************************
  !   Take in entries a Matrix to which the SVD is wanted
  !		E_val is the vector of eigen value in a decreasing order, V the orthonormal matrix built with the eigen vector from M^t M 
  !*************************************************************************

    IMPLICIT NONE
    REAL*8    , DIMENSION(:,:), INTENT(IN)    :: A
    REAL*8    , DIMENSION(:,:), ALLOCATABLE, INTENT(InOUT) :: V , S , U 
    INTEGER					  , INTENT(IN)    :: n, m ! size of the system  which can be rectangular m * n 
    !------------------------------------------------------
	INTEGER :: i, j , k 
	REAL*8 :: eps , lambda , res 
	!------------------------------------------------------
	REAL*8    , DIMENSION(:,:), ALLOCATABLE    :: A_bis, A_inter 
	REAL*8    , DIMENSION(:), ALLOCATABLE      :: NORM_svd 
	!------------------------------------------------------
	
	
	
	! definition de la matrice S diagonale et de la matrice V des vecteurs propre de A^t A 
	CALL Deflation(A,U,S,V,m,n)
	
	!Il nous reste à faire le calcul de U ou de V en fonction de la taille de A 
	
	
	IF (m<=n) THEN  ! on a reconstruit V par deflation 
		! il nous reste donc à reconstruire U 
		U = 0.d0 
		DO k=1, m
			Do i =1, m 
				Do j =1, n 
					U(i,k) = U(i,k) + A(i,j)*V(j,k)
				END DO 
			END DO 
		END DO
	
		!il reste à faire la normalisation pour obtenir un vecteur de norme 1 
		Do k =1, m
			res = 0 
			Do i = 1, m
				res = res + U(i,k)**2
			END DO 
			res = sqrt(res) 
			U(:,k) = U(:,k)/res 
		END DO 

	ELSE ! n < m 
	! On connait la matrice U par déflation on reconstruit donc la matrice V 
		V = 0.d0
		DO k=1, n
			Do i =1, n
				Do j =1, m
					V(i,k) = V(i,k) + A(j,i)*U(j,k)
				END DO 
			END DO 
		END DO
		
		!on normalise les vecteurs formants la matri ce V 
		Do k =1, n
			res = 0 
			Do i = 1 , n
				res = res + V(i,k)**2
			END DO 
			res = sqrt(res)
			V(:,k) = V(:,k)/res 
		END DO 
	END IF 
	! on a donc defini les trois matrice définissant la svd de A 
	
	! on va faire un test pour verifier que l'on a bien reconstruit la bonne factorisation 
	! Calcul de A avec la svd
	ALLOCATE(A_bis(m,n),A_inter(m,n) )
	
	! on commence par le produit matricielle U*S
	A_bis = 0.d0 ; a_inter = 0.d0
	Do i = 1, m 
		DO j = 1, n 
			Do k = 1, m 
				IF(k==j) THEN ! comme S est diagonale on evite des operations inutiles 
					A_inter(i,j) = A_inter (i,j) + U(i,k)*S(k,j)
				END IF 
			END DO 
		END DO 
	END DO 
	
	Do i = 1, m 
		DO j = 1, n 
			Do k = 1, n
				A_bis(i,j) = A_bis(i,j) + A_inter(i,k)*V(j,k) 
			END DO 
		END DO 
	END DO 
	
	ALLOCATE(norm_SVD(n)) ; norm_SVD = 0.d0
	
	DO i = 1 ,n 
		DO j = 1, m
			norm_SVD(i) = norm_SVD(i) + (A(j,i)-A_bis(j,i))**2
		END DO 
		norm_SVD(i) = sqrt(norm_SVD(i)) 
	END DO 
	
	!PRINT*,"residu SVD", norm_SVD
	
	DEALLOCATE(A_inter, A_bis) ; DEALLOCATE(norm_SVD)
	
	
	
	
    END SUBROUTINE SVD_facto
  !==================================!
  
  !routine de calcul de la pseudo inverse 
  
    !==================================!
  SUBROUTINE Pseudo_Inverse(A,Ac,m,n)
  !==================================!
  
  
  !*************************************************************************
  !   Take in entries a Matrix to which the SVD is wanted
  !		E_val is the vector of eigen value in a decreasing order, V the orthonormal matrix built with the eigen vector from M^t M 
  !*************************************************************************
    IMPLICIT NONE
    REAL*8    , DIMENSION(:,:), INTENT(IN)    			 :: A ! matrice dont on veut calculer la pseudo inverse 
    REAL*8    , DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Ac ! pseudo inverse de A 
    INTEGER					  , INTENT(IN)    			 :: n, m ! size of the system  which can be rectangular m * n 
    !------------------------------------------------------
	INTEGER :: i, j , k , minleg
	REAL*8    , DIMENSION(:,:), ALLOCATABLE   :: V , S , U , A_inter 
	REAL*8    , DIMENSION(:,:), ALLOCATABLE   :: St ! pseudo inverse de S 


	!------------------------------------------------------
	
	
	CALL  SVD_facto(A,U,S,V,m,n)
	
	! on a calculer precedement la decomposition en valeur singuliere 
	ALLOCATE(St(n,m)) ; St = 0.d0 
	IF(n<m) THEN 
		minleg = n 
	ELSE 
		minleg = m
	END IF 
	Do i = 1, minleg 
		St(i,i) = 1.d0/S(i,i) 
	END DO
	ALLOCATE(A_inter(n,m)) ; A_inter = 0.d0 
	
	! Operation mat-mat V * S pseudo inverse noté St 
	Do i =1, n 
		DO j = 1,m 
			Do k = 1, n
				A_inter(i,j) = A_inter(i,j) + V(i,k)*St(k,j)
			END DO 
		END DO 
	END DO 
	
	!Operation matrice-matrice entre A_inter ( = V * S^pseudo) et U transposé
	ALLOCATE(Ac(n,m)) ; Ac = 0.d0 
	Do i =1 , n 
		Do j = 1, m 
			do K =1, m 
				Ac(i,j) = Ac(i,j) + A_inter(i,k)*U(j,k)
			end do 
		END DO 
	END DO 
	
	DEALLOCATE(A_inter)
  
    END SUBROUTINE Pseudo_Inverse
  !==================================!
  
  

END MODULE SVD_Mod
