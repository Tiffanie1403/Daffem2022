MODULE Locate

  USE Types
  USE Algebra
  USE Algebra2
  
  IMPLICIT NONE

CONTAINS

  !===============================================!
  SUBROUTINE LocateNode(X, mesh, nele, found)
  !===============================================!
  !partir de l'element ou l'on cherche son projeté 
  !utilise les coordonnées barycentrique
  
    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN)    :: X ! vecteur de dimen 2 du noeud quel'on cherche sur le maillage
    type(MeshType)      , INTENT(IN)    :: mesh !maillage
    INTEGER             , INTENT(INOUT) :: nele ! premiere element sur le quel on recherche le n,oeud
    LOGICAL             , INTENT(OUT)   :: found ! true si le noeud est trouvé
    !------------------------------------------
    type(Element) :: ele
    !--------------------
    REAL*8, DIMENSION(SIZE(X)+1) :: barCoor, B
    !-----------------------------------------------------
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: Checked
    INTEGER :: nNodesPerElement, nDim
    INTEGER :: i, oldNele, nele2
    !-----------------------------------

    nDim = SIZE(X) ; nNodesPerElement = nDim+1
    ALLOCATE(Checked(mesh%ne,nDim+1))
    checked = .FALSE.
    found = .TRUE.

    DO
        ele = mesh%ele(nele)
        CALL computeBarycentricCoordinates(Ele, barCoor, x)
        found = .TRUE.
        oldNele = nele
        DO i = 1, nDim + 1
            IF ( barCoor(i) < 0.d0 .AND. ABS(barCoor(i)) > 1.d-12 ) THEN
                found = .FALSE.
                IF ( mesh%ele(nele)%adj(i) /= -1 ) THEN
                   nele2 = mesh%ele(nele)%adj(i)/ele%nNodesPerelement
					IF ( .NOT. checked(nele2, i) ) THEN
					   nele = mesh%ele(nele)%adj(i)/ele%nNodesPerElement
					   Checked(nele,i) = .TRUE.
					   EXIT
                    END IF
                END IF
            END IF
        END DO
        IF ( found .eqv. .TRUE. ) THEN
          EXIT
        END IF
        IF ( oldNele == nele ) THEN
          found = .FALSE.
          EXIT
        END IF
    END DO
    DEALLOCATE(Checked)

  END SUBROUTINE LocateNode
  !====================================!

  !========================================================!
  SUBROUTINE ComputeBarycentricCoordinates(ele, BarCoor, X)
  !========================================================!

    IMPLICIT NONE
    type(element)       , INTENT(IN)  :: ele
    REAL*8, DIMENSION(:), INTENT(OUT) :: BarCoor
    REAL*8, DIMENSION(:), INTENT(IN)  :: X
    !--------------------------------------------
    REAL*8, DIMENSION(SIZE(X)+1,SIZE(X)+1) :: A, Am1
    !---------------------------------------------------
    REAL*8, DIMENSION(SIZE(X)+1) :: B
    !---------------------------------
    INTEGER :: i, j, nDim
    !---------------------
    
    nDim = SIZE(X)
    DO i = 1, nDim
       DO j = 1, ele%nNodesPerElement
          A(i,j) = ele%Coor(j,i)
       END DO
    END DO
    A(nDim+1,:) = 1
    B(1:nDim) = X
    B(nDim+1) = 1
    
   
    !CALL InverseMatrix_maxPivot(A, Am1)
    CALL Inverse3x3Matrix(A, Am1)
    CALL CheckInversion(A,Am1) !rajouter par moi
    barCoor = MATMUL(Am1,B)
    
  END SUBROUTINE ComputeBarycentricCoordinates
  !========================================================!

 
END MODULE Locate
