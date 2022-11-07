MODULE BCSelection

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !==================================================!
  SUBROUTINE SelectBC(Ele,  BCType, Data, iFace, iBC)
  !==================================================!

    IMPLICIT NONE
    type(element)      , INTENT(INOUT) :: ele
    CHARACTER(len=20)  , INTENT(OUT)   :: BCType
    type(DataStructure), INTENT(IN)    :: Data
    INTEGER            , INTENT(OUT)   :: iBC
    INTEGER            , INTENT(IN)    :: iFace
    !---------------------------------------------
    INTEGER :: i
    !------------

    BCType = "Error"
    DO i = 1, Data%N_BC
       IF ( ele%EdgRef(iFace) == data%BC_ref(i) ) THEN
          BCType = TRIM(ADJUSTL(Data%BC(i)))
          iBC = i
          EXIT
       END IF
    END DO

  END SUBROUTINE SelectBC
  !=================================================!

  !=============================================!
  SUBROUTINE SelectBC2(Edg,  BCType, Data, iBC)
  !=============================================!

    IMPLICIT NONE
    type(Edge)         , INTENT(IN) :: edg
    CHARACTER(len=20)  , INTENT(OUT)   :: BCType ! resort le type de condition associe a une arete de bord
    type(DataStructure), INTENT(IN)    :: Data
    INTEGER            , INTENT(OUT)   :: iBC
    !---------------------------------------------
    INTEGER :: i
    !------------

    BCType = "Error" ; iBC = 5 ! cas ou on a une arete d'interface 
    DO i = 1, Data%N_BC ! Data%N_BC nb de condition sur le bord, ne prend pas en compte l'interface
       IF ( edg%Ref == data%BC_ref(i) ) THEN
          BCType = TRIM(ADJUSTL(Data%BC(i)))
          iBC = i ! numero associe a la condition
          EXIT
       END IF
    END DO

  END SUBROUTINE SelectBC2
  !=================================================!


END MODULE BCSelection
