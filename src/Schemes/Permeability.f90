MODULE Permeability

  USE Types

  IMPLICIT NONE

CONTAINS

  !===========================================================!
  SUBROUTINE SetNodalInverseMatrixPermeability(Lm1, ele, icas)
  !===========================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: Lm1
    INTEGER               , INTENT(IN)    :: icas
    type(element)         , INTENT(IN)    :: ele
    !---------------------------------------------
    INTEGER :: id, jd
    !------------------

    DO id = 1, SIZE(Lm1,1)
       DO jd = 1, SIZE(Lm1,1)
          Lm1(id,jd) = ele%Lm1(id,jd)
       END DO
    END DO

  END SUBROUTINE SetNodalInverseMatrixPermeability
  !=============================================================!

  !==================================================!
  SUBROUTINE SetNodalMatrixPermeability(L, ele, icas)
  !==================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: L
    INTEGER               , INTENT(IN)    :: icas
    type(element)         , INTENT(IN)    :: ele
    !---------------------------------------------
    INTEGER :: id, jd
    !------------------

    DO id = 1, SIZE(L,1)
       DO jd = 1, SIZE(L,1)
          L(id,jd) = ele%L(id,jd)
       END DO
    END DO

  END SUBROUTINE SetNodalMatrixPermeability
  !==============================================!

END MODULE Permeability
