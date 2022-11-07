MODULE Tools

  USE Types

  IMPLICIT NONE

CONTAINS

  !==============================================!
  SUBROUTINE GetFaceNumber(nDim, ele, edg, iFace)
  !==============================================!

  !*************************************************************************
  ! Gives the number corresponding to the edg in ele
  !
  !  Parameters:
  !
  !    Input, INTEGER  nDim, Dimension of the problem
  !    Input, type(Element) ele, Current element
  !    Input, type(Edge) edg, Current edge
  !    Output, INTEGER  iFace, number corresponding to the edg in ele
  !
  !*************************************************************************



    IMPLICIT NONE
    INTEGER      , INTENT(IN)  :: nDim
    type(Element), INTENT(IN)  :: ele
    type(Edge)   , INTENT(IN)  :: edg
    INTEGER      , INTENT(OUT) :: iFace
    !------------------------------------
    INTEGER, DIMENSION(nDim) :: edgNodes, eleNodes
    !-----------------------------------------------
    INTEGER :: i
    !--------------

    DO i = 1, ndim + 1
       IF ( ele%Edg(i) == edg%num ) THEN
          iFace = i
           EXIT
       END IF
    END DO

  END SUBROUTINE GetFaceNumber
  !==============================================!

END MODULE Tools
