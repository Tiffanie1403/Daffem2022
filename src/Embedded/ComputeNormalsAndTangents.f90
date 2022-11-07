MODULE ComputeNormalsAndTangents

  IMPLICIT NONE

CONTAINS

  !===========================================================!
  SUBROUTINE getNormalAndTangents(d, n, t, nt, nnt, tnt, nDim)
  !===========================================================!

  !*************************************************************************
  !  Defines normals and tangents
  !
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(:) d, distance vector
  !    Input, REAL*8, DIMENSION(:)  nt, normal to the true interface
  !    Input,Output, REAL*8, DIMENSION(:) n, normal to the surrogate
  !    Input,Output,REAL*8, DIMENSION(:) tnt, scalar product tangents and normals (surrogate)
  !    Input,Output,REAL*8, DIMENSION(:,:) t, tangent vector
  !    Output,REAL*8 nnt, scalar product between the two normals
  !    Input,INTEGER nDim, dimension of the resolution (2 or3)
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(:)  , INTENT(IN)    :: d, nt
    REAL*8, DIMENSION(:)  , INTENT(INOUT) :: n, tnt
    REAL*8, DIMENSION(:,:), INTENT(INOUT) :: t
    REAL*8                , INTENT(OUT)   :: nnt
    INTEGER               , INTENT(IN)    :: nDim
    !------------------------------------------------
    REAL*8 :: Nrm, nnn
    !--------------------
    INTEGER :: id, jd
    !------------------

    Nrm = 0.0;
    n = d ! on normalise le vecteur
    DO id = 1, nDim
       Nrm =Nrm + d(id)**2
    END DO
    !calcul de la norme au carre du vecteur definis par la distance
    Nrm = sqrt(Nrm)
    n = n/(Nrm + 1.d-12)    ! on normalise la normal surrogate
    nnt = DOT_PRODUCT(n,nt) !(n_qp,Np) ! si pose probleme vient forcement du calcul de n

    IF ( Nrm > 1.d-12 ) THEN
       IF ( nDim == 2 .OR. ( abs(n(1)) > 1d-8 .OR. abs(n(2)) > 1d-8)) THEN

          t(1,1) = -n(2) ; t(2,1) = n(1)
          IF ( nDim == 3 ) THEN
             t(3,1) = 0.d0
             nnn = 0
             DO  id = 1, nDim
                nnn = nnn + t(id,1)*t(id,1);
             END DO
             nnn = sqrt(nnn)
             t = t/nnn
          END IF
       ELSE
          t(1,1) = -n(3) ; t(3,1) = n(1)
          t(2,1) = 0.d0
          nnn = 0.d0
          DO id = 1, nDim
             nnn = nnn + t(id,1)**2
          END DO
          nnn = sqrt(nnn);
          t(:,1) = t(:,1)/nnn
       END IF
       IF ( nDim == 3 ) THEN
          t(1,2) = n(2)*t(3,1) - n(3)*t(2,1);
          t(2,2) = n(3)*t(1,1) - n(1)*t(3,1);
          t(3,2) = n(1)*t(2,1) - n(2)*t(1,1);
       END IF
       DO id = 1, nDim - 1
          tnt(id) = DOT_PRODUCT(t(:,id),nt)
       END DO
    ELSE ! les deux normales sont confondus
       DO id = 1, nDim
          n(id) = nt(id);
          nnt = 1.d0
          DO jd = 1, nDim - 1
             t(id,jd) = 0.d0
             tnt(jd) = 0.d0
          END DO
       END DO
    END IF

  END SUBROUTINE GetNormalAndTangents
  !================================================!

END MODULE ComputeNormalsAndTangents
