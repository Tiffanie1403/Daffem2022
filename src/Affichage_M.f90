MODULE Affichage_M
! display matrix

  USE Types
  USE PublicVar

  IMPLICIT NONE

CONTAINS

  !===============================================!
  SUBROUTINE Affichage_Mat(Mat,Sol,Mesh,Data)
  !===============================================!

  !*************************************************************************
  !   Displays the global matrix of the system Mat%Coeff
  !
  !  Parameters:
  !
  !    Input, type (SparseMatrix) Mat, Structure for the global matrix
  !    Input, type(SolStructure) Sol, Structure for the solution vector
  !    Input, type(MeshType) Mesh, Structure for the mesh from dom.mesh
  !    Input, type(DataStructure) Data, Data structure from dom.data
  !*************************************************************************


    IMPLICIT NONE

    type(SparseMatrix)    , INTENT(IN)    :: Mat
    type(SolStructure)    , INTENT(IN)    :: Sol
    type(MeshType)        , INTENT(IN)    :: Mesh
    type(DataStructure)   , INTENT(IN)    :: Data


    !--------------------------------------------------------

    INTEGER    :: nDim, NVar, NvActive , Ns, NMax, ivar, jvar
    INTEGER    :: test,nb
    REAL*8 , DIMENSION(:,:), ALLOCATABLE :: RealTab
    !--------------------------------------------------------


    IF (Data%Icing .eqv. .TRUE.) THEN
      Ns=Mesh%NB_VertInterf!nombre de noeud sur l'interface
    ELSE
      Ns=0
    END IF

    NvActive = Sol%NvActive ; nDim = Data%nDim
    NVar = nDim + 1 ; test=0 ; nb=0 ; NMax = 70

    ALLOCATE(RealTAb(30,30))
    RealTAb = 0.d0

    !WRITE(*,*) "Tab de la matrice globale"
  !  WRITE(*,*) "Nombre de noeuds", NvActive, Ns

      DO jvar = 1, 10 ! Nmax*Nvar
           DO ivar = 1, 30
               IF (ivar >= 1 .AND. ivar <= 10) THEN
                 IF (jvar == 1) THEN
                   RealTab(ivar,ivar)=Mat%Coeff(ivar,jvar)
                 ELSE
                   IF (Mat%Ind(ivar,jvar) /= 0) THEN
                     RealTab(ivar,Mat%Ind(ivar,jvar))=Mat%Coeff(ivar,jvar)
                   END IF
                 END IF
               ELSE IF ( ivar >= 11 .AND. ivar <= 20) THEN
                 IF (jvar == 1) THEN
                   RealTab(ivar,ivar-10)=Mat%Coeff(ivar,jvar)
                 ELSE
                   IF (Mat%Ind(ivar-10,jvar) /= 0) THEN ! ie que l'on a bien un voisin
                     RealTab(ivar,Mat%Ind(ivar-10,jvar))=Mat%Coeff(ivar,jvar)
                   END IF
                 END IF
               ELSE
                 IF (jvar == 1) THEN
                   RealTab(ivar,ivar-20)=Mat%Coeff(ivar,jvar)
                 ELSE
                   IF (Mat%Ind(ivar-20,jvar) /= 0) THEN
                     RealTab(ivar,Mat%Ind(ivar-20,jvar))=Mat%Coeff(ivar,jvar)
                   END IF
                 END IF
              END IF
           END DO
        END DO

        DO jvar = Nmax+1,Nmax+10 ! Nmax*Nvar
             DO ivar = 1, 30
                 IF (ivar >= 1 .AND. ivar <= 10) THEN
                   IF (jvar-Nmax == 1) THEN
                     RealTab(ivar,ivar+10)=Mat%Coeff(ivar,jvar)
                   ELSE
                     IF (Mat%Ind(ivar,jvar-Nmax) /= 0) THEN
                       RealTab(ivar,Mat%Ind(ivar,jvar-Nmax)+10)=Mat%Coeff(ivar,jvar)
                     END IF
                   END IF
                 ELSE IF ( ivar >= 11 .AND. ivar <= 20) THEN
                   IF (jvar-Nmax == 1) THEN
                     RealTab(ivar,ivar)=Mat%Coeff(ivar,jvar)
                   ELSE
                     IF (Mat%Ind(ivar-10,jvar-Nmax) /= 0) THEN
                       RealTab(ivar,Mat%Ind(ivar-10,jvar-Nmax)+10)=Mat%Coeff(ivar,jvar)
                     END IF
                   END IF
                 ELSE
                   IF (jvar-Nmax == 1) THEN
                     RealTab(ivar,ivar-10)=Mat%Coeff(ivar,jvar)
                   ELSE
                     IF (Mat%Ind(ivar-20,jvar-Nmax) /= 0) THEN
                       RealTab(ivar,Mat%Ind(ivar-20,jvar-Nmax)+10)=Mat%Coeff(ivar,jvar)
                     END IF
                   END IF
                END IF
             END DO
          END DO

          DO jvar = 2*Nmax+1,2*Nmax+10 ! Nmax*Nvar
               DO ivar = 1, 30
                   IF (ivar >= 1 .AND. ivar <= 10) THEN
                     IF (jvar-2*Nmax == 1) THEN
                       RealTab(ivar,ivar+20)=Mat%Coeff(ivar,jvar)
                     ELSE
                       IF (Mat%Ind(ivar,jvar-2*Nmax) /= 0) THEN
                         RealTab(ivar,Mat%Ind(ivar,jvar-2*Nmax)+20)=Mat%Coeff(ivar,jvar)
                       END IF
                     END IF
                   ELSE IF ( ivar >= 11 .AND. ivar <= 20) THEN
                     IF (jvar-2*Nmax == 1) THEN
                       RealTab(ivar,ivar+10)=Mat%Coeff(ivar,jvar)
                     ELSE
                       IF (Mat%Ind(ivar-10,jvar-2*Nmax) /= 0) THEN
                         RealTab(ivar,Mat%Ind(ivar-10,jvar-2*Nmax)+20)=Mat%Coeff(ivar,jvar)
                       END IF
                     END IF
                 ELSE
                   IF (jvar-2*Nmax == 1) THEN
                       RealTab(ivar,ivar)=Mat%Coeff(ivar,jvar)
                     ELSE
                       IF (Mat%Ind(ivar-20,jvar-2*Nmax) /= 0) THEN
                         RealTab(ivar,Mat%Ind(ivar-20,jvar-2*Nmax)+20)=Mat%Coeff(ivar,jvar)
                       END IF
                     END IF
                  END IF
               END DO
            END DO

    DO jvar = 1, 30 ! Nmax*Nvar
      PRINT*, jvar
      DO ivar = 1, 30
             WRITE(*,*) RealTab(ivar,jvar)
      END DO
    END DO

        DEALLOCATE(RealTab)


  END SUBROUTINE Affichage_Mat
  !======================================================================!

  !=====================================================================================!
  SUBROUTINE Affichage_MatAloc(Aloc,taille,ie)
  !=====================================================================================!

  !*************************************************************************
  !   Displays the local matrix associate to the element ie
  !
  !  Parameters:
  !
  !    Input, type REAL*8, DIMENSION(:,:) Aloc, the local matrix
  !    Input, type INTEGER  taille, the size of the local matrix
  !    Input, type INTEGER ie, the element number associate to the local matrix
  !
  !*************************************************************************


    IMPLICIT NONE

    REAL*8, DIMENSION(:,:), INTENT(IN)    :: Aloc
    INTEGER               , INTENT(IN)    :: taille, ie
    !--------------------------------------------------------
    INTEGER                               :: i,j
    !--------------------------------------------------------


    WRITE(*,*) "Matrice locale pour l'element ", ie
    DO i = 1, taille
         DO j = 1, taille
             WRITE(*,*) Aloc(i,j)
         END DO
      END DO



  END SUBROUTINE Affichage_MatAloc
  !======================================================================!

END MODULE Affichage_M
