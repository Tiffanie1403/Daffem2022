MODULE PublicVar

  USE Types

  IMPLICIT NONE
  REAL*8, PUBLIC :: alphaIP1, alphaIP2, alphaS, beta
  REAL*8, PUBLIC :: alphaP, alphaB, betaN, betaT, hMesh
  REAL*8, PUBLIC :: t1_MUMPS_TRansfer, t2_MUMPS_transfer
  REAL*8, PUBLIC :: lxx, lxy, lyx, lyy, Stefan_nu
  REAL*8, PUBLIC :: lxz, lyz, lzx, lzy, lzz,Twall
  !-----------------------------------------
  REAL*8, PARAMETER, PUBLIC :: PI = ACOS(-1.d0)
  REAL*8, PARAMETER, PUBLIC :: Se = 1.d-12
  !------------------------------------------
  INTEGER, PUBLIC :: debug, KPrint
  INTEGER, PUBLIC :: icas
  !-------------------------
  INTEGER, PARAMETER, PUBLIC :: Dirichlet = 1
  INTEGER, PARAMETER, PUBLIC :: Neumann = 2
  !--------------------------------------------
  REAL*8, PUBLIC :: mu_1, mu_2, pD_pub, hN_pub
  !-----------------------------
  INTEGER, DIMENSION(6,2), PUBLIC :: Permut
  !--------------------------------------------
  REAL*8, PUBLIC :: t1_initialization, t2_initialization
  REAL*8, PUBLIC :: t1_assemble, t2_assemble 
  REAL*8, PUBLIC :: t1_matrixTransfer, t2_matrixTransfer
  REAL*8, PUBLIC :: t1_solve, t2_solve
  REAL*8, PUBLIC :: t1_PostProcess, t2_PostProcess
  REAL*8, PUBLIC :: t1_output, t2_output,Stw,Sti
  INTEGER,PUBLIC :: nbChangInterface, tourMAx
  !---------------------------------------------------------
  REAL*8, PUBLIC  :: Icing_L1, Icing_L2, Icing_T1, Icing_T2, IcingZone1_ref, IcingZone2_Ref,IcingZone3_ref, Icing_L1_prec
  INTEGER, PUBLIC :: Icing_NodePhysicalInterface
  REAL*8, PUBLIC  :: RefZoneAChanger, Icing_alpha_l, Icing_alpha_s, Icing_PrevFixPosition
  REAL*8, PUBLIC  :: Icing_Lambda1, Icing_Lambda2, Icing_Lambda3, Icing_L3, Icing_Lm, Icing_Tm, Icing_c1, Icing_c2
  REAL*8, PUBLIC  :: Icing_rho, Icing_Tinit, Icing_Tr, Icing_hte, Icing_Tw, Icing_chi, Icing_rhoAirfoil, Icing_c3 
  !----------------------------------------------------------------------------------------

CONTAINS

  !========================================================!
  SUBROUTINE GetEmbeddedLogical(ele, N_LS, eleSolve, eleID)
  !========================================================!

  !*************************************************************************
  !   Defines which element is used for the resolution
  !
  !  Parameters:
  !
  !    Input, type(element) ele, Current element
  !    Input, INTEGER N_LS, Number of nodes in ele
  !    Output, LOGICAL eleSolve, Boolean to incorporate the element in the resolution
  !    Output, LOGICAL eleID, !!!!!
  !
  !*************************************************************************

    IMPLICIT NONE
    type(element), INTENT(IN)  :: ele
    INTEGER      , INTENT(IN)  :: N_LS
    LOGICAL      , INTENT(OUT) :: eleSolve, eleID
    !---------------------------------------------
    INTEGER :: i, j
    !---------------

    eleSolve = .TRUE.
    eleId = .FALSE.

    DO i = 1, N_LS
       IF ( ele%Status(i) /= -1 ) THEN
          eleSolve = .FALSE.
       ELSE IF ( ele%status(i) == 2 ) THEN
          IF ( N_LS > 1 ) THEN
             DO j = 1, N_LS
                IF ( j /= i ) THEN
                   IF ( ele%status(j) /= 1 ) THEN
                      eleId = .TRUE.
                   END IF
                END IF
             END DO
          ELSE
             eleId = .TRUE.
          END IF
       END IF
    END DO

  END SUBROUTINE GetEmbeddedLogical
  !========================================================!

  !============================!
  SUBROUTINE SetPublicVar(Data)
  !============================!
  !*************************************************************************
  !   Put some variables public
  !
  !  Parameters:
  !
  !    Input, type(DataStructure) Data, Data structure
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(IN) :: Data
    !-----------------------------------------
    
    alphaIP1 = Data%alphaIP1 ; alphaIP2 = Data%alphaIP2
    alphaP   = Data%alphaP   ; alphaB   = Data%alphaB
    alphaS   = Data%alphaS   ; icas     = Data%icas  
    KPRint   = Data%KPrint 
    
    mu_1 = 1.d0 ; mu_2 = 10.d0 ; nbChangInterface = 0

    Permut(1,1) = 2 ; Permut(1,2) = 3
    Permut(2,1) = 3 ; Permut(2,2) = 1
    Permut(3,1) = 1 ; Permut(3,2) = 2
    Permut(4,1) = 1 ; Permut(4,2) = 4
    Permut(5,1) = 2 ; Permut(5,2) = 4
    Permut(6,1) = 3 ; Permut(6,2) = 4
    
    
    !represente le rapport entre la chaleur sensible et la chaleur latente 
	!valeur sans dimension caractérisant le transfert thermique 
    Stw = Icing_c1*(Icing_Tw-Icing_Tm)/Icing_Lm    ! Stefan number for the water/liquid phase 
    Sti = Icing_c2*(Icing_Tm-Icing_Tinit)/Icing_Lm ! Stefan number for the ice/solid phase 
    
	Icing_alpha_l = Icing_Lambda1/(Icing_c1*Icing_rho)
	Icing_alpha_s = Icing_Lambda2/(Icing_c2*Icing_rho)

    IF (Data%Icing) THEN
      RefZoneAChanger =IcingZone2_ref
    END IF
    
  END SUBROUTINE SetPublicVar
  !============================!

  !============================================!
  FUNCTION SF(ind,eval) RESULT(phi)
  !============================================!

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ind, eval
    !-- -- -- -- -- -- -- -- -- -- --
    REAL*8 :: phi
    !---------------

    IF ( ind == 3 ) THEN
       phi = 0.d0
    ELSE
       IF ( ind == eval ) THEN
          phi = 1.d0
       ELSE IF ( eval == 4 ) THEN
          phi = 0.5d0
       ELSE
          phi = 0.d0
       END IF
    END IF

  END FUNCTION SF
  !=======================================!

  !=============================!
  FUNCTION SF2(i,xj) RESULT(phi)
  !=============================!

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, xj
    !----------------------------
    REAL*8 :: phi
    !-------------

    IF ( i == 3 ) THEN
       phi = 0.d0
    ELSE
       IF ( xj == i ) THEN
          phi = 1.d0
       ELSE IF ( xj == 3 ) THEN
          phi = 0.5d0
       ELSE
          phi = 0.d0
       END IF
    END IF

  END FUNCTION SF2
  !=============================!

  SUBROUTINE SFs(SF_i1, SF_i2, SF_i3, i)

    IMPLICIT NONE
    REAL*8 , INTENT(OUT) :: SF_i1, SF_i2, SF_i3
    INTEGER, INTENT(IN)  :: i
    !------------------------------------------

    SF_i1 = SF2(i,1) ; SF_i2 = SF2(i,2) ; SF_i3 = SF2(i,3)

  END SUBROUTINE SFs

  !=======================================================!
  SUBROUTINE ComputeNT(d1, d2, d3, n1, n2, n3, t1, t2, t3)
  !=======================================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(2), INTENT(IN)  :: d1, d2, d3
    REAL*8, DIMENSION(2), INTENT(OUT) :: n1, n2, n3, t1, t2, t3
    !-----------------------------------------------------------
    REAL*8 :: Nrm1, Nrm2, Nrm3
    !--------------------------

    Nrm1 = sqrt(d1(1)**2 + d1(2)**2)
    Nrm2 = sqrt(d2(1)**2 + d2(2)**2)
    Nrm3 = sqrt(d3(1)**2 + d3(2)**2)

    n1 = d1/Nrm1 ; n2 = d2/Nrm2 ; n3 = d3/Nrm3

    t1(1) = -n1(2) ; t1(2) = n1(1)
    t2(1) = -n2(2) ; t2(2) = n2(1)
    t3(1) = -n3(2) ; t3(2) = n3(1)


  END SUBROUTINE ComputeNT
  !=======================================================!

    !=================================!
  SUBROUTINE ComputeNodalNT(d, n, t)
  !=================================!

    IMPLICIT NONE
    REAL*8, DIMENSION(2), INTENT(IN)  :: d
    REAL*8, DIMENSION(2), INTENT(OUT) :: n, t
    !-----------------------------------------
    REAL*8 :: Nrm
    !-------------

    Nrm = sqrt(d(1)**2 + d(2)**2)

    n = d/Nrm
    t(1) = -n(2) ; t(2) = n(1)

  END SUBROUTINE ComputeNodalNT
  !=======================================================!

  !=================================!
  SUBROUTINE PrintError(RoutineName)
  !=================================!

  !*************************************************************************
  !   Prints an error if the scheme asked for the compution isn't available
  !
  !  Parameters:
  !
  !    Input, CHARACTER(len=*) RoutineName, Name of the scheme
  !
  !*************************************************************************


    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: RoutineName
    !-------------------------------------------

    PRINT*, " Error in "//RoutineName//" :: unknown scheme"
    PRINT*, "       Available :: - CG"
    PRINT*, "                    - CG-Primal"
    PRINT*, "                    - DG"
    PRINT*, "                    - DG-Primal"
    PRINT*, "                    - MSEG"
    PRINT*, "                    - MSEG-Stab"
    PRINT*, " "
    PRINT*, " If new implementation, add in the different select case"
    STOP


  END SUBROUTINE PrintError
  !=================================!

  !===================================!
  SUBROUTINE PrintError3D(RoutineName)
  !===================================!

  !*************************************************************************
  !   Prints an error if the scheme asked for the compution isn't available
  !   3D version of PrintError
  !  Parameters:
  !
  !    Input, CHARACTER(len=*) RoutineName, Name of the scheme
  !
  !*************************************************************************



    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: RoutineName
    !-------------------------------------------

    PRINT*, " Error in "//RoutineName//" :: unknown scheme"
    PRINT*, "       Available :: - CG"
    PRINT*, "                    - CG-Primal"
    PRINT*, "                    - DG"
    PRINT*, "                    - DG-Primal"
    PRINT*, " "
    PRINT*, " If new implementation, add in the different select case"
    STOP


  END SUBROUTINE PrintError3D
  !===================================!

  !subroutine qui ne semble pas etre utilise
  !==================================!
  SUBROUTINE IdentifyI(ele, edg, iOK)
  !==================================!

    IMPLICIT NONE
    type(element), INTENT(IN)  :: ele
    type(Edge)   , INTENT(IN)  :: Edg
    INTEGER      , INTENT(OUT) :: iOK
    !----------------------------------
    INTEGER, DIMENSION(3,2) :: Permut
    !----------------------------------
    INTEGER :: i
    INTEGER :: Ne1, Ne2, Nt1, Nt2
    !-----------------------------

    Permut(1,1) = 2 ; Permut(1,2) = 3
    Permut(2,1) = 3 ; Permut(2,2) = 1
    Permut(3,1) = 1 ; Permut(3,2) = 2

    Ne1 = edg%Vertex(1) ; Ne2 = edg%Vertex(2)
    IF ( Ne1 >= Ne2 ) THEN
       PRINT*, "IdentifyI, ordering necessary..."
       STOP
    END IF

    DO i = 1, 3
       Nt1 = ele%Vertex(Permut(i,1)) ; Nt2 = ele%Vertex(Permut(i,2))
       IF ( Nt1 >= Nt2 ) THEN
          Nt1 = ele%Vertex(Permut(i,2)) ; Nt2 = ele%Vertex(Permut(i,1))
       END IF
       IF ( (Nt1 == Ne1) .AND. (Nt2 == Ne2) ) THEN
          iOK = i
          EXIT
       END IF
    END DO

  END SUBROUTINE IdentifyI
  !==================================!

  !=======================!
  SUBROUTINE OrderTable(V)
  !=======================!

 ! on definit les references associé aux noeuds dans l'ordre croissant
 ! la plus petite reference de noeud est définit en premer dans le tableau
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT) :: V
    !-----------------------------------------
    INTEGER :: Ni, Nj, Ntmp, N
    INTEGER :: i, j
    !---------------------------

    N = SIZE(V)

    DO i = 1, N
       Ni = V(i)
       DO j = i + 1, N
          Nj = V(j)
          IF ( Nj < Ni ) THEN
             Ntmp = Nj
             Nj = Ni
             Ni = Ntmp
             V(j) = Nj
             V(i) = Ni
          END IF
       END DO
    END DO

  END SUBROUTINE OrderTable
  !=======================!


END MODULE PublicVar
