MODULE Data_mod

  USE Types
  USE Algebra
  USE PublicVar

  IMPLICIT NONE

CONTAINS


  !=======================!
  SUBROUTINE GetData(Data)
  !=======================!
  !*************************************************************************
  !   Reads the dom.data file and Fills the Data structure
  !
  !  Parameters:
  !
  !    Output, type(DataStructure) Data, Structure for the date from the dom.data file
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(OUT) 	:: Data
    !------------------------------------------
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: ATmp, ATmpm1
    !---------------------------------------------------
    INTEGER :: i, l !l numerotation de la ligne de lecture du fichier
    !-------------
    LOGICAL :: file_exists
    !----------------------
    CHARACTER(len=99) :: tmp
    !---------------------------

    CALL getarg(1,Data%FileName)

    INQUIRE(FILE=TRIM(ADJUSTL(Data%FileName))//'.data',EXIST=file_exists)
    !verification que le fichier dom.data exist bien
    IF ( file_exists .eqv. .FALSE. ) THEN
       PRINT*, " **ERROR, data file "//TRIM(ADJUSTL(Data%FileName))//'.data'//" does not exist..."
       CALL SYSTEM("ls -ltr")
       STOP !sinon arret du programme
    END IF

    OPEN(unit = 12, file = TRIM(ADJUSTL(Data%FileName))//'.data') ! ouverture du fichiert dom.data dans l'unite 12
    ! Reading of Global info
    l = 1
    READ(12,*, ERR=99, END=99)
    l = l+1
    READ(12,*,ERR=99,END=99) Data%Ndim, Data%icas
    ALLOCATE(ATmp(Data%Ndim,Data%Ndim))
    ALLOCATE(ATmpm1(Data%Ndim,Data%Ndim))
    l = l+1
    READ(12,*,ERR = 99, END = 99) Data%Space_Scheme, Data%Space_Type, Data%InterfaceScheme
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%SolverType, data%newtonIterations, data%nNewtonIter,  &
        data%NewtonTolerance,Data%TScheme, Data%T, Data%Tmax
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%InitialTime, Data%FinalTime 
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%DtChoiceCFL, Data%h, Data%CFL
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%SpeedCst, Data%Speed, Data%SetEquation
    l = l+1
    IF ( .NOT. data%newtonIterations ) data%nNewtonIter = 1 !deja imposer dans le fichier data
    READ(12, *, ERR = 99, END = 99) Data%IterMax, Data%NKrilov, Data%TolAbs, Data%TolRel, Data%KPrint, &
         Data%MUMPS_Option
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%alphaS
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%alphaIP1, Data%alphaIP2, Data%alphaP, Data%alphaB
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%PressureEnrichment, data%testFunctionEnrichment
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%Output, Data%ComputeError
    l = l+1
    READ(12, *, ERR = 99, END = 99) Data%Embedded, data%immersed
    IF ( data%Embedded .AND. data%immersed ) THEN
       PRINT*, " "
       PRINT*, "Error in "//TRIM(ADJUSTL(data%fileName))//".data"
       PRINT*, " Cannot have both embedded and immersed"
       PRINT*, " "
       STOP
    END IF
    l = l+1
    READ(12, *, ERR=99, END=99)
    IF ( .NOT. Data%Embedded .AND. .NOT. data%immersed ) Data%N_LS = 0

    ! Reading of (different) permeability matrices
    l = l+1
    READ(12, *, ERR=98, END=98)
    l = l+1
    READ(12, *, ERR=98, END=98) Data%N_dom ! nombre de sous domaine
    ALLOCATE(Data%Dom(Data%N_dom))
    IF ( Data%Ndim == 2 ) THEN
       ALLOCATE(Data%lambda(Data%N_dom,4))
       ALLOCATE(Data%lambdam1(Data%N_dom,4))
    ELSE ! on considere donc que deux cas possible le 2 ou le trois
       ALLOCATE(Data%lambda(Data%N_dom,9))
       ALLOCATE(Data%lambdam1(Data%N_dom,9))
    END IF
    DO i = 1, Data%N_dom
       l = l+1
       READ(12, *, ERR=98, END=98) Data%Dom(i)
       IF ( Data%Ndim == 2 ) THEN
          l = l+1
          READ(12, *, ERR=98, END=98) ATmp(1,1), ATmp(1,2), ATmp(2,1), ATmp(2,2)
          CALL Inverse2x2Matrix(ATmp, ATmpm1) ! ATmp1 inverse de la matrice ATmp
          Data%lambda(i,1) = ATmp(1,1) ; Data%lambda(i,2) = ATmp(1,2)
          Data%lambda(i,3) = ATmp(2,1) ; Data%lambda(i,4) = ATmp(2,2)
          Data%lambdam1(i,1) = ATmpm1(1,1) ; Data%lambdam1(i,2) = ATmpm1(1,2)
          Data%lambdam1(i,3) = ATmpm1(2,1) ; Data%lambdam1(i,4) = ATmpm1(2,2)
       ELSE
          l = l+1
          READ(12, *, ERR=98, END=98) ATmp(1,1), ATmp(1,2), ATmp(1,3), &
               ATmp(2,1), ATmp(2,2), ATmp(2,3), &
               ATmp(3,1), ATmp(3,2), ATmp(3,3)

          CALL InverseMatrix(ATmp, ATmpm1) !inverse d'une matrice 3*3

          Data%lambda(i,1) = ATmp(1,1) ; Data%lambda(i,2) = ATmp(1,2)
          Data%lambda(i,3) = ATmp(1,3)
          Data%lambda(i,4) = ATmp(2,1) ; Data%lambda(i,5) = ATmp(2,2)
          Data%lambda(i,6) = ATmp(2,3)
          Data%lambda(i,7) = ATmp(3,1) ; Data%lambda(i,8) = ATmp(3,2)
          Data%lambda(i,9) = ATmp(3,3)
          Data%lambdam1(i,1) = ATmpm1(1,1) ; Data%lambdam1(i,2) = ATmpm1(1,2)
          Data%lambdam1(i,3) = ATmpm1(1,3)
          Data%lambdam1(i,4) = ATmpm1(2,1) ; Data%lambdam1(i,5) = ATmpm1(2,2)
          Data%lambdam1(i,6) = ATmpm1(2,3)
          Data%lambdam1(i,7) = ATmpm1(3,1) ; Data%lambdam1(i,8) = ATmpm1(3,2)
          Data%lambdam1(i,9) = ATmpm1(3,3)
       END IF
    END DO
    l = l+1
    READ(12, *, ERR=98, END=98)

    ! Reading of boundary conditions
    l = l+1
    READ(12, *, ERR=97, END=97)
    l = l+1
    READ(12, *, ERR=97, END=97) Data%N_BC
    ALLOCATE(Data%BC(Data%N_BC))
    ALLOCATE(Data%BC_ref(Data%N_BC))
    ALLOCATE(Data%BCVal(Data%N_BC))
    ALLOCATE(Data%alphaD(Data%N_BC), Data%betaN(Data%N_BC), Data%betaT(Data%N_BC))
    DO i = 1, Data%N_BC
       l = l+1
       READ(12, *, ERR=97, END=97) Data%BC_ref(i), Data%BC(i)
       SELECT CASE ( TRIM(ADJUSTL(Data%BC(i))) )
       CASE ( "Dirichlet","dirichlet" )
          l = l+1
          READ(12, *, ERR=97, END=97) Data%alphaD(i), Data%betaT(i)
       CASE ( "Neumann","neumann")
          l = l+1
          READ(12, *, ERR=97, END=97) Data%betaN(i), Data%betaT(i)
       CASE ( "NeumannJump","neumannjump","Neumannjump","neumannJump" )
          l = l+1
                    READ(12, *, ERR=97, END=97) Data%betaN(i), Data%betaT(i)
       CASE DEFAULT
          PRINT*, " **ERROR in conformal BC..."
          PRINT*, "   Available :"
          PRINT*, "     + Dirichlet/dirichlet"
          PRINT*, "     + Neumann/neumann"
          PRINT*, "     + NeumannJump/neumannjump/Neumannjump/neumannJump"
          STOP
       END SELECT
       l = l+1
       READ(12, *, ERR=97, END=97) Data%BCVal(i)
    END DO
    CLOSE(12)
	
	CALL CheckSchemeData(Data)
	CALL CheckTimeSchemeData(Data)

    IF ( Data%Embedded .OR. data%immersed ) THEN
       INQUIRE(FILE='Embedded.data',EXIST=file_exists)
       IF ( file_exists .eqv. .FALSE. ) THEN
          PRINT*, " **ERROR, data file Embedded.data does not exist..."
          CALL SYSTEM("ls -ltr")
          STOP
       END IF
       OPEN(unit = 13 , file = 'Embedded.data')
       l = 1
       READ(13, *, ERR=96, END=96) Data%N_LS
       ALLOCATE(Data%LSDef(Data%N_LS,2))
       ALLOCATE(Data%LS_L(Data%N_LS,4))
       ALLOCATE(data%LS_center(data%N_LS,data%nDim))
       ALLOCATE(Data%EmbeddedBCType(Data%N_LS))
       ALLOCATE(Data%Emb_alphaD(Data%N_LS))
       ALLOCATE(Data%Emb_betaN(Data%N_LS))
       ALLOCATE(Data%Emb_betaT(Data%N_LS))
       ALLOCATE(Data%BCEmbVal(Data%N_LS))
       DO i = 1, Data%N_LS
          l = l+1
          READ(13, *, ERR=96, END=96) Data%EmbeddedBCType(i)
          l = l+1
          READ(13, *, ERR=96, END=96) Data%LSDef(i,:)
          SELECT CASE ( TRIM(ADJUSTL(Data%LSDef(i,1))) )
          CASE ( "Circle","circle","Sphere","sphere" )
             l = l+1
             READ(13, *, ERR=96, END=96) Data%LS_L(i,1)
             l = l+1
             READ(13, *, ERR=96, END=96) data%LS_center(i,:)
          CASE ( "Box","box" )
             l = l+1
             READ(13, *, ERR=96, END=96) Data%LS_L(i,1:2)
             l = l+1
             READ(13, *, ERR=96, END=96) data%LS_center(i,:)
          END SELECT
          SELECT CASE ( TRIM(ADJUSTL(Data%EmbeddedBCType(i))) )
          CASE ( "Dirichlet","dirichlet" )
             l = l+1
             READ(13, *, ERR=96, END=96) Data%Emb_alphaD(i), Data%Emb_betaT(i)
             l = l+1
             READ(13, *, ERR=96, END=96) Data%BCEmbVal(i)
          CASE ( "Neumann","neumann" )
             Data%Emb_alphaD(i) = 0.d0
             l = l+1
             READ(13, *, ERR=96, END=96) Data%Emb_betaN(i), Data%Emb_betaT(i)
             l = l+1
             READ(13, *, ERR=96, END=96) Data%BCEmbVal(i)
          CASE ( "xy_dependancy" )
          CASE ( "NeumannJump","neumannjump","Neumannjump","neumannJump" )
             READ(13, *, ERR=96, END=96) Data%BCEmbVal(i)
          CASE DEFAULT
             PRINT*, ' **ERROR in bondary type...'
             PRINT*, '   Available :'
             PRINT*, '     + Dirichlet/dirichlet'
             PRINT*, '     + Neumann/neumann'
             PRINT*, '     + NeumannJump/Neumannjump/neumannjump/neumannJump'
             STOP
          END SELECT

          IF (i /= Data%N_LS) THEN
             l = l+1
             READ(13, *, ERR=96, END=96)
          END IF
       END DO
       CLOSE(13)
    ELSE
       Data%N_LS = 0
    END IF


    data%icing = .FALSE.

    ! Data for Icing stuff
    IF ( data%icas >= 100 .OR. data%icas == 001 ) THEN
       data%icing = .TRUE.
       INQUIRE(FILE='Icing.data',EXIST=file_exists)
       IF ( file_exists .eqv. .FALSE. ) THEN
          PRINT*, " **ERROR, data file Icing.data does not exist..."
          CALL SYSTEM("ls -ltr")
          STOP
       END IF
       OPEN(unit = 13 , file = 'Icing.data')
       READ(13,*, ERR=95, END=95) IcingZone1_ref, IcingZone2_ref, IcingZone3_ref !ref des zones
       READ(13,*, ERR=95, END=95) Icing_Lambda1, Icing_Lambda2, Icing_Lambda3  !valeurs des coeff lambda
       READ(13,*, ERR=95, END=95) Data%DefInterfaceNode, Icing_NodePhysicalInterface
       READ(13,*, ERR=95, END=95) Data%ResisTance, Data%Hl, Data%Ha, Data%Hi, Data%Lh, Data%L
       READ(13,*, ERR=95, END=95) Icing_L1, Icing_L2, Icing_L3 !valeur des coefficients ei, longueur caractéristique en cas de domaine rectangulaire  
       READ(13,*, ERR=95, END=95) Icing_T1, Icing_T2 !valeurs des conditions de dirichlet
       READ(13,*, ERR=95, END=95) Icing_Lm, Icing_Tm ! Lm energie latente Tm temperature de fusion
       READ(13,*, ERR=95, END=95) Icing_c1,Icing_c2 , Icing_c3 !! valeur des coefficient cw (water est la zone 1) et ci ( ice est la zone 2)
       READ(13,*, ERR=95, END=95) Icing_rho, Icing_Tinit, Icing_Tr, Icing_hte , Icing_Tw, Icing_chi  ! rho (pression), Temperature initiale, temperature exterieur, heat transfer coefficient , temperature sur le bord de la source de chaleur 
       READ(13,*, ERR=95, END=95) Icing_rhoAirfoil  ! rho (pression), Temperature initiale, temperature exterieur, heat transfer coefficient , temperature sur le bord de la source de chaleur 

    END IF
	Icing_NodePhysicalInterface = Icing_NodePhysicalInterface +2 !because the two extreme points are going to be fixed points 
	
	!Coordinate of the resistance in the solid airfoil 
	!IF(Data%Resistance) THEN 
	!	ALLOCATE(CoordResistance(2,2))
	!	CoordResistance(1,1) = Data%Hl 
	!	CoordResistance(2,1) = Data%Hl 
	!	CoordResistance(1,2) = Data%L/2. - Data%Lh/2. ! define in the middle of the domain 
	!	CoordResistance(2,2) = Data%L/2. + Data%Lh/2. 
	!END IF 
    !on definit Dt comme valant CFL*h^2/Lmax si demander dans le fichier de données
    Data%tour = 0  ! permet de savoir si on est à l'étape d'initialisation ou à qu'elle tour on en est
    Data%ErrorT = 0.d0 ; Data%ErrorBx = 0.d0 ; Data%ErrorBy = 0.d0

	IF(data%icas /= 001) THEN  
		IF(DATA%TSCHEME == "EULER") THEN
			IF(Data%DtChoiceCFL .AND. Data%Icing) THEN
				IF(Icing_Lambda1>Icing_Lambda2) THEN
					Data%DT = (Data%CFL*Data%h**2)/Icing_Lambda1
				ELSE
					Data%DT = (Data%CFL*Data%h**2)/Icing_Lambda2
				END IF
			END IF
		ELSE IF (DATA%TSCHEME=="NICOLSON" .OR. Data%TSCHEME == "BDF2") THEN !both are second order accurate 
			IF(Data%SpeedCst) THEN 
				Data%Dt = Data%h*Data%CFL/(2*Data%Speed) !general case, if it's not then it's replaced in the main function 
			END IF 
		ELSE
			PRINT*,"Dt not available for this time scheme"
			STOP
		END IF
		
		IF(Data%SpeedCst) THEN 
			IF( Data%Speed > (Data%h*Data%CFL/(2*Data%DT))) THEN
				PRINT*,"You are trying to move the interface following a constant speed"
				STOP "Speed to high for the geometry"
			END IF
		END IF 
		PRINT*,"Value of DT",Data%Dt 
	END IF 
	
	IF (Data%SpeedCst .AND. Data%ResisTance) THEN 
		PRINT*,"Problem in the reading of the output files"
		STOP 
	END IF 

    RETURN ! -- OK
99  CONTINUE
    WRITE(tmp,*) l
    PRINT*, " ** Error in Data.f90"
    PRINT*, "  First part of the datafile must contain:"
    PRINT*, " "
    PRINT*, " empty line (or comment)"
    PRINT*, " data%Ndim (int), data%icas (int)"
    PRINT*, " data%space_Scheme (char)"
    PRINT*, " data%solverType (char), data%newtonIterations (int), data%nNewtonIter(int), data%NewtonTolerance (double), "&
          "data%TScheme, data%DT, Data%T, Data%Tmax"
    PRINT*, " data%iterMax (int), data%NKrilov (int), data%tolAbs (double), data%tolRel (double), "&
         "data%KPrint (double), data%MUMPS_Option(int)"
    PRINT*, " data%alphaS"
    PRINT*, " data%alphaIP1 (double), data%alphaIP2(double), data%alphaP (double), data%alphaB (double)"
    PRINT*, " data%pressureEnrichment (bool), data%testFunctionEnrichment (bool)"
    PRINT*, " data%Output (bool), data%computeError (bool)"
    PRINT*, " data%embedded (bool), data%immersed (bool)"
    PRINT*, " empty line (or comment) "
    PRINT*, " "
    PRINT*, " --> Check line "//TRIM(ADJUSTL(tmp))
    PRINT*, " "
    STOP

98  CONTINUE
    WRITE(tmp,*) l
    PRINT*, " ** Error in Data.f90"
    PRINT*, "  Second part of the datafile must contain:"
    PRINT*, " "
    PRINT*, " empty line (or comment)"
    PRINT*, " data%N_dom (int)"
    PRINT*, " data%dom(i) (int)"
    PRINT*, " data%lambda(i, 1:nDim,1:nDim) (nDim x nDim double)"
    PRINT*, " empty line (or comment)"
    PRINT*, " "
    PRINT*, " --> Check line "//TRIM(ADJUSTL(tmp))
    PRINT*, " "
    STOP

97  CONTINUE
    WRITE(tmp,*) l
    PRINT*, " ** Error in Data.f90"
    PRINT*, "  Third part of the datafile must contain:"
    PRINT*, " "
    PRINT*, " empty line (or comment)"
    PRINT*, " data%N_BC (int)"
    PRINT*, " data%BC_ref(i) (int), Data%BC(i) (char)"
    PRINT*, " data%alphaD(i) (double), data%betaT(i) (double)"
    PRINT*, " data%BCVal(i) (double)"
    PRINT*, " "
    PRINT*, " --> Check line "//TRIM(ADJUSTL(tmp))
    PRINT*, " "
    STOP

96  CONTINUE
    WRITE(tmp,*) l
    PRINT*, " ** Error in Data.f90"
    PRINT*, "  Embedded.data file must contain:"
    PRINT*, " "
    PRINT*, " data%N_LS (int)"
    PRINT*, " data%EmbeddedBCType(i) (char)"
    PRINT*, " data%LSDef(i,:) (char char)"
    PRINT*, " data%LS_L(i,:) (ncharacterisctic size of geom x double)"
    PRINT*, " data%LS_center(i,:) (center of the geometry, double)"
    PRINT*, " data%emb_alphaD(i) (double), data%Emb_betaT(i) (double)"
    PRINT*, " data%BCEmbVal(i) (double)"
    PRINT*, " empty line or comment"
    PRINT*, " "
    PRINT*, " --> Check line "//TRIM(ADJUSTL(tmp))
    PRINT*, " "
    STOP

95  CONTINUE
    PRINT*, " **Error in  Icing.data"
    PRINT*, "  Icing.data must contains"
    PRINT*, "IcingZone1_ref (int), IcingZone2_ref (int)"
    PRINT*, "Icing_Lambda1 (double), Icing_Lambda2 (double)"
    PRINT*, "Icing_L1 (double), Icing_L2 (double), Icing_L3"
    PRINT*, "Icing_T1 (double), Icing_T2 (double)"
    STOP

  END SUBROUTINE GetData
  !=======================!


  !====================================!
  SUBROUTINE Inverse2x2Matrix(Mat, Am1)
  !====================================!

  !*************************************************************************
  !   Inversion of a 2x2 matrix
  !
  !  Parameters:
  !
  !    Input, REAL*8, DIMENSION(2,2) Mat, matrix we want to inverse
  !    Output, REAL*8, DIMENSION(2,2) Am1, inverse of mat
  !
  !*************************************************************************


    IMPLICIT NONE
    REAL*8, DIMENSION(2,2), INTENT(IN)  :: Mat
    REAL*8, DIMENSION(2,2), INTENT(OUT) :: Am1
    !----------------------------------------
    REAL*8, DIMENSION(2,2) :: Tmp
    REAL*8 :: a, b, c, d, det, delta, ev1, ev2
    REAL*8 :: a1, b1, c1
    !------------------------

    a = mat(1,1) ; b = mat(1,2) ; c = mat(2,1) ; d = mat(2,2)

    det = a*d - b*c
    IF ( abs(det) < 0.000001d0 ) THEN
       PRINT*, ' ERROR :: Permeability matrix singular...'
       PRINT*, abs(det)
       STOP
    END IF

    a = Mat(1,1) ; b = Mat(1,2) ; c = Mat(2,2)
    a1 = 1.d0 ; b1 = -(a+c) ; c1 = a*c - b*b
    delta = b1*b1 - 4.d0*a1*c1
    ev1 = (-b1 - sqrt(delta)) / (2.d0*a1)
    ev2 = (-b1 + sqrt(delta)) / (2.d0*a1)
    IF ( (ev1 < 0.d0) .OR. (ev2 < 0.d0) ) THEN
       PRINT*, ' ERROR Permeability matrix not positive... '
       PRINT*, ev1, ev2
       STOP
    END IF


    a = mat(1,1) ; b = mat(1,2) ; c = mat(2,1) ; d = mat(2,2)

    Am1(1,1) = d
    Am1(1,2) = -b
    Am1(2,1) = -c
    Am1(2,2) = a
    Am1 = 1.d0/det * Am1
    Tmp = MATMUL(Am1,Mat)

  END SUBROUTINE Inverse2x2Matrix
  !====================================!


  !===============================!
  SUBROUTINE CheckSchemeData(Data)
  !===============================!

  !*************************************************************************
  ! Checks if the scheme (Cg,Dg,..) defined in dom.data is implemented!
  ! Stops the code if it's not
  !  Parameters:
  !
  !    Input, type(DataStructure) Data, Data structure from dom.data file
  !
  !*************************************************************************


    IMPLICIT NONE
    type(DataStructure), INTENT(IN) :: Data
    !---------------------------------------

    SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
    CASE ( "CG" )
		IF(Data%Space_Type =="EXPLICIT") THEN
			PRINT*, " %%%% Continuous Galerkin Explicit scheme in mixed form %%%% "
			PRINT*, " "
		ELSE IF (Data%Space_Type =="IMPLICIT") THEN
			PRINT*, " %%%% Continuous Galerkin Implicit scheme in mixed form %%%% "
			PRINT*, " "
		END IF 
    CASE ( "CG-Primal" )       
       IF(Data%Space_Type =="EXPLICIT") THEN
			PRINT*, " %%%% Continuous Galerkin Explicit scheme in primal form %%%% "
			PRINT*, " "
		ELSE IF (Data%Space_Type =="IMPLICIT") THEN
			PRINT*, " %%%% Continuous Galerkin Implicit scheme in primal form %%%% "
			PRINT*, " "
		END IF 
		
    CASE ( "DG" )
       PRINT*, " %%%% Discontinuous Galerkin scheme in mixed form %%%% "
       PRINT*, " "
    CASE ( "DG-Primal" )
       PRINT*, " %%%% Discontinuous Galerkin scheme in primal form %%%% "
       PRINT*, " "
    CASE ( "DEFAULT" )
       CALL PrintError("CheckSchemeData in Data.f90")
    END SELECT

  END SUBROUTINE CheckSchemeData
  !===============================!

    !===============================!
    SUBROUTINE  CheckTimeSchemeData(Data)
    !===============================!

    !*************************************************************************
    ! Checks if the scheme for the time discretization defined in dom.data
    ! is implemented, it stops the code if it's not
    !  Parameters:
    !
    !    Input, type(DataStructure) Data, Data structure from dom.data file
    !
    !*************************************************************************


      IMPLICIT NONE
      type(DataStructure), INTENT(INOUT) :: Data
      !---------------------------------------

      SELECT CASE ( TRIM(ADJUSTL(Data%TScheme)) )
      CASE ("EULER","euler","Euler","euLER","eULER","eulER","euleR","EuLeR","EUL","eul","EU")
         PRINT*, " %%%% Backward Euler scheme for time discretization %%%% "
         PRINT*, " "
         Data%TScheme="EULER"
       CASE ("NICOLSON","Nicolson","nicolson","NICOL","NICO","NICOLS","nico","nicol")
         PRINT*, " %%%% Crank-Nicolson scheme for time discretization %%%% "
         PRINT*, " "
         Data%TScheme="NICOLSON"
       CASE ("BDF2","bdf2","Bdf2","bDf2","bdF2","BDf2","BdF2","bDF2")
		 PRINT*, " %%%% Backward Differentiation Scheme Order 2 for time discretization %%%% "
         PRINT*, " "
         Data%TScheme="BDF2"
       CASE DEFAULT
          PRINT*, " **ERROR in temporal Scheme"
          PRINT*, "   Available :"
          PRINT*, "     + EULER"
          PRINT*, "     + NICOLSON"
          PRINT*, "		+ BDF2"
          CALL PrintError("CheckTimeSchemeData in Data.f90")
          STOP
      END SELECT

    END SUBROUTINE CheckTimeSchemeData
    !===============================!

END MODULE Data_mod
