MODULE L2Norm

  USE Types
  USE ExactSolutiont
  USE ExactFunctionJumpt

  IMPLICIT NONE

CONTAINS

  !========================================!
  SUBROUTINE ComputeL2Norm(Mesh, Sol, Data)
  !========================================!

  !*************************************************************************
  !   Displays the residual norms of the different unknowns after resolution
  !
  !  Parameters:
  !
  !    Input, type(DataStructure) Data, Data structure from the dom.data file
  !    Input, type(MeshType) Mesh, Mesh structure from the dom.mesh file
  !    Input, type(SolStructure) Sol, Structure for the solution
  !
  !*************************************************************************



    IMPLICIT NONE
    type(DataStructure), INTENT(IN)    :: Data
    type(MeshType)     , INTENT(IN)    :: Mesh
    type(SolStructure) , INTENT(INOUT) :: Sol
    !----------------------------------------
    type(element) :: ele
    !--------------------
    REAL*8, DIMENSION(:), ALLOCATABLE :: Exa
    !------------------------------------------
    REAL*8 :: NrmL2P, NrmL2B
    REAL*8 :: x, y, z,res,length
    !----------------------------------------------
    Integer :: ie, i, NU_i, Pos_i,Nu_i2, Nu_j,ib
    INTEGER :: Ne, Nv, id, nDim, nNodesPerElement
    !----------------------------------------------
    LOGICAL :: EleSolve, OnOUTBond 
    !-------------------


    Nv = Mesh%Nv ; Ne = Mesh%Ne ; nDim = Data%nDim
    ALLOCATE ( Exa(nDim+1) ) ; res =0.d0
    Exa=0.d0 ; NrmL2P = 0.d0 ; NrmL2B = 0.d0
	!length = interface position + 2
	length = Icing_L1+4*Data%h ! we add 4 lengths of elements 


		SELECT CASE ( TRIM(ADJUSTL(Data%Space_Scheme)) )
		CASE ("CG","CG-Primal")
			DO ie = 1, Ne

				ele = Mesh%Ele(ie)
				nNodesPerElement = ele%nNodesPerElement			
				
				eleSolve = .TRUE.
				IF (.NOT. Data%Embedded .AND. Data%Icing) THEN
					DO i = 1, Data%N_LS
					IF ( ele%Status(i) /= -1 ) eleSolve = .FALSE.
					END DO
				END IF

				IF ( eleSolve ) THEN
					DO i = 1, nNodesPerElement
						NU_i = Ele%Vertex(i)
						!cas icing 
						NU_i2=Ele%SurVertex(i) !pour avoir la bonne numerotation du noeud
						
						! we only take the points in an area aroud the interface position 
						IF(mesh%Vertex(NU_i)%X(1) <= length) THEN 
							! we also remove points on the boundary for the calcul of the error 
							IF(mesh%Vertex(NU_i)%X(1)/=0 .AND. mesh%Vertex(NU_i)%X(2)/=0 .AND. &
							mesh%Vertex(NU_i)%X(1)/=Icing_L2 .AND. mesh%Vertex(NU_i)%X(2)/=Icing_L3) THEN 

								IF ( .NOT. Data%Icing ) THEN
									CALL ExactSolt(ele%L, Exa(:), mesh%Vertex(NU_i)%X ,Data%T, ele, nDim)
								ELSE
									CALL ExactSolJumpt(ele%L, Exa(:), mesh%Vertex(NU_i)%X ,Data%T, ele, nDim)
								END IF

								DO id = 1,2
									NrmL2B = NrmL2B + ele%A/(nDim+1)*(Exa(id)-Sol%Bcg(NU_i2,id))**2
									!PRINT*,"b",id, Exa(id),Sol%Bcg(NU_i2,id)
								END DO

								NrmL2P = NrmL2P + ele%A/(nDim+1)*(Exa(nDim+1)- Sol%Pcg(NU_i2))**2
							
							END IF 
						END IF 
				    END DO
			    END IF
		    END DO
		   
		CASE DEFAULT
			CALL PrintError("ComputeL2Norm")
		END SELECT
		
		NrmL2B  = sqrt(NrmL2B)
		NrmL2P  = sqrt(NrmL2P)
		DEALLOCATE(Exa)
		
		PRINT*, '*****************************************'
		PRINT*, ' '
		PRINT*, ' L2 Error :'
		PRINT*, '   ** T :', NrmL2P
		IF(Data%Space_Scheme=="CG") PRINT*, '   ** B :', NrmL2B
		!if Primal formulation the error on B is the error made by the Green Gauss reconstruction 
		IF(Data%Space_Scheme=="CG-Primal") PRINT*, '   error Green Gauss ** B :', NrmL2B
		PRINT*, ' '
		PRINT*, '*****************************************'
		PRINT*, ' '
	    WRITE(112,*)Data%T,NrmL2P, NrmL2B
	    
  END SUBROUTINE ComputeL2Norm
  !========================================!


END MODULE L2Norm
