MODULE Hash

  USE Types

  IMPLICIT NONE

CONTAINS

  !=======================================!
  SUBROUTINE FillHash2D(Mesh, N_ele, Hash)
  !=======================================!

  !*************************************************************************
  !  Hash table in 2D
  !
  !  Parameters:
  !
  !    Input,Output, type(MeshType) Mesh, mesh structure from the dom.mesh file
  !    Input,Output, type(HashTable2D) Hash, Hash table
  !    Input, INTEGER N_ele, element number

  !
  !*************************************************************************

    IMPLICIT NONE
    type(MeshType)   , INTENT(INOUT) :: Mesh
    type(HashTable2D), INTENT(INOUT) :: Hash
    INTEGER          , INTENT(IN)    :: N_ele
    !------------------------------------------
    INTEGER :: alpha, beta, cle, nd1, nd2, ModCle, pt1, pt2, NewLoc
    INTEGER :: i, j, iter, Next, Npts, k
    !-----------------------------------------------------------------

    alpha = 7
    beta = 11

    DO i = 1,2
       DO j = i+1,3
          pt1 = Mesh%Ele(N_ele)%Vertex(i)
          pt2 = Mesh%Ele(N_ele)%Vertex(j)
          IF (pt1 > pt2) THEN
             pt1 = Mesh%Ele(N_ele)%Vertex(j)
             pt2 = Mesh%Ele(N_ele)%Vertex(i)
          END IF

          cle = alpha*pt1 + beta*pt2
          ModCle = mod(cle,Hash%SizeH) + 1

          DO
             IF( Hash%Item(ModCle)%Tri == -1 ) THEN
                Hash%Item(ModCle)%Tri = 3*N_ele + 5 - (i + j)
                Hash%Item(ModCle)%Vertex(1) = pt1
                Hash%Item(ModCle)%Vertex(2) = pt2
                Mesh%Ned = Mesh%Ned + 1
                Hash%Edg(Mesh%Ned)%Vertex(1) = pt1
                Hash%Edg(Mesh%Ned)%Vertex(2) = pt2
                Hash%Edg(Mesh%Ned)%Tri(1) = N_ele
                Hash%EdgNum(ModCle) = Mesh%Ned

                DO k = 1, 3
                   IF ( Mesh%Ele(N_ele)%Edg(k) == -1 ) THEN
                      Mesh%Ele(N_ele)%Edg(k) = Mesh%Ned
                      EXIT
                   END IF
                END DO
                !IF ( ( pt1 == 2 ) .AND. ( pt2 == 5 ) ) THEN
                !   PRINT*, Hash%Compt, Mesh%Ned, N_ele
                !END IF
                EXIT
             ELSE
                nd1 = Hash%Item(ModCle)%Vertex(1)
                nd2 = Hash%Item(ModCle)%Vertex(2)
                IF ( (pt1 == nd1) .AND. (pt2 == nd2) ) THEN
                   Mesh%Ele(N_ele)%Adj(6 - (i + j)) = Hash%Item(ModCle)%Tri
                   Mesh%Ele(Hash%Item(ModCle)%Tri/3)%Adj(mod(Hash%Item(ModCle)%Tri,3) + 1) = 3*N_ele + 5 - (i + j)
                   Hash%Edg(Hash%EdgNum(ModCle))%Tri(2) = N_ele
                   DO k = 1, 3
                      IF ( Mesh%Ele(N_ele)%Edg(k) == -1 ) THEN
                         Mesh%Ele(N_ele)%Edg(k) = Hash%EdgNum(ModCle)
                         EXIT
                      END IF
                   END DO
                   EXIT
                ELSE
                   IF ( Hash%Item(ModCle)%Nxt == 0 ) THEN
                      Hash%Item(ModCle)%Nxt = Hash%Compt
                      Hash%Compt = Hash%Compt + 1
                      ModCle = Hash%SizeH + Hash%Item(ModCle)%Nxt
                   ELSE
                      ModCle = Hash%SizeH + Hash%Item(ModCle)%Nxt
                   END IF
                END IF
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE FillHash2D
  !=======================================!


  !=======================================!
  SUBROUTINE FillHash3D(Mesh, N_ele, Hash)
  !=======================================!

  !*************************************************************************
  !  Hash table in 3D
  !
  !  Parameters:
  !
  !    Input,Output, type(MeshType) Mesh, mesh structure from the dom.mesh file
  !    Input,Output, type(HashTable2D) Hash, Hash table
  !    Input, INTEGER N_ele, element number

  !
  !*************************************************************************


    IMPLICIT NONE
    type(MeshType)   , INTENT(INOUT) :: Mesh
    type(HashTable3D), INTENT(INOUT) :: Hash
    INTEGER          , INTENT(IN)    :: N_ele
    !------------------------------------------
    type(element) :: TmpEle
    INTEGER, DIMENSION(3) :: node, idx
    INTEGER :: alpha, beta, gamma, cle, nd1, nd2, nd3, ModCle, pt1, pt2, NewLoc
    INTEGER :: i, j, iter, Next, Npts, k, nodeTmp, l, m
    INTEGER :: Nmin, Nmax, imin, imax, Nmid, AdjNum
    !-----------------------------------------------------------------------

    alpha = 7
    beta = 11
    gamma = 13

    DO i = 1,2
       DO j = i+1,3
          DO k = j+1, 4

             node(1) = Mesh%Ele(N_ele)%Vertex(i)
             node(2) = Mesh%Ele(N_ele)%Vertex(j)
             node(3) = Mesh%Ele(N_ele)%Vertex(k)

             IF ( ( node(1) < node(2) ) .AND. ( node(1) < node(3) ) ) THEN
                IF ( node(2) > node(3) ) THEN
                   Nmin = node(1) ; Nmid = node(3) ; Nmax = node(2)
                ELSE
                   Nmin = node(1) ; Nmid = node(2) ; Nmax = node(3)
                END IF
             ELSE IF ( ( node(2) < node(1) ) .AND. ( node(2) < node(3) ) ) THEN
                IF ( node(1) < node(3) ) THEN
                   Nmin = node(2) ; Nmid = node(1) ; Nmax = node(3)
                ELSE
                   Nmin = node(2) ; Nmid = node(3) ; Nmax = node(1)
                END IF
             ELSE
                IF ( node(1) < node(2) ) THEN
                   Nmin = node(3) ; Nmid = node(1) ; Nmax = node(2)
                ELSE
                   Nmin = node(3) ; Nmid = node(2) ; Nmax = node(1)
                END IF
             END IF
             node(1) = Nmin ; node(2) = Nmid ; node(3) = Nmax

             cle = alpha*node(1) + beta*node(2) + gamma*node(3)
             ModCle = mod(cle,Hash%SizeH) + 1

             IF ( (node(1) > node(2)) .OR. (node(1) > node(3)) .OR. (node(2) > node(3)) ) THEN
                PRINT*, 'ERROR in hashtable 3D...'
                PRINT*, Mesh%Ele(N_ele)%Vertex(i), Mesh%Ele(N_ele)%Vertex(j), Mesh%Ele(N_ele)%Vertex(k)
                PRINT*, node
                STOP
             END IF

             DO
                IF( Hash%Item(ModCle)%Tri == -1 ) THEN
                   Hash%Item(ModCle)%Tri = 4*N_ele + 9 - (i + j + k)
                   Hash%Item(ModCle)%Vertex(1) = node(1)
                   Hash%Item(ModCle)%Vertex(2) = node(2)
                   Hash%Item(ModCle)%Vertex(3) = node(3)
                   Mesh%Ned = Mesh%Ned + 1
                   ! Note for Leo, allocate in hash asshole
                   Hash%Edg(Mesh%Ned)%Vertex(1) = node(1)
                   Hash%Edg(Mesh%Ned)%Vertex(2) = node(2)
                   Hash%Edg(Mesh%Ned)%Vertex(3) = node(3)
                   Hash%Edg(Mesh%Ned)%Tri(1) = N_ele
                   Hash%EdgNum(ModCle) = Mesh%Ned
                   DO l = 1, 4
                      !IF ( Mesh%Ele(N_ele)%Edg(l) == -1 ) THEN
                         !PRINT*, Mesh%Ele(N_ele)%Vertex(l)
                      !PRINT*, node
                         IF ( (Mesh%Ele(N_ele)%Vertex(l) /= node(1)) .AND. &
                              (Mesh%Ele(N_ele)%Vertex(l) /= node(2)) .AND. &
                              (Mesh%Ele(N_ele)%Vertex(l) /= node(3)) ) THEN


                            Mesh%Ele(N_ele)%Edg(l) = Mesh%Ned
                            EXIT
                         END IF
                      !END IF
                      END DO
                      !PRINT*, " "
                   EXIT
                ELSE
                   nd1 = Hash%Item(ModCle)%Vertex(1)
                   nd2 = Hash%Item(ModCle)%Vertex(2)
                   nd3 = Hash%Item(ModCle)%Vertex(3)
                   IF ( (node(1) == nd1) .AND. (node(2) == nd2) .AND. &
                        (node(3) == nd3)  ) THEN

                      !Mesh%Ele(N_ele)%Adj(10 - (i + j + k)) = Hash%Item(ModCle)%Tri
                      !Mesh%Ele(Hash%Item(ModCle)%Tri/4)%Adj(&
                      !     mod(Hash%Item(ModCle)%Tri,4) + 1) = &
                      !     4*N_ele + 9 - (i + j + k)

                      Hash%Edg(Hash%EdgNum(ModCle))%Tri(2) = N_ele


                      DO l = 1, 4
                         IF ( (Mesh%Ele(N_ele)%Vertex(l) /= nd1) .AND. &
                              (Mesh%Ele(N_ele)%Vertex(l) /= nd2) .AND. &
                              (Mesh%Ele(N_ele)%Vertex(l) /= nd3) ) THEN
                            Mesh%Ele(N_ele)%Edg(l) = Hash%EdgNum(ModCle)
                            Mesh%Ele(N_ele)%Adj(l) = Hash%Item(ModCle)%Tri
                            EXIT
                         END IF
                      END DO
                      DO l = 1, 4
                         TmpEle = Mesh%Ele(Hash%Item(ModCle)%Tri/4)
                         IF ( (TmpEle%Vertex(l) /= nd1) .AND. &
                              (TmpEle%Vertex(l) /= nd2) .AND. &
                              (TmpEle%Vertex(l) /= nd3) ) THEN
                            Mesh%Ele(Hash%Item(ModCle)%Tri/4)%Adj(l) = &
                                 4*N_ele + 9 - (i + j + k)
                         END IF
                      END DO
                      EXIT
                   ELSE
                      IF ( Hash%Item(ModCle)%Nxt == 0 ) THEN
                         Hash%Item(ModCle)%Nxt = Hash%Compt
                         Hash%Compt = Hash%Compt + 1
                         ModCle = Hash%SizeH + Hash%Item(ModCle)%Nxt
                      ELSE
                         ModCle = Hash%SizeH + Hash%Item(ModCle)%Nxt
                      END IF
                   END IF
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE FillHash3D
  !=======================================!

END MODULE Hash
