MODULE MODULE_CREATE_SURFACE
    
    TYPE SURFACE
        INTEGER :: NG     ! NUMBER OF GRIDS
        INTEGER :: NG1    ! NUMBER OF GRIDS - SPAN DIRECTION
        INTEGER :: NG2    ! NUMBER OF GRIDS - CHORD DIRECTION
        INTEGER :: N, M
        INTEGER :: NP     ! NUMBER OF PANELS
        INTEGER :: ID
        REAL*8, ALLOCATABLE :: GRIDS(:,:)      ! NUMBER OF GRIDS | X,Y,Z,SPAN,CHORD
        INTEGER, ALLOCATABLE :: CN(:,:)        ! NUMBER OF PANELS - CONNECTIVITY
        REAL*8, ALLOCATABLE :: PANELS(:,:)     ! NUMBER OF PANELS - X,Y,Z,AXY,AXZ,AYZ
        REAL*8, ALLOCATABLE :: CDEF(:), SDEF(:)
    END TYPE SURFACE
    
    TYPE COMPONENT
        INTEGER :: NC
        TYPE(SURFACE), ALLOCATABLE :: SUF(:)
        INTEGER :: NG, NE
        INTEGER, ALLOCATABLE :: E(:,:)
        REAL*8, ALLOCATABLE :: G(:,:)
    END TYPE COMPONENT
    
    
    TYPE CONFIG
        REAL*8, ALLOCATABLE :: LE(:,:)  ! NUMBER OF POINTS DEFINING THE LEADING EDGE | X,Y,Z,CHORD,INCIDENCE
        REAL*8, ALLOCATABLE :: TE(:,:)
        INTEGER :: NLE
        INTEGER :: DIRS, MDIR
    END TYPE CONFIG
    
    CONTAINS
    
    !============================================================
    SUBROUTINE WRITE_VTK(COMP)
        TYPE(COMPONENT), INTENT(INOUT) :: COMP
        INTEGER :: I, J, SG, S, OG
        
        OPEN(UNIT=40, FILE='surface.vtk', STATUS='UNKNOWN')
        
        SG = 0
        SE = 0
        OG = 0
        
        DO I=1, COMP%NC
            DO J=1, COMP%SUF(I)%NG
                SG = SG + 1
                COMP%G(SG, :) = COMP%SUF(I)%GRIDS(J, :)
            END DO
            DO J=1, COMP%SUF(I)%NP
                SE = SE + 1
                COMP%E(SE, 1) = COMP%SUF(I)%CN(J, 1) + OG
                COMP%E(SE, 2) = COMP%SUF(I)%CN(J, 2) + OG
                COMP%E(SE, 3) = COMP%SUF(I)%CN(J, 3) + OG
                COMP%E(SE, 4) = COMP%SUF(I)%CN(J, 4) + OG
            END DO
            OG = OG + COMP%SUF(I)%NG
        END DO
        
        
        WRITE(40, '("# vtk DataFile Version 2.0")')
        WRITE(40, '("surface.vtk")')
        WRITE(40, '("ASCII")')
        WRITE(40, '("DATASET UNSTRUCTURED_GRID")')
        !WRITE(40, *)
        
        WRITE(40, '("POINTS",I8," float")')COMP%NG
        
        DO I=1, COMP%NG
            WRITE(40,'(3F12.6)')COMP%G(I,1), COMP%G(I,2), COMP%G(I,3)
        END DO
        
        WRITE(40,'("CELLS ",I8,I8)')COMP%NE, COMP%NE*5
        
        DO I=1, COMP%NE
            WRITE(40,'(I2,4I8)')4, COMP%E(I,1)-1, COMP%E(I,2)-1, COMP%E(I,3)-1, COMP%E(I,4)-1
        END DO

        !WRITE(40, *)
        
        WRITE(40,'("CELL_TYPES ",I8)')COMP%NE
        
        DO I=1, COMP%NE
            WRITE(40,'(I1)')9
        END DO
        
        
        CLOSE(40)
    
    END SUBROUTINE
    
    !============================================================    
    SUBROUTINE DEF_SURFACE(C, S, OFFSET_GRIDS, OFFSET_ELEMENTS)
    
        TYPE(CONFIG), INTENT(INOUT) :: C
        TYPE(SURFACE), INTENT(INOUT) :: S
        INTEGER, INTENT(IN) :: OFFSET_GRIDS, OFFSET_ELEMENTS  
        REAL*8 :: PT(3,2), LCHORD
        INTEGER :: I, J, L
        INTEGER :: K
        
        ALLOCATE(S%PANELS(S%NP,6))
        ALLOCATE(S%CN(S%NP,4))
        ALLOCATE(S%GRIDS(S%NG ,5))
        
        OPEN(UNIT=21, FILE='grids.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        OPEN(UNIT=22, FILE='elements.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        OPEN(UNIT=23, FILE='panels_info.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        OPEN(UNIT=24, FILE='grids_info.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
                
        K = 0
        
        ! Grids
        DO I=1, S%NG2   ! chordwise
            DO J=1, S%NG1  ! spanwise
                K = K + 1
                S%GRIDS(K, 4) = J
                S%GRIDS(K, 5) = I
                IF (C%DIRS == 2) THEN
                    
                    S%GRIDS(K, 2) = C%LE(1,2) + (C%LE(C%NLE, 2) - C%LE(1, 2)) * S%SDEF(J)

                    PT(2,1) = S%GRIDS(K, 2)
                    PT(2,2) = S%GRIDS(K, 2)
                    
                    ! Point at the Leading Edge
                    CALL LININTERP(C%LE(:,2),C%LE(:,1),C%NLE,PT(2,1),PT(1,1))
                    CALL LININTERP(C%LE(:,2),C%LE(:,3),C%NLE,PT(2,1),PT(3,1))
                    ! Point at the Trailing Edge
                    CALL LININTERP(C%TE(:,2),C%TE(:,1),C%NLE,PT(2,2),PT(1,2))
                    CALL LININTERP(C%TE(:,2),C%TE(:,3),C%NLE,PT(2,2),PT(3,2))
                    
                    S%GRIDS(K, 1) = PT(1,1) + (PT(1,2) - PT(1,1)) * S%CDEF(I)

                    CALL LININTERP(PT(1,:),PT(3,:),2,S%GRIDS(K, 1),S%GRIDS(K, 3))
                ELSE IF (C%DIRS == 3) THEN

                    S%GRIDS(K, 3) = C%LE(1,3) + (C%LE(NLE, 3) - C%LE(1, 3)) * S%SDEF(J)

                    PT(3,1) = S%GRIDS(K, 3)
                    PT(3,2) = S%GRIDS(K, 3)
                    ! Point at the Leading Edge
                    CALL LININTERP(C%LE(:,3),C%LE(:,1),C%NLE,PT(3,1),PT(1,1))
                    CALL LININTERP(C%LE(:,3),C%LE(:,2),C%NLE,PT(3,1),PT(2,1))
                    ! Point at the Trailing Edge
                    CALL LININTERP(C%TE(:,3),C%TE(:,1),C%NLE,PT(3,2),PT(1,2))
                    CALL LININTERP(C%TE(:,3),C%TE(:,2),C%NLE,PT(3,2),PT(2,2))

                    S%GRIDS(K, 1) = PT(1,1) + (PT(1,2) - PT(1,1)) * S%CDEF(I)

                    CALL LININTERP(PT(1,:),PT(2,:),2,S%GRIDS(K, 1),S%GRIDS(K, 2))
                END IF
                WRITE(24,'(2I8,3F12.4)')S%ID, K+OFFSET_GRIDS,  S%GRIDS(K, 1), S%GRIDS(K, 2), S%GRIDS(K, 3)
                WRITE(21,'(I8,3F12.4)')K+OFFSET_GRIDS, S%GRIDS(K, 1), S%GRIDS(K, 2), S%GRIDS(K, 3)
            END DO
        END DO
        
        K = 0 ! Spanwise
        J = 1 ! Chordwise
        !Connectivity
        DO I = 1, S%NP
            IF (K == S%N) THEN !Span strip
                K = 1
                J = J + 1
            ElSE
                K = K + 1
            END IF
            IF (C%MDIR == 1) THEN
                S%CN(I,1) = I + (J - 1)
                S%CN(I,2) = S%NG1 + I + (J - 1)
                S%CN(I,3) = S%NG1 + I + 1 + (J - 1)
                S%CN(I,4) = I + 1 + (J - 1)
            ELSE IF (C%MDIR == 0) THEN
                S%CN(I,1) = I + (J - 1)
                S%CN(I,4) = S%NG1 + I + (J - 1)
                S%CN(I,3) = S%NG1 + I + 1 + (J - 1)
                S%CN(I,2) = I + 1 + (J - 1)
            END IF
            WRITE(22, '(4I8)')S%CN(I,1)+OFFSET_GRIDS,S%CN(I,2)+OFFSET_GRIDS,S%CN(I,3)+OFFSET_GRIDS,S%CN(I,4)+OFFSET_GRIDS
        END DO
        
        
        !Panels
        DO I = 1, S%NP
            S%PANELS(I,1) = 0.25*(S%GRIDS(S%CN(I,1),1) + S%GRIDS(S%CN(I,2),1) + S%GRIDS(S%CN(I,3),1) + S%GRIDS(S%CN(I,4),1))
            S%PANELS(I,2) = 0.25*(S%GRIDS(S%CN(I,1),2) + S%GRIDS(S%CN(I,2),2) + S%GRIDS(S%CN(I,3),2) + S%GRIDS(S%CN(I,4),2))
            S%PANELS(I,3) = 0.25*(S%GRIDS(S%CN(I,1),3) + S%GRIDS(S%CN(I,2),3) + S%GRIDS(S%CN(I,3),3) + S%GRIDS(S%CN(I,4),3))
            CALL AREA_QUAD(S%GRIDS(S%CN(I,1),1), S%GRIDS(S%CN(I,2),1), S%GRIDS(S%CN(I,3),1), S%GRIDS(S%CN(I,4),1),  &
                           S%GRIDS(S%CN(I,1),2), S%GRIDS(S%CN(I,2),2), S%GRIDS(S%CN(I,3),2), S%GRIDS(S%CN(I,4),2),  &
                           S%GRIDS(S%CN(I,1),3), S%GRIDS(S%CN(I,2),3), S%GRIDS(S%CN(I,3),3), S%GRIDS(S%CN(I,4),3),  &
                           S%PANELS(I,4), S%PANELS(I,5), S%PANELS(I,6))
            
            WRITE(23,'(3I8,3F12.4,3F12.8)')S%ID, I+OFFSET_ELEMENTS, I, S%PANELS(I,1:6)
            
        END DO
        
        
        CLOSE(21)
        CLOSE(22)
        CLOSE(23)
        CLOSE(24)
    
    END SUBROUTINE
    !============================================================
    

    !============================================================
    SUBROUTINE SET_CONFIG(C, COMP)
    
        TYPE(CONFIG), INTENT(INOUT) :: C
        TYPE(COMPONENT), INTENT(INOUT) :: COMP
        INTEGER :: I, J
        REAL*8 :: DEG2RAD
        
        DEG2RAD = ACOS(-1.0) / 180
        
        OPEN(UNIT=10, FILE='config.txt', STATUS='OLD')

        OPEN(UNIT=21, FILE='grids.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        OPEN(UNIT=22, FILE='elements.txt', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')      
        
        READ(10,*)C%MDIR
        READ(10,*)C%DIRS

        READ(10,*)C%NLE        
        ALLOCATE(C%LE(C%NLE, 5))
        ALLOCATE(C%TE(C%NLE, 3))
        
        DO I=1,C%NLE
            
            READ(10,*)(C%LE(I,J),J=1,5)
            
            C%TE(I,1) = C%LE(I,1) + C%LE(I,4)
            IF (C%DIRS == 2) THEN
                C%TE(I,2) = C%LE(I,2)
                IF (C%LE(I,5) .NE. 0) THEN
                    C%TE(I,3) = C%LE(I,3) - C%LE(I,4) * TAN(C%LE(I,5) * DEG2RAD)
                ELSE
                    C%TE(I,3) = C%LE(I,3)
                END IF
            ELSE IF (C%DIRS == 3) THEN
                C%TE(I,3) = C%LE(I,3)
                IF (C%LE(I,5) .NE. 0) THEN
                    C%TE(I,2) = C%LE(I,2) - C%LE(I,4) * TAN(C%LE(I,5) * DEG2RAD) 
                ELSE
                    C%TE(I,2) = C%LE(I,2)
                END IF
            END IF
                        
        END DO
            
        READ(10, *)COMP%NC
        ALLOCATE(COMP%SUF(COMP%NC))
        
        COMP%NG = 0
        COMP%NE = 0
        
        DO I = 1, COMP%NC
            READ(10, *) COMP%SUF(I)%ID
            READ(10, *) COMP%SUF(I)%N
            READ(10, *) COMP%SUF(I)%M
            
            COMP%SUF(I)%NG2 = COMP%SUF(I)%M + 1
            COMP%SUF(I)%NG1 = COMP%SUF(I)%N + 1
            COMP%SUF(I)%NG = COMP%SUF(I)%NG1 * COMP%SUF(I)%NG2 
            COMP%SUF(I)%NP = COMP%SUF(I)%M * COMP%SUF(I)%N
 
            COMP%NG = COMP%NG + COMP%SUF(I)%NG
            COMP%NE = COMP%NE + COMP%SUF(I)%NP
            
            ALLOCATE(COMP%SUF(I)%CDEF(COMP%SUF(I)%NG2))
            ALLOCATE(COMP%SUF(I)%SDEF(COMP%SUF(I)%NG1))
            
            READ(10, *)(COMP%SUF(I)%CDEF(J), J = 1, COMP%SUF(I)%NG2)
            READ(10, *)(COMP%SUF(I)%SDEF(J), J = 1, COMP%SUF(I)%NG1)
            
        END DO
        
        WRITE(21,* )COMP%NG
        WRITE(22,* )COMP%NE
        
        ALLOCATE(COMP%G(COMP%NG, 3))
        ALLOCATE(COMP%E(COMP%NE, 4))
        
        
        CLOSE(10)
        CLOSE(21)
        CLOSE(22)

    
    END SUBROUTINE
    !============================================================

    
    
    !===================================================
    !
    ! LINEAR INTERPOLATION
    !
    !==================================================   
    SUBROUTINE LININTERP(X,Y,NP,XEX,YS)
         
         INTEGER,INTENT(IN) :: NP
         REAL*8,INTENT(IN) :: X(NP), Y(NP), XEX
         REAL*8,INTENT(OUT) :: YS


         IF(X(2) > X(1)) THEN                                                                                      
            IF(XEX < X(1)) THEN
              !YS = Y(1)
              YS = (Y(2)-Y(1))*(XEX-X(1))/(X(2)-X(1)) + Y(1)                                                         
            ELSE                                                                                                     
               IF(XEX >= X(NP)) THEN 
                  !YS = Y(NP)
                  YS = (Y(NP)-Y(NP-1))*(XEX-X(NP))/(X(NP)-X(NP-1))+Y(NP)                                               
               ELSE                                                                                                   
                  DO I = 2, NP                                                                                         
                     IF(XEX <= X(I)) THEN                                                                               
                        YS = (Y(I)-Y(I-1))*(XEX-X(I-1))/(X(I)-X(I-1))+Y(I-1)                                             
                        EXIT                                                                    
                     ENDIF                                                                                              
                  ENDDO                                                                                                
               ENDIF                                                                                                  
            ENDIF                                                                                                    
         ELSE                                                                                                       
            IF(XEX > X(1)) THEN
                !YS = Y(1)
                YS  = (Y(2)-Y(1))*(XEX-X(1))/(X(2)-X(1)) + Y(1)                                                        
            ELSE                                                                                                     
                IF(XEX <= X(NP)) THEN
                   !YS = Y(NP)
                   YS  = (Y(NP)-Y(NP-1))*(XEX-X(NP))/(X(NP)-X(NP-1))+Y(NP)                                              
                ELSE                                                                                                   
                   DO I = 2, NP                                                                                         
                       IF(XEX >= X(I)) THEN                                                                               
                           YS = (Y(I)-Y(I-1))*(XEX-X(I-1))/(X(I)-X(I-1))+Y(I-1)                                             
                           EXIT                                                                    
                       ENDIF                                                                                              
                   ENDDO                                                                                                
                ENDIF                                                                                                  
            ENDIF                                                                                                    
         ENDIF  

    END SUBROUTINE
    !==================================================
    
    !==================================================
    SUBROUTINE AREA_QUAD(X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, AXY, AXZ, AYZ)
    
        REAL*8, INTENT(IN) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4
        REAL*8, INTENT(INOUT) :: AXY, AXZ, AYZ
        REAL*8 :: A1, A2, A3, B1, B2, B3, VX1, VY1, VZ1
        
    !     THIS ROUTINE CALCULATE THE PANEL AREAS

    !     PRIMARY DIAGONAL (1-3)
          A1 = X3 - X1
          A2 = Y3 - Y1
          A3 = Z3 - Z1

    !      SECUNDARY DIAGONAL (2-4)
          B1 = X4 - X2
          B2 = Y4 - Y2
          B3 = Z4 - Z2

    !     VECTORIAL PRODUCT (A X B)
          VX1 = (A2 * B3 - A3 * B2)
          VY1 = (A3 * B1 - A1 * B3)
          VZ1 = (A1 * B2 - A2 * B1)
          
    !     PROJECTED AREAS AT EACH PLANE
          AXY = 0.5 * VZ1
          AXZ = 0.5 * VY1
          AYZ = 0.5 * VX1
          
    END SUBROUTINE
    !==================================================
    
    
    !==================================================
    SUBROUTINE ROT(PT, SWP, DIE, PROT)
    
        REAL*8, INTENT(INOUT) :: PROT(3)
        REAL*8, INTENT(IN) :: PT(3), DIE, SWP
        REAL*8 :: PI, DEG2RAD
        
        PI = ACOS(-1)
        DEG2RAD = PI / 180
        
        PROT(1) = PT(1) * COS(DEG2RAD * SWP) + PT(2) * SIN(DEG2RAD * SWP) * COS(DEG2RAD * DIE) - PT(3) * SIN(DEG2RAD * DIE) * SIN(DEG2RAD * SWP)
        PROT(2) = - PT(1) * SIN(DEG2RAD * SWP) + PT(2) * COS(DEG2RAD * SWP) * COS(DEG2RAD * DIE) - PT(3) * SIN(DEG2RAD * DIE) * COS(DEG2RAD * SWP)
        PROT(3) = PT(2) * SIN(DEG2RAD * DIE) + PT(3) * COS(DEG2RAD * DIE)
        
    END SUBROUTINE
    !==================================================
    
    
END MODULE