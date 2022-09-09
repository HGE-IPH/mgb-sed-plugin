        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 30 12:09:45 2019
        MODULE INTEP__genmod
          INTERFACE 
            SUBROUTINE INTEP(QAUX,DT1,NT,Q,IFIN)
              INTEGER(KIND=4), INTENT(IN) :: IFIN
              INTEGER(KIND=4), INTENT(IN) :: NT
              REAL(KIND=4), INTENT(IN) :: QAUX(NT)
              REAL(KIND=4), INTENT(IN) :: DT1
              REAL(KIND=4), INTENT(INOUT) :: Q(IFIN)
            END SUBROUTINE INTEP
          END INTERFACE 
        END MODULE INTEP__genmod
