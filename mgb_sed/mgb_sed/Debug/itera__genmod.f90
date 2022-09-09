        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 30 12:10:01 2019
        MODULE ITERA__genmod
          INTERFACE 
            SUBROUTINE ITERA(J1,J2,ARK1,ARK2,DXAUX)
              INTEGER(KIND=4), INTENT(IN) :: J2
              INTEGER(KIND=4), INTENT(IN) :: J1
              REAL(KIND=4), INTENT(IN) :: ARK1(NP((J1)))
              REAL(KIND=4), INTENT(IN) :: ARK2(NP((J2)))
              REAL(KIND=4), INTENT(IN) :: DXAUX
            END SUBROUTINE ITERA
          END INTERFACE 
        END MODULE ITERA__genmod
