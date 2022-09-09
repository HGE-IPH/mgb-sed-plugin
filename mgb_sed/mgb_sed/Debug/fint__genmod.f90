        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 30 12:09:48 2019
        MODULE FINT__genmod
          INTERFACE 
            FUNCTION FINT(X,Y,N,ABC) RESULT(FINT_0)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=4), INTENT(IN) :: X(N)
              REAL(KIND=4), INTENT(IN) :: Y(N)
              REAL(KIND=4), INTENT(IN) :: ABC
              REAL(KIND=4) :: FINT_0
            END FUNCTION FINT
          END INTERFACE 
        END MODULE FINT__genmod
