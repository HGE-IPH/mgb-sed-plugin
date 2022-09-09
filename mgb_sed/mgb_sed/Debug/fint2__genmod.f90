        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 30 12:09:47 2019
        MODULE FINT2__genmod
          INTERFACE 
            FUNCTION FINT2(X,Y,Q,NX,NY,X1,Y1) RESULT(FINT2_0)
              INTEGER(KIND=4), INTENT(IN) :: NY
              INTEGER(KIND=4), INTENT(IN) :: NX
              REAL(KIND=4), INTENT(IN) :: X(NX)
              REAL(KIND=4), INTENT(IN) :: Y(NY)
              REAL(KIND=4), INTENT(IN) :: Q(NX,NY)
              REAL(KIND=4), INTENT(IN) :: X1
              REAL(KIND=4), INTENT(IN) :: Y1
              REAL(KIND=4) :: FINT2_0
            END FUNCTION FINT2
          END INTERFACE 
        END MODULE FINT2__genmod
