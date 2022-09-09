        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 30 12:09:45 2019
        MODULE DSPRSTX__genmod
          INTERFACE 
            SUBROUTINE DSPRSTX(SA,IJA,X,B,N,NMAX)
              INTEGER(KIND=4) :: NMAX
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: SA(NMAX)
              INTEGER(KIND=4) :: IJA(NMAX)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: B(N)
            END SUBROUTINE DSPRSTX
          END INTERFACE 
        END MODULE DSPRSTX__genmod
