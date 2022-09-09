        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 30 12:09:48 2019
        MODULE LINBCG__genmod
          INTERFACE 
            SUBROUTINE LINBCG(N,B,X,ITOL,TOL,ITMAX,ITER,ERR,NMAX)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: B(N)
              REAL(KIND=8) :: X(N)
              INTEGER(KIND=4) :: ITOL
              REAL(KIND=8) :: TOL
              INTEGER(KIND=4) :: ITMAX
              INTEGER(KIND=4) :: ITER
              REAL(KIND=8) :: ERR
              INTEGER(KIND=4) :: NMAX
            END SUBROUTINE LINBCG
          END INTERFACE 
        END MODULE LINBCG__genmod
