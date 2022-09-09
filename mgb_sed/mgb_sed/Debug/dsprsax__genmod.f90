        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 30 12:09:46 2019
        MODULE DSPRSAX__genmod
          INTERFACE 
            SUBROUTINE DSPRSAX(SA,IJA,X,B,N,NMAX)
              INTEGER(KIND=4) :: NMAX
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: SA(NMAX)
              INTEGER(KIND=4) :: IJA(NMAX)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: B(N)
            END SUBROUTINE DSPRSAX
          END INTERFACE 
        END MODULE DSPRSAX__genmod
