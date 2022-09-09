      SUBROUTINE asolve(n,b,x,itrnsp,NMAX)
	  USE SPAMAT ! Modulo com matriz esparsa sa e ponteiro ija
	  IMPLICIT NONE
      INTEGER n,itrnsp,i,NMAX
      !DOUBLE PRECISION x(n),b(n),SA(NMAX)
	  REAL(8) x(n),b(n)

      !PARAMETER (NMAX=1000)
      !COMMON /mat/ sa(NMAX),ija(NMAX)
      do i=1,n
!        x(i)=b(i)/1.0
		x(i)=b(i)/sa(i)
	  enddo
      return
      END