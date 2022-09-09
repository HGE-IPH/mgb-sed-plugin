      SUBROUTINE atimes(n,x,r,itrnsp,NMAX)
	  USE SPAMAT ! Modulo com matriz esparsa sa e ponteiro ija
	  IMPLICIT NONE
      INTEGER n,itrnsp,NMAX
!	  DOUBLE PRECISION x(n),r(n),SA(NMAX)
      REAL(8) x(n),r(n)

 
 

      if (itrnsp.eq.0) then
        call dsprsax(sa,ija,x,r,n,NMAX)
      else
        call dsprstx(sa,ija,x,r,n,NMAX)
      endif
      return
      END