      SUBROUTINE Dsprsax(sa,ija,x,b,n,NMAX)
	  INTEGER NMAX
      INTEGER n,ija(NMAX)
      !DOUBLE PRECISION b(n),sa(NMAX),x(n)
	  REAL(8) b(n),sa(NMAX),x(n)
      INTEGER i,k

      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in sprsax'
      do  i=1,n
        b(i)=sa(i)*x(i)
        do k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
	    enddo
	  enddo
      return
      END