      SUBROUTINE Dsprstx(sa,ija,x,b,n,NMAX)
      INTEGER NMAX
	  INTEGER n,ija(NMAX)
	  !DOUBLE PRECISION b(n),sa(NMAX),x(n)
      REAL(8) b(n),sa(NMAX),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) pause 'mismatched vector and matrix in sprstx'
      do i=1,n
        b(i)=sa(i)*x(i)
      enddo
      do i=1,n
        do k=ija(i),ija(i+1)-1
          j=ija(k)
          b(j)=b(j)+sa(k)*x(i)
        enddo
	  enddo
      return
      END