      FUNCTION snrm(n,sx,itol,NMAX)
	  INTEGER NMAX
	  LOGICAL valor
      INTEGER n,itol,i,isamax
      !DOUBLE PRECISION sx(NMAX),snrm
	  REAL(8) sx(NMAX),snrm
      if (itol.le.3)then
        snrm=0.
        do 11 i=1,n
			snrm=snrm+sx(i)**2
11      continue
		valor = isnan(snrm)
		if (valor == .TRUE.) THEN
		snrm = 0.0
		else
		snrm = SQRT(snrm) !SRQT
		endif
      else
        isamax=1
        do 12 i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12      continue
        snrm=abs(sx(isamax))
      endif
      return
      END