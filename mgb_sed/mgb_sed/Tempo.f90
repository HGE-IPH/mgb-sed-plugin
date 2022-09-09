	SUBROUTINE TEMPO(IBOBO,ISEED)
	!esta subrotina serve para estimar 
	!o tempo de processamento do programa
	USE PORTLIB
	IMPLICIT NONE
	INTEGER IBOBO
	INTEGER ITIME3,ITIME2,ITIME1 !   CONTADOR DE TEMPO
	INTEGER ISEED
	
	TEMPOS: SELECT CASE(IBOBO)
	CASE(0)   !   INICIA A CONTAGEM DE TEMPO
		itime1 = TIME()
		ISEED  = ITIME1
	CASE(1)	  !   TERMINA CONTAGEM DO TEMPO
	    ITIME1 = ISEED  !@ DCB
		itime2 = TIME()
		ITIME3 = ITIME2-ITIME1
		if(itime3<60)                       print *, 'TEMPO TOTAL ', itime3,' SEGUNDOS'
		if(itime3>=60 .AND. itime3<3600)    print *, 'TEMPO TOTAL ', itime3/60.0,' MINUTOS'
		if(itime3>=3600)                    print *, 'TEMPO TOTAL ', itime3/3600.0,' HORAS'
	END SELECT TEMPOS
	RETURN
	END