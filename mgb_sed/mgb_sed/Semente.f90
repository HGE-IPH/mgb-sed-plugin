	SUBROUTINE SEMENTE(ISEED)
	!usa relogio do computador para gerar semente do processo aleatorio
	USE PORTLIB
	INTEGER ISEED
	ISEED=TIME()
	RETURN
	END
