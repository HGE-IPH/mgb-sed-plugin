	SUBROUTINE INTERPMUSK(QMUSKNL,CELMUSKNL,BMUSKNL)
	!ESTA ROTINA FAZ UMA INTERPOLA??O R?PIDA PARA MUSKINGUN CUNGE NAO LINEAR

	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER KH

	REAL QMUSKNL,CELMUSKNL,BMUSKNL

	DO KH=2,NTMUP-1 !PERCORRE TABELA
		IF(QMUSKNL<QMUP(IC,KH))EXIT
	ENDDO
	!SEMPRE ASSUME QUE O PONTO EST? ENTRE KH-1 E KH (MESMO NOS EXTREMOS)
	!WRITE(*,*)IT,IC,KH,QMUSKNL
	CELMUSKNL=CMUP(IC,KH-1)+(CMUP(IC,KH)-CMUP(IC,KH-1))*(QMUSKNL-QMUP(IC,KH-1))/(QMUP(IC,KH)-QMUP(IC,KH-1))
	BMUSKNL=BMUP(IC,KH-1)+(BMUP(IC,KH)-BMUP(IC,KH-1))*(QMUSKNL-QMUP(IC,KH-1))/(QMUP(IC,KH)-QMUP(IC,KH-1))
	




	RETURN
	END