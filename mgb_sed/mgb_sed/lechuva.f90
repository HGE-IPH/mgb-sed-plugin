	SUBROUTINE LECHUVA
	!esta subrotina le os arquivos chuvabin.hig
	!onde estao os dados de chuva 
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER KC


	!LE CHUVA DA CELULA DO INTERVALO DE TEMPO ATUAL
	READ(FILPLU,REC=ITCHUVA)(P(KC),KC=1,NC)
    !write(*,*)p(204)


	RETURN
	END