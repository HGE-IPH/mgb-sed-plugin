	SUBROUTINE LEPREV
	!esta subrotina le os arquivos DE PREVISÃO DE CHUVA
	!onde estao os dados de chuva para cada célula do modelo
	
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER KC
	!write(*,*) itchuva


	!LE CHUVA DA CELULA DO INTERVALO DE TEMPO ATUAL
	READ(FILPREV,REC=ITCHUVA)(P(KC),KC=1,NC)


	RETURN
	END