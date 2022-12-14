	SUBROUTINE INTEP(QAUX,DT1,NT,Q,IFIN)
	! Passa serie temporal de DT para DT1 sendo DT1<DT
	! Interpola dados de series temporais para refinar discretiza??o temporal
	!-----------------------------------------------------------------------
	! Descri??o das vari?veis de entrada:
	!
	! QAUX(.) = intervalos de tempo dos dados em segundos (NT x 1)
	! DT1 = intervalo de tempo de calculo (s) 
	! NT = numero de intervalos de tempo dos dados
	! IFIN = numero de intervalos de tempo de calculo
	! 
	! Descri??o das vari?veis de sa?da:
	!
	! Q(.) = serie temporal que entra com intervalo de tempo dos dados
	!		 e sai com intevalo de tempo de calculo (IFIN x 1)
	!
	! Descri??o das vari?veis locais:
	!
	! QAUX2(.) = armazena serie temporal Q no intervalo de tempo original (NT x 1)
	! ACON = intervalo de tempo em segundos
	! K,IA = auxiliares
	! FINT = auxiliar rotina de interpola??o
	!----------------------------------------------------------------------
	!
	IMPLICIT NONE
	! Declara??o de vari?veis:
	! Vari?veis de entrada:
	integer,intent(in)::NT,IFIN
	real, intent(in):: QAUX(NT),DT1

	! Vari?veis de sa?da:
	real, intent(inout)::Q(IFIN)

	! Vari?veis locais de c?lculo:
	real:: QAUX2(NT),ACON,FINT
	integer::K,IA
	!----------------------------------------------------------------------------

	! Armazena serie temporal Q em QAUX2
	QAUX2=Q(1:NT)

	ACON=0

	DO IA=2,IFIN
		! Calcula intervalo de tempo de calculo em segundos:
		ACON=(IA-1)*DT1
		! Interpola Q correspondente
		Q(IA)=FINT(QAUX,QAUX2,NT,ACON)

	ENDDO


	RETURN
	END
