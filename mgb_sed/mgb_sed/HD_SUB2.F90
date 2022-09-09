	SUBROUTINE SUB2(K1,K2,K3,I)
	! Define posição dos coeficientes não nulos de uma linha (reorganiza linhas)
	!-----------------------------------------------------------------------
	! Descrição das variáveis da entrada:
	!  
	! K1,K2,K3 = indices das variáveis de estado correspondentes a eq. das confluencias 
	!                    (vazão ou profundidade de determinada secao)
	! I = número da confluencia
	!
	! Descrição das variáveis locais:
	!
	! K(.) = vetor 3 x 1 que armazena K1,K2 e K3
	! KJ = numero da linha na organizacao inicial de ICOL
	! KK = numero da linha na organizacao final de ICOL
	! MH = numero de coeficientes nao nulos da equação + 1
	! M,J = contadores
	!---------------------------------

	! Declaração de variáveis:
	USE PAR1_MOD
	USE MAT_MOD

	IMPLICIT NONE

	! Variáveis de entrada:
	integer,intent(in):: K1,K2,K3,I

	! Variáveis locais de cálculo:
	integer:: K(3),KJ,KK,MH,M,J 
	!----------------------------------------------------------------------------

	! Inicializa vars.:
	K(1)=K1
	K(2)=K2
	K(3)=K3
	! Linha na organizacao inicial
	KJ=NBOUN+NREAC*2+(I-1)*3
	DO J=1,3
		KJ=KJ+1 ! Linha da equação da confluencia:

		! Identifica linha na organizacao final (igual a coluna do coeficinte K(J)):
		KK=K(J)
		JLIN(KJ)=KK
		MH=ICOL(KJ,1) ! Tipo de equação (continuidade ou energia) e numero de coeficientes + 1
		
		! Copia posição (coluna) dos coeficientes em ICOL para variável ICAUX
		DO M=1,MH
			ICAUX(KK,M)=ICOL(KJ,M)
		ENDDO
	ENDDO

	RETURN
	END
