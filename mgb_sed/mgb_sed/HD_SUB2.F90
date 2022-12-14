	SUBROUTINE SUB2(K1,K2,K3,I)
	! Define posi??o dos coeficientes n?o nulos de uma linha (reorganiza linhas)
	!-----------------------------------------------------------------------
	! Descri??o das vari?veis da entrada:
	!  
	! K1,K2,K3 = indices das vari?veis de estado correspondentes a eq. das confluencias 
	!                    (vaz?o ou profundidade de determinada secao)
	! I = n?mero da confluencia
	!
	! Descri??o das vari?veis locais:
	!
	! K(.) = vetor 3 x 1 que armazena K1,K2 e K3
	! KJ = numero da linha na organizacao inicial de ICOL
	! KK = numero da linha na organizacao final de ICOL
	! MH = numero de coeficientes nao nulos da equa??o + 1
	! M,J = contadores
	!---------------------------------

	! Declara??o de vari?veis:
	USE PAR1_MOD
	USE MAT_MOD

	IMPLICIT NONE

	! Vari?veis de entrada:
	integer,intent(in):: K1,K2,K3,I

	! Vari?veis locais de c?lculo:
	integer:: K(3),KJ,KK,MH,M,J 
	!----------------------------------------------------------------------------

	! Inicializa vars.:
	K(1)=K1
	K(2)=K2
	K(3)=K3
	! Linha na organizacao inicial
	KJ=NBOUN+NREAC*2+(I-1)*3
	DO J=1,3
		KJ=KJ+1 ! Linha da equa??o da confluencia:

		! Identifica linha na organizacao final (igual a coluna do coeficinte K(J)):
		KK=K(J)
		JLIN(KJ)=KK
		MH=ICOL(KJ,1) ! Tipo de equa??o (continuidade ou energia) e numero de coeficientes + 1
		
		! Copia posi??o (coluna) dos coeficientes em ICOL para vari?vel ICAUX
		DO M=1,MH
			ICAUX(KK,M)=ICOL(KJ,M)
		ENDDO
	ENDDO

	RETURN
	END
