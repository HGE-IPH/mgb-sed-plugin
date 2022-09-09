	SUBROUTINE SUB(N1,N2,J)
	! Define coluna dos coeficientes das equações dinamica e continuidade 
	! para trechos e energia para confluencias
	!-----------------------------------------------------------------------
	! Descrição das variáveis da entrada:
	!  
	! N1,N2 = indices das secoes de montante e jusante
	! J = indice da linha da matriz AA
	! 
	! Descrição das variáveis locais:
	!
	! M,K,NUP = variáveis auxiliares
	!---------------------------------

	! Declaração de variáveis:
	USE MAT_MOD

	IMPLICIT NONE

	! Variáveis de entrada:
	integer,intent(in):: N1,N2
	integer, intent(inout):: J

	! Variáveis locais de cálculo:
	integer:: M,K,NUP 
	!----------------------------------------------------------------------------

	!Linha da equação
	J=J+1
	JLIN(J)=J
	! Identificador do tipo de equação (continuidade, dinamica ou energia em confluencia)
	ICOL(J,1)=5
	M=2
	NUP=N1

	! Define coluna dos coeficientes:
	DO 
		K=(NUP-1)*2+1
		ICOL(J,M)=K
		ICOL(J,M+1)=K+1
		IF(M==4)EXIT
		NUP=N2
		M=M+2
	ENDDO


	RETURN
	END
