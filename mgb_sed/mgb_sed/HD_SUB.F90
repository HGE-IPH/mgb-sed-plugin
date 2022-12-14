	SUBROUTINE SUB(N1,N2,J)
	! Define coluna dos coeficientes das equa??es dinamica e continuidade 
	! para trechos e energia para confluencias
	!-----------------------------------------------------------------------
	! Descri??o das vari?veis da entrada:
	!  
	! N1,N2 = indices das secoes de montante e jusante
	! J = indice da linha da matriz AA
	! 
	! Descri??o das vari?veis locais:
	!
	! M,K,NUP = vari?veis auxiliares
	!---------------------------------

	! Declara??o de vari?veis:
	USE MAT_MOD

	IMPLICIT NONE

	! Vari?veis de entrada:
	integer,intent(in):: N1,N2
	integer, intent(inout):: J

	! Vari?veis locais de c?lculo:
	integer:: M,K,NUP 
	!----------------------------------------------------------------------------

	!Linha da equa??o
	J=J+1
	JLIN(J)=J
	! Identificador do tipo de equa??o (continuidade, dinamica ou energia em confluencia)
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
