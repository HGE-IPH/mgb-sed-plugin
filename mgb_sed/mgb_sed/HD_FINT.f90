	FUNCTION FINT(X,Y,N,ABC)
	! Subrotina de interpolação
	!-----------------------------------------------------------------------
	! Descrição das variáveis de entrada:
	!
	! X(.) = vetor N x 1 com valores da variável x em ordem crescente
	! Y(.) = vetor N x 1 com valores da variável y
	! N = número de pontos na tabela
	! ABC = valor de x
	! 
	! Descrição das variáveis de saída:
	!
	! FINT = valor interpolado de y
	!
	! Descrição das variáveis locais:
	!
	! I,J = variáveis auxiliares
	!
	!----------------------------------------------------------------------
	!
	IMPLICIT NONE
	! Declaração de variáveis:
	! Variáveis de entrada:
	integer,intent(in)::N
	real, intent(in):: X(N),Y(N),ABC
	! Variáveis de saída:
	real :: FINT
	! Variáveis locais de cálculo:
	integer:: I, J
	!----------------------------------------------------------------------------
	
	
	DO I=2,N
		IF (ABC<X(I)) THEN
			! Interpolação ou extrapolação quando X(1)>ABC
			J=I-1
			FINT=Y(J)+(Y(I)-Y(J))*(ABC-X(J))/(X(I)-X(J))
			EXIT
		ELSEIF (ABC==X(I)) THEN
			FINT=Y(I)		
			EXIT
		ELSEIF (I==N) THEN
			! Extrapolação:
			J=I-1
			FINT=Y(J)+(Y(I)-Y(J))*(ABC-X(J))/(X(I)-X(J))
		ENDIF
	ENDDO

	END FUNCTION FINT

