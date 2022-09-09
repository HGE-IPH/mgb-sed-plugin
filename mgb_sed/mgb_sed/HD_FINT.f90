	FUNCTION FINT(X,Y,N,ABC)
	! Subrotina de interpola��o
	!-----------------------------------------------------------------------
	! Descri��o das vari�veis de entrada:
	!
	! X(.) = vetor N x 1 com valores da vari�vel x em ordem crescente
	! Y(.) = vetor N x 1 com valores da vari�vel y
	! N = n�mero de pontos na tabela
	! ABC = valor de x
	! 
	! Descri��o das vari�veis de sa�da:
	!
	! FINT = valor interpolado de y
	!
	! Descri��o das vari�veis locais:
	!
	! I,J = vari�veis auxiliares
	!
	!----------------------------------------------------------------------
	!
	IMPLICIT NONE
	! Declara��o de vari�veis:
	! Vari�veis de entrada:
	integer,intent(in)::N
	real, intent(in):: X(N),Y(N),ABC
	! Vari�veis de sa�da:
	real :: FINT
	! Vari�veis locais de c�lculo:
	integer:: I, J
	!----------------------------------------------------------------------------
	
	
	DO I=2,N
		IF (ABC<X(I)) THEN
			! Interpola��o ou extrapola��o quando X(1)>ABC
			J=I-1
			FINT=Y(J)+(Y(I)-Y(J))*(ABC-X(J))/(X(I)-X(J))
			EXIT
		ELSEIF (ABC==X(I)) THEN
			FINT=Y(I)		
			EXIT
		ELSEIF (I==N) THEN
			! Extrapola��o:
			J=I-1
			FINT=Y(J)+(Y(I)-Y(J))*(ABC-X(J))/(X(I)-X(J))
		ENDIF
	ENDDO

	END FUNCTION FINT

