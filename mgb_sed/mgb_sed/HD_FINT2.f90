	FUNCTION FINT2 (X,Y,Q,NX,NY,X1,Y1)
	! Subrotina de interpola��o
	!-----------------------------------------------------------------------
	! Descri��o das vari�veis de entrada:
	!
	! X(.) = vetor NX x 1 com valores da vari�vel x em ordem crescente
	! Y(.) = vetor NY x 1 com valores da vari�vel y em ordem crescente
	! Q(.,.) = matriz NX x NY da vari�vel dependente Q
	! NX = n�mero de pontos no vetor X e linhas da tabela Q
	! NY = n�mero de pontos no vetor Y e colunas da tabela Q
	! X1 = valor procurado de x
	! Y1 = valor procurado de y
	! 
	! Descri��o das vari�veis de sa�da:
	!
	! FINT = valor interpolado de Q
	!
	! Descri��o das vari�veis locais:
	!
	! I,L,M = vari�veis auxiliares
	! QA,QB,AUX = vari�veis auxiliares para interpola��o 
	!----------------------------------------------------------------------
	!
	IMPLICIT NONE
	! Declara��o de vari�veis:
	! Vari�veis de entrada:
	integer,intent(in)::NX,NY
	real, intent(in):: X(NX),Y(NY),Q(NX,NY),X1,Y1
	! Vari�veis de sa�da:
	real :: FINT2
	! Vari�veis locais de c�lculo:
	integer:: I,L,M
	real:: QA,QB,AUX
	!----------------------------------------------------------------------------
	
	! Verifica localiza��o de X1:
	!DO I=2,NX-1
	DO I=2,NX
		IF(X1<=X(I))EXIT
	ENDDO 
	L=I
	! Verifica localiza��o de Y1:
	!DO I=2,NY-1
	DO I=2,NY
		IF(Y1<=Y(I))EXIT
	ENDDO
	M=I

	! Interpola valores:
	AUX=(X1-X(L-1))/(X(L)-X(L-1))
	QA=Q(L-1,M-1)+AUX*(Q(L,M-1)-Q(L-1,M-1))
	QB=Q(L-1,M)+AUX*(Q(L,M)-Q(L-1,M))

	FINT2=QA+(Y1-Y(M-1))*(QB-QA)/(Y(M)-Y(M-1))

	END FUNCTION FINT2