	FUNCTION FINT2 (X,Y,Q,NX,NY,X1,Y1)
	! Subrotina de interpolação
	!-----------------------------------------------------------------------
	! Descrição das variáveis de entrada:
	!
	! X(.) = vetor NX x 1 com valores da variável x em ordem crescente
	! Y(.) = vetor NY x 1 com valores da variável y em ordem crescente
	! Q(.,.) = matriz NX x NY da variável dependente Q
	! NX = número de pontos no vetor X e linhas da tabela Q
	! NY = número de pontos no vetor Y e colunas da tabela Q
	! X1 = valor procurado de x
	! Y1 = valor procurado de y
	! 
	! Descrição das variáveis de saída:
	!
	! FINT = valor interpolado de Q
	!
	! Descrição das variáveis locais:
	!
	! I,L,M = variáveis auxiliares
	! QA,QB,AUX = variáveis auxiliares para interpolação 
	!----------------------------------------------------------------------
	!
	IMPLICIT NONE
	! Declaração de variáveis:
	! Variáveis de entrada:
	integer,intent(in)::NX,NY
	real, intent(in):: X(NX),Y(NY),Q(NX,NY),X1,Y1
	! Variáveis de saída:
	real :: FINT2
	! Variáveis locais de cálculo:
	integer:: I,L,M
	real:: QA,QB,AUX
	!----------------------------------------------------------------------------
	
	! Verifica localização de X1:
	!DO I=2,NX-1
	DO I=2,NX
		IF(X1<=X(I))EXIT
	ENDDO 
	L=I
	! Verifica localização de Y1:
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