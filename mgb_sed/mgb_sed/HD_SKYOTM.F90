	SUBROUTINE SKYOTM
	! Rotina skyline otimizada. Utiliza ponteiros para operacoes e elimina decisões no algoritmo. Computa somente
	! operacoes algebricas pre definidas na rotina SKYLINE0.
	!
	!
	!
	! Resolve sistema de equações lineares pela eliminação de Gauss utilizando 
	! esquema especial de armazenamento, onde:
	! 
	! AA(.) = vetor que armazena coeficientes do sistema de esquações lineares
	!	   Armazena de forma sequencial:
	!			 - Elemento i da diagonal principal:
	!			 - Coeficientes da coluna acima do elemento de baixo para cima
	!            - Elementos na linha i a esquerda da diagonal principal ordenados da direita a esquerda
	! IHIGH(.) = numero de elementos acima da diagonal principal na coluna i, incluindo o elemento
	!				 A(i,i) da diagonal
	! IR(.) = número de elementos a esquerda da diagonal principal na linha i, excluindo o elemento	
	!				 A(i,i) da diagonal
	! NUM =  numero de equações (NBOUN + 2*NREAC + 3*NCONF)
	! BB(.) = vetor com coeficientes independentes do sistema de equações linear
	! XI(.) = vetor de variáveis de estado (profundidade e vazão em cada seção)
	! IDIAG(.) = posição do elemento i da diagonal no vetor AA
	! IOP(.) = armazena posicoes no vetor AA dos coeficientes que devem ser operados.
	!		- coluna j sendo zerada. J
	!       - linha i sendo operada. JJ
	!		- posicoes JAUX e L1 no vetor AA dos elementos A(i,j) e A(j,j)
	!		- numero de elementos operados na linha.
	!		- Para cada operacao na linha i zerando a coluna j:
	!				A(i,k)=A(i,k)-A(j,k)*COEF
	!		, armazena:
	!		- Posicao J2 em AA do elemento operado A(i,k) na linha k
	!		- Posicao J3 em AA do elemento localizado na linha j que opera A(i,k)
 	!	
	!-------------------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	!
	! NN1 = total de equações menos 1
	! J = contador (linhas da matriz A)
	! L1 = posição do elemento da diagonal principal no vetor AA
	! M = contador (linha a partir do qual se zeram coef. da coluna J -  no primeiro loop)
	! JJ = contador - linha a ser operada para zerar coeficiente A(JJ,J)
	! JAUX = auxiliar
	! COEF = fator utilizado na operação entre linhas (Ex: A(j,:)=A(j,:)-A(i,:)*COEF)
	! JM = auxiliar
	! JB = auxiliar
	! J1,J2,J3 = auxiliares
	! L = posição do ultimo elemento da diagonal 
	! K,K1 = auxiliares
	!--------------------------------------------------------------------------------
	! Declaração de variáveis:
	USE MAT_MOD

	IMPLICIT NONE

	! Variáveis locais de cálculo:

	integer:: i,NN1,J,L1,M,JJ,JAUX
	real(8):: COEF
	integer:: JM,JB,J1,J2,J3,L,K,K1
	integer:: cont,cont2,N0,niOP
	real(8):: aux
	!----------------------------------------------------------------------------

	NN1=NUM-1
	AA2=dble(AA)
	BB2=dble(BB)
	XI2=dble(XI)
	
	! Loop de operacoes algebricas para zerar elementos inferiores a diagonal principal de A:
	! Utiliza ponteiro IOP para indicar como fazer operacoes.
	i=0
	do
		i=i+1
		J=IOP(i)	!coluna j sendo zerada
		i=i+1
		JJ=IOP(i)   !linha i sendo operada
		i=i+1
		JAUX=IOP(i) ! posicao JAUX no vetor AA do elemento A(i,j)
		i=i+1
		L1=IOP(i) ! posicao L1 no vetor AA do elemento A(j,j)
		i=i+1
		niOP=IOP(i)

		! Calcula coef. multiplicador para operar linhas.
		! Linha(JJ)=Linha(JJ)-COEF*Linha(J)
		COEF=AA2(JAUX)/AA2(L1)
		! Opera termo independente da linha i (JJ)
		BB2(JJ)=BB2(JJ)-BB2(J)*COEF
		
		if (niOP==0) cycle

		! Opera termos dependentes (matriz A) da linha i (JJ):
		! Linha(JJ)=Linha(JJ)-COEF*Linha(J)
		! Loop dos termos a serem operados:
		do k=1,niOP
			i=i+1
			J2=IOP(i) ! Posicao J2 em AA do elemento operado A(i,k) na linha k
			i=i+1
			J3=IOP(i) ! Posicao J3 em AA do elemento localizado na linha j que opera A(i,k)
			! Opera coef. A(JJ,JM)
			AA2(J2)=AA2(J2)-AA2(J3)*COEF
		enddo
		if (i==NOP) exit ! Verifica se chega ao final das operacoes algebricas.	
	enddo

	! Depois de zerar elementos abaixo da diagonal principal, resolve sistema por substituição:

	! Posição do coeficiente A(NUM,NUM):
	L=IDIAG(NUM)
	XI2(NUM)=BB2(NUM)/AA2(L)
 
	! Calcula solução para XI(I) de I=NUM-1 até 1
	DO J=2,NUM
		K=NUM-J+1		
		i=i+1
		N0=IOP(i) !numero de elementos operados na linha i zerando a coluna j.
		if (N0>0) then
			do l=1,N0
				i=i+1
				M=IOP(i)
				i=i+1
				J1=IOP(i)
				BB2(K)=BB2(K)-AA2(J1)*XI2(M) ! Isola termo da diagonal principal
			enddo		
		endif
		i=i+1
		J1=IOP(i) ! Posicao J1 em AA do coef A(i,i)
		XI2(K)=BB2(K)/AA2(J1)
	ENDDO


	! Depois de zerar elementos abaixo da diagonal principal, resolve sistema por substituição:
	! Posição do coeficiente A(NUM,NUM):
!	L=IDIAG(NUM)
!	XI2(NUM)=BB2(NUM)/AA2(L)
	! Calcula solução para XI(I) de I=NUM-1 até 1
!	DO J=2,NUM
!		K=NUM-J+1
!		K1=K+1
!		! Computa XI por substituição: 
!		DO M=K1,NUM
!			JAUX=M-K+1
!			IF(JAUX>IHIGH(M))CYCLE ! Verifica se coef. é não nulo:

			! Posição do coef. no vetor AA
!			J1=IDIAG(M)+JAUX-1
!			if (AA2(J1)==0.0) cycle

!			BB2(K)=BB2(K)-AA2(J1)*XI2(M)
!		ENDDO
!		J1=IDIAG(K)
!		XI2(K)=BB2(K)/AA2(J1)
!	ENDDO



	XI=real(XI2)
	XI3=XI


	RETURN
	END
