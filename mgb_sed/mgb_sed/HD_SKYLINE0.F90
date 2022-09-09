	SUBROUTINE SKYLINE0
	! Rotina auxiliar para preparar ponteiros das operacoes do algoritmo de eliminacao de gauss com 
	! armazenamento skyline.
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
	write(*,*) 'Conta numero de ponteiros'



	! Marca todos os elementos nao nulos:
	AA2=0.0
	do i=1,NUM
		K=ICOL(i,1)
		do j=2,K
			L=ICOL(i,j)
			AA2(L)=1.0
		enddo
	enddo


	! Primeiro conta o numero de elementos necessários em IOP(.)
	! Loop das colunas a serem zeradas abaixo da diagonal principal:
	NOP=0
	DO J=1,NN1
		! Itentifica posição em do elemento da diagonal principal da linha J em AA
		L1=IDIAG(J)
		M=J+1
		
		! Loop das linhas a serem operadas para zerar elementos A(JJ,J): 
		DO JJ=M,NUM
			! Numero de elementos entre coef a ser zerado A(JJ,J) e diagonal principal A(JJ,JJ) 
			JAUX=JJ-J
			! Verifica se elemento a ser zerado é não nulo:
			IF(IR(JJ)<JAUX)CYCLE 
			! Computa posição do elemento a ser zerado no vetor AA:
			JAUX=IDIAG(JJ)+IHIGH(JJ)-1+JAUX
			! Verifica se elemento a ser zerado já é igual a zero:
			IF(AA2(JAUX)==0.0)CYCLE      
			! Se A(JJ,J) é não nulo, linha JJ deve ser operada para zerar A(JJ,J)

			! Opera linha JJ:
			! Linha(JJ)=Linha(JJ)-COEF*Linha(J)

			NOP=NOP+1 !coluna j sendo zerada
			NOP=NOP+1 !linha i sendo operada
			NOP=NOP+1 ! posicao JAUX no vetor AA do elemento A(i,j)
			NOP=NOP+1 ! posicao L1 no vetor AA do elemento A(j,j)
			NOP=NOP+1 !numero de elementos operados na linha i zerando a coluna j.

!			COEF=AA2(JAUX)/AA2(L1)

!			BB2(JJ)=BB2(JJ)-BB2(J)*COEF


			! Opera elementos JM da linha JJ:
			DO JM=M,NUM
				JB=JM-J+1
				! Verifica se elemento JM da linha J é não nulo
				IF(JB>IHIGH(JM))CYCLE
				
				! Posição de A(JM,J) no vetor AA:
				J3=IDIAG(JM)+JB-1

				if (AA2(J3)==0.0) cycle		!!RP		

				IF(JM>=JJ)THEN
					! Elemento A(JJ,JM) a esquerda da diagonal principal:	
					J1=JM-JJ+1
					J2=IDIAG(JM)+J1-1 ! Posição de A(JJ,JM) no vetor AA
					! Computa coef. A(JJ,JM)
!					AA2(J2)=AA2(J2)-AA2(J3)*COEF
					AA2(J2)=1.0		
					NOP=NOP+1 ! Posicao J2 em AA do elemento operado A(i,k) na linha k
					NOP=NOP+1 ! Posicao J3 em AA do elemento localizado na linha j que opera A(i,k)
				ELSE
					! Elemento A(JJ,JM) a direita da diagonal principal:
					J1=JJ-JM
					J2=IDIAG(JJ)+IHIGH(JJ)-1+J1 ! Posição de A(JJ,JM) no vetor AA
					! Computa coef. A(JJ,JM)
!					AA2(J2)=AA2(J2)-AA2(J3)*COEF
					AA2(J2)=1.0
					NOP=NOP+1 ! Posicao J2 em AA do elemento operado A(i,k) na linha k
					NOP=NOP+1 ! Posicao J3 em AA do elemento localizado na linha j que opera A(i,k)
				ENDIF
			ENDDO
		ENDDO
	ENDDO



	! Depois de zerar elementos abaixo da diagonal principal, resolve sistema por substituição:

	! Posição do coeficiente A(NUM,NUM):
!	L=IDIAG(NUM)
	
	NOP2=NOP

!	XI2(NUM)=BB2(NUM)/AA2(L)

	! Calcula solução para XI(I) de I=NUM-1 até 1
	DO J=2,NUM
		
		K=NUM-J+1
		K1=K+1
		
		NOP2=NOP2+1 ! Numero de coeficientes computados no loop.:

		! Computa XI por substituição: 
		DO M=K1,NUM
			JAUX=M-K+1
			IF(JAUX>IHIGH(M))CYCLE ! Verifica se coef. é não nulo:
			! Posição do coef. no vetor AA
			J1=IDIAG(M)+JAUX-1

			if (AA2(J1)==0.0) cycle
!			BB2(K)=BB2(K)-AA2(J1)*XI2(M)

			NOP2=NOP2+1 ! Coluna M do coef. computado em XI(.)
			NOP2=NOP2+1 ! Posicao J1 em AA do coef A(.,.)

		ENDDO
		
		NOP2=NOP2+1 ! Posicao J1 em AA do coef A(i,i)
!		J1=IDIAG(K)
!		XI2(K)=BB2(K)/AA2(J1)
	ENDDO



!@ DCB set/2012	write(*,*) NOP,NOP2
!	read(*,*)

	allocate(IOP(NOP2))



	! Marca todos os elementos nao nulos:
	AA2=0.0
	do i=1,NUM
		K=ICOL(i,1)
		do j=2,K
			L=ICOL(i,j)
			AA2(L)=1.0
		enddo
	enddo

	! Armazena posicoes no vetor AA dos coeficientes que devem ser operados.
	! Loop das colunas a serem zeradas abaixo da diagonal principal:
	NOP=0
	DO J=1,NN1
		! Itentifica posição em do elemento da diagonal principal da linha J em AA
		L1=IDIAG(J)
		M=J+1
		
		! Loop das linhas a serem operadas para zerar elementos A(JJ,J): 
		DO JJ=M,NUM
			! Numero de elementos entre coef a ser zerado A(JJ,J) e diagonal principal A(JJ,JJ) 
			JAUX=JJ-J
			! Verifica se elemento a ser zerado é não nulo:
			IF(IR(JJ)<JAUX)CYCLE 
			! Computa posição do elemento a ser zerado no vetor AA:
			JAUX=IDIAG(JJ)+IHIGH(JJ)-1+JAUX
			! Verifica se elemento a ser zerado já é igual a zero:
			IF(AA2(JAUX)==0.0)CYCLE      
			! Se A(JJ,J) é não nulo, linha JJ deve ser operada para zerar A(JJ,J)

			! Opera linha JJ:
			! Linha(JJ)=Linha(JJ)-COEF*Linha(J)

			NOP=NOP+1 !coluna j sendo zerada
			IOP(NOP)=J
			NOP=NOP+1 !linha i sendo operada
			IOP(NOP)=JJ
			NOP=NOP+1 ! posicao JAUX no vetor AA do elemento A(i,j)
			IOP(NOP)=JAUX
			NOP=NOP+1 ! posicao L1 no vetor AA do elemento A(j,j)
			IOP(NOP)=L1
			NOP=NOP+1 !numero de elementos operados na linha i zerando a coluna j.
			N0=NOP			
			niOP=0
!			COEF=AA2(JAUX)/AA2(L1)

!			BB2(JJ)=BB2(JJ)-BB2(J)*COEF


			! Opera elementos JM da linha JJ:
			DO JM=M,NUM
				JB=JM-J+1
				! Verifica se elemento JM da linha J é não nulo
				IF(JB>IHIGH(JM))CYCLE
				
				! Posição de A(JM,J) no vetor AA:
				J3=IDIAG(JM)+JB-1

				if (AA2(J3)==0.0) cycle		!!RP		

				IF(JM>=JJ)THEN
					! Elemento A(JJ,JM) a esquerda da diagonal principal:	
					J1=JM-JJ+1
					J2=IDIAG(JM)+J1-1 ! Posição de A(JJ,JM) no vetor AA
					! Computa coef. A(JJ,JM)
!					AA2(J2)=AA2(J2)-AA2(J3)*COEF
					AA2(J2)=1.0		
					niOP=niOP+1
					NOP=NOP+1 ! Posicao J2 em AA do elemento operado A(i,k) na linha k
					IOP(NOP)=J2
					NOP=NOP+1 ! Posicao J3 em AA do elemento localizado na linha j que opera A(i,k)
					IOP(NOP)=J3
				ELSE
					! Elemento A(JJ,JM) a direita da diagonal principal:
					J1=JJ-JM
					J2=IDIAG(JJ)+IHIGH(JJ)-1+J1 ! Posição de A(JJ,JM) no vetor AA
					! Computa coef. A(JJ,JM)
!					AA2(J2)=AA2(J2)-AA2(J3)*COEF
					AA2(J2)=1.0		
					niOP=niOP+1
					NOP=NOP+1 ! Posicao J2 em AA do elemento operado A(i,k) na linha k
					IOP(NOP)=J2
					NOP=NOP+1 ! Posicao J3 em AA do elemento localizado na linha j que opera A(i,k)
					IOP(NOP)=J3
				ENDIF
			ENDDO
			IOP(N0)=niOP ! numero de elementos operados na linha i para zerar coluna j.			
		ENDDO
	ENDDO






	! Depois de zerar elementos abaixo da diagonal principal, resolve sistema por substituição:
	! Posição do coeficiente A(NUM,NUM):
	L=IDIAG(NUM)
	
	NOP2=NOP

!	XI2(NUM)=BB2(NUM)/AA2(L)

	! Calcula solução para XI(I) de I=NUM-1 até 1
	DO J=2,NUM
		
		K=NUM-J+1
		K1=K+1
		NOP2=NOP2+1 !numero de elementos operados na linha i zerando a coluna j.
		N0=NOP2			
		niOP=0
		! Computa XI por substituição: 
		DO M=K1,NUM
			JAUX=M-K+1
			IF(JAUX>IHIGH(M))CYCLE ! Verifica se coef. é não nulo:
			! Posição do coef. no vetor AA
			J1=IDIAG(M)+JAUX-1

			if (AA2(J1)==0.0) cycle
!			BB2(K)=BB2(K)-AA2(J1)*XI2(M)
			niOP=niOP+1
			NOP2=NOP2+1 ! Coluna M do coef. computado em XI(.)
			IOP(NOP2)=M
			NOP2=NOP2+1 ! Posicao J1 em AA do coef A(.,.)
			IOP(NOP2)=J1
		ENDDO

		IOP(N0)=niOP		

		J1=IDIAG(K)
		NOP2=NOP2+1 ! Posicao J1 em AA do coef A(i,i)
		IOP(NOP2)=J1
!		XI2(K)=BB2(K)/AA2(J1)
	ENDDO


!@ DCB set/2012	write(*,*) NOP,NOP2
!	read(*,*)



	RETURN
	END
