	SUBROUTINE SKYLIN
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

	integer:: NN1,J,L1,M,JJ,JAUX
	real(8):: COEF
	integer:: JM,JB,J1,J2,J3,L,K,K1
	integer:: cont,cont2
	real(8):: aux
	!----------------------------------------------------------------------------
cont=0
!cont2=0
	NN1=NUM-1

	AA2=dble(AA)
	BB2=dble(BB)
	XI2=dble(XI)

	
	! Loop das colunas a serem zeradas abaixo da diagonal principal:
!if (1==2) then
	DO J=1,NN1
		! Itentifica posição em do elemento da diagonal principal da linha J em AA
		L1=IDIAG(J)
		M=J+1
		
		! Pode paralelizar aqui: RP:
		!!$OMP PARALLEL NUM_THREADS(4)
		!!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(JAUX,COEF,JM,JB,J3,J1,J2)

		! Loop das linhas a serem operadas para zerar elementos A(JJ,J): 
		DO JJ=M,NUM
			! Numero de elementos entre coef a ser zerado A(JJ,J) e diagonal principal A(JJ,JJ) 
			JAUX=JJ-J
			! Verifica se elemento a ser zerado é não nulo:
			IF(IR(JJ)<JAUX)CYCLE 
			! Computa posição do elemento a ser zerado no vetor AA:
			JAUX=IDIAG(JJ)+IHIGH(JJ)-1+JAUX
			! Verifica se elemento a ser zerado já é igual a zero:
			IF(AA2(JAUX)==0)CYCLE      
			! Se A(JJ,J) é não nulo, linha JJ deve ser operada para zerar A(JJ,J)
!cont2=cont2+1			
			! Opera linha JJ:
			! Linha(JJ)=Linha(JJ)-COEF*Linha(J)
			if (AA2(L1)==0.0) then
!@ DCB ago/2012				write(7000,*) 'zero no denominador'
				write(*,*) 'zero no denominador'
			endif
			COEF=AA2(JAUX)/AA2(L1)
			BB2(JJ)=BB2(JJ)-BB2(J)*COEF
			! Pode paralelizar aqui tambem: RP:
			!!$OMP PARALLEL NUM_THREADS(4)
			!!$OMP DO PRIVATE(JB,J3,J1,J2)

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


!cycle
!cont=cont+1
					AA2(J2)=AA2(J2)-AA2(J3)*COEF
				ELSE

					! Elemento A(JJ,JM) a direita da diagonal principal:
					J1=JJ-JM
					J2=IDIAG(JJ)+IHIGH(JJ)-1+J1 ! Posição de A(JJ,JM) no vetor AA
					! Computa coef. A(JJ,JM)

!cycle
!cont=cont+1
					AA2(J2)=AA2(J2)-AA2(J3)*COEF
				ENDIF
				
!				if (isNaN(AA2(J2))) then
!					write(*,*) JM,JJ,AA2(J2),COEF,AA2(JAUX),AA2(L1),BB2(JJ),AA2(J3)
!					read(*,*)
!				endif


			ENDDO
			!!$OMP END DO
			!!$OMP END PARALLEL
		ENDDO
		!!$OMP END DO
		!!$OMP END PARALLEL
	ENDDO
!endif
!AA=1.0
!BB=1.0
!read(*,*)
!cont=0
	! Depois de zerar elementos abaixo da diagonal principal, resolve sistema por substituição:

	! Posição do coeficiente A(NUM,NUM):
	L=IDIAG(NUM)

	if (AA2(L)==0.0) then
			write(*,*) 'zero no denominador 2'
	endif

	XI2(NUM)=BB2(NUM)/AA2(L)
	! Calcula solução para XI(I) de I=NUM-1 até 1

	DO J=2,NUM
		K=NUM-J+1
		K1=K+1
!cycle		
		! Computa XI por substituição: 
		DO M=K1,NUM
			JAUX=M-K+1
			IF(JAUX>IHIGH(M))CYCLE ! Verifica se coef. é não nulo:
			! Posição do coef. no vetor AA

!cont=cont+1
			J1=IDIAG(M)+JAUX-1
			BB2(K)=BB2(K)-AA2(J1)*XI2(M)
		ENDDO
!cycle
		J1=IDIAG(K)

		if (AA2(J1)==0.0) then
				write(*,*) 'zero no denominador 2'
!@ DCB ago/2012				write(7000,*) 'zero no denominador 2'
		endif
	
		XI2(K)=BB2(K)/AA2(J1)
!		if (isNaN(XI2(K)).or.XI2(K)>100000000.0.or.XI2(K)<-100000000.0) then
!			write(*,*) K,AA2(J1),XI2(K)
!			read(*,*)
!		endif
	ENDDO

!XI=1.0
!write(*,*) cont
!read(*,*)

	XI=real(XI2)



	RETURN
	END
