	MODULE MAT_MOD
	! Declaração de variáveis globais relativas as matrizes de coeficientes do modelo hidrodinâmico e 
	! armazenamento especial da matriz esparsa.
	! 
	!
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descrição das variáveis:  ***********************************
	!
	!
	! AA(.) = vetor que armazena coeficientes do sistema de esquações lineares (N1 x 1)
	!			N1 é o numero de elementos não nulos mais os elementos nulos eventualmente armazenados
	!         Armazena de forma sequencial:
	!			 - Elemento i da diagonal principal:
	!			 - Coeficientes da coluna acima do elemento de baixo para cima
	!            - Elementos na linha i a esquerda da diagonal principal ordenados da direita a esquerda
	! IHIGH(.) = numero de elementos acima da diagonal principal na coluna i, incluindo o elemento
	!				 A(i,i) da diagonal (NUM x 1)
	! IR(.) = número de elementos a esquerda da diagonal principal na linha i, excluindo o elemento	
	!				 A(i,i) da diagonal (NUM x 1)
	! NUM =  numero de equações (NBOUN + 2*NREAC + 3*NCONF)
	! BB(.) = vetor com coeficientes independentes do sistema de equações linear (NUM x 1)
	! XI(.) = vetor de variáveis de estado (profundidade e vazão em cada seção), (NX*2 x 1) armazenados conforme a seguir:
	!         XI=[HO(1) QO(1) HO(2) QO(2) ... HO(NX) QO(NX)
	! IDIAG(.) = posição do elemento i da diagonal no vetor AA (NUM X 1)
	! ICOL(.,.) = identifica localização dos elementos não nulos da matriz de coeficientes dependentes e depois no vetor AA (NUM x 5).
	!			ICOL(I,1) identifica tipo de linha: 
	!					ICOL(I,1)=2 : eq. de condição de contorno de Q ou H (1 coef. não nulo)
	!					ICOL(I,1)=3 eq. de condição de contorno de curva-chave (2 coef. não nulo)
	!					ICOL(I,1)=4 eq. de conservação de massa de confluencia (3 coef. não nulo)
	!					ICOL(I,1)=5 eq. dinamica ou continuidade em trecho ou energia em conflu. (4 coef. não nulo)
	!			ICOL(I,J), J=2,5: inicialmente identifica coluna do coeficiente na matriz A e depois identifica 
	!								posição do coef. no vetor AA 
	! ICAUX(.,.) = idem ao ICOL, auxiliar para reordenamento das linhas da matriz de coeficientes dependentes.
	! JLIN(.) = identifica localização de determinada linha (equação) da matriz A original na matriz 
	!					A reordenada (NUM x 1)
	!------------------------------------------------------------------------------------------------

	IMPLICIT NONE
	SAVE

	REAL,ALLOCATABLE:: AA(:)
	INTEGER,ALLOCATABLE:: IHIGH(:),IR(:) 
	INTEGER NUM
	REAL,ALLOCATABLE:: BB(:),XI(:)
	INTEGER,ALLOCATABLE:: IDIAG(:),ICOL(:,:),ICAUX(:,:),JLIN(:)


	REAL,ALLOCATABLE:: AA2(:)
	REAL,ALLOCATABLE:: BB2(:),XI2(:)

	integer:: ndimAA




	integer,allocatable:: IOP(:)
	integer:: NOP,NOP2
	real,allocatable:: XI3(:)

	END MODULE