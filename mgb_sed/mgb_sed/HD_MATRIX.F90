	SUBROUTINE MATRIX
	! Define arranjo da matriz de coeficientes e os vetores IDIAG, IR e IHIGH:
	!-------------------------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	!
	! INC(.) = identifica se existe condicao de contorno em secao de confluencia (NCONF x 1)
	! ISU = auxiliar. ISU=1 se cond. contorno de H e ISU=2 se cond. contorno de Q
	! I,J,K,M,N1,N2,N3,JK = auxiliares
	! MH,KX,K1,K2,K3,LM,L,JJ = auxiliares
	! L1 = auxiliar
	! JDIF = aux
	! JAU = aux
	! KT = aux
	!--------------------------------------------------------------------------------------
	! Declaração de variáveis:
	USE PAR1_MOD
	USE MAT_MOD
	USE TIME_MOD
	USE AUX_MOD

	IMPLICIT NONE

	! Variáveis locais de cálculo:

	integer,allocatable:: INC(:)
	integer:: ISU
	integer:: I,J,K,M,N1,N2,N3,JK
	integer:: MH,KX,K1,K2,K3,LM,L,JJ
	integer:: L1,JDIF,JAU,KT
	!-------------------------------------------------------------------------------------


	! Alloca variáveis: 
	NUM=NBOUN+2*NREAC+3*NCONF
	ALLOCATE(XI(2*NX))
	ALLOCATE(IHIGH(NUM),IR(NUM),BB(NUM),IDIAG(NUM),ICOL(NUM,5),ICAUX(NUM,5),JLIN(NUM))

	IF(NCONF>0)THEN
		ALLOCATE (INC(NCONF))
	ENDIF



	! Inicialmente define posição das equações e colunas dos coeficientes de cada uma delas.
	! As equações são inicialmente ordenadas da seguinte forma: eq. confluencias, eq. trechos e eq. confluencias

	! (A) Define posição dos coeficientes das eq. de condição de contorno:
	DO J=1,NBOUN
		! Define linha da equação:
		JLIN(J)=J
		! Verifica tipo de CC
		IF(bounFLAG(J)<=4) THEN
			! Série de vazão ou nível(profundidade) d'água:

			! Verifica tipo de série temporal: 
			ISU=1 ! Nivel
			IF (bounFLAG(J)<=3) ISU=2 ! Vazão

			! Calcula coluna do coef. da equação:
			K=(NBO(J)-1)*2+ISU
			ICOL(J,1)=2  ! Identificador do tipo de eq. (CC de Série de Vazão ou nível)
			ICOL(J,2)=K  ! Armazena coluna do coef. da equacao.
		ELSE
			! CC do tipo curva de descarga
			ICOL(J,1)=3  ! ! Identificador do tipo de Eq (CC de curva de descarga)
			! Calcula coluna do coef. da equação:
			K=(NBO(J)-1)*2
			! Armazena colunas dos coeficientes
			ICOL(J,2)=K+1
			ICOL(J,3)=K+2
		ENDIF
	ENDDO

	J=NBOUN


	! (B) Define posição dos coeficientes das eq. dos trechos (Continuidade e Dinamica)
	DO I=1,NREAC
		L1=NST(I,1) ! Seção de montante

		! Define linha e coluna dos coeficientes da eq. dinamica e continuidade
		CALL SUB(L1,NST(I,2),J)
		CALL SUB(L1,NST(I,2),J)
	ENDDO


	! (C) Define posição dos coeficientes das eq. das confluencias (Continuidade e Energia)

	IF(NCONF/=0)THEN
		DO I=1,NCONF
			! Linha da primeira eq. da confluencia
			J=J+1
			JLIN(J)=J
			ICOL(J,1)=4  ! Identificador do tipo de equação (continuidade confluencia)
			! Define coluna dos coeficientes da equação da continudade: 
			DO M=1,3
				K=(IABS(NCC(I,M))-1)*2+2
				ICOL(J,M+1)=K
			ENDDO
			
			! Define coluna dos coeficientes das equações de energia: 
			
			! Identifica seções da confluencia:
			N1=NCC(I,1)
			N2=NCC(I,2)
			N3=IABS(NCC(I,3))

			! Verifica tipo de confluencia:
			IF(NCC(I,3)>=0)THEN
				! Divergente:
				! Define linha e coluna dos coeficientes das eq. de energia
				CALL SUB(N3,N1,J)
				CALL SUB(N3,N2,J)
			ELSE
				! Convergente:
				! Define linha e coluna dos coeficientes das eq. de energia
				CALL SUB(N1,N3,J)
				CALL SUB(N2,N3,J)
			ENDIF
		ENDDO
	ENDIF

	! Linhas das equações e colunas dos coeficientes conhecidos:
	! Reordenar linhas e colunas para obter diagonal principal somente com elementos não nulos


	! Inicialização de variáveis:
	DO I=1,J
		ICAUX(I,1)=1
	ENDDO
	! Numero total de equações NUM=NBOUN + 2*NREAC + 3*NCONF 
	NUM=J    
	IF (NCONF>0) INC=0



!@ DCB ago/2012	do i=1,num
!@ DCB ago/2012		write(6787678,*) i,(icol(i,j),j=1,icol(i,1))
!@ DCB ago/2012	enddo


	!******************************************************************************
	!ORGANIZA AS LINHAS PARA TER ELEMENTOS DISTINTOS DE ZERO NA DIAGONAL PRINCIPAL
	!******************************************************************************

	!(A) Ordena linhas das equações de condição de contorno:
	DO J=1,NBOUN
		! Coluna do primeiro coeficiente da eq.
		K=ICOL(J,2)
		! Posiciona eq. de forma a manter coef. na diagonal principal	
		JLIN(J)=K   ! Linha = Coluna primeiro coeficiente
		
		! Copia posição (coluna) dos coeficientes em ICOL para variável ICAUX
		MH=ICOL(J,1)
		DO M=1,MH
			ICAUX(K,M)=ICOL(J,M)
		ENDDO

		! Seção da CC
		KX=NBO(J)

		! Verifica se CC se localiza em seção de confluencia:
		! Neste caso define linhas das equações da confluencia:
		IF(NCONF/=0)THEN
			DO I=1,NCONF
				! Posição do nivel HO de cada secao da confluencia do vetor de var. estado XI
				K1=(NCC(I,1)-1)*2+1
				K2=(NCC(I,2)-1)*2+1
				K3=(IABS(NCC(I,3))-1)*2+1
				! Em principio as equacoes das confluencias sao posicionadas de forma a 
				!manter dos coef. de HO na diagonal principal.
				
				! Equacoes de confluencias deverão ser posicionadas com coeficientes 
				! XI(K1),XI(K2) e XI(K3) na diagonal principal:

				! Cond. contorno na secao 1 da confluencia:
				IF(KX==NCC(I,1))THEN 
					! Verifica se confluencia é divergente:
					IF(NCC(I,3)>0)THEN
						! Se é divergente, QO da secao 3 fica na diagonal principal
						K3=K3+1
					ENDIF
					! Testa tipo de CC
					IF(bounFLAG(J)==4)THEN
						K1=K1+1 ! Se é tipo nivel dagua, QO da secao 1 fica na diagonal principal
					ENDIF
					! QO da secao 2 fica na diag. principal na eq. da continuidade
					! define linhas das equações da confluencia:
					CALL SUB2(K2+1,K1,K3,I)
					! Identifica que equacoes da confluencia I tem posicao definida
					INC(I)=1
				ENDIF

				! Cond. contorno na secao 2 da confluencia:
				IF(KX==NCC(I,2))THEN  !CONDIC DE CONTORNO NO OUTRO BRACO DA CONFLUENCIA
					! Verifica se confluencia é divergente:
					IF(NCC(I,3)>0)THEN
						! Se é divergente, QO da secao 3 fica na diagonal principal
						K3=K3+1
					ENDIF
					! Testa tipo de CC
					IF(bounFLAG(J)==4)THEN
						K2=K2+1 ! Se é tipo nivel dagua, QO da secao 2 fica na diagonal principal
					ENDIF
					
					! QO da secao 1 fica na diag. principal na eq. da continuidade
					
					! Define linhas das equações da confluencia:
				
					CALL SUB2(K1+1,K3,K2,I) !IDENTIFICA LINHA
					! Identifica que equacoes da confluencia I tem posicao definida
					INC(I)=1
				ENDIF

				! Cond. contorno na secao 3 da confluencia:
				!IF(KX.NE.IABS(NCC(I,3)))EXIT ! Acho que deve ser Cycle 
				IF(KX.NE.IABS(NCC(I,3)))CYCLE ! Acho que deve ser Cycle 
			
				! Testa tipo de CC
				IF(bounFLAG(J)==4)THEN
					K3=K3+1 ! Se é tipo nivel dagua, QO da secao 3 fica na diagonal principal
				ENDIF

				! QO da secao 2 fica na diag. principal na eq. da continuidade
				
				! Define linhas das equações da confluencia:
				CALL SUB2(K2+1,K1+1,K3,I)
				! Identifica que equacoes da confluencia I tem posicao definida
				INC(I)=1
			ENDDO
		ENDIF
	ENDDO

	J=NBOUN

	!(B) Ordena linhas dos trechos:
	DO I=1,NREAC
		! Posição da vazão QO nas secoes do trecho no vetor de var. estado XI
		K1=(NST(I,1)-1)*2+2 ! Montante
		K2=(NST(I,2)-1)*2+1       ! Jusante
		J=J+1

		! Em principio as equacoes dos trechos sao ordenadas de forma a manter 
		! QO(NST(I,1)) e HO(NST(I,2)) na diagonal principal:

		! Testa se linha ja possui equação de CC
		IF(ICAUX(K1,1)>1)THEN
			! HO da seção de montante na diagonal principal:
			K1=K1-1
		ENDIF
		! Armazena posição da linha da equação:
		JLIN(J)=K1
		! Copia posição (coluna) dos coeficientes em ICOL para variável ICAUX
		DO M=1,5
			ICAUX(K1,M)=ICOL(J,M)  !PEGA AS COLUNAS DOS COEFICIENTES
		ENDDO
		
		J=J+1
		
		! Testa se linha ja possui equação de CC
		IF(ICAUX(K2,1)>1)THEN
			! QO da seção de jusante na diagonal principal:
			K2=K2+1
		ENDIF
		! Armazena posição da linha da equação:
		JLIN(J)=K2
		! Copia posição (coluna) dos coeficientes em ICOL para variável ICAUX
!		write(*,*) I,K2,J,NUM,NBOUN,NREAC,NCONF,NBOUN+2*NREAC+3*NCONF,NST(I,1),NST(I,2)
		DO M=1,5
			ICAUX(K2,M)=ICOL(J,M)
		ENDDO
	ENDDO

	!J=NBOUN + 2*NREAC

	!(C) Ordena linhas das confluencias

	IF(NCONF>0)THEN
		DO I=1,NCONF
			! Verifica se posicao das eq. da confluencia I ja não esta definida:
			IF(INC(I)<=0)THEN
				! Posição do nivel HO de cada secao da confluencia do vetor de var. estado XI
				K1=(NCC(I,1)-1)*2+1
				K2=(NCC(I,2)-1)*2+1
				K3=(IABS(NCC(I,3))-1)*2+1
				! Em principio as equacoes das confluencias sao posicionadas de forma a 
				!manter dos coef. de HO na diagonal principal.

				! Verifica tipo de confluencia:
				IF(NCC(I,3)>=0)THEN
					! Divergente:
					CALL SUB2(K3+1,K1,K2,I)
					! HO secao 1 na diagonal principal na equação da energia
					! HO secao 2 na diagonal principal na equação da energia
					! QO secao 3 na diagonal principal na equacao da continuidade
				ELSE
					! Convergente:
					CALL SUB2(K1+1,K3,K2+1,I)
					! QO secao 1 na diagonal principal na equação da continuidade
					! QO secao 2 na diagonal principal na equação da energia
					! HO secao 3 na diagonal principal na equacao da energia
				ENDIF
			ENDIF
		ENDDO
	ENDIF


	! Posição das equações e respectivos coeficientes definida e armazenada em ICAUX 
!@ DCB ago/2012	do i=1,num
!@ DCB ago/2012		write(6787679,*) i,(icaux(i,j),j=1,icaux(i,1))
!@ DCB ago/2012	enddo


	! Copia matriz auxiliar ICAUX na variavel ICOL
	DO I=1,NUM
		L1=ICAUX(I,1)
		DO M=1,L1
			ICOL(I,M)=ICAUX(I,M)
		ENDDO
	ENDDO




	! Define ponteiros do esquema de armazenamento especial da matriz esparsa IHIGH,IDIAG,IR.

	!(A) Inicializa vetores:
	IDIAG(1)=1
	DO K=1,NUM
		IR(K)=0
		IHIGH(K)=1
	ENDDO

	!(B) Define vetores IDIAG,IR,IHIGH:
	N1=0
	N2=0
	! Loop das equacoes:
	DO K=1,NUM
		! K identifica a linha
		
		! Tipo de equação:
		M=ICOL(K,1)
		N2=N2+M-1
		! Posição do elemento da diagonal principal da linha K no vetor AA
		N1=N1+1
		IDIAG(K)=N1
		
		LM=0
		
		! Verifica coeficientes da linha K:
		DO L=2,M
			! Coluna do coeficiente:
			K1=ICOL(K,L)
			! Verifica posição em relação a diagonal principal:
			L1=K1-K
			IF(L1<0)THEN
				! Coef. a esquerda da diagonal principal
				! Verifica se é o coeficiente mais afastado (a esquerda) da diagonal principal:
				IF(LM>L1)THEN
					LM=L1
				ENDIF
				! Armazena numero de elementos a esquerda da diagonal na linha K:
				IR(K)=IABS(LM)
			ELSEIF(L1>0)THEN
				! Coef. a direita da diagonal principal:
				
				! Verifica se é o coeficiente mais afastado (acima) da diagonal principal:
				IF(IHIGH(K1)<L1+1)THEN
					! Armazena numero de elementos acima da diagonal na linha K1 (inclui diagonal):
					IHIGH(K1)=L1+1
				ENDIF
			ENDIF
		ENDDO
		! Posição do ultimo elemento armazenado em AA (numero de elementos armazenados de A(1:K,1:K))
		N1=N1+IHIGH(K)-1+IR(K)



	ENDDO


	! N1 é o numero de elementos nulos e nao nulos armazenados.
	! N2 é o numero de elementos nao nulos
	write(*,*) 'Numero de elementos armazenados em AA:',N1
	write(*,*) 'Numero de elementos nao nulos:',N2
	!read(*,*)

	! Aloca vetor AA:
	ALLOCATE (AA(N1)) 


	ALLOCATE (AA2(N1),BB2(NUM),XI2(NUM),XI3(NUM)) 
	ndimAA=N1
	!(C) Define posição de cada coeficiente no vetor AA e armazena em ICOL:
	! Loop das equações:
	DO K=1,NUM
		! Tipo de equação:
		L=ICOL(K,1)
		! Loop dos coef.:
		DO JJ=2,L
			! Distancia entre coef. e diagonal principal:
			JDIF=ICOL(K,JJ)-K
			! Verifica posição em relação a diagonal principal: 
			IF(JDIF<0)THEN
				! Coef. a esquerda da diagonal principal
				! Posicao do coef. em AA:
				JAU=IDIAG(K)+IHIGH(K)-1+IABS(JDIF)
				! Obs.: Posicao do elemento da diagonal + numero de elementos acima da diagonal + distancia da diagonal
			ELSE
				! Coef. acima ou na diagonal principal
				! Calcula coluna do coef.: 
				KT=JDIF+K
				! Posicao do coef. em AA:
				JAU =IDIAG(KT)+JDIF
				! Obs.: Posicao do elemento da diagonal + distancia da diagonal
			ENDIF
			! Armazena posicao do coef. em AA na variavel ICOL:
			ICOL(K,JJ)=JAU
		ENDDO
	ENDDO
	! A partir de agora ICOL armazena posicao dos coeficientes no vetor AA 

    IF(NCONF>0)THEN
		DEALLOCATE (INC)
	ENDIF


	RETURN
	END