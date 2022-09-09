	SUBROUTINE MATRIX2
	! Define arranjo da matriz de coeficientes e os vetores IDIAG, IR e IHIGH:
	! Rotina adaptada para esquema de armazenamento de dados utilizado no método dos gradientes conjugados.
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
	USE SPAMAT

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




	! Linhas das equações e colunas dos coeficientes conhecidos:
	! Reordenar linhas e colunas para obter diagonal principal somente com elementos não nulos


	! Inicialização de variáveis:
	DO I=1,J
		ICAUX(I,1)=1
	ENDDO
	! Numero total de equações NUM=NBOUN + 2*NREAC + 3*NCONF 
	NUM=J    
	IF (NCONF>0) INC=0


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


	! Copia matriz auxiliar ICAUX na variavel ICOL
	DO I=1,NUM
		L1=ICAUX(I,1)
		DO M=1,L1
			ICOL(I,M)=ICAUX(I,M)
		ENDDO
	ENDDO







	!******************************************************************************
	!DEFINE PONTEIROS DO ESQUEMA DE ARMAZENAMENTO DA MATRIZ ESPARSA DO MÉT. GRAD. BICONJUGADOS
	!******************************************************************************
	! Define vetor ija
	! Define posição de cada coeficiente no vetor AA e armazena em ICOL:

	allocate(IJA(NAOZEROMAX),IJA2(NAOZEROMAX))

	! Inicializa vars.:
	ija(1)=NUM+2
	k=NUM+1
	
	N2=0

	! Loop equacoes (linhas):
	do 	i=1,NUM
		L=ICAUX(i,1)
		! Loop colunas:
		do j=2,L
			if (ICAUX(i,j)==i) then
				! Elementos da diagonal principal:
				
				! Posicao do coeficiente no vetor AA:
				ICOL(i,j)=i
		
				N2=N2+1
			else
				! Elementos fora da diagonal principal:
				! Contador:
				k=k+1
				! Coluna do coeficiente:
				ija(k)=ICAUX(i,j)
				! Posicao do coeficiente no vetor AA:
				ICOL(i,j)=k
			endif
		enddo
		ija(i+1)=k+1
	enddo



	N1=k
	N2=NUM-N2+1
	N2=N1-N2

	! Posição das equações e respectivos coeficientes definida e armazenada em ICAUX 
	! Posicao dos coeficientes no vetor AA armazenada em ICOL


	


	! N1 é o numero de elementos nulos e nao nulos armazenados.
	! N2 é o numero de elementos nao nulos
	write(*,*) 'Numero de secoes,confluencias,cond.contorno e subtrechos:',nX,nConf,nboun,nREAC
	write(*,*) 'Numero de equacoes:',NUM
	write(*,*) 'Numero de elementos armazenados em AA:',N1
	write(*,*) 'Numero de elementos nao nulos:',N2
!	read(*,*)

	! Aloca vetor AA:
	dimSA=N1
	ALLOCATE (AA(dimSA)) 
	ALLOCATE (sa(dimSA)) 
	ija2=ija
	deallocate(ija)
	allocate(ija(dimSA))
	ija=ija2(1:dimSA)
	deallocate(ija2)
	allocate(xxi(NUM),bbb(NUM))
	
	! Inicializa vetor AA
	AA=0.0
	sa=0.0

	RETURN
	END