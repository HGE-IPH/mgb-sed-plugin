	SUBROUTINE INPUT2
	! Subrotina para leitura e organização de dados de entrada:
	!-------------------------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	!
	! REF(.) = fator de correção das condicoes de contorno de nivel (NBOUN x 1)
	!			! Ex.: correcao de zero da régua

	! QAUX(.) = auxiliar. Acumula DT dos dados (NThd x 1)
	! I,L,J,K,M,IK = auxiliares
	! FINT = auxiliar rotina de interpolacao
	! FLAG = auxiliar
	! NBOUNTS1 = numero de condicoes de contorno tipo serie temporal
	! NBOUNTS2 = numero de condicoes de contorno tipo serie temporal lidas em arquivo externo
	! NBOUNRC = numero de condicoes de contorno tipo curva de descarga
	! COUNT,COUNT2 = auxiliar
	! H1,Q1 = nivel e vazao auxuliar
	! TRASH = aux.
	! BBaux(.),HHaux(.) = largura e profundidade máxima da secao transversal auxiliares de leitura (nX x 1)
	! nTCC = numero de intervalo de tempo dos dados de condicao de contorno.
	! dTCC = intervalo de tempo dos dados de cond. contorno.
	!--------------------------------------------------------------------------------------
	! Declaração de variáveis:
	USE PAR1_MOD
	USE CHEIA_MOD
	USE MAT_MOD
	USE TIME_MOD
	USE BARR_MOD
	USE PP_MOD
	USE AUX_MOD
	USE VARS_MAIN

	IMPLICIT NONE

	! Variáveis locais de cálculo:

	real,allocatable:: REF(:),QAUX(:)
	integer:: I,L,J,K,M,IK
	real:: FINT
	integer:: FLAG,NBOUNTS1,NBOUNTS2,NBOUNRC,COUNT,COUNT2
	real:: H1,Q1
	character(5):: TRASH
	real:: aux(10)
	real,allocatable:: BBaux(:),HHaux(:)
	integer:: nTCC
	real:: dTCC
	real,allocatable:: FAFaux(:)
	integer:: i1,i2,i3
	real:: aux2
	real:: dV
	!-------------------------------------------------------------------------------------




	! O que fazer para simplificar codigo:



	
	
	! Colocar na rotina coef1 teste se dx<valor limite e utilizar equacao simplificada: 	
	! Simplificar rotina coef2. Tirar parte de energia cinetica e perda de carga localizada
	! VARRER CODIGO VERIFICANDO INCONSISTENCIAS.

	! VERIFICAR UNIDADES DE COMPRIMENTOS GERADOS NO PREPRO E UTILIZADOS NO HIDRODINAMICO


	! VERIFICAR UNIDADES DA CONTRIBUIÇÂO LATERAL

	
	write(*,*) ' Lê dados do arquivo MainInfo.txt'
	! Lê dados do arquivo MainInfo.txt:
	OPEN(301,FILE='.\input\HD_MainInfo.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	! Numero de trechos:
	READ(301,*) nTr
	READ(301,*) 
	! Numero de secoes transversais:
	READ(301,*) NX
	READ(301,*) 
	! Numero de subtrechos:
	READ(301,*) NREAC
	READ(301,*) 
	! Numero de confluencias:
	READ(301,*) NCONF
	READ(301,*) 
	! Numero de condicoes de contorno:
	READ(301,*) NBOUN
	READ(301,*) 
	! Flag que indica se dados estao em termos de profundidade ou nivel dagua:
	READ(301,*) LD
	READ(301,*) 
	! Intervalo de tempo dos dados de condicao de contorno:
	READ(301,*) dTCC
	READ(301,*)
	! Intervalo de tempo de calculo:
	READ(301,*) dThd
	READ(301,*)
	! Flag calcula condicoes iniciais:
	READ(301,*) INIC
	CLOSE(301)


	!Aloca variáveis:
	ALLOCATE (NST(NREAC,2),NCC(NCONF,3),NBO(NBOUN),ZO(NX),DXT(NREAC),DXC(NCONF*2),ALFA(NCONF*2),NP(NX))
	ALLOCATE (BETA(NX))
	ALLOCATE (THETA(NREAC))



	ALLOCATE (NPF(NX),N(NX))
	!Aloca variáveis:
	ALLOCATE (HO(NX),QO(NX))
	ALLOCATE (HOprev(NX),QOprev(NX))

	ALLOCATE (HO1(NX),QO1(NX))
	ALLOCATE (HminFLAG(NX),Qmin(NX))
	ALLOCATE (HminFLAGprev(NX),Qminprev(NX))
	HminFLAG=0
	Qmin=-1.0

	ALLOCATE (QL2(NREAC,2))

	! Aloca variáveis:
	ALLOCATE (T(NX),A(NX),R(NX),CK(NX),CKY(NX),SFF(NX))

	ALLOCATE (IDhd(nTr),TrCELL(nTr),TrNST(nTr,2),nXSec(nTr),TrL(nTr),TrMan(nTr),TrXSec(nX))


!	ALLOCATE(hdFLAG(nC)) ! RP tirar depois ! So teste

	! Verificar

	! Numero de intervalos de tempo:
	! Hidrodinamico dentro de 1 intervalo de tempo do MGB
	nThd=DTP/DThd+1
	! Total do hidrodinamico:	
	nThd2=NT*DTP/DThd+1
	! Total dados condicao de contorno:
	nTCC=NT*DTP/DTCC+1
	! Verificar isto:

	dtHD0=dtHD

	ALLOCATE (HQB(nThd,NBOUN))
	ALLOCATE(REF(NBOUN))
	ALLOCATE(QAUX(nThd2))
	ALLOCATE (bounID(nBoun),bounFLAG(nBoun))
	ALLOCATE (reacFLAG(nReac))

	! Inicializa:
	DXT=0.0
	ICONF=1
	NP=0

	reacFLAG=1

	! Lê dados do arquivo HD_Trecho.txt:
	write(*,*) ' Lê dados do arquivo HD_Trecho.txt'
	OPEN(301,FILE='.\input\HD_Trecho.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	DO iTr=1,nTr
		! Trecho, subbacia, bacia, secao montante, secao jusante, numero de secoes, comprimento trecho, coef. rug. manning:
		READ(301,*) IDhd(iTr),TrCELL(iTr),I,TrNST(iTr,1),TrNST(iTr,2),nXSec(iTr),TrL(iTr),TrMan(iTr)
		! Tirar depois, somente enquanto cell.hig esta errado
!		hdFLAG(TrCELL(iTr))=1 ! Teste
	ENDDO	
	CLOSE(301)

	! Lê dados do arquivo HD_SubTrechos.txt:
	write(*,*) ' Lê dados do arquivo HD_SubTrechos.txt'
	OPEN(301,FILE='.\input\HD_SubTrecho.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	DO I=1,NREAC 
		! Subtrecho, secoes de montante e jusante e dx:
		READ(301,*) J,NST(I,1),NST(I,2),DXT(I)
		! Testa se subtrecho é curto:
		IF (DXT(I)<=100.0) reacFLAG(I)=2 
		

! Eliminar trechos muito curtos: ! RP
		IF (DXT(I)>100.0.and.DXT(I)<5000.0) DXT(I)=5000.0 !*********

	ENDDO
	CLOSE(301)



	! Lê dados do arquivo HD_Conflu.txt:
	write(*,*) ' Lê dados do arquivo HD_Conflu.txt'
	OPEN(301,FILE='.\input\HD_Conflu.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	IF (NCONF>0) THEN
		DO I=1,NCONF
			! Conflu, secoes 1, 2 e 3, trechos,sentido do escoamento:
			READ(301,*) J, NCC(I,1),NCC(I,2),NCC(I,3), (aux(k),k=1,3), FLAG
			IF (FLAG>=0) NCC(I,3)=-NCC(I,3) ! Muda sinal da secao 3 para ficar de acordo com dados de entrada do IPH-IV
		ENDDO
	ENDIF
	CLOSE(301)
	! Define distancia entre secoes da confluencia igual a zero:
	DXC=0.0
	! Define coef. ALFA igual a 1:
	ALFA=1.0

	! Define coef. BETA igual a 1:
	BETA=1.0	

	! Lê dados do arquivo HD_XSecInfo.txt:
	write(*,*) ' Lê dados do arquivo HD_XSecInfo.txt'
	ALLOCATE(BBaux(NX),HHaux(NX))
	OPEN(301,FILE='.\input\HD_XSecInfo.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	DO iX=1,nX
		! Secao, trecho de rio, distancia a inicio do trecho, largura, prof. maxima, nivel do fundo,
		! numero de pontos das tabelas cota x (area, largura e raio hidráulico), cota x Manning, cota x largura planicie:
		READ(301,*) I,TrXSec(iX),aux(1),BBaux(iX),HHaux(iX),ZO(iX),NP(iX),N(iX),NPF(iX)
	ENDDO
	CLOSE(301)

	! Aloca variáveis:
	COUNT=MAXVAL(NP)
	IF (COUNT/=0) THEN
		ALLOCATE(HA(COUNT,NX),AR(COUNT,NX),RR(COUNT,NX),BA(COUNT,NX))
	ELSE
		ALLOCATE(HA(1,NX),AR(1,NX),RR(1,NX),BA(1,NX))
	ENDIF

	COUNT=MAXVAL(N)
	IF (COUNT/=0) THEN
		ALLOCATE(F(COUNT,NX))
		ALLOCATE(HUX(COUNT,NX))
	ELSE
		ALLOCATE(F(1,NX))
	ENDIF

	COUNT=MAXVAL(NPF)
	! OBS.: Verificar se deixaremos centralizado na secao ou no trecho:????????????????????????????????????????
	IF (COUNT>0)  THEN
		ALLOCATE(HF(COUNT,NX),FAF(COUNT,NX)) ! (HF(COUNT,nTr),FAF(COUNT,nTr))
		ALLOCATE(AFI(NX))
!		ALLOCATE(Afl(COUNT,nTr),Zfl(COUNT,nTr),nPfl(nTr))
		ALLOCATE(Afl(COUNT,nX),Zfl(COUNT,nX),nPfl(nX))
		ALLOCATE(Vfl(COUNT,nX))
	ENDIF




	! Lê dados do arquivo HD_XSecData.txt:
	write(*,*) ' Lê dados do arquivo HD_XSecData.txt'
	OPEN(301,FILE='.\input\HD_XSecData.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	READ(301,*)

	! Lê tabela prof.ou nivel x (area, largura e raio hidráulico):
	DO iX=1,NX
		IF (NP(iX)==0) THEN
			! Usa secao retangular:
			HA(1,iX)=HHaux(iX)
			BA(1,iX)=BBaux(iX)
		ELSE
			! Usa batimetria detalhada:
			DO J=1,NP(iX)
				READ(301,*) FLAG,HA(J,iX),AR(J,iX),RR(J,iX),BA(J,iX)
			ENDDO
		ENDIF
	ENDDO

	! Lê tabela prof.ou nivel x Manning:
	READ(301,*)
	READ(301,*)
	DO iX=1,NX
		IF (N(iX)==0) THEN
			! Usa Manning dos trechos:
			iTr=TrXsec(iX)
			F(1,iX)=TrMan(iTr)
		ELSE
			! Usa curva nivel coef de Manning:
			DO J=1,N(I)
				READ(301,*) FLAG, HUX(J,I),F(J,I)
			ENDDO
		ENDIF
	ENDDO
	CLOSE(301)


	Afl=0.0
	Zfl=0.0
	! Lê tabela prof.ou nivel x area alagada:
	write(*,*) ' Lê dados do arquivo HD_FPlainData.txt'
	OPEN(301,FILE='.\input\HD_FPlainData.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	READ(301,*)
	nPfl=0
	DO WHILE(.NOT.EOF(301))
		READ(301,*) iX,aux(1),Zfl(nPfl(iX)+1,iX),Afl(nPfl(iX)+1,iX)
		nPfl(iX)=nPfl(iX)+1
	ENDDO
	CLOSE(301)

	! Calcula largura da planicie de cada seção
	
	! Calculo considerando area alagada por minibacia:
!	DO iX=1,nX
!		iTr=TrXSec(iX)
!		! Niveis iguais a curva nivel x area alagada
!		HF(:,iX)=Zfl(:,iTr)
!			
!		! Largura planicie igual a area alagada (km) /comprimento rio (km) - largura rio (m)
!		IF (NP(iX)==0) THEN
!			! Usa secao retangular:
!			FAF(:,iX)=1000.0*Afl(:,iTr)/TrL(iTr)-BA(1,iX)
!		ELSE
!			! Usa batimetria detalhada:
!			FAF(:,iX)=1000.0*Afl(:,iTr)/TrL(iTr)-BA(NP(iX),iX)
!		ENDIF
!		FAF(:,iX)=max(FAF(:,iX),0.0)
!	ENDDO

	FAF=0.0
	HF=0.0

	!*************
	! Para mudar coef. de rugosidade de manning:
!	F=0.031
	
	
	
	
	! Para corrigir largura no rio Solimoes:
!	do ix=1,nx
!		if (BA(1,ix)>1000.0) BA(1,ix)=3000.0
!	enddo

	! Refinar depois forma de calcular largura: !RP
	DO i=1,nReac
		iX=NST(i,2)
		iTr=TrXSec(iX)

		if (DXT(i)==0.0) cycle ! Pula trecho curto

		! Niveis iguais a curva nivel x area alagada
		HF(:,iX)=Zfl(:,iX)	
		
		!*************
		! Para mudar referencia de nivel da planicie:
		HF(:,iX)=HF(:,iX) !5.0
		!**************


		! Largura planicie igual a area alagada (km) /comprimento rio (km) - largura rio (m)
		IF (NP(iX)==0) THEN
			! Usa secao retangular:
			FAF(:,iX)=1000000.0*Afl(:,iX)/DXT(i)-BA(1,iX)
		ELSE
			! Usa batimetria detalhada:
			FAF(:,iX)=1000000.0*Afl(:,iX)/DXT(i)-BA(NP(iX),iX)
		ENDIF

		FAF(:,iX)=max(FAF(:,iX),0.0)
		
!		do j=1,NPF(iX)
!			if (HF(j,iX)<ZO(iX)) then
!				FAF(j,iX)=0.0
!			elseif (HF(j,iX)<ZO(iX)+HA(1,iX)) then
!				FAF(j,iX)=FAF(j,iX)*(HF(j,iX)-ZO(iX))/HA(1,iX)
!			endif
!		enddo


	ENDDO

	! Considera planicie nas secoes de inicio de trecho iguais as das secoes seguintes:
	do iTr=1,nTr
		iX=TrNST(iTr,1)
		HF(:,iX)=HF(:,iX+1)
		FAF(:,iX)=FAF(:,iX+1)
		NPF(iX)=NPF(iX+1)
	enddo


!*************************************************************
!	HF=HF+10.0 ! Para Mudar referencia da planicie de inundacao

!	FAF=0.0 ! Para desativar planicie inundacao:


    ! Testes final do amazonas a partir da confluencia com o madeira
    
 !   do iX=1,nX
 !       iTr=TrXSec(iX)
 !       iC=TrCELL(iTr)
 !       iB=IBAC(iC)
        
 !       if (iC<6823) cycle
  !      if (iB<179) cycle
       
  !      BA(1,iX)=0.8*BA(1,iX)
        
      !  if (iB<179) cycle
        !write(*,*) iB,iTr,iC,iX
        !read(*,*)
     !   FAF(:,iX)=2.0*FAF(:,iX)
        
  !  enddo

	
		

	count=MAXVAL(NPF)
	allocate(FAFaux(count))







!@ DCB ago/2012	write(8002,*) 'iX  ZO(iX)  BA HA  HF  HF  FAF'
	do iX=1,nX
		do i=1,NPF(iX)
			if (HF(i,iX)>=ZO(iX)) then
!@ DCB ago/2012				write(8002,'(I7,6F10.2)') iX,ZO(iX),BA(1,iX),HA(1,iX),HF(1,iX),HF(i,iX),FAF(i,iX)
				exit
			endif
		enddo
	enddo
!@ DCB ago/2012	close(8002)
	
!@ DCB ago/2012	write(9001,*) 'iX  Z  FAF'
	do iX=1,nX
		do i=1,NPF(iX)
!@ DCB ago/2012				write(9001,*) iX,HF(i,iX),FAF(i,iX)
		enddo
	enddo
!@ DCB ago/2012	close(9001)
		

!************************************************************************************
	! Lê dados do arquivo HD_CCInfo.txt:
	write(*,*) ' Lê dados do arquivo HD_CCInfo.txt'
	NBOUNTS1=0
	NBOUNTS2=0
	NBOUNRC=0
	REF=0.0
	OPEN(301,FILE='.\input\HD_CCInfo.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	DO 
		! Condição de contorno, seção transversal, tipo, fator de correção da régua:
		READ(301,*) I,NBO(I),aux(1),bounFLAG(I),REF(I)
		! Compatibiliza com variáveis de entrada do IPH-IV:
		SELECTCASE (bounFLAG(I))
		CASE (1,2) !Condicao de contorno de vazões (resultado do MGB de série de vazao):
			NBOUNTS1=NBOUNTS1+1
		CASE (3) !Condicao de contorno de vazões lido em arquivo externo:
			NBOUNTS1=NBOUNTS1+1
			NBOUNTS2=NBOUNTS2+1
		CASE (4) !Condicao de contorno de níveis d'água lido em arquivo externo:
			NBOUNTS1=NBOUNTS1+1
			NBOUNTS2=NBOUNTS2+1
		CASE (5) !Curva chave:	
			NBOUNRC=NBOUNRC+1
		CASE (6) ! Equação de Manning:
		ENDSELECT
		IF (I==NBOUN) EXIT
	ENDDO
	CLOSE(301)




	! Lê dados do arquivo HD_CC_TSData.txt:
	write(*,*) ' Lê dados do arquivo HD_CC_TSData.txt'
	IF (NBOUNTS2/=0) THEN
		ALLOCATE (bounTS(nThd2,NBOUNTS2))

		OPEN(301,FILE='.\input\HD_CC_TSData.txt',STATUS='OLD',ACTION='READ')


		bounTS=0.0
		DO I=1,NTCC
			! Le dados de condicao de contorno em serie temporal:

			READ(301,*) (bounTS(I,J),J=1,NBOUNTS2) 
			
		ENDDO
		CLOSE(301)
	ENDIF


	! Define vetor bounID:
	bounID=0
	J=0
	K=0
	DO I=1,NBOUN
		IF(bounFLAG(I)==4.OR.bounFLAG(I)==5) THEN
			J=J+1
			bounID(I)=J			
		ELSEIF (bounFLAG(I)==5) THEN
			K=K+1
			bounID(I)=K			
		ENDIF
	ENDDO



	! Talvez tirar essa correcao ou pelo menos melhorar codigo:
	!************** OBS************************************8

	! Correção dos dados:
!	REF=0.0
	DO I=1,NBOUN
		IF(bounFLAG(I)==4) THEN
			! Fator de correção:
			J=bounID(I)
			bounTS(:,J)=bounTS(:,J)+REF(I)
!			bounTS(:,1)=10.00 !Teste
			! Corrige niveis d'água:
			IF (LD/=0) THEN
				bounTS(:,J)=bounTS(:,J)-ZO(NBO(I))
			ENDIF
		ENDIF
	ENDDO



					
	IF (NBOUNRC/=0) THEN	
		! Lê dados do arquivo HD_CC_RCData.txt:
		write(*,*) ' Lê dados do arquivo HD_CC_RCData.txt'
		OPEN(301,FILE='.\input\HD_CC_RCData.txt',STATUS='OLD',ACTION='READ')
		READ(301,*)
		
		! Alloca vars.:
		ALLOCATE (NPX(NBOUNRC))
		NPX=0
		J=0
		K=0
		COUNT=0
		! Conta o numero máximo de pontos de curva chave:
		DO WHILE (.NOT.EOF(301))
			! Condição de contorno, número da curva chave
			READ(301,*) I

			IF (I==J) THEN
				NPX(K)=NPX(K)+1
			ELSE
				K=K+1
				NPX(K)=NPX(K)+1
				J=I
			ENDIF		
			COUNT=MAX(NPX(K),COUNT)
		ENDDO

		! COUNT é o numero maximo de pontos na tabela de curva chave
		ALLOCATE(QT(COUNT,NBOUNRC),HT(COUNT,NBOUNRC))
		! Comeca no inicio do arquivo.
		REWIND(301) 
		! Lê dados do arquivo HD_CCRCData.txt:
		READ(301,*)
		DO I=1,NBOUNRC
			DO J=1,NPX(I)
				! Condição de contorno, nível, vazão
				READ(301,*) K,HT(J,I),QT(J,I)
			ENDDO
		ENDDO
		CLOSE(301)
	ENDIF





!***********************************************************************************************************
	! Lê dados das barragens, arquivo HD_StrHdData:
	write(*,*) ' Lê dados do arquivo HD_StrHdData.txt'
	OPEN(301,FILE='.\input\HD_StrHdData.txt',STATUS='OLD',ACTION='READ')
	READ(301,*)
	READ(301,*) LTOT
	IF (LTOT/=0) THEN
		READ(301,*)
		READ(301,*) NCMAX  !Numero máximo de curvas por barragem:
		READ(301,*)
		READ(301,*) NPMMAX,NPJMAX  !Máximo numero de pontos das curvas H1xH2xQ referentas 
								   ! a cotas de montante e jusante:
		! Aloca variáveis:
		ALLOCATE(NIT(LTOT),NPB(LTOT,NCMAX),NPB2(LTOT,NCMAX),JIT(LTOT,NCMAX),H2BR(NPJMAX,LTOT,NCMAX),HBR(NPMMAX,LTOT,NCMAX),QBR(NPMMAX,NPJMAX,LTOT,NCMAX))
		ALLOCATE(barrReac(LTOT))

		HBR=0.0
		H2BR=0.0
		QBR=0.0
		READ(301,*)
		! Lê dados das barragens:
		DO I=1,LTOT
			READ(301,*)
			READ(301,*)
			READ(301,*)
			READ(301,*) J !Subtrecho
			K=NST(J,1)
			barrReac(I)=J ! Armazena trecho da barragem
			reacFLAG(J)=3 	!NST(J,1)=-NST(J,1) ! Flag que indica que tem barragem no subtrecho
			READ(301,*)
			READ(301,*) NIT(I) ! Numero de curvas
			READ(301,*)
			READ(301,*) FLAG   ! Flag se existe efeito de jusante
			DO J=1,NIT(I)
				READ(301,*)
				READ(301,*)
				IF (FLAG==1) THEN
					READ(301,*) NPB(I,J),NPB2(I,J)
					READ(301,*)
					READ(301,*) JIT(I,J)
					READ(301,*)
					READ(301,*) TRASH, (H2BR(L,I,J), L=1,NPB2(I,J))
					DO L=1,NPB(I,J)
						READ(301,*) HBR(L,I,J), (QBR(L,M,I,J),M=1,NPB2(I,J))
					ENDDO
					NPB(I,J)=-NPB(I,J)
				ELSE
					READ(301,*) NPB(I,J)
					READ(301,*)
					READ(301,*) JIT(I,J)
					READ(301,*)
					READ(301,*)
					DO L=1,NPB(I,J)
						READ(301,*) HBR(L,I,J),QBR(L,1,I,J)
					ENDDO	
				ENDIF
			ENDDO
		ENDDO	    
	ENDIF
	CLOSE(301)


!************************************************************************************************************************************



	! PENSAR EM TROCAR PARA HOTSTART:

	! Condições Iniciais:
	IF (INIC==0) THEN
		! Lê dados do arquivo CondInicial.txt
		OPEN(301,FILE='.\input\HD_CondInic.txt',STATUS='OLD',ACTION='READ')
		READ(301,*)
		DO I=1,NX
			! Condicao inicial de nivel ou profundidade e vazão em cada secao:
			READ(301,*) J,HO(I),QO(I) 
		ENDDO
		CLOSE(301)
		! Passa dados de níveis d'água para profundidade
		IF(LD/=0)THEN
!			HO=HO-ZO	! RP Verificar depois.  
		ENDIF
	ENDIF


	! Se DT dos dados é diferente que DT de calculo, 
	! interpola séries de condição de contorno e contribuição lateral:
	IF(dThd/=dTCC.and.NBOUNTS2/=0)THEN
		QAUX(1)=0.0
		! Acumula DT dos dados:
		DO IK=2,NTCC
			QAUX(IK)=QAUX(IK-1)+DTCC
		ENDDO

		!Interpola condições de contorno:
		DO J=1,NBOUNTS2
			CALL INTEP(QAUX,dThd,nTCC,bounTS(:,J),nThd2)
		ENDDO


		! Colocar depois para dados observados tambem REVER

		! Novos numero e tamanho de intervalos de tempo (igual ao de calculo):
		
	ENDIF
	

	ALLOCATE(TP2(nThd2))
	do i=1,nThd2
		TP2(i)=(i-1)*dthd
	enddo



! Correção do Codigo: RP fev/2011
! Motivo, intervalo de tempo de QCONTORM é sempre 1 h. Estava causando erros quando hidrodinamico rodava com dt/= 1 h
!	ALLOCATE(TP1(nThd))
!	do i=1,nthd
!		TP1(i)=(i-1)*dthd
!	enddo
	ALLOCATE(TP1(25))
	do i=1,25
		TP1(i)=(i-1)*DTP/(25-1)
	enddo
!**************************************8



!	bounTS(:,bounID(nBoun))=bounTS(1,bounID(nBoun))

!	Qaux(1)=0.0
!	do ik=2,nthd2
!		QAUX(IK)=dthd*(ik-1)
!		write(120123,*) Qaux(ik)/60./60./24.,bounTS(ik,1)
!	enddo





	! Correcao de niveis de fundo muito diferentes em confluencias
	IF (NCONF>0) THEN
		DO I=1,NCONF
			! Conflu, secoes 1, 2 e 3, trechos,sentido do escoamento:
			
!			WRITE(*,*) I, NCC(I,1),NCC(I,2),NCC(I,3), (ZO(ABS(NCC(I,J))),J=1,3),MAX(ZO(ABS(NCC(I,1))),ZO(ABS(NCC(I,2))),ZO(ABS(NCC(I,3))))-MIN(ZO(ABS(NCC(I,1))),ZO(ABS(NCC(I,2))),ZO(ABS(NCC(I,3))))
			i1=ABS(NCC(I,1))
			i2=ABS(NCC(I,2))
			i3=ABS(NCC(I,3))
			IF (MAX(ZO(i1),ZO(i2),ZO(i3))-MIN(ZO(i1),ZO(i2),ZO(i3))>0.5) THEN
				write(*,*) 'Erro: diferença nos niveis de fundo em secoes de confluencia'
				PAUSE
			ENDIF
			aux2=(ZO(i1)+ZO(i2)+ZO(i3))/3.0
			ZO(i1)=aux2
			ZO(i2)=aux2
			ZO(i3)=aux2
			aux2=(HO(i1)+HO(i2)+HO(i3))/3.0
			HO(i1)=aux2
			HO(i2)=aux2
			HO(i3)=aux2

!			DO J=1,3
!				ZO(ABS(NCC(I,J)))=MIN(ZO(ABS(NCC(I,1))),ZO(ABS(NCC(I,2))),ZO(ABS(NCC(I,3))))			
!			ENDDO	
		ENDDO
	ENDIF
	! Corrige niveis de fundo em trechos curtos
	do i=1,nreac
		if (reacFLAG(i)/=2) cycle
		i1=NST(i,1)
		i2=NST(i,2)
		aux2=(ZO(i1)+ZO(i2))/2.0
		ZO(i1)=aux2
		ZO(i2)=aux2
		aux2=(HO(i1)+HO(i2))/2.0
		HO(i1)=aux2
		HO(i2)=aux2


	enddo


	! 
	Vfl=0.0
	! Calcula relação entre nível d'água e volume armazenado na planície de inundação jun/2011
	DO i=1,nReac
		iX=NST(i,2)
		iTr=TrXSec(iX)

		if (DXT(i)==0.0) cycle ! Pula trecho curto

		! Niveis iguais a curva nivel x area alagada
		HF(:,iX)=Zfl(:,iX)	
		do j=2,nPf(iX)
			if (HF(j,iX)<=ZO(iX)) cycle
!			dV=1000000.0*0.5*(Afl(j,iX)+Afl(j-1,iX))*(HF(j,iX)-HF(j-1,iX))
!			dV=dV-DXT(i)*(HF(j,iX)-ZO(iX))*(BA(1,iX)+BA(1,iX-1))*0.5 ! Tira volume no rio. Considera somente seção retangular
			dV=0.5*(FAF(j,iX)+FAF(j-1,iX))*DXT(i)*(HF(j,iX)-HF(j-1,iX)) !@ DCB_HD_Sed
			dV=max(dV,0.0)
			Vfl(j,iX)=Vfl(j-1,iX)+dV
		enddo		

	ENDDO
	



	! Abre arquivos de saida:
!@ DCB ago/2012	open(700,FILE='.\output\Qout.txt')
!@ DCB ago/2012	open(701,FILE='.\output\Yout.txt')
!@ DCB ago/2012	open(702,FILE='.\output\Zout.txt')

	RETURN
	END
