	subroutine CondInic
	! Subrotina para gerar condicoes iniciais do modelo
	! Estima condicoes iniciais e aquece o modelo por alguns anos utilizando os primeiros anos de dados
	!-------------------------------------------------------------------------------------
	!
	! Descri��o das vari�veis locais:
	!
	!--------------------------------------------------------------------------------------
	! Declara��o de vari�veis:

	USE AUX_MOD
	USE VARS_MAIN

	implicit none

	! Vari�veis locais de c�lculo:
	integer JUS,i,NN,TK
	!-------------------------------------------------------------------------------------

	
	!ZERANDO VARIAVEIS
	PM=0.0
	PM2=0.0
	SIM=0.0 !GUARDA LAMINA INTERCEPTADA
	EM=0.0 !GUARDA LAMINA EVAPOTRANSPIRA��O
	WBM=0.0 !LAMINA ARMAZENADA MEDIA
	DSUM=0.0 !ESCOMENTO SUPERFICIAL DA BACIA
	DINM=0.0 !ESCOMENTO INTERMEDI�RIO DA BACIA
	DBAM=0.0 !ESCOMENTO SUBTERR�NEO DA BACIA
			
	
	!****************************CONDI��ES INICIAIS*****************************************
	
	! Volume interceptacao:
	SI=0.0
	! Umidade do solo:
	DO IC=1,NC
		IB=IBAC(IC)
		DO IU=1,NU
			W(IC,IU)=WM(IB,IU)*0.40 !CONSIDERA QUE SOLO TEM 40% AGUA NO IN�CIO DA SIMULA��O
		ENDDO
	ENDDO
	! Volumes nos reservatorios:
	DO IC=1,NC
		IB=IBAC(IC)
		VBAS(IC)=QESP(IB)*ACEL(IC)*CB(IB)*3600.
		VINT(IC)=0.0
		VSUP(IC)=0.0

	ENDDO
	! Vaz�o nos reservatorios:
	DO IC=1,NC
		! Bacia
		IB=IBAC(IC)
		
		! Vaz�o subterranea:
		TK=CB(IB)*3600.
		QBAS(IC)=VBAS(IC)/TK
		! Vaz�o subsuperficial:
		TK=CI(IB)*TIND(IC)
		QINT(IC)=VINT(IC)/TK
		! Vaz�o superficial:
		TK=CS(IB)*TIND(IC)
		QSUP(IC)=VSUP(IC)/TK
	ENDDO
	
	!Temperatura no dia anterior:
	TA=20.0

	!****************************FIM DAS CONDI��ES INICIAIS*********************************

	
	! Condicoes iniciais de vazao nos rios
	!Condicoes inicias passam para a rotina condinic
!********************************************************************************
	!CONDI��ES INICIAIS DA VAZAO NO RIO
	! OBS RP Melhorar condicoes iniciais do MGB, usar QREF no lugar de vazao zero nos rios ********************

	QM2=0.0
	PMB2=0.0
	PMI2=0.0
	PMS2=0.0
	QM1=0.0
	QJ1=0.0
	QJ2=0.0
	QCEL1=0.0
	QCEL2=0.0
	QRIOINI=0.0
	QCONTORM=0.0
	QCONTORJ=0.0
	
	
	DO IC=1,NC
		QCEL2(IC)=QBAS(IC)+QINT(IC)+QSUP(IC) !SOMA VAZOES GERADAS NA CELULA
		QCEL1(IC)=QCEL2(IC)

!*****************************************************************
		! Testa se propaga��o na c�lula � computada com modelo hidrodin�mico:
		IF (hdFLAG(IC)==1) CYCLE
!*****************************************************************

		IF(NSUBT(IC).GT.0)THEN	!CELULA  COM RIO

			!ACUMULA VAZAO
			QM2(IC)=QM2(IC)+QCEL2(IC)
			
			!VAZAO CONSTANTE NO TRECHO
			QJ2(IC)=QM2(IC)	

			!O QUE SAI DO RIO DE UMA CEL. VAI PRA OUTRA
			JUS=CELJUS(IC)
			if (JUS>0) QM2(JUS)=QM2(JUS)+QJ2(IC) 

		ELSE !CELULA FONTE
			QJ2(IC)=QCEL2(IC)
			QM2(IC)=QCEL2(IC)

			JUS=CELJUS(IC)
			
			if (JUS>0) QM2(JUS)=QM2(JUS)+QCEL2(IC)
		ENDIF


		! Verifica se tem substitui��o de vazao:
		do i=1,NUMSUBST 
			if (iC==ISUBST(i)) then
				QM2(JUS)=QLIDO(i,1)
				exit			
			endif
		enddo
	ENDDO
	QJ1=QJ2
	QM1=QM2

	!CONDI��O DE CONTORNO MUSKINGUM-CUNGE
	do iC=1,nC	
		QRIOINI(IC,:)=QM2(IC)
		QCONTORM(IC,:)=QM2(IC)
		QCONTORJ(IC,:)=QJ2(IC)
	enddo


	! Verifica minibacias que estao a montante de subtituicao de vazoes:
	sbtFLAG=-1
!	do i=1,NUMSUBST 
!		iC=ISUBST(i)
!		sbtFLAG(iC)=1
!	enddo
	
!	do iC=1,nC
!		! Desce at� encontrar 
!		JUS=iC
!		do
!			if (JUS<0) then
!				sbtFLAG(iC)=0
!				exit
!			endif
!			if (sbtFLAG(JUS)/=-1) then
!				sbtFLAG(iC)=sbtFLAG(JUS)
!				exit
!			endif
!			JUS=CELJUS(JUS)
!		enddo
!	enddo
!	do i=1,NUMSUBST 
!		iC=ISUBST(i)
!		sbtFLAG(iC)=0
!	enddo

	


	!*****************************************************************
		! Modelo hidrodin�mico:
		IF (hdFLAG0>0) then
			 CALL HD_CondInic
		ENDIF
	!*****************************************************************


	! Aquecimento do modelo MGB
	! Simula bacia nos NN primeiros intervalos de tempo:
	NN=365




	return
	end subroutine
