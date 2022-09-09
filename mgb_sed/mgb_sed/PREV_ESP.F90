	SUBROUTINE PREV_ESP
	! Esta rotina faz hindcast de Ensemble Streamflow Prediction (ESP) e Reverse Ensemble Streamflow Prediction (revESP)
	! com a finalidade de analise das fontes de incertezas na previsão com o modelo.
	! Detalhes sobre as técnicas utilizadas são encontrados em Wood and Lettenmaier (2008)
	!
	! 
	! Autoria: Rodrigo Paiva dez/2010
	!------------------------------------------------------------------------------------------------
	! ********************************  Descrição dos arquivos lidos:  ***********************************
	!	
	! Par_revESP.txt - arquivo com informações
	! StateVars.txt - arquivo com todas as variáveis de estado do modelo em todos os intervalo de tempo diários
	! Obs.: As vars de estado são armazenadas na seguite ordem:
	! 
	! W,VBAS,VINT,VSUP,TA,QM2,QJ2,SI,QJ1,QRIOINI,QO,HO,Qmin,HminFLAG Vamos armazenar mais que o necessario.
	! ********************************  Descrição das variáveis principais:  ***********************
	!
	! ESPrevESP_flag = Tipo de simulação:
	! 1 - Simulação referencia
	! 2 - Hindcast of Ensemble Streamflow Prediction
	! 3 - Hindcast of Reverse Ensemble Streamflow Prediction
	! KK = Numero de previsões por ano (K):
	! Data_ESP(.,.) = Dia e mes das K previsões (K,2)
	! lead = Horizonte previsão (lead time)
	! ESP_flag = Tipo de ensemble de forcantes:
	! 1 - Precipitação + Vars Meteorologicas
	! 2 - Precipitação
	! 3 - Vars Meteorologicas
	! revESP_flag = Tipo de ensemble de vars de estado:
	! 1 - Tudo
	! 2 - Águas superficiais (Vsup, Q, z + volume planicie)
	! 3 - Umidade do solo (W)
	! 4 - Águas subterrâneas (Vbas, Vint)
	! nTano = numero de anos
	! nStatVar=numero de variaveis de estado
	!------------------------------------------------------------------------------------------------
	USE VARS_MAIN
	USE AUX_MOD
	USE TIME_MOD
	USE PAR1_MOD
	IMPLICIT NONE

	integer:: ESPrevESP_flag,KK,ESP_flag,revESP_flag,iKK,nTano,iS,iM,lead,nTfim,nStatVar,flag1,iAnos
	integer,allocatable:: Data_ESP(:,:)
	real,allocatable:: stateVec(:)	 
	integer:: iTchuva2,iTclim,iTtemp2,iTvar,iTfor,iTprev
	integer julday	
	integer IMESkk,IDIAkk,iKK0
	integer ilead,iprev


	INTEGER ITTEMP
	INTEGER K,KHID !CONTADORES
	INTEGER KB,JCB,MWP 	!CONTADORES
	INTEGER IANTE,ITP,IHORIZ,INIPREV,IFIPREV,LINIARQ
	CHARACTER (15) ARQPREV !NOME DO ARQUIVO DE CHUVA PREVISTA MODELO GLOBAL CPTEC
	!CHARACTER (17) ARQPREV !NOME DO ARQUIVO DE CHUVA PREVISTA MODELO ETA
	INTEGER ITEMP
	REAL ERROANT1,ERROANT2,ERROCORR
	REAL RIANTE,RIHORIZ
	INTEGER IOPCHU
	CHARACTER (4) ABOBO
	INTEGER IDIAPREV,IANOPREV,IMESPREV,IP

	REAL QPMENS(29)
	INTEGER MDIAS,CONTADOR
	INTEGER ITINI,IBOBO
	CHARACTER (1) BARRA
	INTEGER NTRECH,NTMUSK,LT,I,J
	INTEGER IANTEMES,LDIA(NT),LMES(NT),LANO(NT)





!*************************************************************************************
	! CODIGO:

	OPEN(FILHID,FILE='.\output\VAZAO_revESP.HIG',STATUS='UNKNOWN')

	OPEN(991,FILE='.\output\Qout_ens.txt',STATUS='UNKNOWN')
	

	iprev=0
	! Leitura do arquivo de parametos Par_ESP
	OPEN(555,FILE='.\input\Par_revESP.mgb',STATUS='OLD',ACTION='READ')
	do i=1,4
		read(555,*)
	enddo
	read(555,*) ESPrevESP_flag
	read(555,*)
	read(555,*) KK
	read(555,*)
	allocate(Data_ESP(KK,2))
	do iKK=1,KK
		read(555,*) Data_ESP(iKK,1),Data_ESP(iKK,2)
	enddo
	read(555,*)
	read(555,*) lead
	do i=1,4
		read(555,*)
	enddo
	read(555,*) ESP_flag	
	do i=1,5
		read(555,*)
	enddo
	read(555,*) revESP_flag	
	read(555,*) 
	read(555,*) nTano	

	close(555)

	
	! Inicio simulação:
	nStatVar=sum(NSUBT)+nC
	nStatVar=nStatVar+nC*nU+6*nC+nC*nU+nC+4*nX

	allocate(stateVec(nStatVar))
	! Condições iniciais:
	if (ESPrevESP_flag==1) open(316,FILE='.\input\StateVars.bin',STATUS='REPLACE',RECL=nStatVar,FORM='UNFORMATTED',ACCESS='DIRECT')
	if (ESPrevESP_flag==1) open(317,FILE='.\output\Qcel.txt',STATUS='UNKNOWN')

	if (ESPrevESP_flag==3) open(316,FILE='.\input\StateVars.bin',STATUS='old',RECL=nStatVar,FORM='UNFORMATTED',ACCESS='DIRECT')
	ITINI=1
	IT=0
	nTfim=nT
	do while (IT<nTfim)
		IT=IT+1
		write(*,*) "iT= ", iT
		! Passa variáveis de estado:
		IF(IT==ITINI)THEN !PRIMEIRO INTERVALO DE TEMPO, CONDIÇÕES INICIAIS ESTIMADAS
			call CondInic
		ELSE !OUTROS INTERVALOS DE TEMPO, CONDIÇÕES INICIAIS DADAS PELOS VALORES GUARDADOS
			!DO SOLO
			W=WPREV
			!DOS RESERVATORIOS SUBTERRANEO, SUBSUPERFICIAL E SUPERFICIAL e vazao nos rios Muskingum Cunge
			VBAS=VBASPREV
			VINT=VINTPREV
			VSUP=VSUPPREV
			TA=TAPREV
			QM2=QM2PREV
			QJ2=QJ2PREV
			SI=SIPREV
			QRIOINI=QRIOINIPREV
			QCEL2=QCEL2PREV
			! Incluir variaveis de estado do modelo hidrodinamico !RP
			! Interceptação
			QO=QOprev
			HO=HOprev
			Qmin=Qminprev
			HminFLAG=HminFLAGprev
		ENDIF

		! Modelo:
		JDIA=IDINI+IT-1 !VERIFICA QUAL É O DIA DO CALENDÁRIO 
		CALL CALDAT(JDIA,IMES,IDIA,IANO)
		!SUBROTINA DE LEITURA E PREPARACAO DA CHUVA
		ITCHUVA=IT
		CALL LECHUVA
		!SUBROTINA DE PREPARAÇÃO DOS DADOS CLIMATOL.
		TONTEM=TA
		CALL LECLIMA
		!SUBROTINA DA CELULA
		CALL CELULA
		!SUBROTINA DA REDE DE DRENAGEM
		CALL REDE
		! Modelo hidrodinâmico:
		IF (hdFLAG0>0) CALL Hidrodinamico2

		! Guarda variaveis de estado antes da previsão:
		!DO SOLO
		Wprev=W
		!DOS RESERVATORIOS SUBTERRANEO, SUBSUPERFICIAL E SUPERFICIAL e vazao nos rios Muskingum Cunge
		VBASprev=VBAS
		VINTprev=VINT
		VSUPprev=VSUP
		TAprev=TA
		QM2prev=QM2
		QJ2prev=QJ2
		SIprev=SI
		QRIOINIprev=QRIOINI
		QCEL2prev=QCEL2
		! Incluir variaveis de estado do modelo hidrodinamico !RP
		QOprev=QO
		HOprev=HO
		Qminprev=Qmin
		HminFLAGprev=HminFLAG

		!******************************************************************************
		! Escreve variáveis de estado:
		i=0
		do iC=1,nC
			do iU=1,nU
				i=i+1
				StateVec(i)=W(iC,iU)
			enddo
		enddo
		do iC=1,nC
				i=i+1
				StateVec(i)=VBAS(iC)
		enddo
		do iC=1,nC
				i=i+1
				StateVec(i)=VINT(iC)
		enddo
		do iC=1,nC
				i=i+1
				StateVec(i)=VSUP(iC)
		enddo
		do iC=1,nC
				i=i+1
				StateVec(i)=TA(iC)
		enddo
		do iC=1,nC
				i=i+1
				StateVec(i)=QM2(iC)
		enddo
		do iC=1,nC
				i=i+1
				StateVec(i)=QJ2(iC)
		enddo
		do iC=1,nC
			do iU=1,NU
				i=i+1
				StateVec(i)=SI(iC,iU)
			enddo
		enddo
		do iC=1,nC
			do j=1,NSUBT(IC)+1
				i=i+1
				StateVec(i)=QRIOINI(iC,j)
			enddo
		enddo
		do iC=1,nC
				i=i+1
				StateVec(i)=QCEL2(iC)
		enddo
		do iX=1,nX
				i=i+1
				StateVec(i)=QO(iX)
		enddo
		do iX=1,nX
				i=i+1
				StateVec(i)=HO(iX)
		enddo
		do iX=1,nX
				i=i+1
				StateVec(i)=Qmin(iX)
		enddo
		do iX=1,nX
				i=i+1
				StateVec(i)=HminFLAG(iX)
		enddo
		if (ESPrevESP_flag==1) write(316,REC=iT) (StateVec(i),i=1,nStatVar)
		if (ESPrevESP_flag==1) write(317,778) (QJ2(iC),iC=1,nC)
		
		!**********************************************************************


		ITTEMP=IT
		iTprev=iT	!!!

		! Verifica se é dia de previsão:
		flag1=0
		do iKK=1,KK
			if (IDIA==Data_ESP(iKK,1).and.IMES==Data_ESP(iKK,2)) then
				IMESkk=IMES
				IDIAkk=IDIA
				iKK0=iKK
				flag1=1
			endif						 
		enddo

		if (flag1==1.and.ESPrevESP_flag/=1) then
		
			iprev=iprev+1
			! Loop do ensemble
			do iAnos=1,nTAno
				! Recupera vars de estado:
				W=WPREV
				VBAS=VBASPREV
				VINT=VINTPREV
				VSUP=VSUPPREV
				TA=TAPREV
				QM2=QM2PREV
				QJ2=QJ2PREV
				SI=SIPREV
				QRIOINI=QRIOINIPREV
				QCEL2=QCEL2PREV
				QO=QOprev
				HO=HOprev
				Qmin=Qminprev
				HminFLAG=HminFLAGprev

				iTchuva2=iT
				iTclim=iT
				iTvar=iT
				iTfor=iT
					
				! Prepara ensemble:
				selectcase (ESPrevESP_flag)
				case(2)
					! 2 - Hindcast of Ensemble Streamflow Prediction
					! Corrige o iT:
					CALL CALDAT(IDINI,i,j,IANO)
					IANO=IANO+iAnos-1
					JDIA=iT+IDINI-1
					CALL CALDAT(JDIA,IMES,IDIA,j)

					JDIA=JULDAY(IMES,IDIA,IANO)
					! ESP_flag = Tipo de ensemble de forcantes:
					! 1 - Precipitação + Vars Meteorologicas
					! 2 - Precipitação
					! 3 - Vars Meteorologicas
					iTchuva2=iT
					iTclim=iT
					if (ESP_flag==1.or.ESP_flag==2) iTchuva2=JDIA-IDINI+1
					if (ESP_flag==1.or.ESP_flag==3) iTclim=JDIA-IDINI+1
					iTfor=JDIA-IDINI+1
				case(3)
					! 3 - Hindcast of Reverse Ensemble Streamflow Prediction
					! Corrige o iT:
					CALL CALDAT(IDINI,i,j,IANO)
					IANO=IANO+iAnos-1

					JDIA=iT+IDINI-1
					CALL CALDAT(JDIA,IMES,IDIA,j)
					JDIA=JULDAY(IMES,IDIA,IANO)
					iTvar=JDIA-IDINI+1
					
					! Le variáveis de estado:
					read(316,REC=iTvar) (StateVec(i),i=1,nStatVar)
					
					! revESP_flag = Tipo de ensemble de vars de estado:
					! 1 - Tudo
					! 2 - Águas superficiais (Vsup, Q, z + volume planicie)
					! 3 - Umidade do solo (W)
					! 4 - Águas subterrâneas (Vbas, Vint)
					i=0
					do iC=1,nC
						do iU=1,nU
							i=i+1
							if (revESP_flag==1.or.revESP_flag==3)	W(iC,iU)=StateVec(i)
						enddo
					enddo
					
					do iC=1,nC
							i=i+1
							if (revESP_flag==1.or.revESP_flag==4)	VBAS(iC)=StateVec(i)
					enddo
					
					do iC=1,nC
							i=i+1
							if (revESP_flag==1.or.revESP_flag==4)	VINT(iC)=StateVec(i)
					enddo
					

					

					do iC=1,nC
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) VSUP(iC)=StateVec(i)
					enddo
					
					
					do iC=1,nC
							i=i+1
							TA(iC)=StateVec(i)
					enddo
	
					
					do iC=1,nC
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) QM2(iC)=StateVec(i)
					enddo
					do iC=1,nC
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) QJ2(iC)=StateVec(i)
					enddo
					do iC=1,nC
						do iU=1,NU
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) SI(iC,iU)=StateVec(i)
						enddo
					enddo
					do iC=1,nC
						do j=1,NSUBT(IC)+1
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) QRIOINI(iC,j)=StateVec(i)
						enddo
					enddo
					do iC=1,nC
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) QCEL2(iC)=StateVec(i)
					enddo
					do iX=1,nX
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) QO(iX)=StateVec(i)
					enddo
					do iX=1,nX
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) HO(iX)=StateVec(i)
					enddo
					do iX=1,nX
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) Qmin(iX)=StateVec(i)
					enddo
					do iX=1,nX
							i=i+1
							if (revESP_flag==1.or.revESP_flag==2) HminFLAG(iX)=StateVec(i)
					enddo
					
				endselect

			
				! Loop do lead time:
				ilead=0				
				do while (iTprev<iTTEMP+lead)
					iTprev=iTprev+1	!!!
					iT=iT+1
					if (iT>nT) iT=1	!!!
					iTchuva2=iTchuva2+1
					if (iTchuva2>nT) iTchuva2=1 ! Pega dados do primeiro ano
					iTclim=iTclim+1
					if (iTclim>nT) iTclim=1 ! Pega dados do primeiro ano
					iTfor=iTfor+1
					if (iTfor>nT) iTfor=1 ! Pega dados do primeiro ano
					iTvar=iTvar+1
					if (iTvar>nT) iTvar=1 ! Pega dados do primeiro ano

					ilead=ilead+1
					write(*,*) iprev,iKK0,IDIAkk,IMESkk,ITTEMP,iAnos,IDIA,IMES,IANO,ilead,iTprev,iTfor,iTvar,itchuva2,itclim	!!!
					! Modelo:
					JDIA=IDINI+iTchuva2-1 !VERIFICA QUAL É O DIA DO CALENDÁRIO 
					CALL CALDAT(JDIA,IMES,IDIA,IANO)
					!SUBROTINA DE LEITURA E PREPARACAO DA CHUVA
					ITCHUVA=iTchuva2
					CALL LECHUVA
					!SUBROTINA DE PREPARAÇÃO DOS DADOS CLIMATOL.
					iTtemp2=iT
					iT=iTclim
					JDIA=IDINI+iT-1 !VERIFICA QUAL É O DIA DO CALENDÁRIO 
					CALL CALDAT(JDIA,IMES,IDIA,IANO)
					TONTEM=TA
					CALL LECLIMA
					
					iT=iTtemp2
					JDIA=IDINI+iT-1 !VERIFICA QUAL É O DIA DO CALENDÁRIO 
					CALL CALDAT(JDIA,IMES,IDIA,IANO)
					!SUBROTINA DA CELULA
					CALL CELULA
					!SUBROTINA DA REDE DE DRENAGEM
					CALL REDE
					! Modelo hidrodinâmico:
					IF (hdFLAG0>0) CALL Hidrodinamico2


					write(991,777) iprev,iKK0,IDIAkk,IMESkk,ITTEMP,iAnos,IDIA,IMES,IANO,ilead,iTprev,iT,iTfor,iTvar,(QJ2(iC),iC=1,nC) 
					write(992,*) iprev,iKK0,IDIAkk,IMESkk,ITTEMP,iAnos,IDIA,IMES,IANO,ilead,iTprev,iT,iTfor,iTvar
				enddo
				IT=ITTEMP
				iTprev=ITTEMP	
			enddo
			
		endif



	enddo
		



!*****************************************************************************************8	




71	FORMAT(13I6,<NUMHIDG>F12.3)
75	FORMAT(I6,24F7.0)
777	FORMAT(14I7,<nC>F12.3)
778	FORMAT(<nC>F12.3)

     
	RETURN
	END

