	SUBROUTINE CELULA	       	                 
	!esta subrotina comanda o loop das celulas
	USE VARS_MAIN
	use aux_mod
!	
!************ DIOGO BUARQUE ********************
	USE SED_VARS
!************ DIOGO BUARQUE ********************
!
	IMPLICIT NONE
	!RETARDOS E TEMPO DE CONCENTRAÇÃO
	REAL TKB(NC),TKI(NC),TKS(NC)
	!LAMINAS DRENADAS
	REAL DSUP,DINT,DBAS
	REAL DSB,DIB,DBB
	!VARIAVEIS AUXILIARES
	!	REAL WX,WMX,PAX,TAX,VVX,URX,SOLX,BX,KIX,KBX,RIAFX,ALBX,ZX,EX,PX,RLX,SIX
	!	REAL XL,T1,T2,RSX,CAPX,ETX,WCX,VBX,VIX,VSX,GLAT,PPU,WXB,SIB,ETB,DCAP
	!OUTRAS
	!	REAL CLATE,ED,ES
	!INTEGER IB,IU
	INTEGER JULDAY

	!VERIFICA A QUE MES CORRESPONDE O DIA JULIANO_____
	CALL CALDAT(JDIA,IMES,IDIA,IANO)

    !mudanças climáticas
!    SELECT CASE(IMES)
!        CASE (1) !JANEIRO
!            TA=TA+2.7
!            P=P*1.168
!        CASE (2) !FEVEREIRO
!            TA=TA+2.5
!            P=P*1.115
!        CASE (3) !MARÇO
!            TA=TA+2.7
!            P=P*1.101
!        CASE (4) !ABRIL
!            TA=TA+2.8
!            P=P*1.053
!        CASE (5) !MAIO
!            TA=TA+3.1
!            P=P*1.261
!        CASE (6) !JUNHO
!            TA=TA+3.1
!            P=P*1.266
!        CASE (7) !JULHO
!            TA=TA+3.0
!            P=P*1.140
!        CASE (8) !AGOSTO
!            TA=TA+2.8
!            P=P*1.150
!        CASE (9) !SETEMBRO
!            TA=TA+2.7
!            P=P*0.987
!        CASE (10) !OUTUBRO
!            TA=TA+2.8
!            P=P*1.166
!        CASE (11) !NOVEMBRO
!            TA=TA+2.6
!            P=P*1.115
!        CASE (12) !DEZEMBRO
!            TA=TA+2.7
!            P=P*1.053
!    END SELECT

	
    JDIA=JDIA-JULDAY(1,1,IANO)+1

	TWS=0.0 ! RP
	TWS1=0.0
	TWS2=0.0
	TWS3=0.0
	TWS4=0.0
	TWS5=0.0
	TWS6=0.0
	TWS7=0.0
	Ebacia=0.0
	Pbacia=0.0	


	DO IC=1,NC	 	!LOOP DAS CELULAS

		if (sbtFLAG(iC)==1) cycle ! RP nao computa celulas a montante de celula com substituicao de dados

		IB=IBAC(IC)
!@ DCB_sed ###########################################################################
	    !@ Calcula apenas sub-bacia de interesse.
		!IF(IB>SUBfim.OR.IB<SUBini)CYCLE
!@ DCB_sed ###########################################################################
	

		EVQ(IC)=0.0 !EVAPORAÇÃO DIRETA DA SUPERFÍCIE LÍQUIDA
		DO IU=1,NU !LOOP DOS USOS

			IF(PUSO(IC,IU).LT.0.0001)THEN !NAO TEM ESTE USO NESTA CELULA
				CYCLE !PASSA PARA O PROXIMO USO
			ENDIF
 
 			PX=P(IC) !A CHUVA INICIA IGUAL, MAS É INTERCEPTADA
			WX=W(IC,IU)			  !VARIAVEL AUXILIAR
			SIX=SI(IC,IU)		  !VARIAVEL AUXILIAR
			WMX=WM(IB,IU)		  !VARIAVEL AUXILIAR
			PAX=PA(IC)            !VARIAVEL AUXILIAR
			TAX=TA(IC)			  !VARIAVEL AUXILIAR
			VVX=VV(IC)			  !VARIAVEL AUXILIAR
			URX=UR(IC)			  !VARIAVEL AUXILIAR
			SOLX=SOL(IC)		  !VARIAVEL AUXILIAR
			XL=PLAM(IB,IU)		  !VARIAVEL AUXILIAR
			T1=TONTEM(IC)         !TEMP. DO DIA ANTERIOR 
			T2=TA(IC)             !TEMPERATURA DO DIA ATUAL
			RIAFX=RIAF(IU,IMES)   !VARIAVEL AUXILIAR
			ALBX=ALB(IU,IMES)     !VARIAVEL AUXILIAR
			ZX=Z(IU,IMES)		  !VARIÁVEL AUXILIAR
			RSX=RS(IU,IMES)		  !VARIÁVEL AUXILIAR
			BX=B(IB,IU)	  !VARIAVEL AUXILIAR
			KIX=KINS(IB,IU)  !VARIAVEL AUXILIAR
			KBX=KBAS(IB,IU)  !VARIAVEL AUXILIAR
			CAPX=CAP(IB,IU)	 !VARIAVEL AUXILIAR
			WCX=WC(IB,IU)	 !VARIAVEL AUXILIAR

		!	IF(PAX.GT.130.0.OR.PAX.LT.70.0)STOP 'ERRO PRESSAO SUB CELULA'
			IF(TAX.LT.-50.0.OR.TAX.GT.50.0)STOP 'ERRO TEMP SUB CELULA'
			IF(VVX.LT.0.0.OR.VVX.GT.100.)STOP 'ERRO VEL VENTO SUB CELULA'
			IF(URX.LT.0.0.OR.URX.GT.100.01)STOP 'ERRO UMIDADE SUB CELULA'
		!	IF(SOLX.LT.0.0.OR.SOLX.GT.24.)STOP 'ERRO INSOLACAO SUB CELULA'

			!	CALCULA A RADIAÇÃO SOLAR LÍQUIDA
			!define a latitude da célula
			GLAT=Y(IC)
			CALL RADIACAO  !RLX É A RADIAÇÃO LIQUIDA


			!CALCULA EVAPOTRANSPIRAÇÀO E INTERCEPTA A CHUVA
			CALL EVAPO

!			ET(IC,IU)=ETX

!Linha para testar modelo com evapo zero:
!			EIX=0.0
!			EX=0.0
!			ET(IC,IU)=0.0
!			PX=P(IC)



			!FAZ O BALANÇO HÍDRICO DO SOLO, SE FOR ÁREA COBERTA POR ÁGUA (WM=0.0) PASSA DIRETO
			IF(WMX.GT.0.001) THEN
				CALL SOLO(PX,EX,WX,WMX,BX,KIX,KBX,XL,DSUP,DINT,DBAS,CAPX,WCX,DTP)
				ETX=EX+REIX
			ELSE
				IF (HDFLAG(iC)==0) THEN
					DSUP=PX
					DINT=0.0
					DBAS=0.0
					CAPX=0.0
					EVQ(IC)=(EIX*1000.*ACEL(IC)*(PUSO(IC,IU)/100.))/DTP !EVAPORAÇÃO DIRETA DAS SUPERFÍCIES LÍQUIDAS EM M3/S
				ELSE	
					!*******
					! Consideracao sobre evaporacao em planicie de inundacao e água: !RP
					DSUP=PX-EIX 
					DINT=0.0
					DBAS=0.0
					CAPX=0.0
					!*******
				ENDIF
				ETX=EIX
			ENDIF

			
			ET(IC,IU)=ETX

			W(IC,IU)=WX
			SI(IC,IU)=SIX
			CAF(IC,IU)=CAPX

			PPU=PUSO(IC,IU)/100. !FRAÇÃO DE USO DO SOLO
			WXB=WXB+WX*PPU !ARMAZENAMENTO MEDIO NA BACIA
			SIB=SIB+(P(IC)-PX)*PPU !INTERCEPTAÇÃO MÉDIA NA BACIA
			ETB=ETB+ETX*PPU !EVAPOTRANSPIRAÇÃO MÉDIA NA BACIA
!			
!
!
!************ DIOGO BUARQUE ********************
!@ ARMAZENA LÂMINA DO ESCOAMENTO SUPERFICIAL (MM) DE CADA USO DE CADA CÉLULA
!@ PARA UTILIZAÇÃO NA SUBROTINA MUSLE
			QSAUX(IC,IU+2) = DSUP !@ (mm)
!************ DIOGO BUARQUE ********************
!
!
!
			!CORRIGE AS UNIDADES; MULTIPLICA PELA ÁREA; MULTIPLICA PELA PROPORÇÃO DE USO
			!OS VALORES DE DSUP, DINT E DBAS ESTAO EM MM/DTP - CONVERTE PARA M3/S
			DSUP=((DSUP*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DSUP=((DSUP*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			DINT=((DINT*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DINT=((DINT*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			DBAS=((DBAS*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(DTP)	!DBAS=((DBAS*1000.)*ACEL(IC)*(PUSO(IC,IU)/100.))/(3600.*24.)
			!FLUXO CAPILAR TOTAL - CONVERTE MM/DIA PARA M3/S
			DCAP=DCAP+(PUSO(IC,IU)/100.)*(CAPX*1000.*ACEL(IC))/(DTP) !DCAP=DCAP+(PUSO(IC,IU)/100.)*(CAPX*1000.*ACEL(IC))/(3600.*24.) 
			!ATUALIZA VOLUMES (EM M3)
			VBAS(IC)=VBAS(IC)+DBAS*DTP   !VBAS(IC)=VBAS(IC)+DBAS*3600.*24.
			VINT(IC)=VINT(IC)+DINT*DTP   !VINT(IC)=VINT(IC)+DINT*3600.*24.
			VSUP(IC)=VSUP(IC)+DSUP*DTP   !VSUP(IC)=VSUP(IC)+DSUP*3600.*24.

			TWS(iC)=TWS(iC)+(SI(IC,IU)+W(IC,IU))*(PUSO(IC,IU)/100.) ! Cálculo em mm

			TWS1(iC)=TWS1(iC)+(W(IC,IU))*(PUSO(IC,IU)/100.) ! Cálculo em mm


		


		ENDDO !FIM DO LOOP DOS USOS

		!RETIRA FLUXO CAPILAR DO RESERVATORIO SUBTERRANEO
!@ DCB 26-05-2011 **************************************************************************
!		VBAS(IC)=VBAS(IC)-DCAP*DTP !VBAS(IC)=VBAS(IC)-DCAP*3600.*24. !ATUALIZA VOLUME
		VBAS(IC)=max(VBAS(IC)-DCAP*DTP,0.0) !VBAS(IC)=VBAS(IC)-DCAP*3600.*24. !ATUALIZA VOLUME
!@ DCB 26-05-2011 **************************************************************************
		DCAP=0.0
		!CALCULA VAZOES DAS CELULAS
		TKB(IC)=CB(IB)*3600. !CB DEFINIDO EM HORAS NO ARQUIVO DE PARAMETROS
		TKI(IC)=CI(IB)*TIND(IC) !CI DEFINIDO NO ARQUIVO DE PARAMETROS
		TKS(IC)=CS(IB)*TIND(IC) !CS DEFINIDO NO ARQUIVO DE PARAMETROS
		QBAS(IC)=VBAS(IC)/TKB(IC)
		QINT(IC)=VINT(IC)/TKI(IC)
		QSUP(IC)=VSUP(IC)/TKS(IC)


!
!
!************ DIOGO BUARQUE ********************
!@ ARMAZENA PARÂMETRO DE RETARDO DO ESCOAMENTO SUPERFICIAL E O VOLUME
!@ TOTAL GERADO NA CÉLULA PARA UTILIZAÇÃO NA SUBROTINA MUSLE
		QSAUX(IC,1) = TKS(IC)
		QSAUX(IC,2) = VSUP(IC)
!************ DIOGO BUARQUE ********************
!
!
!
		!ATUALIZA VOLUMES (EM M3)
		VBX=VBAS(IC)-QBAS(IC)*DTP       !VBX=VBAS(IC)-QBAS(IC)*3600.*24.
		VIX=VINT(IC)-QINT(IC)*DTP       !VIX=VINT(IC)-QINT(IC)*3600.*24.
		VSX=VSUP(IC)-QSUP(IC)*DTP       !VSX=VSUP(IC)-QSUP(IC)*3600.*24.
		
		
		!VERIFICA SE NAO SECOU RESERVATORIO
		
	
		
		! Considera que reservatorio superficial pode ter volume negativo para considerar evaporacao na planicie de inundacao: !RP
		! Codigo original MGB:
!		IF(VSX.LT.0.0)THEN 		!SUPERFICIAL
!			QSUP(IC)=VSUP(IC)/(DTP)            !QSUP(IC)=VSUP(IC)/(3600.*24.)
!			VSUP(IC)=0.0
!		ELSE
!			VSUP(IC)=VSX
!		ENDIF
		
		!SUPERFICIAL:

		
		if (hdFLAG(iC)==0) then
			IF(VSX.LT.0.0)THEN 		!SUPERFICIAL
				QSUP(IC)=VSUP(IC)/(DTP)            !QSUP(IC)=VSUP(IC)/(3600.*24.)
				VSUP(IC)=0.0
			ELSE
				VSUP(IC)=VSX
			ENDIF
		else

			IF(VSUP(IC)<0.0)THEN 		
				! Vazao negativa é retirada instantaneamente no proximo intervalo de tempo, sem defasagem do reservatorio linear:
				QSUP(IC)=VSUP(IC)/(DTP)  !QSUP(IC)=VSUP(IC)/(3600.*24.)
				VSUP(IC)=0.0
			ELSEIF(VSX<0.0)THEN
				QSUP(IC)=VSUP(IC)/(DTP)            !QSUP(IC)=VSUP(IC)/(3600.*24.)
				VSUP(IC)=0.0
			ELSE
				VSUP(IC)=VSX
			ENDIF		
		
		endif
	

		IF(VIX.LT.0.0)THEN !SUB-SUPERFICIAL
			QINT(IC)=VINT(IC)/(DTP)             !QINT(IC)=VINT(IC)/(3600.*24.)
			VINT(IC)=0.0
		ELSE
			VINT(IC)=VIX
		ENDIF
		IF(VBX.LT.0.0)THEN !SUBTERRANEO
			QBAS(IC)=VBAS(IC)/(DTP)             !QBAS(IC)=VBAS(IC)/(3600.*24.)
			VBAS(IC)=0.0
		ELSE
			VBAS(IC)=VBX
		ENDIF
		!FIM DA VERIFICAÇÃO DE RESERVATORIO SECO

		!IB É O NUMERO DA BACIA A QUAL PERTENCE A CELULA IC

		PM2(IB,IT)=PM2(IB,IT)+P(IC)/KCB(IB) !GUARDA PRECIPITAÇÃO
		DSB=(QSUP(IC)*3600.*24.)/(1000.*ACEL(IC))
		DIB=(QINT(IC)*3600.*24.)/(1000.*ACEL(IC))
		DBB=(QBAS(IC)*3600.*24.)/(1000.*ACEL(IC))
		
		SIM(IB,IT)=SIM(IB,IT)+SIB/KCB(IB) !GUARDA LAMINA INTERCEPTADA
		EM(IB,IT)=EM(IB,IT)+ETB/KCB(IB) !GUARDA LAMINA EVAPOTRANSPIRAÇÃO
		WBM(IB,IT)=WBM(IB,IT)+WXB/KCB(IB) !LAMINA ARMAZENADA MEDIA
		DSUM(IB,IT)=DSUM(IB,IT)+DSB/KCB(IB) !ESCOMENTO SUPERFICIAL DA BACIA
		DINM(IB,IT)=DINM(IB,IT)+DIB/KCB(IB) !ESCOMENTO INTERMEDIÁRIO DA BACIA
		DBAM(IB,IT)=DBAM(IB,IT)+DBB/KCB(IB) !ESCOMENTO SUBTERRÂNEO DA BACIA
		


		! Evapo média:
		if (IB<=97.or.IB>=169) Ebacia=Ebacia+ETB*ACEL(iC)/maxval(ACUR)
		! Pmedia:
		if (IB<=97.or.IB>=169) Pbacia=Pbacia+P(IC)*ACEL(iC)/maxval(ACUR)
		
		E0media(iC)=ETB			

		SIB=0.0 !ZERA VARIAVEIS TEMPORARIAS
		ETB=0.0
		WXB=0.0
		DSB=0.0
		DIB=0.0
		DBB=0.0 !ZERA VARIAVEIS TEMPORARIAS


	

		! Calculo de indicadores para evapotranspiração:
		! Evaporação de superficies liquidas:
		E0agua(iC)=EIX
		! Evaporação utilizando toda energia que chega no topo da atmosfera (calculado pela constante solar)
		E0TOPO(iC)=1000.*STO/(MESP*CLATE)
		! Evaporação utilizando radiação de onda curta que chega na superficie terrestre (dados de entrada)
		E0SUP(iC)=1000.*SSUP/(MESP*CLATE)




		TWS(iC)=TWS(iC)+(VBAS(IC)+VINT(IC)+VSUP(IC))/(ACEL(IC)*1000.0) ! Cálculo em mm

		TWS2(iC)=TWS2(iC)+(VINT(IC))/(ACEL(IC)*1000.0) ! Cálculo em mm
		TWS3(iC)=TWS3(iC)+(VBAS(IC))/(ACEL(IC)*1000.0) ! Cálculo em mm
		TWS4(iC)=TWS4(iC)+(VSUP(IC))/(ACEL(IC)*1000.0) ! Cálculo em mm

	ENDDO !FIM DO LOOP DAS CÉLULAS


		

79	FORMAT(I5,12F6.1)
	RETURN
	END