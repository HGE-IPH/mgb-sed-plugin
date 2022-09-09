!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abril de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SED_REDE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - SUBROTINA PARA PROPAGAR A CERGA DE SEDIMENTO AO LONGO DA DEDE DE DRENAGEM
!@
!@ ---------------------------------------------------------------------------------------------
!@  VARIÁVEIS
!@ Descrição das variáveis locais:
!@
!@ RUGMAN   = Rugosidade de Manning para MC
!@ QmTREC   = Vazão média no trecho (m^3/s)
!@ HmTREC   = Profundidade média no trecho por Manning, assumindo Raio = HmTREC (m)
!@ VmTREC   = Velocidade média no trecho (m/s)
!@ RHm      =
!@ Uat      =
!@ SfTREC   =
!@ VISC     =
!@ FatM     =
!@ FatN     =
!@ FatP     =
!@ DCar     =
!@ REX      =
!@ UcWs     =
!@ FMY      =
!@ FNY      =
!@ LOGCTS   =
!@ HK       =
!@ ERROH    =
!@ ERROQ    =
!@ TOL      =
!@ QK       =
!@ QK2      =
!@ DEPaux   =
!@ EROSaux  =
!@ CARGTREC = carga de sedimento final no trecho (ton)
!@ ---------------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

SUBROUTINE SED_REDE

!@ **********************************
USE AUX_MOD
USE VARS_MAIN
USE SED_VARS
!@ **********************************

IMPLICIT NONE

!@ **********************************
!@ Declaração de variáveis:
REAL HmJ, VmJ, RHmJ, SfTREC, UatJ
REAL LIMe, DeT, LIMd
REAL PONDt
REAL EQ1, EQ2, EQ3, EQ4
REAL DEPaux, DEPaux2(3), EROSaux, EROSaux2(3)
INTEGER iTREC, i, JUS
INTEGER VERIF1, VERIF2, VERIF3

REAL LIMeaux(3), LIMdaux(3)
REAL iENTRAaux(3), iSAIaux(3)

REAL QM1aux, QM2aux, QJ1aux, QJ2aux !@ DCB CONTINUIDADE 

ALLOCATE(CSJ2aux(nSEDmini,3))
!@ **********************************


!@ *****************************************************************************************
!@ EM CADA INTERVALO DE TEMPO O QUE ERA i+1 VIRA i
!@ -----------------------------------------------------------------------------------------
DO IC=1,NC
	CSM1(IC,:)      = CSM2(IC,:)
	CSJ1(IC,:)      = CSJ2(IC,:)
	VolTREC1(IC)    = VolTREC2(IC)
ENDDO
!@ *****************************************************************************************


!@ INICIALIZANDO VARIÁVEIS
PONDt   = 1.0   !@ Ponderador temporal da discretização por diferenças finitas
iTREC   = 0     !@ Código do trecho HD da minibacia
iSEDaux = 0     !@ Contador de minibacia com cálculo de sedimentos
CARGM   = 0.0   !@ Carga a montante do trecho
CSJ2    = 0.0


DO IC = 1,NC    !@ INICIO DO LOOP DAS MINIBACIAS

    IB=IBAC(IC)
    
    !@ CALCULA APENAS SUB-BACIAS SELECIONADAS
!	IF(IB>SUBfim .OR. IB<SUBini)CYCLE
	IF(IB>82 .OR. IB<57)CYCLE

	iSEDaux = iSEDaux + 1 !@ Contador de minibacia com cálculo de sedimentos
	
	QJ1aux = 10.0           !@ DCB CONTINUIDADE 
	QM1aux = 10.0           !@ DCB CONTINUIDADE 
	QJ2aux = 10.0           !@ DCB CONTINUIDADE 
	QM2aux = 10.0           !@ DCB CONTINUIDADE 
	PSC(IC,1) = 3.*50.0     !@ DCB CONTINUIDADE 
	PSC(IC,2:4) = 50.0      !@ DCB CONTINUIDADE 


!@ #########################################################################################
! ------------------------------------------------------------------------------------------
!@ ##########################      TRECHOS SIMULADOS COM MC       ##########################
! ------------------------------------------------------------------------------------------
!@ #########################################################################################
	IF (hdFLAG(IC) == 0) THEN

        !@ *****************************************************************************************
	    !@ CARACTERÍSTICAS HIDRÁULICAS DA SEÇÃO DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
        !@ Declividade de atrito no trecho (mínima de 0.01 m/km para evitar problemas com o
        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vazões maiores
        !@ que 5000 m3/s)
		    SfTREC = max(DECL(IC),0.00001)
    	
        !@ Profundidade média na seção de jusante por Manning assumindo Raio = HmTREC (m)
!		    HmJ = ((RUGMAN*QJ2(IC))/(BRIO(IC)*SfTREC**0.5))**(3.0/5.0)
		    HmJ = ((RUGMAN*QJ2aux)/(BRIO(IC)*SfTREC**0.5))**(3.0/5.0)  !@ DCB CONTINUIDADE 

        !@ Volume médio de água no trecho
            VolTREC2(IC) = HmJ*BRIO(IC)*SRIO(IC)

        !@ Velocidade média na seção de jusante (m/s)
	        VmJ = 0.0
!	        IF (QJ2(IC) > 0.0) THEN
	        IF (QJ2aux > 0.0) THEN    !@ DCB CONTINUIDADE 
!		        VmJ = QJ2(IC)/(BRIO(IC)*HmJ)
		        VmJ = QJ2aux/(BRIO(IC)*HmJ)    !@ DCB CONTINUIDADE 
            ENDIF

        !@ Raio Hidráulico médio da seção de jusante
		    RHmJ = HmJ

        !@ Velocidade de atrito na seção de jusante (m^2/s)
		    UatJ = sqrt(9.81*RHmJ*SfTREC)
	    !@ *****************************************************************************************


    !@ #########################################################################################
    !@ -------------------    CELULA COM RIO    ------------------------------------------------
    !@ #########################################################################################
	    IF(NSUBT(IC).GT.0)THEN
	
	    !@ *****************************************************************************************
	    !@ CONCENTRAÇÃO DE SEDIMENTOS NA SEÇÃO DE JUSANTE AO FINAL DO PASSO DE TEMPO(TON/M3)
	    !@ -----------------------------------------------------------------------------------------
	    !@  Equação de transporte: d(AC)/dt + d(AuC)/dx = ql
		    DO i= 1,3
		    !@ Concentração na seção de montante (ton/m3)
		        CSM2(IC,i) = 0.0
!                IF (QM2(IC) > 0.0) THEN  !@ Só tem concentração a montante se tem vazão
!                    CSM2(IC,i) = CARGM(IC,i)/(QM2(IC)*DTP)
                IF (QM2aux > 0.0) THEN  !@ DCB CONTINUIDADE 
                    CSM2(IC,i) = CARGM(IC,i)/(QM2aux*DTP)    !@ DCB CONTINUIDADE 
                ENDIF
                
             !@ Concentração na seção de jusante (ton/m3)
                CSJ2(IC,i) = 0.0
                !@ DCB CONTINUIDADE 
                IF (QJ2aux > 0.0) THEN  !@ Só tem concentração a jusante se tem vazão
                    EQ1 = PONDt*QM2aux*CSM2(IC,i)
                    EQ2 = (1.0-PONDt)*(QJ1aux*CSJ1(IC,i)-QM1aux*CSM1(IC,i))
                    EQ3 = VolTREC1(IC)*CSJ1(IC,i)/DTP
                    EQ4 = VolTREC2(IC)/DTP + PONDt*QJ2(IC)
                    CSJ2(IC,i) = (EQ1 - EQ2 + EQ3 + PSC(IC,i+1)/DTP)/EQ4    !@ TON/M3
                ENDIF
                !@ DCB CONTINUIDADE 
!                IF (QJ2(IC) > 0.0) THEN  !@ Só tem concentração a jusante se tem vazão
!                    EQ1 = PONDt*QM2(IC)*CSM2(IC,i)
!                    EQ2 = (1.0-PONDt)*(QJ1(IC)*CSJ1(IC,i)-QM1(IC)*CSM1(IC,i))
!                    EQ3 = VolTREC1(IC)*CSJ1(IC,i)/DTP
!                    EQ4 = VolTREC2(IC)/DTP + PONDt*QJ2(IC)
!                    CSJ2(IC,i) = (EQ1 - EQ2 + EQ3 + PSC(IC,i+1)/DTP)/EQ4    !@ TON/M3
!                ENDIF
            ENDDO
        !@ *****************************************************************************************     


        !@ *****************************************************************************************
        !@ CAPACIDADE DE TRANSPORTE DO TRECHO DE RIO EM TON/M3 (FORMULA DE YANG)
        !@ -----------------------------------------------------------------------------------------
            FracS(IC,:) = 0.0
            !@ DCB CONTINUIDADE CTS(IC,:) = 0.0
            !@ DCB CONTINUIDADE IF (QJ2(IC) > 0.0) THEN !@ Só tem capacidade a jusante se tem vazão
                
                IF (CSJ2(IC,1) > 0.0) THEN
    	            FracS(IC,1) = CSJ2(IC,1)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de areia
	                FracS(IC,2) = CSJ2(IC,2)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de silte
	                FracS(IC,3) = CSJ2(IC,3)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de argila
                ENDIF
                !@ Capacidade de transporte
!                CALL YANG(UatJ, VmJ, SfTREC, IC, IT)
                CTS(IC,:) = CSJ2(IC,:)  !@ DCB CONTINUIDADE 
            !@ DCB CONTINUIDADE ENDIF
        !@ *****************************************************************************************


        !@ *****************************************************************************************
        !@ VERIFICAÇÃO DE TENDÊNCIA A EROSÃO OU DEPOSIÇÃO
        !@ -----------------------------------------------------------------------------------------
            EROSaux2    = 0.0   !@ IMPRESSÃO - Erosão no rio no passo de tempo (ton/dia)
            DEPaux2     = 0.0   !@ IMPRESSÃO - Deposição no rio no passo de tempo (ton/dia)
            LIMdaux     = 0.0   !@ IMPRESSÃO - Limitador da deposição no rio no passo de tempo
            LIMeaux     = 0.0   !@ IMPRESSÃO - Limitador da erosão no rio no passo de tempo
            iENTRAaux   = 0.0   !@ IMPRESSÃO - Carga que entra no rio no passo de tempo (ton/dia)
            iSAIaux     = 0.0   !@ IMPRESSÃO - Carga que sai do rio no passo de tempo (ton/dia)
            
            DO i= 1,3   !@ INICIO DO LOOP DAS PARTICULAS
            
                DEPaux  = 0.0
                EROSaux = 0.0
                LIMd    = 0.0
		        LIMe    = 0.0
!		        IF (QJ2(IC) > 0.0) THEN   !@ Só há passagem de sedimentos se há vazão de jusante!
		        IF (QJ2aux > 0.0) THEN   !@ DCB CONTINUIDADE 
		        
		        !@ DEPOSIÇÃO ---------------------------------------------------------------------
		            IF (CSJ2(IC,i) >= CTS(IC,i)) THEN
		            !@ limitador temporal da deposição (igual ao HEC-RAS 4.0) - PARA D > 0.062 mm
		                IF (DMP(i) <= 0.000062) THEN        !@ limite máximo silte (Wu, 2008)
		                    CSJ2(IC,i) = CSJ2(IC,i)         !@ Toda carga de finos passa para jusante
!		                    CSJ2(IC,i) = CTS(IC,i)          !@ Toda carga de finos passa para jusante
		                ELSEIF (DMP(i) <= 0.000125) THEN    !@ limite máximo areia muito fina (Wu, 2008)
		                    DeT = Hmj        ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2aux*DTP   !@ DCB CONTINUIDADE *QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP - DEPaux)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE *QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ELSEIF (DMP(i) <= 0.00025) THEN     !@ limite máximo areia fina (Wu, 2008)
		                    DeT = HmJ/2.5    ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2aux*DTP   !@ DCB CONTINUIDADE *QJ2(IC)*DTP*LIMd  !@ (TON)
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP - DEPaux)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE *QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)  !@ (TON/M3)
		                ELSE
		                    DeT = HmJ/11.24  ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2aux*DTP   !@ DCB CONTINUIDADE *QJ2(IC)*DTP*LIMd  !@ (TON)
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP - DEPaux)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE *QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)  !@ (TON/M3)
		                ENDIF
                   LIMdaux(i) = LIMd
                   DEPaux2(i) = DEPaux
		        
		        !@ EROSÃO -----------------------------------------------------------------------
		            ELSE
		            !@ limitador temporal da erosão (igual ao HEC-RAS 4.0) - igual para todas as partículas
		                IF (DMP(i) <= 0.000062) THEN    !@ silte médio (Wu, 2008)
		                CSJ2(IC,i) = CSJ2(IC,i)
		                ELSE
		                    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)/(30.*HmJ)), 1.0)
		                    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*QJ2aux*DTP  !@ DCB CONTINUIDADE *QJ2(IC)*DTP*LIMe
    		                EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux
    		                CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP + EROSaux)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE *QJ2(IC)*DTP + EROSaux)/(QJ2(IC)*DTP)
		                ENDIF
                   LIMeaux(i) = LIMe
                   EROSaux2(i) = EROSaux
		            ENDIF
		        ELSE
                    DEPaux = PSC(IC,i+1) + CSM2(IC,i)*QM2aux*DTP   !@ DCB CONTINUIDADE *QM2(IC)*DTP    !@ Tudo vira deposição se não há vazão de saída!
                    DEPaux2(i) = DEPaux
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
                    CSJ2(IC,i) = 0.
		        ENDIF
        !@ *******************************************************************************


            !@ *******************************************************************************
		    !@ BALANÇO DE SEDIMENTO NO TRECHO DE RIO
		    !@ -------------------------------------------------------------------------------
		    !@ Descarga sólida na seção de jusante (ton/s)
                QSJ2(IC,i) = CSJ2(IC,i)*QJ2aux !@ DCB CONTINUIDADE *QJ2(IC)
    	        JUS=CELJUS(IC)
		        IF (JUS>0) THEN
		           CARGM(JUS,i) = CARGM(JUS,i) + QSJ2(IC,i)*DTP !@ (ton/dia)
		        ENDIF
		    !@ *******************************************************************************


            !@ armazena carga de entrada na minibacia (ton/dia)
                iENTRAaux(i)      = CARGM(IC,i) + PSC(IC,i+1) + EROSaux - DEPaux
            !@ acumula carga de entrada na minibacia (ton)
                iENTRA(iSEDaux,i) = iENTRA(iSEDaux,i) + iENTRAaux(i)
            !@ armazena carga de saída na minibacia (ton/dia)
                iSAIaux(i)        = CSJ2(IC,i)*QJ2aux*DTP  !@ DCB CONTINUIDADE *QJ2(IC)*DTP
            !@ acumula carga de saída na minibacia (ton)
                iSAI(iSEDaux,i)   = iSAI(iSEDaux,i) + iSAIaux(i)
            !@ calcula erro no intervalo de tempo
                iERRO(IC,i)       = 0.
                IF (iENTRAaux(i) > 0.0)   iERRO(IC,i) = (iENTRAaux(i) - iSAIaux(i))*100./(iENTRAaux(i))
  		
		    ENDDO   !@ FIM DO LOOP DAS PARTICULAS

            !@ Salva acumulado para cada minibacia com calculo de sedimento
            DEPT(iSEDaux,:) = DEPTREC(IC,:)
            EROT(iSEDaux,:) = EROSTREC(IC,:)
!            DEPT(iSEDaux,:) = DEPT(iSEDaux,:) + DEPTREC(IC,:)
!            EROT(iSEDaux,:) = EROT(iSEDaux,:) + EROSTREC(IC,:)


    !@ #########################################################################################
    !@ -------------------    CELULA SEM RIO    ------------------------------------------------
    !@ #########################################################################################       
        ELSE
	    
	    !@ *****************************************************************************************
	    !@ CONCENTRAÇÃO DE SEDIMENTOS NA SEÇÃO DE JUSANTE AO FINAL DO PASSO DE TEMPO(TON/M3)
	    !@ -----------------------------------------------------------------------------------------
	        DO i = 1,3
                !@ A carga da minibacia vai compor a concentração de jusante
                CSJ2(IC,i) = 0.0
!                IF (QJ2(IC) > 0.0) CSJ2(IC,i) = PSC(IC,i+1)/(QJ2(IC)*DTP)    ! ton/m3
                IF (QJ2aux > 0.0) CSJ2(IC,i) = PSC(IC,i+1)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE 
            ENDDO
	    !@ *****************************************************************************************
    
	    
	    !@ *****************************************************************************************
        !@ CAPACIDADE DE TRANSPORTE DO TRECHO DE RIO EM TON/M3 (FORMULA DE YANG)
        !-------------------------------------------------------------------------------------------
            FracS(IC,:) = 0.0
            !@ DCB CONTINUIDADE CTS(IC,:) = 0.0
            !@ DCB CONTINUIDADE IF (QJ2(IC) > 0.0) THEN !@ Só tem capacidade a jusante se tem vazão
                
                IF (CSJ2(IC,1) > 0.0) THEN
    	            FracS(IC,1) = CSJ2(IC,1)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de areia
	                FracS(IC,2) = CSJ2(IC,2)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de silte
	                FracS(IC,3) = CSJ2(IC,3)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de argila
                ENDIF
                !@ Capacidade de transporte
                !@ DCB CONTINUIDADE CALL YANG(UatJ, VmJ, SfTREC, IC, IT)
                CTS(IC,:) = CSJ2(IC,:)  !@ DCB CONTINUIDADE
            !@ DCB CONTINUIDADE ENDIF
        !@ *****************************************************************************************
        
        
        !@ *****************************************************************************************
        !@ VERIFICAÇÃO DE TENDÊNCIA A EROSÃO OU DEPOSIÇÃO
        !@ -----------------------------------------------------------------------------------------
            EROSaux2    = 0.0   !@ IMPRESSÃO
            DEPaux2     = 0.0   !@ IMPRESSÃO
            LIMdaux     = 0.0   !@ IMPRESSÃO
            LIMeaux     = 0.0   !@ IMPRESSÃO
            iENTRAaux   = 0.0   !@ IMPRESSÃO
            iSAIaux     = 0.0   !@ IMPRESSÃO
            
            DO i= 1,3   !@ INICIO DO LOOP DAS PARTICULAS
            
                DEPaux  = 0.0
                EROSaux = 0.0
                LIMd    = 0.0
		        LIMe    = 0.0
!		        IF (QJ2(IC) > 0.0) THEN   !@ Só há passagem de sedimentos se há vazão de jusante!
		        IF (QJ2aux > 0.0) THEN   !@ DCB CONTINUIDADE 
		        
		        !@ DEPOSIÇÃO ---------------------------------------------------------------------
		            IF (CSJ2(IC,i) >= CTS(IC,i)) THEN
		            !@ limitador temporal da deposição (igual ao HEC-RAS 4.0) - PARA D > 0.062 mm
		                IF (DMP(i) <= 0.000062) THEN        !@ limite máximo silte (Wu, 2008)
		                    CSJ2(IC,i) = CSJ2(IC,i)         !@ Toda carga de finos passa para jusante
!                            CSJ2(IC,i) = CTS(IC,i)          !@ Toda carga de finos passa para jusante
		                ELSEIF (DMP(i) <= 0.000125) THEN    !@ limite máximo areia muito fina (Wu, 2008)
		                    DeT = Hmj        ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2aux*DTP  !@ DCB CONTINUIDADE QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP - DEPaux)/(10.0*DTP)    !@ DCB CONTINUIDADE QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ELSEIF (DMP(i) <= 0.00025) THEN     !@ limite máximo areia fina (Wu, 2008)
		                    DeT = HmJ/2.5    ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2aux*DTP  !@ DCB CONTINUIDADE *QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP - DEPaux)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE *QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ELSE
		                    DeT = HmJ/11.24  ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2aux*DTP   !@ DCB CONTINUIDADE *QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP - DEPaux)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE *QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ENDIF
                   LIMdaux(i) = LIMd
                   DEPaux2(i) = DEPaux
		        
		        !@ EROSÃO -----------------------------------------------------------------------
		            ELSE
		            !@ limitador temporal da erosão (igual ao HEC-RAS 4.0) - igual para todas as partículas
		                IF (DMP(i) <= 0.000062) THEN    !@ silte médio (Wu, 2008)
		                CSJ2(IC,i) = CSJ2(IC,i)
		                ELSE
		                    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)/(30.*HmJ)), 1.0)
		                    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*QJ2aux*DTP  !@ DCB CONTINUIDADE *QJ2(IC)*DTP*LIMe
    		                EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux
    		                CSJ2(IC,i) = (CSJ2(IC,i)*QJ2aux*DTP + EROSaux)/(QJ2aux*DTP)    !@ DCB CONTINUIDADE *QJ2(IC)*DTP + EROSaux)/(QJ2(IC)*DTP)
		                ENDIF
                   LIMeaux(i) = LIMe
                   EROSaux2(i) = EROSaux
		            ENDIF
		        ELSE
		        !@ PODE HAVER UM ERRO AQUI, POIS SE A CONCENTRAÇÃO DE JUSANTE FOI NULA DEVIDO À VAZÃO
		        !@ NULA DE JUSANTE, NÃO APENAS O QUE ENTRA DA BACIA É DEPOSITADO, MAS O QUE VINHA DE
		        !@ MONTANTE TAMBÉM DEPOSITA NO PASSO DE TEMPO, VISTO QUE NÃO SAI NADA...
!		            DEPaux = PSC(IC,i+1)
                    DEPaux = PSC(IC,i+1) + CSM2(IC,i)*QJ2aux*DTP  !@ DCB CONTINUIDADE QM2(IC)*DTP    !@ Tudo vira deposição se não há vazão de saída!
                    DEPaux2(i) = DEPaux
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
                    CSJ2(IC,i) = 0.
		        ENDIF
        !@ *******************************************************************************    
                                   
            
            !@ *******************************************************************************
		    !@ BALANÇO DE SEDIMENTO NO TRECHO DE RIO
		    !@ -------------------------------------------------------------------------------
		    !@ Descarga sólida na seção de jusante (ton/s)
                QSJ2(IC,i) = CSJ2(IC,i)*QJ2aux    !@ DCB CONTINUIDADE QJ2(IC)
    	        JUS=CELJUS(IC)
		        IF (JUS>0) THEN
		           CARGM(JUS,i) = CARGM(JUS,i) + QSJ2(IC,i)*DTP !@ (ton/dia)
		        ENDIF
		    !@ *******************************************************************************
            
            
            !@ armazena carga de entrada na minibacia (ton/dia)
                iENTRAaux(i)      = PSC(IC,i+1) + EROSaux - DEPaux
            !@ acumula carga de entrada na minibacia (ton)
                iENTRA(iSEDaux,i) = iENTRA(iSEDaux,i) + iENTRAaux(i)
            !@ armazena carga de saída na minibacia (ton/dia)
                iSAIaux(i)        = CSJ2(IC,i)*QJ2aux*DTP  !@ DCB CONTINUIDADE *QJ2(IC)*DTP
            !@ acumula carga de saída na minibacia (ton)
                iSAI(iSEDaux,i)   = iSAI(iSEDaux,i) + iSAIaux(i)
            !@ calcula erro no intervalo de tempo
                iERRO(IC,i)       = 0.
                IF (iENTRAaux(i) > 0.0)   iERRO(IC,i) = (iENTRAaux(i) - iSAIaux(i))*100./(iENTRAaux(i))
  		
		    ENDDO   !@ FIM DO LOOP DAS PARTICULAS
		    
		    !@ Salva acumulado para cada minibacia com calculo de sedimento
            DEPT(iSEDaux,:) = DEPTREC(IC,:)
            EROT(iSEDaux,:) = EROSTREC(IC,:)
            
            
        ENDIF
        

!@ #########################################################################################
! ------------------------------------------------------------------------------------------
!@ ##########################      TRECHOS SIMULADOS COM HD       ##########################
! ------------------------------------------------------------------------------------------
!@ #########################################################################################
	ELSE
		
		!@ *****************************************************************************************
	    !@ CARACTERÍSTICAS HIDRÁULICAS DA SEÇÃO DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
	        !@ Declividade de atrito no trecho (mínima de 0.01 m/km para evitar problemas com o
	        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vazões maiores
	        !@ que 5000 m3/s)

		    SfTREC = max(DECL(IC),0.00001)
    	
	        !@ Profundidade média na seção de jusante por Manning assumindo Raio = HmTREC (m)
		    HmJ = ((RUGMAN*abs(QJ2(IC)))/(BRIO(IC)*SfTREC**0.5))**(3.0/5.0)

            !@ Volume médio de água no trecho
            VolTREC2(IC) = HmJ*BRIO(IC)*SRIO(IC)
!                if (isNaN(VolTREC2(IC))) then
!                    write(*,*) 'IT, IC, i = ', IT, IC, i
!                    write(*,*) 'HmJ, B    = ', HmJ, BRIO(IC)
!                    write(*,*) 'Rug, Sf   = ', RUGMAN, SfTREC
!                    write(*,*) 'QJ2, Srio = ', QJ2(IC), SRIO(IC)
!                    pause
!                endif
            

	        !@ Velocidade média na seção de jusante (m/s)
	        VmJ = 0.0
	        IF (QJ2(IC) > 0.0) THEN
		        VmJ = QJ2(IC)/(BRIO(IC)*HmJ)
            ENDIF

	        !@ Raio Hidráulico médio da seção de jusante
		    RHmJ = HmJ

	        !@ Velocidade de atrito na seção de jusante (m^2/s)
		    UatJ = sqrt(9.81*RHmJ*SfTREC)
	    !@ *****************************************************************************************


    !@ #########################################################################################
    !@ -------------------    CELULA COM RIO    ------------------------------------------------
    !@ #########################################################################################		
	    IF(NSUBT(IC).GT.0)THEN
	
	    !@ *****************************************************************************************
	    !@ CONCENTRAÇÃO DE SEDIMENTOS NA SEÇÃO DE JUSANTE AO FINAL DO PASSO DE TEMPO(TON/M3)
	    !@ -----------------------------------------------------------------------------------------
	    !@  Equação de transporte: d(AC)/dt + d(AuC)/dx = ql
		    DO i= 1,3
		    !@ Concentração na seção de montante (ton/m3)
		        CSM2(IC,i) = 0.0
                IF (QM2(IC) > 0.0) THEN  !@ Só tem concentração a montante se tem vazão
                    CSM2(IC,i) = CARGM(IC,i)/(QM2(IC)*DTP)
                ENDIF
                
             !@ Concentração na seção de jusante (ton/m3)
                CSJ2(IC,i) = 0.0
                IF (QJ2(IC) > 0.0) THEN  !@ Só tem concentração a jusante se tem vazão
                    EQ1 = PONDt*QM2(IC)*CSM2(IC,i)
                    EQ2 = (1.0-PONDt)*(QJ1(IC)*CSJ1(IC,i)-QM1(IC)*CSM1(IC,i))
                    EQ3 = VolTREC1(IC)*CSJ1(IC,i)/DTP
                    EQ4 = VolTREC2(IC)/DTP + PONDt*QJ2(IC)
                    CSJ2(IC,i) = (EQ1 - EQ2 + EQ3 + PSC(IC,i+1)/DTP)/EQ4    !@ TON/M3
                ENDIF
                
!!                if (isNaN(CSJ2(IC,i))) then
!                if (isNaN(VolTREC2(IC))) then
!                    write(*,*) 'IT, IC, i = ', IT, IC, i
!                    write(*,*) 'QM1, QM2  = ', QM1(IC), QM2(IC)
!                    write(*,*) 'QJ1, QJ2  = ', QJ1(IC), QJ2(IC)
!                    write(*,*) 'CM1, CM2  = ', CSM1(IC,i), CSM2(IC,i)
!                    write(*,*) 'CJ1, PSC  = ', CSJ1(IC,i), PSC(IC,i+1)
!                    write(*,*) 'V1, V2    = ', VolTREC1(IC), VolTREC2(IC)
!                    pause
!                endif

            ENDDO
        !@ *****************************************************************************************
!        write(*,*) 'IT, IC = ', IT, IC
!        write(*,*) 'CM1    = ', CSM1(IC,1), CSM1(IC,2), CSM1(IC,3)
!        write(*,*) 'CJ1    = ', CSJ1(IC,1), CSJ1(IC,2), CSJ1(IC,3)
!        write(*,*) 'CM2    = ', CSM2(IC,1), CSM2(IC,2), CSM2(IC,3)
!        write(*,*) 'CJ2    = ', CSJ2(IC,1), CSJ2(IC,2), CSJ2(IC,3)


        !@ *****************************************************************************************
        !@ CAPACIDADE DE TRANSPORTE DO TRECHO DE RIO EM TON/M3 (FORMULA DE YANG)
        !@ -----------------------------------------------------------------------------------------
            FracS(IC,:) = 0.0
            CTS(IC,:) = 0.0
            IF (QJ2(IC) > 0.0) THEN !@ Só tem capacidade a jusante se tem vazão
                
                IF (CSJ2(IC,1) > 0.0) THEN
    	            FracS(IC,1) = CSJ2(IC,1)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de areia
	                FracS(IC,2) = CSJ2(IC,2)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de silte
	                FracS(IC,3) = CSJ2(IC,3)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de argila
                ENDIF
                !@ Capacidade de transporte
                CALL YANG(UatJ, VmJ, SfTREC, IC, IT)
            ENDIF
        !@ *****************************************************************************************


        !@ *****************************************************************************************
        !@ VERIFICAÇÃO DE TENDÊNCIA A EROSÃO OU DEPOSIÇÃO
        !@ -----------------------------------------------------------------------------------------
            EROSaux2    = 0.0   !@ IMPRESSÃO - Erosão no rio no passo de tempo (ton/dia)
            DEPaux2     = 0.0   !@ IMPRESSÃO - Deposição no rio no passo de tempo (ton/dia)
            LIMdaux     = 0.0   !@ IMPRESSÃO - Limitador da deposição no rio no passo de tempo
            LIMeaux     = 0.0   !@ IMPRESSÃO - Limitador da erosão no rio no passo de tempo
            iENTRAaux   = 0.0   !@ IMPRESSÃO - Carga que entra no rio no passo de tempo (ton/dia)
            iSAIaux     = 0.0   !@ IMPRESSÃO - Carga que sai do rio no passo de tempo (ton/dia)
            
            DO i= 1,3   !@ INICIO DO LOOP DAS PARTICULAS
            
                DEPaux  = 0.0
                EROSaux = 0.0
                LIMd    = 0.0
		        LIMe    = 0.0
		        IF (QJ2(IC) > 0.0) THEN   !@ Só há passagem de sedimentos se há vazão de jusante!
		        
		        !@ DEPOSIÇÃO ---------------------------------------------------------------------
		            IF (CSJ2(IC,i) >= CTS(IC,i)) THEN
		            !@ limitador temporal da deposição (igual ao HEC-RAS 4.0) - PARA D > 0.062 mm
		                IF (DMP(i) <= 0.000062) THEN        !@ limite máximo silte (Wu, 2008)
		                    CSJ2(IC,i) = CSJ2(IC,i)         !@ Toda carga de finos passa para jusante
!		                    CSJ2(IC,i) = CTS(IC,i)          !@ Toda carga de finos passa para jusante
		                ELSEIF (DMP(i) <= 0.000125) THEN    !@ limite máximo areia muito fina (Wu, 2008)
		                    DeT = Hmj        ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ELSEIF (DMP(i) <= 0.00025) THEN     !@ limite máximo areia fina (Wu, 2008)
		                    DeT = HmJ/2.5    ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd  !@ (TON)
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)  !@ (TON/M3)
		                ELSE
		                    DeT = HmJ/11.24  ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd  !@ (TON)
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux  !@ (TON)
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)  !@ (TON/M3)
		                ENDIF
                   LIMdaux(i) = LIMd
                   DEPaux2(i) = DEPaux
		        
		        !@ EROSÃO -----------------------------------------------------------------------
		            ELSE
		            !@ limitador temporal da erosão (igual ao HEC-RAS 4.0) - igual para todas as partículas
		                IF (DMP(i) <= 0.000062) THEN    !@ silte médio (Wu, 2008)
		                CSJ2(IC,i) = CSJ2(IC,i)
		                ELSE
		                    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)/(30.*HmJ)), 1.0)
		                    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*QJ2(IC)*DTP*LIMe
    		                EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux
    		                CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP + EROSaux)/(QJ2(IC)*DTP)
		                ENDIF
                   LIMeaux(i) = LIMe
                   EROSaux2(i) = EROSaux
		            ENDIF
		        ELSE
                    DEPaux = PSC(IC,i+1) + CSM2(IC,i)*QM2(IC)*DTP    !@ Tudo vira deposição se não há vazão de saída!
                    DEPaux2(i) = DEPaux
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
                    CSJ2(IC,i) = 0.
		        ENDIF
        !@ *******************************************************************************

! descomentando este balanço dá pau, com valores NaN!!!!
            !@ *******************************************************************************
		    !@ BALANÇO DE SEDIMENTO NO TRECHO DE RIO
		    !@ -------------------------------------------------------------------------------
		    !@ Descarga sólida na seção de jusante (ton/s)
!		        write(*,*) 'IT, IC, part   = ', IT, IC, i
!		        write(*,*) 'QJ2, CJ2, QSJ2 = ', QJ2(IC), CSJ2(IC,i), QSJ2(IC,i)
                QSJ2(IC,i) = CSJ2(IC,i)*QJ2(IC)
    	        JUS=CELJUS(IC)
		        IF (JUS>0) THEN
		           CARGM(JUS,i) = CARGM(JUS,i) + QSJ2(IC,i)*DTP !@ (ton/dia)
		        ENDIF
		        
!		        if (isNaN(CARGM(JUS,i))) then
!		        write(*,*) 'IT, IC, part   = ', IT, IC, i
!		        write(*,*) 'QJ2, CJ2, QSJ2 = ', QJ2(IC), CSJ2(IC,i), QSJ2(IC,i)
!		        pause
!		        endif
		    !@ *******************************************************************************


            !@ armazena carga de entrada na minibacia (ton/dia)
                iENTRAaux(i)      = CARGM(IC,i) + PSC(IC,i+1) + EROSaux - DEPaux
            !@ acumula carga de entrada na minibacia (ton)
                iENTRA(iSEDaux,i) = iENTRA(iSEDaux,i) + iENTRAaux(i)
            !@ armazena carga de saída na minibacia (ton/dia)
                iSAIaux(i)        = CSJ2(IC,i)*QJ2(IC)*DTP
            !@ acumula carga de saída na minibacia (ton)
                iSAI(iSEDaux,i)   = iSAI(iSEDaux,i) + iSAIaux(i)
            !@ calcula erro no intervalo de tempo
                iERRO(IC,i)       = 0.
                IF (iENTRAaux(i) > 0.0)   iERRO(IC,i) = (iENTRAaux(i) - iSAIaux(i))*100./(iENTRAaux(i))
  		
		    ENDDO   !@ FIM DO LOOP DAS PARTICULAS

            !@ Salva acumulado para cada minibacia com calculo de sedimento
            DEPT(iSEDaux,:) = DEPTREC(IC,:)
            EROT(iSEDaux,:) = EROSTREC(IC,:)
!            DEPT(iSEDaux,:) = DEPT(iSEDaux,:) + DEPTREC(IC,:)
!            EROT(iSEDaux,:) = EROT(iSEDaux,:) + EROSTREC(IC,:)


    !@ #########################################################################################
    !@ -------------------    CELULA SEM RIO    ------------------------------------------------
    !@ #########################################################################################        
        ELSE
	    
	    !@ *****************************************************************************************
	    !@ CONCENTRAÇÃO DE SEDIMENTOS NA SEÇÃO DE JUSANTE AO FINAL DO PASSO DE TEMPO(TON/M3)
	    !@ -----------------------------------------------------------------------------------------
	        DO i = 1,3
                !@ A carga da minibacia vai compor a concentração de jusante
                CSJ2(IC,i) = 0.0
                IF (QJ2(IC) > 0.0) CSJ2(IC,i) = PSC(IC,i+1)/(QJ2(IC)*DTP)    ! ton/m3
            ENDDO	    
	    !@ *****************************************************************************************

	    
	    !@ *****************************************************************************************
        !@ CAPACIDADE DE TRANSPORTE DO TRECHO DE RIO EM TON/M3 (FORMULA DE YANG)
        !-------------------------------------------------------------------------------------------
            FracS(IC,:) = 0.0
            CTS(IC,:) = 0.0
            IF (QJ2(IC) > 0.0) THEN !@ Só tem capacidade a jusante se tem vazão
                
                IF (CSJ2(IC,1) > 0.0) THEN
    	            FracS(IC,1) = CSJ2(IC,1)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de areia
	                FracS(IC,2) = CSJ2(IC,2)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de silte
	                FracS(IC,3) = CSJ2(IC,3)/sum(CSJ2(IC,:)) !@ porcentagem da concentração de argila
                ENDIF
                !@ Capacidade de transporte
                CALL YANG(UatJ, VmJ, SfTREC, IC, IT)
            ENDIF
        !@ *****************************************************************************************
        
        
        !@ *****************************************************************************************
        !@ VERIFICAÇÃO DE TENDÊNCIA A EROSÃO OU DEPOSIÇÃO
        !@ -----------------------------------------------------------------------------------------
            EROSaux2    = 0.0   !@ IMPRESSÃO
            DEPaux2     = 0.0   !@ IMPRESSÃO
            LIMdaux     = 0.0   !@ IMPRESSÃO
            LIMeaux     = 0.0   !@ IMPRESSÃO
            iENTRAaux   = 0.0   !@ IMPRESSÃO
            iSAIaux     = 0.0   !@ IMPRESSÃO
            
            DO i= 1,3   !@ INICIO DO LOOP DAS PARTICULAS
            
                DEPaux  = 0.0
                EROSaux = 0.0
                LIMd    = 0.0
		        LIMe    = 0.0
		        IF (QJ2(IC) > 0.0) THEN   !@ Só há passagem de sedimentos se há vazão de jusante!
		        
		        !@ DEPOSIÇÃO ---------------------------------------------------------------------
		            IF (CSJ2(IC,i) >= CTS(IC,i)) THEN
		            !@ limitador temporal da deposição (igual ao HEC-RAS 4.0) - PARA D > 0.062 mm
		                IF (DMP(i) <= 0.000062) THEN        !@ limite máximo silte (Wu, 2008)
		                    CSJ2(IC,i) = CSJ2(IC,i)         !@ Toda carga de finos passa para jusante
!                            CSJ2(IC,i) = CTS(IC,i)          !@ Toda carga de finos passa para jusante
		                ELSEIF (DMP(i) <= 0.000125) THEN    !@ limite máximo areia muito fina (Wu, 2008)
		                    DeT = Hmj        ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ELSEIF (DMP(i) <= 0.00025) THEN     !@ limite máximo areia fina (Wu, 2008)
		                    DeT = HmJ/2.5    ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ELSE
		                    DeT = HmJ/11.24  ! Distância efetiva
		                    LIMd = min(WSP(i)*DTP/DeT, 1.0)
		                    DEPaux = (CSJ2(IC,i) - CTS(IC,i))*QJ2(IC)*DTP*LIMd
		                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
		                    CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP - DEPaux)/(QJ2(IC)*DTP)
		                ENDIF
                   LIMdaux(i) = LIMd
                   DEPaux2(i) = DEPaux
		        
		        !@ EROSÃO -----------------------------------------------------------------------
		            ELSE
		            !@ limitador temporal da erosão (igual ao HEC-RAS 4.0) - igual para todas as partículas
		                IF (DMP(i) <= 0.000062) THEN    !@ silte médio (Wu, 2008)
		                CSJ2(IC,i) = CSJ2(IC,i)
		                ELSE
		                    LIMe = min( (1.0 + exp(-1.0)) - exp(-SRIO(IC)/(30.*HmJ)), 1.0)
		                    EROSaux = (CTS(IC,i) - CSJ2(IC,i))*QJ2(IC)*DTP*LIMe
    		                EROSTREC(IC,i) = EROSTREC(IC,i) + EROSaux
    		                CSJ2(IC,i) = (CSJ2(IC,i)*QJ2(IC)*DTP + EROSaux)/(QJ2(IC)*DTP)
		                ENDIF
                   LIMeaux(i) = LIMe
                   EROSaux2(i) = EROSaux
		            ENDIF
		        ELSE
		        !@ PODE HAVER UM ERRO AQUI, POIS SE A CONCENTRAÇÃO DE JUSANTE FOI NULA DEVIDO À VAZÃO
		        !@ NULA DE JUSANTE, NÃO APENAS O QUE ENTRA DA BACIA É DEPOSITADO, MAS O QUE VINHA DE
		        !@ MONTANTE TAMBÉM DEPOSITA NO PASSO DE TEMPO, VISTO QUE NÃO SAI NADA...
!		            DEPaux = PSC(IC,i+1)
                    DEPaux = PSC(IC,i+1) + CSM2(IC,i)*QM2(IC)*DTP    !@ Tudo vira deposição se não há vazão de saída!
                    DEPaux2(i) = DEPaux
                    DEPTREC(IC,i) = DEPTREC(IC,i) + DEPaux
                    CSJ2(IC,i) = 0.
		        ENDIF
        !@ *******************************************************************************    
                                   
            
            !@ *******************************************************************************
		    !@ BALANÇO DE SEDIMENTO NO TRECHO DE RIO
		    !@ -------------------------------------------------------------------------------
		    !@ Descarga sólida na seção de jusante (ton/s)
                QSJ2(IC,i) = CSJ2(IC,i)*QJ2(IC)
    	        JUS=CELJUS(IC)
		        IF (JUS>0) THEN
		           CARGM(JUS,i) = CARGM(JUS,i) + QSJ2(IC,i)*DTP !@ (ton/dia)
		        ENDIF
		    !@ *******************************************************************************
            
            
            !@ armazena carga de entrada na minibacia (ton/dia)
                iENTRAaux(i)      = PSC(IC,i+1) + EROSaux - DEPaux
            !@ acumula carga de entrada na minibacia (ton)
                iENTRA(iSEDaux,i) = iENTRA(iSEDaux,i) + iENTRAaux(i)
            !@ armazena carga de saída na minibacia (ton/dia)
                iSAIaux(i)        = CSJ2(IC,i)*QJ2(IC)*DTP
            !@ acumula carga de saída na minibacia (ton)
                iSAI(iSEDaux,i)   = iSAI(iSEDaux,i) + iSAIaux(i)
            !@ calcula erro no intervalo de tempo
                iERRO(IC,i)       = 0.
                IF (iENTRAaux(i) > 0.0)   iERRO(IC,i) = (iENTRAaux(i) - iSAIaux(i))*100./(iENTRAaux(i))
  		
		    ENDDO   !@ FIM DO LOOP DAS PARTICULAS
		    
		    !@ Salva acumulado para cada minibacia com calculo de sedimento
            DEPT(iSEDaux,:) = DEPTREC(IC,:)
            EROT(iSEDaux,:) = EROSTREC(IC,:)
            
            
        ENDIF		
		
		
		
		
		
		
		
		iTREC = iTREC + 1	!@ Código do trecho HD da minibacia
	ENDIF

    !@ #########################################################################################


!@ *************************************
	!@ ORDENANDO MINIBACIAS DAS SUBBACIAS DE INTERESSE
    CSJ2aux(iSEDaux,1) = CSJ2(IC,1)
    CSJ2aux(iSEDaux,2) = CSJ2(IC,2)
    CSJ2aux(iSEDaux,3) = CSJ2(IC,3)
!@ *************************************


ENDDO   !@ FIM DO LOOP DAS MINIBACIAS


!@ GRAVA CARGAS DE SEDIMENTO DAS MINIBACIAS QUE CHEGAM A REDE DE DRENAGEM
WRITE(RIOCAR,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,1)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_areia.txt  (mg/L)
WRITE(RIOCSI,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,2)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_silte.txt  (mg/L)
WRITE(RIOCCL,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,3)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_argila.txt (mg/L)


DEALLOCATE(CSJ2aux)


RETURN
END SUBROUTINE