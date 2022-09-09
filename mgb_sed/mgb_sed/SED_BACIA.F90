!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA MUSLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA DETERMINA A PERDA DE SOLO PARA CADA USO DA MINIBACIA (MINIBACIA)
!@
!@	EQUAÇÃO DA MUSLE (Williams, 1975)
!@
!@	Sed = 11.8*[(Qsup*qpico*Aphru)^(0.56)]*K*C*P*LS*Rgros
!@ sendo:
!@	    Sed    = ton - perda de solo
!@	    Qsup   = mmH2O - volume de escoamento superficial
!@	    qpico  = m3/s - taxa de pico do escoamento superficial
!@	    Aphru  = km2 - área do bloco, ou área média do pixel do bloco
!@	    K      = 0.013*(ton*m2*h)/(m3*ton*cm) - fator de erodibilidade do solo
!@	    C,P,LS = adimensional - fatores da MUSLE
!@	    Rgros  = exp(-0.053*Mrocha) - fator de rocha
!@	    Mrocha = porcentagem de rocha na primeira camada de solo
!@
!@-------------------------------------------------------------------------------------
!@  VARIÁVEIS
!@  TKS     = retardos do escoamento superficial (s)
!@  DSUP    = lâmina do escoamento superficial (mm)
!@  VSUP    = volume do escoamento supeficial na minibacia (m3)
!@  sbtFLAG = flag para indicar minibacia com substituição de dados
!@  IB,IBAC = código da subbacia da minibacia
!@  IU,NU   = usos das minibacias
!@  NPXU    = número de pixels de cada uso em cada minibacia
!@  PUSO    = porcentagem de uso em cada minibacia
!@  SLC     = perda de solo da minibacia (ton)
!@  SLU     = perda de solo de cada bloco da minibacia (ton)
!@  SLX     = carga e sedimentos remanescente na minibacia no passo de tempo
!@  SLXU    = carga e sedimentos remanescente por bloco da minibacia no passo de tempo
!@  PSCaux  = aporte das frações de sedimentos que chegam à drenagem das minibacias das subbacias selecionadas
!@-------------------------------------------------------------------------------------
!@
!@-------------------------------------------------------------------------------------
!@  VARIÁVEIS LOCAIS
!@  QSUM    = lâmina total do escoamento superficial na minibacia (mm)
!@  QSAUX   = armazena TKS e VSUP da minibacia e DSUP de cada bloco da minibacia
!@  QSUPaux = vazão superficial que chega à drenagem apenas das minibacias das subbacias selecionadas 2013_02_06 (antes QCEaux)
!@  QCELaux = vazão que chega à drenagem apenas das minibacias das subbacias selecionadas 2013_02_06
!@  FRAC    = porcentagem de sedimento em cada bloco da minibacia
!@-------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

Subroutine SED_BACIA

USE SED_VARS
USE VARS_MAIN
USE aux_mod

IMPLICIT NONE

!@ RETARDOS DO ESCOAMENTO SUPERFICIAL (s)
REAL TKS
!@ LAMINAS DO ESCOAMENTO SUPERFICIAL (mm)
REAL DSUP, QSUM

REAL QSUPaux(NC),FRAC(NU),QPICaux(NC)
REAL QCELaux(NC) !@ 2013_02_06
INTEGER iaux, iUSO,i


ALLOCATE(PSCaux(NC,4))
PSC = 0.0
PSCaux = 0.0

        !@ QSAUX armazena
            !@      TKS da minibacia IC
            !@	    VSUP da minibacia IC
            !@	    DSUP bloco 1 da minibacia IC
            !@	    DSUP bloco 2 da minibacia IC
            !@	    DSUP bloco 3 da minibacia IC
            !@	    	*
            !@	    	*
            !@	    	*
            !@	    DSUP bloco N da minibacia IC
!----------------------
iaux = 0    !@ contador de minibacia com cálculo de sedimentos
!----------------------


!VERIFICA A DATA CORRESPONDE O DIA JULIANO
CALL CALDAT(IDINI + IT - 1,IMES,IDIA,IANO)

QPICaux = 0.0
DO IC=1,NC	!@ LOOP DAS MINIBACIAS
	
	IF (sbtFLAG(iC)==1) CYCLE !@ Nao computa minibacias a montante de minibacias com substituicao de dados
	!@ Neste caso, para sedimentos, ainda terá que ser fornecida uma concentração de sedimentos

    IB=IBAC(IC)
    
!@ DCB 30/04/1011 ###########################################################################
	!IF(IB>SUBfim.OR.IB<SUBini)CYCLE
    !IF(IB>82 .OR. IB<57)CYCLE
!@ DCB 30/04/1011 ###########################################################################

	iaux = iaux + 1 !@ contador de minibacia com cálculo de sedimentos

	QSUM = 0.0  !@ INICIALIZA LÂMINA TOTAL DO ESCOAMENTO SUPERFICIAL NA MINIBACIA (mm)
	
	DO IU=1,NU	!@ LOOP DOS USOS

		APIX = 0.0
		FCTE = 0.0
		QPIC = 0.0

		IF (NPXU(IC,IU) .LE. 0.0) THEN  !@ NAO TEM PIXEL(S) COM ESTE USO NESTA MINIBACIA
			CYCLE !@ PASSA PARA O PROXIMO USO
		ENDIF

		IF(PUSO(IC,IU).LT.0.0001)THEN !@ NAO TEM ESTE USO NESTA MINIBACIA
			CYCLE !@ PASSA PARA O PROXIMO USO
		ENDIF
		

		DSUP = QSAUX(IC,IU+2)
		QSUM = QSUM + DSUP !@ LÂMINA TOTAL DO ESCOAMENTO SUPERFICIAL NA MINIBACIA (mm)

		
		!@ ÁREA MÉDIA DOS PIXELS DE CADA USO DA MINIBACIA (KM^2)
		APIX = (ACEL(IC)*(PUSO(IC,IU)/100.))/NPXU(IC,IU) !@ (KM^2)
		APIX = APIX*100. !@ (Ha = hectares)


		!@ FATOR CONSTANTE DA MUSLE PARA CADA USO DA MINIBACIA
		FCTE = 11.8*(APIX**0.56)*Kusle(IB,IU)*Cusle(IB,IU)*Pusle(IB,IU)*Rgros(IB,IU)*LSAcu(IC,IU)
		APIX = APIX/100. !@ (KM^2)


		IF (DSUP < 0.0) DSUP = 0.0 !@ CORREÇÃO DE DSUP NEGATIVO VINDO DA ROTINA "CELULA"

		!@ TAXA DE PICO DO ESCOAMENTO SUPERFICIAL PELA MÉTODO RACIONAL, CONSIDERANDO
        !@ A INTENSIDADE DA PRECIPITAÇÃO DENTRO DAS 24 HORAS (DIA) (M^3/S)
        IF (P(IC) > 0.0) QPIC = ((DSUP/P(IC))*(P(IC)/24.)*APIX)/3.6 !@ (M^3/S)
!		IF (P(IC) > 0.0) QPIC = ((DSUP/P(IC))*(P(IC)/(TIND(IC)/3600.))*APIX)/3.6 !@ (M^3/S)
		QPICaux(iaux) = QPICaux(iaux) + QPIC !@ acumula taxa de pico, apenas para gravação
!		write(*,*) 'IC, P, Tc = ', IC, P(IC), TIND(IC)/3600.
!		pause
		
		!@ PERDA DE SOLO DA MINIBACIA (TON)
		!SLC(IC) = SLC(IC) + FCTE*(DSUP**0.56)*(QPIC**0.56) !@ (TON)
		
		
		!@ PERDA DE SOLO PARA CADA USO DA MINIBACIA (TON)
		!@ VERIFICAR: Se DSUP precisa ser dividido pela área da bacia
		!@ RESPOSTA: Acho que não, pois o mm é pela área do bloco!!
		SLU(IC,IU) = SLU(IC,IU) + FCTE*(DSUP**0.56)*(QPIC**0.56) !@ (TON)
!        SLU(IC,IU) = SLU(IC,IU) + FCTE
!        SLU(IC,IU) = SLU(IC,IU) + FCTE*(DSUP**0.56)
!        SLU(IC,IU) = SLU(IC,IU) + FCTE*(QPIC**0.56) !@ (TON)
		
	ENDDO	! FIM DO LOOP DOS USOS


!    !@ MULTIPLICAR O SLU PELO SDR, ANTES DE FAZER A PROPAGAÇÃO NA REDE. O SDR AQUI É POR MINIBACIA, ENTÃO
!    !@ É CONSTANTE NOS BLOCOS. COM A OPÇÃO DE VARIAR POR BLOCO É SÓ MULTIPLICAR O SDRi DO BLODO i PELO SLUi DO
!    !@ BLOCO.
!    DO IU = 1,NU
!        !@ O SDR com o parametro de calibração não pode ser maior que 1!
!        SLU(IC,IU) = SLU(IC,IU)*min(Ksdr(IB,IU)*SDR(IC,NU+1),1.)
!    ENDDO


    FRAC = 0.0  !@ porcentagem de sedimento de cada bloco da minibacia
    if (sum(SLU(IC,:)) /= 0.0) then
        FRAC(:) = SLU(IC,:)/sum(SLU(IC,:))
    endif


	!@ PARÂMETRO DE RETARDO DO ESCOAMENTO SUPERFICIAL
	TKS = QSAUX(IC,1)


	!@ DESCARGA SÓLIDA DAS MINIBACIAS (TON/S) NO TEMPO IT
	!QSC(IC)   = SLC(IC)/TKS
	QSU(IC,:) = SLU(IC,:)/TKS	!@ 21/02/2011


	!@ ATUALIZA PERDA DE SOLO DAS MINIBACIAS (TON) NO TEMPO IT
	!SLX     = SLC(IC)   - QSC(IC)*DTP !@ verifica o que sobra de sedimentos na minibacia
	SLXU(:) = SLU(IC,:) - QSU(IC,:)*DTP !@ 21/02/2011


!@ **********************************************************************
!@ *** 21/02/2011 *******************************************************
!@ **********************************************************************
	DO IU=1,NU
		IF(SLXU(IU) < 0.0)THEN				!@ TKS < DTP
			QSU(IC,IU) = SLU(IC,IU)/(DTP)
			SLU(IC,IU) = 0.0
		ELSE
			SLU(IC,IU) = SLXU(IU)
		ENDIF
	ENDDO
!@ **********************************************************************
!@ *** 21/02/2011 *******************************************************
!@ **********************************************************************



    !@ APORTE DE SEDIMENTOS DAS BLOCOS DAS MINIBACIAS PARA O RIO (TON)
	PSU(IC,:) = QSU(IC,:)*DTP	!@ 21/02/2011



	!@ APORTE DE SEDIMENTOS DAS MINIBACIAS PARA O RIO (TON)
	PSC(IC,1) = sum(QSU(IC,:))*DTP
	DO iUSO=1,NU-1
	    PSC(IC,2) = PSC(IC,2) + PSC(IC,1)*FRAC(iUSO)*Mareia(IB,iUSO)/100	!@ APORTE DE AREIA QUE CHEGA À DRENAGEM
	    PSC(IC,3) = PSC(IC,3) + PSC(IC,1)*FRAC(iUSO)*Msilte(IB,iUSO)/100	!@ APORTE DE SILTE QUE CHEGA À DRENAGEM
	    PSC(IC,4) = PSC(IC,4) + PSC(IC,1)*FRAC(iUSO)*Margila(IB,iUSO)/100	!@ APORTE DE ARGILA QUE CHEGA À DRENAGEM
    ENDDO
    !PSC(IC,3) = PSC(IC,3)*1.5
    !PSC(IC,4) = PSC(IC,4)*1.5

!@ VERIFICA O QUE OCORRE SE TUDO POR CONSIDERADO EM SUSPENSÃO
!PSC(IC,3) = PSC(IC,3) + PSC(IC,2)/2.
!PSC(IC,4) = PSC(IC,4) + PSC(IC,2)/2.
!PSC(IC,2) = 0.
!@ *************************************
	!@ ORDENANDO MINIBACIAS DAS SUBBACIAS DE INTERESSE
	PSCaux(iaux,1) = PSC(IC,1)	!@ ARMAZENA O APORTE DE SEDIMENTOS QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	QSUPaux(iaux)   = QSUP(IC)	!@ ARMAZENA A VAZÃO SUPERFICIAL DA MINIBACIA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	QCELaux(iaux)   = QCEL2(IC)	!@ ARMAZENA A VAZÃO DA MINIBACIA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS 2013_02_06
	
	DO iUSO=1,NU-1
	    PSCaux(iaux,2) = PSCaux(iaux,2) + PSCaux(iaux,1)*FRAC(iUSO)*Mareia(IB,iUSO)/100     !@ ARMAZENA O APORTE DE AREIA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	    PSCaux(iaux,3) = PSCaux(iaux,3) + PSCaux(iaux,1)*FRAC(iUSO)*Msilte(IB,iUSO)/100     !@ ARMAZENA O APORTE DE SILTE QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
	    PSCaux(iaux,4) = PSCaux(iaux,4) + PSCaux(iaux,1)*FRAC(iUSO)*Margila(IB,iUSO)/100    !@ ARMAZENA O APORTE DE ARGILA QUE CHEGA À DRENAGEM APENAS NAS MINIBACIAS DAS SUBBACIAS SELECIONADAS
    ENDDO
!@ *************************************

ENDDO	!@ FIM DO LOOP DAS MINIBACIAS


!@ GRAVA CARGAS DE SEDIMENTO DAS MINIBACIAS QUE CHEGAM A REDE DE DRENAGEM
WRITE(FILSED, '(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,1),IC=1,iaux)  !@ SED_MINI.txt - total (TON/dia)
WRITE(FILSEDa,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,2),IC=1,iaux)  !@ SEDa_MINI.txt - silte (TON/dia)
WRITE(FILSEDs,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,3),IC=1,iaux)  !@ SEDs_MINI.txt - argila (TON/dia)
WRITE(FILSEDc,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (PSCaux(IC,4),IC=1,iaux)  !@ SEDc_MINI.txt - areia (TON/dia)

!@ GRAVA A TAXA DE PICO ACUMULADA NAS MINIBACIAS
!WRITE(TAXPIC,'(4I10,<iaux>F15.3)')  IT, IDIA, IMES, IANO, (QPICaux(IC),IC=1,iaux)  !@ TAXpico_MINI.txt - acumulado (m3/s)

!@ GRAVA VAZÃO SUPERFICIAL DAS MINIBACIAS QUE CHEGA A REDE DE DRENAGEM
WRITE(FILQCE,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (QSUPaux(IC),IC=1,iaux) !@ QSUP_MINI.txt (m3/s)

!@ GRAVA VAZÃO DAS MINIBACIAS QUE CHEGA A REDE DE DRENAGEM 2013_02_06
WRITE(FILCEL,'(4I10,<iaux>F15.3)') IT, IDIA, IMES, IANO, (QCELaux(IC),IC=1,iaux) !@ QCEL_MINI.txt (m3/s) 2013_02_06


DEALLOCATE(PSCaux)



RETURN

END