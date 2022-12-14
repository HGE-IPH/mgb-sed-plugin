!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abril de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SED_REDE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - SUBROTINA PARA PREPARA??O DAS VARI?VEIS DO ESCOAMENTO E SEDIMENTOS PARA PROPAGA??O DA CERGA
!@   DE SEDIMENTOS AO LONGO DA DEDE DE DRENAGEM
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
!@ Declara??o de vari?veis:
REAL HmJx, VmJx, RHmJx, SfTRECx, UatJx, Qlimx
INTEGER iTREC
REAL iENTRAaux(3), iSAIaux(3)
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

!@ INICIALIZANDO VARI?VEIS
iTREC   = 0     !@ C?digo do trecho HD da minibacia
iSEDaux = 0     !@ Contador de minibacia com c?lculo de sedimentos
CARGM   = 0.0   !@ Carga a montante do trecho
CSJ2    = 0.0
Qlimx   = 0.01

DO IC = 1,NC    !@ INICIO DO LOOP DAS MINIBACIAS
    IB=IBAC(IC)
    
    !@ CALCULA APENAS SUB-BACIAS SELECIONADAS
!	IF(IB>SUBfim .OR. IB<SUBini)CYCLE
	IF(IB>82 .OR. IB<57)CYCLE

	iSEDaux = iSEDaux + 1 !@ Contador de minibacia com c?lculo de sedimentos


    !@ #########################################################################################
    ! ------------------------------------------------------------------------------------------
    !@ ##########################      TRECHOS SIMULADOS COM MC       ##########################
    ! ------------------------------------------------------------------------------------------
    !@ #########################################################################################
	IF (hdFLAG(IC) == 0) THEN

        !@ *****************************************************************************************
	    !@ CARACTER?STICAS HIDR?ULICAS DA SE??O DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
        !@ Declividade de atrito no trecho (m?nima de 0.01 m/km para evitar problemas com o
        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vaz?es maiores
        !@ que 5000 m3/s)
	    SfTRECx = max(DECL(IC),0.00001) 
    	
        !@ Profundidade m?dia na se??o de jusante por Manning assumindo Raio = HmTREC (m)
        HmJx = ((RUGMAN*QJ2(IC))/(BRIO(IC)*SfTRECx**0.5))**(3.0/5.0)

        !@ Volume m?dio de ?gua no trecho
        VolTREC2(IC) = HmJx*BRIO(IC)*SRIO(IC)*1000.

        !@ Velocidade m?dia na se??o de jusante (m/s)
        VmJx = 0.0
        IF (QJ2(IC) > Qlimx) VmJx = QJ2(IC)/(BRIO(IC)*HmJx)

        !@ Raio Hidr?ulico m?dio da se??o de jusante
        RHmJx = HmJx

        !@ Velocidade de atrito na se??o de jusante (m^2/s)
        UatJx = sqrt(9.81*RHmJx*SfTRECx)
	    !@ *****************************************************************************************

        CALL SED_PROPAG(HmJx, VmJx, RHmJx, SfTRECx, UatJx, Qlimx)
      

    !@ #########################################################################################
    ! ------------------------------------------------------------------------------------------
    !@ ##########################      TRECHOS SIMULADOS COM HD       ##########################
    ! ------------------------------------------------------------------------------------------
    !@ #########################################################################################
	ELSE

        !@ *****************************************************************************************
	    !@ CARACTER?STICAS HIDR?ULICAS DA SE??O DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
        !@ Declividade de atrito no trecho (m?nima de 0.01 m/km para evitar problemas com o
        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vaz?es maiores
        !@ que 5000 m3/s)
	    SfTRECx = max(DECL(IC),0.00001)
    	
        !@ Profundidade m?dia na se??o de jusante por Manning assumindo Raio = HmTREC (m)
        HmJx = ((RUGMAN*abs(QJ2(IC)))/(BRIO(IC)*SfTRECx**0.5))**(3.0/5.0)

        !@ Volume m?dio de ?gua no trecho
        VolTREC2(IC) = HmJx*BRIO(IC)*SRIO(IC)*1000.

        !@ Velocidade m?dia na se??o de jusante (m/s)
        VmJx = 0.0
        IF (QJ2(IC) > Qlimx) VmJx = QJ2(IC)/(BRIO(IC)*HmJx)

        !@ Raio Hidr?ulico m?dio da se??o de jusante
        RHmJx = HmJx

        !@ Velocidade de atrito na se??o de jusante (m^2/s)
        UatJx = sqrt(9.81*RHmJx*SfTRECx)
	    !@ *****************************************************************************************

        CALL SED_PROPAG(HmJx, VmJx, RHmJx, SfTRECx, UatJx, Qlimx)
        
        
		iTREC = iTREC + 1	!@ C?digo do trecho HD da minibacia
		
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