!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abril de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SED_REDE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - SUBROTINA PARA PREPARAÇÃO DAS VARIÁVEIS DO ESCOAMENTO E SEDIMENTOS PARA PROPAGAÇÃO DA CERGA
!@   DE SEDIMENTOS AO LONGO DA DEDE DE DRENAGEM
!@
!@ *********************************************************************************************

SUBROUTINE SED_REDE

!@ **********************************
USE AUX_MOD
USE VARS_MAIN
USE SED_VARS
USE PAR1_MOD    !@ DCB_HD_Sed
USE TIME_MOD    !@ DCB_HD_Sed
USE PP_MOD      !@ DCB_HD_Sed
!@ **********************************

IMPLICIT NONE

!@ **********************************
!@ Declaração de variáveis:
REAL HmJx, VmJx, RHmJx, SfTRECx, UatJx, Qlimx
REAL iENTRAaux(3), iSAIaux(3)
INTEGER SEDTR, i, iX1, iX2      !@ DCB_HD_Sed
REAL QtFL, AtFL, H1             !@ DCB_HD_Sed
REAL FINT, HFL                  !@ DCB_HD_Sed
REAL VFL_SED                    !@ DCB_HD_Sed
!REAL QTROCaux(nSEDmini), QTROC(NC)    !@ DCB_HD_Sed

!ALLOCATE(CSJ2aux(nSEDmini,3))
REAL CSJ2aux(nSEDmini,3)
!@ **********************************


!@ *****************************************************************************************
!@ EM CADA INTERVALO DE TEMPO O QUE ERA i+1 VIRA i
!@ -----------------------------------------------------------------------------------------
DO IC=1,NC
	CSM1(IC,:)      = CSM2(IC,:)
	CSJ1(IC,:)      = CSJ2(IC,:)
	VolTREC1(IC)    = VolTREC2(IC)
ENDDO
HTR1 = HTR2 !@ DCB_HD_Sed
HTR2 = 0.0  !@ DCB_HD_Sed
!@ *****************************************************************************************

!@ INICIALIZANDO VARIÁVEIS
iSEDaux = 0     !@ Contador de minibacia com cálculo de sedimentos
CARGM   = 0.0   !@ Carga a montante do trecho
CSJ2    = 0.0
Qlimx   = 0.01
!QTROCaux    = 0.0   !@ DCB_HD_Sed
!QTROC       = 0.0   !@ DCB_HD_Sed

DO IC = 1,NC    !@ INICIO DO LOOP DAS MINIBACIAS
    IB=IBAC(IC)
    
    !@ CALCULA APENAS SUB-BACIAS SELECIONADAS
!	IF(IB>SUBfim .OR. IB<SUBini)CYCLE
	IF(IB>82 .OR. IB<57)CYCLE

	iSEDaux = iSEDaux + 1 !@ Contador de minibacia com cálculo de sedimentos

    VFL_SED = 0.0   !@ DCB_HD_Sed Volume total de água na planície da minibacia
    H1      = 0.0   !@ DCB_HD_Sed Profundidade média no trecho de rio

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
	    SfTRECx = max(DECL(IC),0.00001) 
    	
        !@ Profundidade média na seção de jusante por Manning assumindo Raio = HmTREC (m)
        HmJx = ((RUGMAN*QJ2(IC))/(BRIO(IC)*SfTRECx**0.5))**(3.0/5.0)

        !@ Volume médio de água no trecho
        VolTREC2(IC) = HmJx*BRIO(IC)*SRIO(IC)*1000.

        !@ Velocidade média na seção de jusante (m/s)
        VmJx = 0.0
        IF (QJ2(IC) > Qlimx) VmJx = QJ2(IC)/(BRIO(IC)*HmJx)

        !@ Raio Hidráulico médio da seção de jusante
        RHmJx = HmJx

        !@ Velocidade de atrito na seção de jusante (m^2/s)
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
	    !@ CARACTERÍSTICAS HIDRÁULICAS DA SEÇÃO DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
        !@ Declividade de atrito no trecho (mínima de 0.01 m/km para evitar problemas com o
        !@ calculo da CT por YANG, pois valores menores geram CT=0 para vazões maiores
        !@ que 5000 m3/s)
        SEDTR = CellTr(IC)                  !@ Código do trecho de minibacia
        iX1  = TrNST(SEDTR,1)               !@ Código da seção de montante do trecho
        iX2  = TrNST(SEDTR,2)               !@ Código da seção de jusante do trecho
	    SfTRECx = max(SFF(iX2),0.00001)     !@ Declividade da linha de energia na seção de jusante
        
        !@ Profundidade d'água na seção de jusante obtida pelo HD
        HmJx = HO(iX2)

        !@ Volume médio de água no trecho
        VolTREC2(IC) = HmJx*BA(1,iX2)*TrL(SEDTR)*1000.

        !@ Velocidade média na seção de jusante (m/s)
        VmJx = 0.0
        IF (QJ2(IC) > Qlimx) VmJx = QJ2(IC)/(BA(1,iX2)*HmJx)

        !@ Raio Hidráulico médio da seção de jusante
        RHmJx = 2.0*HmJx + BA(1,iX2)

        !@ Velocidade de atrito na seção de jusante (m^2/s)
        UatJx = sqrt(9.81*RHmJx*SfTRECx)
        !@ *****************************************************************************************

        !@ Calculo da vazao de troca rio-planicie:
		QtFL = 0.0   !@ Vazão de troca total na minibacia
		AtFL = 0.0   !@ Área total da planície de inundação na minibacia
        HTR2 = HO
		!@ Soma das vazões de troca rio-planície e áreas alagadas para cada secao do trecho:
		DO i = iX1+1,iX2
			!@ Largura alagada: AFI(iX)
    		!@ Vazao de troca rio planicie em um subtrecho
			QtFL = QtFL + ( AFI(i)*TrL(SEDTR)*1000./(nXSec(SEDTR)-1) )*( HTR2(i) - HTR1(i) + HTR2(i-1) - HTR1(i-1) )/(2.*DTP) !@ (m3/s)
			AtFL = AtFL +   AFI(i)*TrL(SEDTR)*1000./(nXSec(SEDTR)-1)   !@ (m2)
		ENDDO
        
        !@ Calculo do volume na planicie:
        !@ Soma o volume de cada subtrecho:
		DO i = iX1+1,iX2
            ! Calcula volume na planicie de inundação e no rio:
    		! Nivel d'água:
		    H1 = 0.5*(HO(i) + HO(i-1))
		    IF (LD/=0) H1 = H1 + ZO(i)
		    ! Interpola volume:	
		    IF (Zfl(1,i)<=H1 .and. nPfl(i)>0) THEN
			    ! Interpola area alagada:
			    ! Obs.: Considera minimo entre volume interpolado e volume no nivel máximo para evitar erro na extrapolacao.
			    VFL_SED = VFL_SED + min(FINT(Zfl(1:nPfl(i),i),Vfl(1:nPfl(i),i),nPfl(i),H1),Vfl(nPfl(i),i))
		    ENDIF
        ENDDO

!        QTROC(IC) = QtFL
        
        !@ Cálculo da profundidade média da planície
        HFL = VFL_SED/AtFL

        CALL SED_PROPAG(HmJx, VmJx, RHmJx, SfTRECx, UatJx, Qlimx)
		
	ENDIF
    !@ #########################################################################################


    !@ *************************************
	!@ ORDENANDO MINIBACIAS DAS SUBBACIAS DE INTERESSE
    CSJ2aux(iSEDaux,1) = CSJ2(IC,1)
    CSJ2aux(iSEDaux,2) = CSJ2(IC,2)
    CSJ2aux(iSEDaux,3) = CSJ2(IC,3)
!    QTROCaux(iSEDaux)  = QTROC(IC)
    !@ *************************************


ENDDO   !@ FIM DO LOOP DAS MINIBACIAS

!@ GRAVA CARGAS DE SEDIMENTO DAS MINIBACIAS QUE CHEGAM A REDE DE DRENAGEM
WRITE(RIOCAR,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,1)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_areia.txt  (mg/L)
WRITE(RIOCSI,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,2)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_silte.txt  (mg/L)
WRITE(RIOCCL,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (CSJ2aux(IC,3)*(10.**(6.)),IC=1,iSEDaux)   !@ CONC_RIO_argila.txt (mg/L)
!WRITE(VAZFLP,'(3I10,<iSEDaux>F15.3)')  IDIA, IMES, IANO, (QTROCaux(IC),IC=1,iSEDaux)  !@ DCB_HD_Sed

!DEALLOCATE(CSJ2aux)


RETURN
END SUBROUTINE