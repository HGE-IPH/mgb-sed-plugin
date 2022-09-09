!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abr de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SED_INICIAL >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - SUBROTINA PARA GERAR CONDIÇÕES INICIAIS DE SEDIMENTOS NOS TRECHOS DE RIO
!@
!@ OBS.: Verificar se pode ser zero!!!
!@
!@-------------------------------------------------------------------------------------
!@ Descrição das variáveis locais:
!@
!@--------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

SUBROUTINE SED_INICIAL

!@ **********************************
USE AUX_MOD
USE VARS_MAIN
USE SED_VARS
!@ **********************************

IMPLICIT NONE

!@ **********************************
!@ Declaração de variáveis:
REAL HmJ, VmJ, RHmJ, SfTREC, UatJ
INTEGER iTREC, i
INTEGER CONTICaux(NC)

!@ **********************************
!@ Alocação de variáveis:
ALLOCATE(WSP(3), CTS(NC,3), Fform(3), DMP(3))
ALLOCATE(CSM1(NC,3), CSM2(NC,3), CSJ1(NC,3), CSJ2(NC,3), VolTREC1(NC), VolTREC2(NC))
ALLOCATE(QSJ2(NC,3))
ALLOCATE(FracS(NC,3))
ALLOCATE(DEPTREC(NC,3), DEPT(NC,3), EROSTREC(NC,3), EROT(NC,3), CARGM(NC,3))
ALLOCATE(iENTRA(NC,3), iSAI(NC,3),iERRO(NC,3))
ALLOCATE(CellTr(NC)) !@ DCB_HD_Sed


!@ ***************************************************************************************
!@ DEFININDO PARÂMETROS
DMP(1) = 0.0008		!@ Diâmetro característico para as partículas de areia (m)
DMP(2) = 0.000016	!@ Diâmetro característico para as partículas de silte (m)
DMP(3) = 0.000001	!@ Diâmetro característico para as partículas de argila (m)

Fform(1) = 0.7	    !@ fator de forma para as partículas de areia
Fform(2) = 0.7	    !@ fator de forma para as partículas de silte
Fform(3) = 0.7	    !@ fator de forma para as partículas de argila

VISC   = 10.**(-6)  !@ viscosidade cinemática da água
RUGMAN = 0.030	    !@ constante para todo trecho, como em PARCUNGE
!@ ***************************************************************************************


!@ ***************************************************************************************
!@ INICIALIZANDO VARIÁVEIS
iTREC    = 0
VolTREC1 = 0.
VolTREC2 = 0.
DEPTREC  = 0. !@ carga acumulada de sedimento depositado no trecho (TON)
DEPT     = 0.
EROSTREC = 0. !@ carga acumulada de sedimento erodido no trecho (TON)
EROT     = 0.
iENTRA   = 0. !@ acumulado da carga que entra nas minibacias (TON/dia)
iSAI     = 0. !@ acumulado da carga que sai das minibacias (TON/dia)
iERRO    = 0. !@ erro entre o que entra e o que sai das minibacias (%)
!@ ***************************************************************************************


!@ ***************************************************************************************
!@ VELOCIDADE DE QUEDA DAS PARTÍCULAS (Wu e Wang, 2006) - M/S
CALL FALLVELOCITY
!@ ***************************************************************************************

FracS = 1./3. !@ Inicialmente considera igual porcentagem de partículas no trecho de rio

iSEDaux = 0

!@ DETERMINA TRECHOS DAS MINIBACIAS !@ DCB_HD_Sed
CellTr = 0.0                        !@ DCB_HD_Sed
DO IC = 1, nTr                      !@ DCB_HD_Sed
    CellTr(TrCell(IC)) = IC         !@ DCB_HD_Sed
ENDDO                               !@ DCB_HD_Sed

DO IC = 1,NC

    IB=IBAC(IC)
    
    !@ DCB 30/04/1011 ##############################################################################
!	IF(IB>SUBfim .OR. IB<SUBini)CYCLE !@ Para calcular apenas sub-bacia de interesse.
	IF(IB>82 .OR. IB<57)CYCLE
    !@ DCB 30/04/1011 ##############################################################################


    !@ DCB sed 20_09_2012 ##############################################################################
    if (IB>56.AND.IB<83) then
	iSEDaux = iSEDaux + 1 !@ contador de minibacia com cálculo de sedimentos
	CONTICaux(iSEDaux) = IC !@ armazena temporariamente os codigos das minibacias simuladas
	endif
    !@ DCB sed 20_09_2012 ##############################################################################


    ! ------------------------------------------------------------------------------------------
    !@ ##########################      TRECHOS SIMULADOS COM MC       ##########################
    ! ------------------------------------------------------------------------------------------
	IF (hdFLAG(IC) == 0) THEN

        !@ *****************************************************************************************
	    !@ CARACTERÍSTICAS HIDRÁULICAS DA SEÇÃO DE JUSANTE
	    !@ -----------------------------------------------------------------------------------------
	    !@ Declividade de atrito no trecho (mínima de 0.01 m/km para evitar problemas com o
	    !@ calculo da CT por YANG, pois valores menores geram CT=0 para vazões maiores
	    !@ que 5000 m3/s)
		SfTREC = max(DECL(IC),0.00001)
		
	    !@ Profundidade média na seção de jusante por Manning assumindo Raio = HmTREC (m)
		HmJ = ((RUGMAN*QJ2(IC))/(BRIO(IC)*SfTREC**0.5))**(3.0/5.0)

        !@ Volume médio de água no trecho
        VolTREC2(IC) = HmJ*BRIO(IC)*SRIO(IC)*1000.

	    !@ Velocidade média na seção de jusante (m/s)
		VmJ = 0.0
	    IF (QJ2(IC) > 0.01) VmJ = QJ2(IC)/(BRIO(IC)*HmJ)

	    !@ Raio Hidráulico médio da seção de jusante
		RHmJ = HmJ

	    !@ Velocidade de atrito na seção de jusante (m^2/s)
		UatJ = sqrt(9.81*RHmJ*SfTREC)
        !@ ***************************************************************************************

        !@ *****************************************************************************************
        !@ CAPACIDADE DE TRANSPORTE DO TRECHO DE RIO EM TON/M3 (FORMULA DE YANG) 
!        CALL YANG(UatJ, VmJ, SfTREC, IC, IT)
        !@ *****************************************************************************************
            
        !@ *****************************************************************************************
        !@ CONCENTRAÇÃO INICIAL DE SEDIMENTOS NA SEÇÃO DE JUSANTE - 50% da CTS 
        CSJ2(IC,:) = 0.             !(TON/M3)
        CSJ1 = CSJ2                 !(TON/M3)
        CSM2 = CSJ2                 !(TON/M3)
        CSM1 = CSM2                 !(TON/M3)
        !@ *****************************************************************************************
	    FracS(IC,1) = 0.    !@ porcentagem da concentração de areia
	    FracS(IC,2) = 0.    !@ porcentagem da concentração de silte
	    FracS(IC,3) = 0.    !@ porcentagem da concentração de argila

    ! ------------------------------------------------------------------------------------------
    !@ ##########################      TRECHOS SIMULADOS COM HD       ##########################
    ! ------------------------------------------------------------------------------------------
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
        VolTREC2(IC) = HmJ*BRIO(IC)*SRIO(IC)*1000.

	    !@ Velocidade média na seção de jusante (m/s)
		VmJ = 0.0
	    IF (QJ2(IC) > 0.01) VmJ = QJ2(IC)/(BRIO(IC)*HmJ)

	    !@ Raio Hidráulico médio da seção de jusante
		RHmJ = HmJ

	    !@ Velocidade de atrito na seção de jusante (m^2/s)
		UatJ = sqrt(9.81*RHmJ*SfTREC)
        !@ ***************************************************************************************

        !@ *****************************************************************************************
        !@ CAPACIDADE DE TRANSPORTE DO TRECHO DE RIO EM TON/M3 (FORMULA DE YANG) 
!        CALL YANG(UatJ, VmJ, SfTREC, IC, IT)
        !@ *****************************************************************************************
            
        !@ *****************************************************************************************
        !@ CONCENTRAÇÃO INICIAL DE SEDIMENTOS NA SEÇÃO DE JUSANTE - 50% da CTS 
        CSJ2(IC,:) = 0.             !(TON/M3)
        CSJ1 = CSJ2                 !(TON/M3)
        CSM2 = CSJ2                 !(TON/M3)
        CSM1 = CSM2                 !(TON/M3)
        !@ *****************************************************************************************
	    FracS(IC,1) = 0.    !@ porcentagem da concentração de areia
	    FracS(IC,2) = 0.    !@ porcentagem da concentração de silte
	    FracS(IC,3) = 0.    !@ porcentagem da concentração de argila

		iTREC = iTREC + 1	!@ Código do trecho HD da minibacia

	ENDIF


ENDDO

nSEDmini = iSEDaux  !@ contador de minibacias simuladas
allocate(CONTIC(nSEDmini))
allocate(Qmini(nSEDmini,NT))    !@ DCB set/2012
CONTIC = CONTICaux(1:nSEDmini)  !@ armazena o codigo das minibacias simuladas

RETURN
END SUBROUTINE