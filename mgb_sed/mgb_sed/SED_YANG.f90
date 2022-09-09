!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abr de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA YANG >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA CALCULA A CAPACIDADE DE TRANSPORTE DAS PARTÍCULAS DE SEDIMENTO NO TRECHO DE RIO
!@   PELA FORMULAÇÃO DE YANG
!@
!@ *********************************************************************************************

SUBROUTINE YANG(Uat, VmTREC, SfTREC, IC, IT)

USE SED_VARS

IMPLICIT NONE

INTEGER i, IC, IT
REAL Uat, VmTREC, SfTREC
REAL REX(3), UcWs(3)
REAL FMY(3), FNY(3), LOGCTS, LOGCTSaux(3)
REAL aux

DO i= 1,3
    !@ Reynolds de atrito - mínimo igual a 1.2 (VERIFICAR OPÇÕES MELHORES)
    REX(i) = max(Uat*DMP(i)/VISC,1.2)

    !@ Estimativa da Velocidade Crítica sobre Velocidade de Queda (Uc/Ws)
	IF (REX(i) < 70.) THEN
		UcWs(i) = 0.66 + 2.5/(log10(REX(i)) - 0.06)
	ELSE
		UcWs(i) = 2.05
	ENDIF

    !@ Parâmetros M e N da equação de capacidade de transporte de YANG
	IF (DMP(i) <= 0.002) THEN
		FMY(i) = 5.435 - 0.286*log10(WSP(i)*DMP(i)/VISC) - 0.457*log10(Uat/WSP(i))
		FNY(i) = 1.799 - 0.409*log10(WSP(i)*DMP(i)/VISC) - 0.314*log10(Uat/WSP(i))
	ELSE
		FMY(i) = 6.681 - 0.633*log10(WSP(i)*DMP(i)/VISC) - 4.816*log10(Uat/WSP(i))
		FNY(i) = 2.874 - 0.305*log10(WSP(i)*DMP(i)/VISC) - 0.282*log10(Uat/WSP(i))
	ENDIF

    aux = VmTREC*SfTREC/WSP(i) - UcWs(i)*SfTREC
    LOGCTS = 0.0
    CTS(IC,i) = 0.0
    IF (aux > 0.) THEN
        LOGCTS = FMY(i) + FNY(i)*log10( aux )
        LOGCTSaux(i) = LOGCTS
        !@ Potencial de transporte do trecho de rio 
	    CTS(IC,i) = 10.**(LOGCTS) !@ (PPM)
        !@ Potencial de transporte do trecho de rio
	    CTS(IC,i) = CTS(IC,i)/(10.**(6.)) !@ (TON/M3)
        !@ Capacidade de transporte do trecho de rio
	    CTS(IC,i) = CTS(IC,i)*FracS(IC,i) !@ (TON/M3)
    ENDIF
  
ENDDO

!IF (IC == 3834) THEN
!write(*,*) 'IT, IC, Sf  = ', IT, IC, SfTREC
!write(*,*) 'DMP, Visc   = ', DMP(1), VISC
!write(*,*) 'REX, WcWs   = ', REX(1), UcWs(1)
!write(*,*) 'FMY, FNY    = ', FMY(1), FNY(1)
!write(*,*) 'VERIF, Log  = ', aux, LOGCTSaux(1)
!write(*,*) 'CTS*, Frac  = ', CTS(IC,1), FracS(IC,1)
!pause
!ENDIF


RETURN
END SUBROUTINE