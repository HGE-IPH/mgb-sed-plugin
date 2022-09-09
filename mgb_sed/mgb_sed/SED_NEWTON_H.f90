!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abr de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA NEWTON_H >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA CALCULA PROFUNDIDADE D'ÁGUA NO TRECHO UTILIZANDO A EQUAÇÃO DE MANNING E O
!@ MÉTODO ITERATIVO DE NEWTON-RAPHSON
!@
!@ *********************************************************************************************

SUBROUTINE NEWTON_H(QmTRECx, HmTRECx, SfTRECx)

USE VARS_MAIN
USE SED_VARS

IMPLICIT NONE

REAL QmTRECx, HmTRECx, SfTRECx
REAL HK, ERROH, ERROQ, TOL, QK, QK2
INTEGER itera

HK    = 0.01*BRIO(IC) !@ Estimativa inicial da profundidade do trecho
ERROH = 9999.
itera = 0
TOL   = 0.000001
DO WHILE (erroH > TOL)
    itera   = itera + 1
    QK      = ((BRIO(IC)*HK)**(5./3.))*(SfTRECx**(1./2.))/(((BRIO(IC)+2.*HK)**(2./3.))*RUGMAN(IC)) !FMF em 06/05/2019 - leitura do novo MINI.GTP
    HmTRECx = HK - (1.-QmTRECx/QK)*(3.*HK*(BRIO(IC)+2.*HK))/(5.*BRIO(IC)+6.*HK)
    QK2     = ((BRIO(IC)*HmTRECx)**(5./3.))*(SfTRECx**(1./2.))/(((BRIO(IC)+2.*HmTRECx)**(2./3.))*RUGMAN(IC)) !FMF em 06/05/2019 - leitura do novo MINI.GTP
    ERROH   = sqrt((HmTRECx - HK)**2.)
    HK      = HmTRECx
    ERROQ   = sqrt((QmTRECx - QK2)**2.)

    IF (MOD(itera,30)==0) then
        TOL = TOL*10
    ENDIF

ENDDO


RETURN
END SUBROUTINE