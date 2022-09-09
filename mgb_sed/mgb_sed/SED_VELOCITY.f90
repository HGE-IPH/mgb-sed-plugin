!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Abr de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA FALLVELOCITY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA CALCULA A VELOCIDADE DE QUEDA DAS PART�CULAS DE SEDIMENTOS
!@
!@ *********************************************************************************************

SUBROUTINE FALLVELOCITY

USE SED_VARS

IMPLICIT NONE

REAL FatM, FatN, FatP, DCar


!@ ***************************************************************************************
!@ VELOCIDADE DE QUEDA DAS PART�CULAS (Wu e Wang, 2006) - M/S
!@ OBS.: Colocar em uma subrotina pr�pria, permitindo outras equa��es
	
!@ AREIA
	FatM = 53.5*exp(-0.65*Fform(1))	!@ Fator M da equa��o
	FatN = 5.65*exp(-2.5*Fform(1))	!@ Fator N da equa��o
	FatP = 0.7 + 0.9*Fform(1)		!@ Fator n da equa��o
	DCar = DMP(1)*((2.65 - 1.)*9.81/(VISC**2))**(1./3.)	!@ Di�metro caracter�stico
	WSP(1) = (sqrt(0.25+(4.*FatN*(DCar**3)/(3.*(FatM**2)))**(1./FatP)) - 0.5)**FatP
	WSP(1) = ((FatM*VISC)/(FatN*DMP(1)))*WSP(1)	!@ Velocidade de queda da part�cula de areia (m/s)

!@ SILTE
	FatM = 53.5*exp(-0.65*Fform(2))
	FatN = 5.65*exp(-2.5*Fform(2))
	FatP = 0.7 + 0.9*Fform(2)
	DCar = DMP(2)*((2.65 - 1.)*9.81/(VISC**2))**(1./3.)
	WSP(2) = (sqrt(0.25+(4.*FatN*(DCar**3)/(3.*(FatM**2)))**(1./FatP)) - 0.5)**FatP
	WSP(2) = ((FatM*VISC)/(FatN*DMP(2)))*WSP(2)	!@ m/s

!@ ARGILA
	FatM = 53.5*exp(-0.65*Fform(3))
	FatN = 5.65*exp(-2.5*Fform(3))
	FatP = 0.7 + 0.9*Fform(3)
	DCar = DMP(3)*((2.65 - 1.)*9.81/(VISC**2))**(1./3.)
	WSP(3) = (sqrt(0.25+(4.*FatN*(DCar**3)/(3.*(FatM**2)))**(1./FatP)) - 0.5)**FatP
	WSP(3) = ((FatM*VISC)/(FatN*DMP(3)))*WSP(3)	!@ m/s

!@ ***************************************************************************************


RETURN
END SUBROUTINE