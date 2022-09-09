!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ Atualizado: Junho 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA PAR_MUSLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA VAI LER OS ARQUIVOS DOS PARÂMETROS DA MUSLE DE CADA BLOCO E O ARQUIVO COM
!@ A TEXTURA DO SOLO DE CADA BLOCO
!@ 
!@ - LÊ E/OU CALCULA O FATOR DE ERODIBILIDADE K PELA EQUAÇÃO DE WILLIAMS (1995)
!@
!@ - LÊ E/OU CALCULA O FATOR DE FRAGMENTOS GROSSEIROS
!@
!@ - LER OS VALORES DO SDR PARA CADA BLOCO DE CADA MINIBACIA OU PARA CADA MINIBACIA
!@
!@-------------------------------------------------------------------------------------
!@  VARIÁVEIS
!@ Fgros   = fator para areia grossa da equação de cálculo da erodibilidade
!@ Fargila = fator para argila da equação de cálculo da erodibilidade
!@ Forg    = fator para carbono orgânico da equação de cálculo da erodibilidade
!@ Fareia  = fator para areia da equação de cálculo da erodibilidade
!@
!@ Kusle   = is the USLE soil erodibility factor [(0.013*ton*m2*hr)/(m3*ton*cm)]
!@ Culse   = is the USLE cover and management factor
!@ Pusle   = is the USLE support practice factor
!@ Rgros   = is the coarse fragment factor
!@
!@ Mareia  = is the percent sand content (0.05-2.00 mm diameter particles)
!@ Msilte  = is the percent silt content (0.002-0.05 mm diameter particles)
!@ Margila = is the percent clay content (< 0.002 mm diameter particles)
!@ Morg    = is the percent of organic mater in the soil
!@ Mrocha  = is the percent rock in the first soil layer
!@
!@ LSAcu   = somatório do fator LS de cada uso em cada minibacia
!@ PSC     = aporte de sedimentos da minibacia para o rio (ton)
!@ SLC     = perda de solo total da minibacia (ton)
!@ SLU     = perda de solo total do bloco da minibacia (ton)
!@ QSC     = descarga sólida das minibacias (ton/s)
!@ QSAUX   = variável que armazena, para o escoamento superficial de cada minibacia, o tempo de retardo (s),
!@           o volume total (m3) e as lâminas (mm) para cada bloco
!@ NPXU    = número de pixels de cada uso em cada minibacia
!@ SDR     = valores do SDR de cada bloco de cada minibacia e de cada minibacia
!@ Ksdr    = parâmetro de calibração do SDR
!@-------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

Subroutine SED_PARAM

USE SED_VARS
USE VARS_MAIN

IMPLICIT NONE

INTEGER iUSO,iCEL,k

REAL Fagros, Fargila, Forg, Fareia

ALLOCATE(Kusle(NU), Cusle(NU), Pusle(NU), Rgros(NU), Ksdr(NU)) !@ Ksdr(NU) DCB 02-07-2011 
ALLOCATE(Mareia(NU), Msilte(NU), Margila(NU), Morg(NU), Mrocha(NU))
ALLOCATE(LSAcu(NC,NU), NPXU(NC,NU), SDR(NC,NU+1)) !@ SDR(NC,NU+1) DCB 02-07-2011 
ALLOCATE(PSC(NC,4), SLC(NC), QSC(NC), QSAUX(NC,NU+2))
ALLOCATE(SLU(NC,IU))
ALLOCATE(QSU(NC,NU), PSU(NC,NU), SLXU(NU))
ALLOCATE(CAREIA(NUMHIDG,NT), CSILTE(NUMHIDG,NT), CARGIL(NUMHIDG,NT)) !CONCENTRAÇÕES PARA GRAVAÇÃO


Kusle = 0.0
Cusle = 0.0
Pusle = 0.0
Rgros = 0.0
Ksdr  = 0.0 !@ DCB 02-07-2011 
PSC   = 0.0
SLC   = 0.0
QSC   = 0.0
QSAUX = 0.0
SLU   = 0.0

QSU  = 0.0
PSU  = 0.0
SLXU = 0.0


!@  QSAUX armazena
!@	TKS da minibacia IC
!@	VSUP da minibacia IC
!@	DSUP bloco 1 da minibacia IC
!@	DSUP bloco 2 da minibacia IC
!@	DSUP bloco 3 da minibacia IC
!@		*
!@		*
!@		*
!@	DSUP bloco N da minibacia IC



!@ *********************************************************************************************
WRITE(*,*)
WRITE(*,*) "DETERMINANDO OS PARAMETROS DA MUSLE..."
WRITE(*,*)


OPEN(FILPAR,FILE='.\input\PARUSO_MUSLE.TXT',STATUS='OLD')
OPEN(FILTEX,FILE='.\input\PARTEXT_MUSLE.TXT',STATUS='OLD')

READ(FILPAR,71)(CABE(k),k=1,4)
READ(FILTEX,72)(CABE(k),k=1,6)

DO iUSO=1,NU-1
	READ(FILPAR,73)Nuso(iUSO),Kusle(iUSO),Cusle(iUSO),Pusle(iUSO),Rgros(iUSO),Ksdr(iUSO)
	READ(FILTEX,74)Nuso(iUSO),Mareia(iUSO),Msilte(iUSO),Margila(iUSO),Morg(iUSO),Mrocha(iUSO)
	
!	!@ Converte % de carbono organico para % de matéria organica
!	Morg(iUSO) = 1.72*Morg(iUSO)
	
	!@ Verifica se o valor do parâmetro K precisa ser calculado
	IF (Kusle(iUSO)==-1.) THEN

		!@		factor que fornece baixo fator de erodibilidade para solos com alto teor de areia grossa e
		!@ alta erodibilidade para solos com areias pequenas
		Fagros = 0.2 + 0.3*exp(-0.256*Mareia(iUSO)*(1. - Msilte(iUSO)/100.))

		!@ Fator que fornece baixo fator de erodibilidade para solos com alta relação argila/silte
		Fargila = (Msilte(iUSO)/(Margila(iUSO) + Msilte(iUSO)))**0.3

		!@ Fator que reduz a erodibilidade do solo para solos com alto conteúdo de carbono orgânico
		Forg = 1. - (0.25*Morg(iUSO)/(Morg(iUSO) + exp(3.72 - 2.95*Morg(iUSO))))

		!@ Fator que reduz a erodibilidade do solo para solos com grandes quantidade de areia
		Fareia = 1. - (0.7*(1. - Mareia(iUSO)/100.)/((1. - Mareia(iUSO)/100.) + exp(-5.51 + 22.9*(1. - Mareia(iUSO)/100.))))

		!@ Fator de erodibilidade do solo de cada bloco
		!@ VERIFICAR: Necessidade de multiplicar por 0.013
		Kusle(iUSO) = Fagros*Fargila*Forg*Fareia

	ENDIF

	!@ Verifica se o valor do parâmetro Rgross precisa ser calculado
	IF (Rgros(iUSO)==-1) THEN
		Rgros(iUSO) = exp(-0.053*Mrocha(iUSO))
	ENDIF

ENDDO

CLOSE(FILPAR)
CLOSE(FILTEX)
!@ *********************************************************************************************

!write(1996,*) (Kusle(iUSO),iUSO=1,NU)

!@ *********************************************************************************************
WRITE(*,*)
WRITE(*,*) "LENDO ARQUIVOS DO PRE_MUSLE..."
WRITE(*,*)

OPEN(FILLSM,FILE='.\input\SED_LSm.txt',STATUS='OLD')
OPEN(FILHRU,FILE='.\input\SED_HRU.txt',STATUS='OLD')
OPEN(FILSDR,FILE='.\input\SED_SDR.txt',STATUS='OLD') !@ DCB 02-07-2011

DO iCEL=1,NC
	READ(FILLSM,'(<NU>F12.4)') (LSAcu(iCEL,iUSO),iUSO=1,NU)
	READ(FILHRU,'(<NU>I8)') (NPXU(iCEL,iUSO),iUSO=1,NU)
	READ(FILSDR,'(<NU>F12.4, F12.4)') (SDR(iCEL,iUSO),iUSO=1,NU+1) !@ DCB 02-07-2011
ENDDO

CLOSE(FILLSM)
CLOSE(FILHRU)
CLOSE(FILSDR) !@ DCB 02-07-2011
!@ *********************************************************************************************


!!@ *********************************************************************************************
!!@ VERIFICA SE LEU CORRETAMENTE OS ARQUIVOS ACIMA
!OPEN(88,FILE='.\output\SED_LSm2.txt')
!OPEN(89,FILE='.\output\SED_hru2.txt')
!DO iCEL=1,NC
!	WRITE(88,'(<NU>F12.4)'), (LSAcu(iCEL,iUSO),iUSO=1,NU)
!	WRITE(89,'(<NU>I8)'), (NPXU(iCEL,iUSO),iUSO=1,NU)
!ENDDO
!CLOSE(88)
!CLOSE(89)
!@ *********************************************************************************************


WRITE(*,*)
WRITE(*,*) 'CALCULO DOS PARAMETROS DA MUSLE ENCERRADO!!!'
WRITE(*,*)


71	FORMAT(5A10)
72	FORMAT(6A10)
73	FORMAT(A10,5F10.5) !@ FORMAT(A10,4F10.5) DCB 02-07-2011 
74	FORMAT(A10,5F10.3)


RETURN

END