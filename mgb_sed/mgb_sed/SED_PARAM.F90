!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ Atualizado: Junho 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA PAR_MUSLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA VAI LER OS ARQUIVOS DOS PAR?METROS DA MUSLE DE CADA BLOCO E O ARQUIVO COM
!@ A TEXTURA DO SOLO DE CADA BLOCO
!@ 
!@ - L? E/OU CALCULA O FATOR DE ERODIBILIDADE K PELA EQUA??O DE WILLIAMS (1995)
!@
!@ - L? E/OU CALCULA O FATOR DE FRAGMENTOS GROSSEIROS
!@
!@ - LER OS VALORES DO SDR PARA CADA BLOCO DE CADA MINIBACIA OU PARA CADA MINIBACIA
!@
!@-------------------------------------------------------------------------------------
!@  VARI?VEIS
!@ Fgros   = fator para areia grossa da equa??o de c?lculo da erodibilidade
!@ Fargila = fator para argila da equa??o de c?lculo da erodibilidade
!@ Forg    = fator para carbono org?nico da equa??o de c?lculo da erodibilidade
!@ Fareia  = fator para areia da equa??o de c?lculo da erodibilidade
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
!@ LSAcu   = somat?rio do fator LS de cada uso em cada minibacia
!@ PSC     = aporte de sedimentos da minibacia para o rio (ton)
!@ SLC     = perda de solo total da minibacia (ton)
!@ SLU     = perda de solo total do bloco da minibacia (ton)
!@ QSC     = descarga s?lida das minibacias (ton/s)
!@ QSAUX   = vari?vel que armazena, para o escoamento superficial de cada minibacia, o tempo de retardo (s),
!@           o volume total (m3) e as l?minas (mm) para cada bloco
!@ NPXU    = n?mero de pixels de cada uso em cada minibacia
!@ SDR     = valores do SDR de cada bloco de cada minibacia e de cada minibacia
!@ Ksdr    = par?metro de calibra??o do SDR
!@-------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

Subroutine SED_PARAM

USE SED_VARS
USE VARS_MAIN

IMPLICIT NONE

INTEGER iUSO,iCEL,k

REAL Fagros, Fargila, Forg, Fareia

ALLOCATE(Kusle(NB,NU), Cusle(NB,NU), Pusle(NB,NU), Rgros(NB,NU), Ksdr(NB,NU)) !@ Ksdr(NU) DCB 02-07-2011 
ALLOCATE(Mareia(NB,NU), Msilte(NB,NU), Margila(NB,NU), Morg(NB,NU), Mrocha(NB,NU))
ALLOCATE(LSAcu(NC,NU), NPXU(NC,NU), SDR(NC,NU+1)) !@ SDR(NC,NU+1) DCB 02-07-2011 
ALLOCATE(PSC(NC,4), SLC(NC), QSC(NC), QSAUX(NC,NU+2))
ALLOCATE(SLU(NC,IU))
ALLOCATE(QSU(NC,NU), PSU(NC,NU), SLXU(NU))
ALLOCATE(CAREIA(NUMHIDG,NT), CSILTE(NUMHIDG,NT), CARGIL(NUMHIDG,NT)) !CONCENTRA??ES PARA GRAVA??O


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

!READ(FILPAR,71)(CABE(k),k=1,4)
!READ(FILTEX,72)(CABE(k),k=1,6)

DO IB=1,NB
    READ(FILPAR,71)(CABE(1))
    READ(FILPAR,71)(CABE(k),k=1,4)
    
    READ(FILTEX,72)(CABE(1))
    READ(FILTEX,72)(CABE(k),k=1,6)

    DO iUSO=1,NU
	    READ(FILPAR,73)Nuso(iUSO),Kusle(IB,iUSO),Cusle(IB,iUSO),Pusle(IB,iUSO),Rgros(IB,iUSO),Ksdr(IB,iUSO)
	    READ(FILTEX,74)Nuso(iUSO),Mareia(IB,iUSO),Msilte(IB,iUSO),Margila(IB,iUSO),Morg(IB,iUSO),Mrocha(IB,iUSO)
	        	
    	write(*,*) 'IC, IB, iUSO = ', IC, IB, iUSO
    !	!@ Converte % de carbono organico para % de mat?ria organica
    !	Morg(iUSO) = 1.72*Morg(iUSO)
    	
	    !@ Verifica se o valor do par?metro K precisa ser calculado
	    IF (Kusle(IB,iUSO)==-1.) THEN

		    !@		factor que fornece baixo fator de erodibilidade para solos com alto teor de areia grossa e
		    !@ alta erodibilidade para solos com areias pequenas
		    Fagros = 0.2 + 0.3*exp(-0.256*Mareia(IB,iUSO)*(1. - Msilte(IB,iUSO)/100.))

		    !@ Fator que fornece baixo fator de erodibilidade para solos com alta rela??o argila/silte
		    Fargila = (Msilte(IB,iUSO)/(Margila(IB,iUSO) + Msilte(IB,iUSO)))**0.3

		    !@ Fator que reduz a erodibilidade do solo para solos com alto conte?do de carbono org?nico
		    Forg = 1. - (0.25*Morg(IB,iUSO)/(Morg(IB,iUSO) + exp(3.72 - 2.95*Morg(IB,iUSO))))

		    !@ Fator que reduz a erodibilidade do solo para solos com grandes quantidade de areia
		    Fareia = 1. - (0.7*(1. - Mareia(IB,iUSO)/100.)/((1. - Mareia(IB,iUSO)/100.) + exp(-5.51 + 22.9*(1. - Mareia(IB,iUSO)/100.))))

		    !@ Fator de erodibilidade do solo de cada bloco
		    !@ VERIFICAR: Necessidade de multiplicar por 0.013
		    Kusle(IB,iUSO) = Fagros*Fargila*Forg*Fareia !/0.013
	    ENDIF

	    !@ Verifica se o valor do par?metro Rgross precisa ser calculado
	    IF (Rgros(IB,iUSO)==-1) THEN
		    Rgros(IB,iUSO) = exp(-0.053*Mrocha(IB,iUSO))
	    ENDIF

    ENDDO
ENDDO

CLOSE(FILPAR)
CLOSE(FILTEX)

!Kusle = 1.5*Kusle
!Cusle = 1.5*Cusle
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