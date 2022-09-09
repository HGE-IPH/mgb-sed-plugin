!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ ATUALIZAÇÃO: Jun 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA PARAMETROS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ - ESTA ROTINA VAI LER OS ARQUIVOS DOS PARÂMETROS DA MUSLE DE CADA BLOCO E O ARQUIVO COM
!@ A TEXTURA DO SOLO DE CADA BLOCO
!@ 
!@ - LÊ E/OU CALCULA O FATOR DE ERODIBILIDADE K PELA EQUAÇÃO DE WILLIAMS (1995)
!@
!@ - LÊ E/OU CALCULA O FATOR DE FRAGMENTOS GROSSEIROS
!@
!@ Kusle   = is the USLE soil erodibility factor [(0.013*ton*m2*hr)/(m3*ton*cm)]
!@ Culse   = is the USLE cover and management factor
!@ Pusle   = is the USLE support practice factor
!@ Rgros   = is the coarse fragment factor
!@
!@ Mareia  = is the percent sand content (0.05-2.00 mm diameter particles)
!@ Msilte  = is the percent silt content (0.002-0.05 mm diameter particles)
!@ Margila = is the percent clay content (< 0.002 mm diameter particles)
!@ Mrocha  = is the percent rock in the first soil layer
!@
!@ *********************************************************************************************

Subroutine PARAMETROS

use VAR

implicit none

integer(2) iB,k

real(KIND=8) Fagros, Fargila, Forg, Fareia

allocate(Kusle(numBLOCO),Cusle(numBLOCO),Pusle(numBLOCO),Rgros(numBLOCO))
allocate(Mareia(numBLOCO),Msilte(numBLOCO),Margila(numBLOCO),Morg(numBLOCO),Mrocha(numBLOCO))



write(*,*)
write(*,*) "6. DETERMINANDO OS PARAMETROS DA MUSLE..."
write(*,*)


OPEN(FILUSO,FILE='.\input\PARUSO_MUSLE.TXT',STATUS='OLD')
OPEN(FILTEX,FILE='.\input\PARTEXT_MUSLE.TXT',STATUS='OLD')

READ(FILUSO,71)(CABE(k),k=1,4)
READ(FILTEX,72)(CABE(k),k=1,6)

DO iB=1,numBLOCO
	READ(FILUSO,73)Nuso(iB),Kusle(iB),Cusle(iB),Pusle(iB),Rgros(iB)
	READ(FILTEX,74)Nuso(iB),Mareia(iB),Msilte(iB),Margila(iB),Morg(iB),Mrocha(iB)
		
	!@ Verifica se o valor do parâmetro K precisa ser calculado
	if (Kusle(iB)==-1) then

		!@		factor que fornece baixo fator de erodibilidade para solos com alto teor de areia grossa e
		!@ alta erodibilidade para solos com areias pequenas
		Fagros = 0.2 + 0.3*exp(-0.256*Mareia(iB)*(1 - Msilte(iB)/100.))

		!@ Fator que fornece baixo fator de erodibilidade para solos com alta relação argila/silte
		Fargila = (Msilte(iB)/(Margila(iB) + Msilte(iB)))**0.3

		!@ Fator que reduz a erodibilidade do solo para solos com alto conteúdo de carbono orgânico
		Forg = 1 - (0.25*Morg(iB)/(Morg(iB) + exp(3.72 - 2.95*Morg(iB))))

		!@ Fator que reduz a erodibilidade do solo para solos com grandes quantidade de areia
		Fareia = 1 - (0.7*(1 - Mareia(iB)/100.)/((1 - Mareia(iB)/100.) + exp(-5.51 + 22.9*(1 - Mareia(iB)/100.))))

		!@ Fator de erodibilidade do solo de cada bloco
		Kusle(iB) = Fagros*Fargila*Forg*Fareia

	endif

	!@ Verifica se o valor do parâmetro Rgross precisa ser calculado
	if (Rgros(iB)==-1) then
		Rgros(iB) = exp(-0.053*Mrocha(iB))
	endif


!write(*,*) "BLOCO = ", iB
!write(*,*) "Nuso  = ", Nuso(iB)
!write(*,*) "Kusle = ", Kusle(iB)
!write(*,*) "Cusle = ", Cusle(iB)
!write(*,*) "Pusle = ", Pusle(iB)
!write(*,*) "Rgros = ", Rgros(iB)

ENDDO

CLOSE(FILUSO)
CLOSE(FILTEX)


write(*,*)
write(*,*) 'CALCULO DOS PARAMETROS DA MUSLE ENCERRADO!!!'
write(*,*)


71	FORMAT(5A10)
72	FORMAT(6A10)
73	FORMAT(A10,4F10.5)
74	FORMAT(A10,5F10.3)


return

end