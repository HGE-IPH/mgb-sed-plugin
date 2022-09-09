!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ ATUALIZAÇÃO: Jun 2011 - alterada rotina de cálculo da área acumulada
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PROGRAMA PRE_SED >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ ROTINA PARA CÁLCULO DOS PARÂMETROS DA MUSLE E DETERMINAÇÃO DO FATOR CONSTANTE PARA CADA
!@ BLOCO DE CADA MINI-BACIA
!@
!@	EQUAÇÃO DA MUSLE (Williams, 1975)
!@
!@	Sed = 11.8*[(Qsup*qpico*Aphru)^(0.56)]*K*C*P*LS*Rgros
!@
!@	Sed    = ton - perda de solo
!@	Qsup   = mmH2O - volume de escoamento superficial
!@	qpico  = m3/s - taxa de pico do escoamento superficial
!@	Aphru  = ha - área do bloco, ou área média do pixel do bloco
!@	K      = 0.013*(ton*m2*h)/(m3*ton*cm) - fator de erodibilidade do solo
!@	C,P,LS = adimensional - fatores da MUSLE
!@	Rgros  = exp(-0.053*Mrocha) - fator de rocha
!@	Mrocha = % de rocha na primeira camada de solo
!@
!@
!@ CÁLCULOS:
!@ - Lê os layers de entrada
!@ - Determina ou ler o layer de Áreas Acumuladas
!@ - Calcula o Fator LS
!@ - Lê ou calcula os demais parâmetros da MUSLE
!@ - Gera o termo constante da equação para entrada no MGB-IPH
!@ 
!@ Sed_cte = 11.8*[Ahru^(0.56)]*K*C*P*LS*Rgros
!@
!@ - 
!@ - 
!@ - 
!@
!@
!@ *********************************************************************************************


Program PRE_SED

Use VAR

implicit none

real(KIND=8) Zx, Zy, Sper, Sgrau, AreaM, SLOPE, COEF, CompL, fatorX
integer linSeg, colSeg, dirSeg



!@ ************************************** CABECALHO ********************************************

write(*,*)
write(*,*)
write(*,*) "------------------------ PROGRAMA PRE-SED 2.0 ------------------------"
write(*,*) "-                                                                    -"
write(*,*) "- OBJETIVO: Calcular parametros da MUSLE e determinar um fator       -"
write(*,*) "-         constante para utilizacao no MGB_SED.                      -"
write(*,*) "-                                                                    -"
write(*,*) "- ARQUIVOS DE ENTRADA:                                               -"
write(*,*) "-  1. MNT.txt        (Matriz do MNT bruto - ascii ArcGIS)            -"
write(*,*) "-  2. DIR.txt        (Matriz de Direcoes de Fluxo - ascii ArcGIS)    -"
write(*,*) "-  3. MINI.txt       (Matriz de Minibacias - ascii ArcGIS)           -"
write(*,*) "-  4. BLOCOS.txt     (Matriz de Blocos - ascii ArcGIS)               -"
write(*,*) "-  5. AREA_ACU.txt   (Matriz de Areas Acumuladas - ascii ArcGIS)     -"
write(*,*) "-                                                                    -"
write(*,*) "-                                                                    -"
write(*,*) "- ARQUIVO DE SAIDA:                                                  -"
write(*,*) "-  1. LS_2D.txt      (Matriz com o fator LS medio para cada HRU)     -"
write(*,*) "-  2. AREA_ACU.txt*  (Matriz de Areas Acumuladas - ascii ArcGIS)     -"
write(*,*) "-                                                                    -"
write(*,*) "-                                                                    -"
write(*,*) "- CONTATO: Diogo Costa Buarque  - diogo.buarque@gmail.com            -"
write(*,*) "-                Instituto de Pesquisas Hidraulicas (IPH/UFRGS)      -"
write(*,*) "-                                                                    -"
write(*,*) "- * Opcional e gerado apenas quando ainda nao existe!                -"
write(*,*) "-                                                                    -"
write(*,*) "- VERSAO: Jun/2011                                                   -"
write(*,*) "----------------------------------------------------------------------"
write(*,*) "(tecle enter)"
read(*,*)
write(*,*)



!@ >>>>>>>>>>>>>>>>> ATENCAO PARA O VALOR NUMERICO DAS DIRECOES >>>>>>>>>>>>>>>>>>>

!@   F  G  H          ArcView:  32 64 128    IDRISI:   64  128  1 
!@   E  *  A                    16  *  1               32   *   2
!@   D  C  B                     8  4  2               16   8   4

!@ definicao da numeracao das direcoes
A=1   
B=2   
C=4   
D=8   
E=16  
F=32  
G=64  
H=128 

!@ definicao do vetor de direcoes daux
ddaux(1)=A
ddaux(2)=B
ddaux(3)=C
ddaux(4)=D
ddaux(5)=E
ddaux(6)=F
ddaux(7)=G
ddaux(8)=H

!@ definicao da posicao relativa dos pixels vizinhos
dlin(A)=0
dcol(A)=1
dlin(B)=1
dcol(B)=1
dlin(C)=1
dcol(C)=0
dlin(D)=1
dcol(D)=-1
dlin(E)=0
dcol(E)=-1
dlin(F)=-1
dcol(F)=-1
dlin(G)=-1
dcol(G)=0
dlin(H)=-1
dcol(H)=1



!@ >>>>>>>>>>>>>>>>>>>>>>> LEITURA DOS ARQUIVOS DE ENTRADA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

write(*,*) "1. INFORMACOES DE ENTRADA..."
	call Entrada



!@ >>>>>>>>>>> CÁLCULO OU LEITURA DO ARQUIVO DE AREA ACUMULADA - AREA.txt >>>>>>>>>>>>>>>>

write(*,*)
write(*,*) "3. OBTENDO AREAS ACUMULADAS..."
write(*,*)

	if (resparea==1 .AND. respls ==1) then

		!@ Calcula áreas acumuladas
		call AREA_ACUM

	elseif (respls==1) then

		!@ Ler arquivo existente de áreas acumuladas
		write(*,*) 
		write(*,*) "3.1 LENDO ARQUIVO DE AREA ACUMULADA..."
		write(*,*)
		write(*,*) "	- ARQUIVO AREA_ACU.txt"
		write(*,*)
		write(*,*)

		open(10,FILE='.\input\AREA_ACU.txt',STATUS='UNKNOWN')
	
		read(10,*) TRASH,nCaux
		read(10,*) TRASH,nLaux
		read(10,*) TRASH,Xminaux
		read(10,*) TRASH,Yminaux
		read(10,*) TRASH,cellsizeaux
		read(10,*) TRASH,noDataFLAG

		if ((nC/=nCaux) .or. (nL/=nLaux) .or. (abs(Xmin-Xminaux)>TOL) .or. (abs(Ymin-Yminaux)>TOL) .or. (abs(cellsize-cellsizeaux)>TOL)) then

			write(*,*) "******************************************************************"
			write(*,*)
			write(*,*) "ATENCAO!!"
			write(*,*)
			write(*,*) "Os parametros abaixo devem ser iguais ao do MNT:"
			write(*,*) "	- numero de colunas;"
			write(*,*) "	- numero de linhas;"
			write(*,*) "	- Xmin;"
			write(*,*) "	- Ymin;"
			write(*,*) "	- tamanho da celula;"
			write(*,*)
			write(*,*) " REVEJA SEUS DADOS DE ENTRADA!!!"
			write(*,*)
			write(*,*) "******************************************************************"
			write(*,*)
			write(*,*) "(tecle enter)"
			read(*,*)
			stop
		endif

        !@ ALOCA VARIÁVEIS PARA A AREA ACUMULADA
        allocate(AREAAcu(nL,nC))

		!@ Le dados:
		do i=1,nL
			read(10,*) (AREAAcu(i,j),j=1,nC)
		enddo
		close(10)

	endif
write(*,*)
write(*,*)



!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> CÁLCULO DO FATOR LS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

write(*,*)
write(*,*) "4. OBTENDO O FATOR LS-2D..."
write(*,*)

    if (respls==1) then

	    !@ Calcula o fator LS_2D
	    call LS_2D_DG1996

    else
    
        !@ Ler arquivo existente com o fator LS-2D
        write(*,*)
        write(*,*) "4.1 LENDO ARQUIVO COM O FATOR LS-2D..."
        write(*,*)
        write(*,*) "	- ARQUIVO LS_2D_DG1996.txt"
	    write(*,*)
	    write(*,*)

	    open(10,FILE='.\input\LS_2D_DG1996.txt',STATUS='UNKNOWN')
	
	    read(10,*) TRASH,nCaux
	    read(10,*) TRASH,nLaux
	    read(10,*) TRASH,Xminaux
	    read(10,*) TRASH,Yminaux
	    read(10,*) TRASH,cellsizeaux
	    read(10,*) TRASH,noDataFLAG

	    if ((nC/=nCaux) .or. (nL/=nLaux) .or. (abs(Xmin-Xminaux)>TOL) .or. (abs(Ymin-Yminaux)>TOL) .or. (abs(cellsize-cellsizeaux)>TOL)) then

		    write(*,*) "******************************************************************"
		    write(*,*)
		    write(*,*) "ATENCAO!!"
		    write(*,*)
		    write(*,*) "Os parametros abaixo devem ser iguais ao do MNT:"
		    write(*,*) "	- numero de colunas;"
		    write(*,*) "	- numero de linhas;"
		    write(*,*) "	- Xmin;"
		    write(*,*) "	- Ymin;"
		    write(*,*) "	- tamanho da celula;"
		    write(*,*)
		    write(*,*) " REVEJA SEUS DADOS DE ENTRADA!!!"
		    write(*,*)
		    write(*,*) "******************************************************************"
		    write(*,*)
		    write(*,*) "(tecle enter)"
		    read(*,*)
		    stop
        endif

        !@ ALOCA VARIÁVEIS DO LS-2D:
        allocate(LS2d(nL,nC))   !,WEIGHT(nL,nC))

	    !@ Le dados:
	    do i=1,nL
		    read(10,*) (LS2d(i,j),j=1,nC)
	    enddo
	    close(10)
	    
	    
	    
	    !@ Calculando LS médio
	    CONT = 0.
	    do lin = 1,nL
	        do col = 1,nC
		
		        !@ Exclindo pixel com valor NoData dos cálculos
		        if (MINI(lin,col)==noDataFLAG .OR. MNT(lin,col) == noDataFLAG) cycle

		        codBAC = MINI(lin,col)		!@ obtem código da minibacia
		        codBLC = BLOCOS(lin,col)	!@ obtem código do bloco
                           
                !@ ------------	ACUMULANDO O FATOR LS DO BLOCO DE CADA MINIBACIA
		        LS_acu(codBAC,codBLC) = LS_acu(codBAC,codBLC) + LS2d(lin,col)
	            !@ **********************************************************************

			    CONT(codBAC,codBLC) = CONT(codBAC,codBLC) + 1
        	enddo
        enddo
        
        deallocate(LS2d)    !,MNT)
	
    endif
write(*,*)
write(*,*)



!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> CÁLCULO DO FATOR K >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
call FATOR_SDR



!@ >>>>>>>>>>>>>>>>>>>>> ESCREVENDO PARÂMETROS PARA A MUSLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write(*,*)
write(*,*) "6. ESCREVENDO ARQUIVOS DE PARAMETROS PARA A MUSLE"
write(*,*)

open(87,file='.\output\SED_SDR.txt')
open(88,file='.\output\SED_LSm.txt')
open(89,file='.\output\SED_HRU.txt')

do lin=1,numMINI
	write(87,'(<numBLOCO>F12.4,F12.4)'), (SDR(lin,j),j=1,numBLOCO),SDRaux(lin)
	write(88,'(<numBLOCO>F12.4)'), (LS_acu(lin,j),j=1,numBLOCO)
	write(89,'(<numBLOCO>I8)'), (CONT(lin,j),j=1,numBLOCO)
enddo

close(87)
close(88)
close(89)





!!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> CÁLCULO DO FATOR K >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!call PARAMETROS    !@ ROTINA INSERIDA NO MODELO
!
!
!
!!@ >>>>>> CÁLCULO DO TERMO CONSTANTE DA MUSLE PARA CADA BLOCO DE CADA MINI-BACIA>>>>>>>>>>
!write(*,*)
!write(*,*) "7. DETERMINANDO TERMOS CONSTANTE DA MUSLE"
!write(*,*)
!
!open(87,file='.\output\SED_cte.txt')
!open(88,file='.\output\SED_LSm.txt')
!open(89,file='.\output\SED_hru.txt')
!do lin=1,numMINI
!	do col=1,numBLOCO
!
!!@ - -------- OBTENDO ÁREA MÉDIA DOS PIXELS DE CADA BLOCO DE CADA MINIBACIA (ha) ---------
!		if (CONT(lin,col) /= 0) then
!			Apixel(lin,col) = Apixel(lin,col)*100/CONT(lin,col)
!		else
!			Apixel(lin,col) = 0.
!		endif
!
!!@ ------------------- OBTENDO O TERMO CONSTANTE DA EQUAÇÃO DA MUSLE ---------------------
!		SED(lin,col) = 11.8*(Apixel(lin,col)**0.56)*Kusle(col)*Cusle(col)*Pusle(col)*Rgros(col)*LS_acu(lin,col)
!
!	enddo
!
!	write(87,'(<numBLOCO>F12.4)'), (SED(lin,j),j=1,numBLOCO)
!	write(88,'(<numBLOCO>F12.4)'), (LS_acu(lin,j),j=1,numBLOCO)
!	write(89,'(<numBLOCO>I8)'), (CONT(lin,j),j=1,numBLOCO)
!
!enddo
!close(87)
!close(88)
!close(89)



write(*,*)
write(*,*) 'PARABENS!!!'
write(*,*) 'A ROTINA PRE_MUSLE ENCERROU!!!'
write(*,*)
write(*,*)
READ(*,*)
STOP

end