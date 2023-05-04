!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Fevereiro de 2010
!@
!@ ATUALIZAÇÃO: Jun 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA SDR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ CALCULA CALCULA O SDR DAS MINIBACIAS PELA METODOLOGIA DE LENHART ET AL (2005)
!@
!@
!@ *********************************************************************************************

Subroutine FATOR_SDR

use VAR

implicit none

integer(2) caminho, linAT, colAT, linSEG, colSEG, codDIR
real(KIND=8) distaux
real(KIND=8) Zx, Zy, Sper
real(KIND=8) Wiaux(numMINI), WiDiaux(numMINI)

!@ ALOCA VARIÁVEIS DO SDR
allocate(Wi(numMINI,numBLOCO), WiDi(numMINI,numBLOCO), SDR(numMINI,numBLOCO), SDRaux(numMINI))
allocate(DIST(nL,nC), WEIGHT(nL,nC))

!@ INICIALIZANDO VARIÁVEIS
DIST = noDataFLAG
WEIGHT = noDataFLAG
Wi = 0.
Wiaux = 0.
WiDi = 0.
WiDiaux = 0.


write(*,*)
write(*,*) "5. CALCULANDO O SDR DAS MINIBACIAS..."
write(*,*)

do lin=1,nL

    !@ Escreve valor da linha a cada 10 linha calculadas
    if (mod(lin,10)==0) then
		write(*,*) "lin: ",lin
	end if

    do col=1,nC
    
        !@ Exclui pixels fora da bacia
        if (MINI(lin,col) /= noDataFLAG) then
        
            !@ Verifica se o pixel analisado pertence à drenagem
            if (DRENAGEM(lin,col) == 1) then
                WEIGHT(lin,col) = 1.    !@ se pertence à drenagem o peso é 1
                DIST(lin,col) = 1.      !@ se pertence à drenagem a distância é nula
                !@ acumula pesos para cada bloco de cada minibacia
                Wi(MINI(lin,col),BLOCOS(lin,col)) = Wi(MINI(lin,col),BLOCOS(lin,col)) + WEIGHT(lin,col)
                !@ acumula pesos para cada minibacia
                Wiaux(MINI(lin,col)) = Wiaux(MINI(lin,col)) + WEIGHT(lin,col)
                !@ acumula peso_x_distancia para cada bloco de cada minibacia
                WiDi(MINI(lin,col),BLOCOS(lin,col)) = WiDi(MINI(lin,col),BLOCOS(lin,col)) + WEIGHT(lin,col)*DIST(lin,col)
                !@ acumula peso_x_distancia para cada minibacia
                WiDiaux(MINI(lin,col)) = WiDiaux(MINI(lin,col)) + WEIGHT(lin,col)*DIST(lin,col)
                cycle !@ PASSA PARA O PROXIMO PIXEL
            endif
        
        
	        !@ **********************************************************************            
            !@ ------------	CALCULANDO COORDENADAS DOS VERTICES DO PIXEL ANALISADO
	            Xesq=Xmin+(col-1)*cellsize
	            Xdir=Xesq+cellsize
	            Yinf=Ymin+(nL-lin)*cellsize
	            Ysup=Yinf+cellsize
!OK	            write(*,*) 'Xesq, Xdir = ', Xesq, Xdir
!OK	            write(*,*) 'Yinf, Ysup = ', Yinf, Ysup
    	        
	        !@ ------------	CALCULA AS DIMENSÕES DO PIXEL ANALISADO
	            DeltaX  = 0.0
			    DeltaY  = 0.0
			    DeltaXY = 0.0
			    CELLx   = 0.0
			    CELLy   = 0.0
			    call Project_COMP(Xesq,Xdir,Ysup,Yinf,DeltaX,DeltaY,DeltaXY,CELLx,CELLy,lado,diag)
	        
	        !@ ------------	OBTENDO DERIVADAS PARCIAIS POR DIFERENÇAS FINITAS
	        !@ FONTE: Wilson e Gallant 2000. Terrain Analysis: Principles and Applications.
	            Zx = 0.
	            Zy = 0.
		        call DERIVA_WG2000(Zx,Zy)

	        !@ ------------	OBTENDO O PONDERADOR PELO DECLIVE DO PIXEL
	        !@ FONTE: Wilson e Gallant 2000. Terrain Analysis: Principles and Applications.
	            Sper = sqrt(Zx**2. + Zy**2.)*100.0			!@ Declive em percentagem
	            WEIGHT(lin,col) = 0.
	            if (Sper /= 0.) WEIGHT(lin,col) = 1./Sper   !@ Matriz de peso para calculo do SDR
	        !@ **********************************************************************
        
        
            DIST(lin,col) = 0.
            caminho = 0.
            linAT = lin !@ linha do pixel analisado
            colAT = col !@ coluna do pixel analisado
            codBAC = MINI(lin,col) !@ código da minibacia do pixel analisado
            codBLC = BLOCOS(lin,col) !@ codigo do bloco do pixel analisado
            do while (caminho==0)
                codDIR = DIR(linAT,colAT)
                linSEG = linAT + dlin(codDIR) !@ linha do pixel seguinte
                colSEG = colAT + dcol(codDIR) !@ coluna do pixel seguinte
                
                !@ Calcula as coordenadas dos vertices do pixel atual
	            Xesq=Xmin+(colAT-1)*cellsize
	            Xdir=Xesq+cellsize
	            Yinf=Ymin+(nL-linAT)*cellsize
	            Ysup=Yinf+cellsize
	                	        
	            !@ Chama rotina para cálculo das dimensões do pixel
	            DeltaX  = 0.0
			    DeltaY  = 0.0
			    DeltaXY = 0.0
			    CELLx   = 0.0
			    CELLy   = 0.0
			    call Project_COMP(Xesq,Xdir,Ysup,Yinf,DeltaX,DeltaY,DeltaXY,CELLx,CELLy,lado,diag)
    			
			    !@ determina posicao relativa ao pixel seguinte
			    distaux = 0.0
	            if ((linAT==linSEG).OR.(colAT==colSEG)) then
		            if (linAT==linSEG) then
			            distaux = DeltaX    !@ comprimento para drenagem na horizontal (km)
		            else
			            distaux = DeltaY    !@ comprimento para drenagem na vertical (km)
		            end if
	            else
		            distaux = DeltaXY       !@ comprimento para drenagem na diagonal (km)
	            end if
    	        
	            DIST(lin,col) = DIST(lin,col) + distaux !@ KILOMETROS

                if (DRENAGEM(linSEG,colSEG)==1 .OR. MINI(linSEG,colSEG)==noDataFLAG) caminho = 1 !@ para o looping quando chega na rede
                
                !@ Caso o pixel seguinte tenha seu comprimento calculado, utiliza-o para ganhar tempo!
                if (DIST(linSEG,colSEG) > 0.) then
                    DIST(lin,col) = DIST(lin,col) + DIST(linSEG,colSEG)
                    caminho = 1
                endif
        
                linAT = linSEG !@ atualiza linha do pixel atual
                colAT = colSEG !@ atualiza coluna do pixel atual
            enddo
            
            !@ acumula pesos para cada bloco de cada minibacia
            Wi(codBAC,codBLC) = Wi(codBAC,codBLC) + WEIGHT(lin,col)
            !@ acumula pesos para cada minibacia
            Wiaux(codBAC) = Wiaux(codBAC) + WEIGHT(lin,col)
            !@ acumula peso_x_distancia para cada bloco de cada minibacia
            WiDi(codBAC,codBLC) = WiDi(codBAC,codBLC) + WEIGHT(lin,col)*DIST(lin,col)
            !@ acumula peso_x_distancia para cada minibacia
            WiDiaux(codBAC) = WiDiaux(codBAC) + WEIGHT(lin,col)*DIST(lin,col)
        endif
    enddo
enddo


!@ ----	CALCULO DO FATOR SDR POR BLOCO DE CADA MINIBACIA E POR MINIBACIA
SDR = 0.
SDRaux = 0.
do lin=1,numMINI
    do col=1,numBLOCO
        !@ Fator SDR por bloco
        if (WiDi(lin,col)/=0) SDR(lin,col) = min(Wi(lin,col)/WiDi(lin,col),1.)
    enddo
    !@ Fator SDR por minibacia
    if (WiDiaux(lin)/=0) SDRaux(lin) = min(Wiaux(lin)/WiDiaux(lin),1.)
enddo
!@ **********************************************************************


!write(*,*)
!write(*,*) "5.1 GRAVANDO ARQUIVO DIST_MAP.txt..."
!write(*,*)
!
!open(80,file='.\output\DIST_MAP.txt')
!	write(80,'(A5,A9,I5)') 'ncols','',nC
!	write(80,'(A5,A9,I5)') 'nrows','',nL
!	write(80,'(A9,A5,F16.12)') 'xllcorner','',Xmin
!	write(80,'(A9,A5,F16.12)') 'yllcorner','',Ymin
!	write(80,'(A8,A6,F18.16)') 'cellsize','',cellsize
!	write(80,'(A12,A2,I5)') 'NODATA_value','',noDataFLAG
!
!
!write(*,*)
!write(*,*) "5.2 GRAVANDO ARQUIVO DECLIVE.txt..."
!write(*,*)
!
!open(81,file='.\output\DECLIVE.txt')
!	write(81,'(A5,A9,I5)') 'ncols','',nC
!	write(81,'(A5,A9,I5)') 'nrows','',nL
!	write(81,'(A9,A5,F16.12)') 'xllcorner','',Xmin
!	write(81,'(A9,A5,F16.12)') 'yllcorner','',Ymin
!	write(81,'(A8,A6,F18.16)') 'cellsize','',cellsize
!	write(81,'(A12,A2,I5)') 'NODATA_value','',noDataFLAG
!	do lin=1,nL
!	    write(80,*) (DIST(lin,col),col=1,nC)
!		write(81,*) (WEIGHT(lin,col),col=1,nC)
!	end do
!close(80)
!close(81)


write(*,*)
write(*,*) 'CALCULO DO SDR POR LENHART ET AL. (2005) ENCERRADO!!!'
write(*,*)




deallocate(Wi,WiDi,DIST)


return

end