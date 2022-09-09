!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ ATUALIZAÇÃO: Jun 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA LS_2D >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ CALCULA O FATOR LS BIDIMENSIONAL DE ACORDO COM DESMET E GOVERS (1996)
!@
!@
!@ *********************************************************************************************

Subroutine LS_2D_DG1996

use VAR

implicit none

real(KIND=8) Zx, Zy, Sper, Sgrau, AreaM, SLOPE, COEF, CompL, fatorX
real(KIND=8) lado1,diag1
integer(2) linSeg, colSeg, dirSeg

!@ ALOCA VARIÁVEIS DO LS-2D:
allocate(LS2d(nL,nC))   !,WEIGHT(nL,nC))
LS2d = noDataFLAG
!WEIGHT = noDataFLAG


!!@ 02-10-2011
!allocate(L2d(nL,nC))
!L2d = noDataFLAG
!allocate(S2d(nL,nC))
!S2d = noDataFLAG
!!@ 02-10-2011


write(*,*)
write(*,*) "4.1 CALCULANDO O FATOR LS-2D..."
write(*,*)


!Apixel = 0.0
lado1 = 1.
diag1 = sqrt(2.)

CONT = 0.
do lin = 1,nL

	if (mod(lin,10)==0) then
		write(*,*) "lin: ",lin
	end if

	do col = 1,nC
		
		!@ Exclindo pixel com valor NoData dos cálculos
		if (MINI(lin,col)==noDataFLAG .OR. MNT(lin,col) == noDataFLAG) cycle    !then

			codBAC = MINI(lin,col)		!@ obtem código da minibacia
			codBLC = BLOCOS(lin,col)	!@ obtem código do bloco

		!@ ------------	OBTENDO DIMENSÕES DO PIXEL
			!@ Calcula as coordenadas dos vertices do pixel atual
			Xesq=Xmin+(col-1)*cellsize
			Xdir=Xesq+cellsize
			Yinf=Ymin+(nL-lin)*cellsize
			Ysup=Yinf+cellsize

			DeltaX  = 0.0
			DeltaY  = 0.0
			DeltaXY = 0.0
			CELLx   = 0.0
			CELLy   = 0.0
			CELLm   = 0.0
			CELL    = 0.0

			!@ Chama rotina para cálculo das dimensões do pixel
			call Project_COMP(Xesq,Xdir,Ysup,Yinf,DeltaX,DeltaY,DeltaXY,CELLx,CELLy,lado1,diag1)

			CELLm = (CELLx + CELLy)/2.	!@ Largura média do pixel (m)
		!@ **********************************************************************


		!@ ------------	OBTENDO ÁREA DO PIXEL
		!@ chama rotina para calcular area do pixel
			AreaCel = 0.0	!@ (km^2)
			
			call Project_AREA(Xesq,Xdir,Ysup,Yinf,AreaCel)

!			!@ Acumulando área de cada bloco de cada mini-bacia
!			!if (lin >= 5120) then
!			!write(*,*) "lin      = ", lin
!			!write(*,*) 'BAC, BLC = ', codBAC, codBLC
!			!endif
!			Apixel(codBAC,codBLC) = Apixel(codBAC,codBLC) + AreaCel	!@ (km^2)
		!@ **********************************************************************


		!@ ------------	OBTENDO ÁREA A MONTANTE DO PIXEL
			AreaM = (AREAAcu(lin,col) - AreaCel)*1000000.	!@ (m^2)
			if (AreaM < 0.) AreaM = 0.
		!@ **********************************************************************


		!@ ------------	OBTENDO DERIVADAS PARCIAIS POR DIFERENÇAS FINITAS
		!@ FONTE: Wilson e Gallant 2000. Terrain Analysis: Principles and Applications.
			call DERIVA_WG2000(Zx,Zy)

		!@ FONTE: ArcGIS 9.2
			!call DERIVA_ArcGIS(Zx,Zy)
		!@ **********************************************************************


		!@ ------------	OBTENDO O DECLIVE DO PIXEL
		!@ FONTE: Wilson e Gallant 2000. Terrain Analysis: Principles and Applications.
			Sper = sqrt(Zx**2. + Zy**2.)*100.0			!@ Declive em percentagem
			Sgrau = (ATAN(Sper/100.0))*180.0/3.141592	!@ Declive em graus
!			WEIGHT(lin,col) = 1./Sper                   !@ Matriz de peso para calculo do SDR
		!@ **********************************************************************


		!@ ------------	OBTENDO O FATOR DE DECLIVIDADE
		!@ FONTE: Wischmeier & Smith (1978)
			SLOPE = 65.4*(sin(Sgrau*3.141592/180.0)**2.) + 4.56*sin(Sgrau*3.141592/180.0) + 0.0654
		!@ **********************************************************************


		!@ ------------	OBTENDO O EXPOENTE DO COMPRIMENTO
			if (Sper >= 5.) then
				COEF = 0.5
			elseif (Sper < 5. .AND. Sper >= 3.) then
				COEF = 0.4
			elseif (Sper < 3. .AND. Sper >= 1.) then
				COEF = 0.3
			else
				COEF = 0.2
			endif
		!@ **********************************************************************


		!@ ------------	OBTENDO O TIPO DE DRENAGEM DO PIXEL
			dirSeg=dir(lin,col)
			linSeg=lin+dlin(dirSeg) !@ identifica linha do pixel de jusante
			colSeg=col+dcol(dirSeg) !@ identifica coluna do pixel de jusante
			
			!@ determina posicao relativa ao pixel seguinte
			if ((lin==linSeg).OR.(col==colSeg)) then
				if (lin==linSeg) then
					fatorX = 1.			!@ drenagem na horizontal
					CELL = CELLx*1000.	!@ dimensão horizontal da célula (m)
				else
					fatorX = 1.			!@ drenagem na vertical
					CELL = CELLy*1000.	!@ dimensão vertical da célula (m)
				end if
			else
				fatorX = sqrt(2.)		!@ drenagem na diagonal
				CELL = CELLm*1000.		!@ dimensão diagonal da célula (m)
			end if
		!@ **********************************************************************


		!@ ------------	OBTENDO O FATOR L DO PIXEL
			CompL = ((AREAAcu(lin,col)*1000000.)**(COEF+1.) - AreaM**(COEF+1.))
			CompL = CompL / ((CELL**(COEF+2.))*(fatorX**(COEF))*(22.13**(COEF)))
		!@ **********************************************************************


		!@ ------------	OBTENDO O FATOR LS DO PIXEL
!		    !@ 02-10-2011
!		    L2d(lin,col) = min(CompL,350.)
!		    S2d(lin,col) = SLOPE
!		    !@ 02-10-2011
		    
!			LS2d(lin,col) = min(CompL,500.)*SLOPE   !@ LIMITA O FATOR L A 500M
            LS2d(lin,col) = CompL*SLOPE
		!@ **********************************************************************


		!@ ------------	ACUMULANDO O FATOR LS DO BLOCO DE CADA MINIBACIA
			LS_acu(codBAC,codBLC) = LS_acu(codBAC,codBLC) + LS2d(lin,col)
!			!@ 02-10-2011
!			LS_acu(codBAC,codBLC) = LS_acu(codBAC,codBLC) + S2d(lin,col)
!			!@ 02-10-2011
		!@ **********************************************************************

			CONT(codBAC,codBLC) = CONT(codBAC,codBLC) + 1
	enddo
enddo


!write(*,*)
!write(*,*) "4.2 GRAVANDO ARQUIVO LS_2D.txt..."
!write(*,*)
!
!open(80,file='.\output\LS_2D_DG1996 - 350.txt')
!	write(80,'(A5,A9,I5)') 'ncols','',nC
!	write(80,'(A5,A9,I5)') 'nrows','',nL
!	write(80,'(A9,A5,F16.12)') 'xllcorner','',Xmin
!	write(80,'(A9,A5,F16.12)') 'yllcorner','',Ymin
!	write(80,'(A8,A6,F18.16)') 'cellsize','',cellsize
!	write(80,'(A12,A2,I5)') 'NODATA_value','',noDataFLAG
!	do lin=1,nL
!		write(80,*) (LS2d(lin,col),col=1,nC)
!	end do
!close(80)

!!@ 02-10-2011
!open(81,file='.\output\L_2D_DG1996.txt')
!	write(81,'(A5,A9,I5)') 'ncols','',nC
!	write(81,'(A5,A9,I5)') 'nrows','',nL
!	write(81,'(A9,A5,F16.12)') 'xllcorner','',Xmin
!	write(81,'(A9,A5,F16.12)') 'yllcorner','',Ymin
!	write(81,'(A8,A6,F18.16)') 'cellsize','',cellsize
!	write(81,'(A12,A2,I5)') 'NODATA_value','',noDataFLAG
!	do lin=1,nL
!    	write(81,*) (L2d(lin,col),col=1,nC)
!	end do
!close(81)
!
!open(82,file='.\output\S_2D_DG1996.txt')
!	write(82,'(A5,A9,I5)') 'ncols','',nC
!	write(82,'(A5,A9,I5)') 'nrows','',nL
!	write(82,'(A9,A5,F16.12)') 'xllcorner','',Xmin
!	write(82,'(A9,A5,F16.12)') 'yllcorner','',Ymin
!	write(82,'(A8,A6,F18.16)') 'cellsize','',cellsize
!	write(82,'(A12,A2,I5)') 'NODATA_value','',noDataFLAG
!	do lin=1,nL
!		write(82,*) (S2d(lin,col),col=1,nC)
!	end do
!close(82)
!!@ 02-10-2011


write(*,*)
write(*,*) 'CALCULO DO LS-2D (Desmet e Govers, 1996) ENCERRADO!!!'
write(*,*)


deallocate(LS2d,AREAAcu)
!deallocate(L2d,AREAAcu)
!deallocate(S2d,AREAAcu)


return

end