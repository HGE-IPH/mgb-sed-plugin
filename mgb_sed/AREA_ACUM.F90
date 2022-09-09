!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ ATUALIZAÇÃO: Jun 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA AREA_ACUM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ CALCULA A ÁREA ACUMULADA A MONTANTE DE CADA PIXEL
!@
!@
!@ *********************************************************************************************

Subroutine AREA_ACUM

use VAR

implicit none

integer codDIR, caminho, linaux, colaux, lin1, col1

	!@ ALOCA VARIÁVEIS PARA A AREA ACUMULADA
	allocate(AREAAcu(nL,nC))
	allocate(CONEXAO(nL,nC))


!@************************************************************************************
!@ CALCULA, PARA CADA PIXEL INTERNO A BACIA, AREA ACUMULADA A MONTANTE
!@************************************************************************************
write(*,*) "3.1 OBTENDO MATRIZ DE CONEXOES..."

CONEXAO = 0
AREAAcu = noDataFLAG
DO lin = 1,nL

	if (mod(lin,10)==0) then
		write(*,*) "lin = ",lin
	end if
	

	DO col = 1,nC
		codDIR = DIR(lin,col)
		IF (codDIR /= noDataFLAG) then
			linaux = lin + dlin(codDIR)
			colaux = col + dcol(codDIR)
			IF (linaux >= 1   .AND.   colaux >= 1   .AND.   linaux <= nL   .AND.   colaux <= nC) THEN
				CONEXAO(linaux,colaux) = CONEXAO(linaux,colaux) + 1
			ENDIF
			Xesq=Xmin+(col-1)*cellsize
			Xdir=Xesq+cellsize
			Yinf=Ymin+(nL-lin)*cellsize
			Ysup=Yinf+cellsize
			CALL Project_AREA(Xesq,Xdir,Ysup,Yinf,AreaCel)
			AREAAcu(lin,col) = AreaCel
		ENDIF
	ENDDO
ENDDO


write(*,*)
write(*,*) "3.3 OBTENDO AREAS ACUMULADAS..."
write(*,*)

DO lin = 1,nL
	if (mod(lin,10)==0) then
		write(*,*) "lin = ",lin
	end if	
	DO col = 1,nC
		codDIR = DIR(lin,col)
		IF (codDIR /= noDataFLAG   .AND.   CONEXAO(lin,col) == 0) THEN
			caminho = 0
			lin1 = lin
			col1 = col
			DO WHILE (caminho == 0)
				linaux = lin1 + dlin(codDIR)
				colaux = col1 + dcol(codDIR)
				IF (linaux < 1   .OR.   colaux < 1   .OR.   linaux > nL   .OR.   colaux > nC) THEN
					caminho = 1
				ELSE
					codDIR = DIR(linaux,colaux)
					IF (codDIR == noDataFLAG) THEN
					 caminho = 1
					ELSEIF (CONEXAO(linaux,colaux) > 1) THEN
						AREAAcu(linaux,colaux) = AREAAcu(linaux,colaux) + AREAAcu(lin1,col1)
						CONEXAO(linaux,colaux) = CONEXAO(linaux,colaux) - 1
						caminho = 1
					ELSE
						AREAAcu(linaux,colaux) = AREAAcu(linaux,colaux) + AREAAcu(lin1,col1)
						lin1 = linaux
						col1 = colaux
						caminho = 0
					ENDIF
				ENDIF
			ENDDO
		ENDIF
	ENDDO
ENDDO



write(*,*)
write(*,*) "3.3 GRAVANDO ARQUIVO AREA_ACU.txt..."
write(*,*)

open(81,file='.\output\AREA_ACU.txt')
	!@ -------- Cabeçalho
	write(81,'(A5,A9,I10)') 'ncols','',nC
	write(81,'(A5,A9,I10)') 'nrows','',nL
	write(81,'(A9,A5,F26.12)') 'xllcorner','',Xmin
	write(81,'(A9,A5,F26.12)') 'yllcorner','',Ymin
	write(81,'(A8,A6,F28.16)') 'cellsize','',cellsize
	write(81,'(A12,A2,I10)') 'NODATA_value','',noDataFLAG
	do lin=1,nL
		!@ -------- Escrevendo layer de area acumulada
		write(81,*) (AREAAcu(lin,col),col=1,nC) !@ (km2)
	end do
close(81)

write(*,*)
write(*,*) 'CALCULO AREA ACUMULADA ENCERRADO!!!'
write(*,*)


deallocate(CONEXAO)


return

end