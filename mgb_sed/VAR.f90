!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ ATUALIZA��O: Jun 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MODULO VAR >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ DEFINE AS VARI�VEIS DA ROTINA
!@
!@ A,B,C,D,E,F,G e H	= identifica��o da numeracao das direcoes
!@ ddaux				= vetor de direcoes auxiliar
!@ dlin e dcol			= vetores com a posicao relativa dos pixels vizinhos
!@ lin					= linha
!@ col					= coluna
!@ nL					= n�mero de linhas da matriz
!@ nC					= n�mero de colunas da matriz
!@ Xmin					= coordenada m�nima em x
!@ Ymin					= coordenada m�nima em y
!@ noDataFLAG			= var�vel com valor -9999 (no data do ArcGIS)
!@ MNT					= matriz com o MNT bruto
!@ AREAAcu	        	= matriz com as �reas acumuladas
!@ MINI					= matriz com as mini-bacia
!@ BLOCOS				= matriz com os blocos
!@ DIR					= matriz com as dire��es de fluxo do ArcGIS
!@    ____________________
!@	!      !      !      !
!@	!  32  !  64  ! 128  !
!@  !______!______!______!
!@	!      !      !      !
!@	!  16  !   X  !  1   !
!@	!______!______!______!
!@	!      !      !      !
!@	!   8  !   4  !  2   !
!@  !______!______!______!
!@



!@ REDE = matriz com a rede de drenagem (1=pixels da rede; 0=pixels externos a rede)
! Zp = cota do pixel
! Zr = cota do pixel da rede
! Zpaux = cota do pixel auxiliar
! nLaux = n�mero de linhas da matriz auxiliar
! nCaux = n�mero de colunas da matriz auxiliar
! Xminaux = coordenada m�nima em x auxiliar
! Yminaux = coordenada m�nima em y auxiliar
! cellsizeaux = tamanho do pixel (resolu��o em graus) auxiliar
!
! HANDval = matriz com as diferencas de cotas entre um pixel e o pixel da rede para o qual ele drena



!@
!@ *********************************************************************************************


Module VAR



integer(2) A, B, C, D, E, F, G, H

integer(2) :: dlin(128), dcol(128), ddaux(8)

character (len=14) :: TRASH

real(KIND=8) Xmin, Ymin, cellsize, TOL

integer(2) nL, nC, noDataFLAG, i, j, lin, col

integer(2) nLaux, nCaux

real(KIND=8) Xminaux, Yminaux, cellsizeaux

integer(2),allocatable :: MNT(:,:), DIR(:,:), MINI(:,:), BLOCOS(:,:), DRENAGEM(:,:), CONEXAO(:,:)
integer(4),allocatable :: CONT(:,:)

real(KIND=8),allocatable :: AREAAcu(:,:), LS2d(:,:), LS_acu(:,:)  !, SED(:,:)    !, Apixel(:,:)
!@ 02-10-2011
!real(KIND=8),allocatable :: L2d(:,:), S2d(:,:)
!@ 02-10-2011
real(KIND=8),allocatable :: WEIGHT(:,:), DIST(:,:), Wi(:,:), WiDi(:,:), SDR(:,:), SDRaux(:)

integer codBAC, codBLC, codSTR

real(KIND=8) Xesq, Xdir, Yinf, Ysup, AreaCel

integer respdist, resparea, respls

real(KIND=8) lado, diag, DeltaX, DeltaY, DeltaXY, CELLx, CELLy, CELL, CELLm

integer numMINI, numBLOCO

INTEGER,PARAMETER:: FILUSO=1	!ARQUIVO DOS PARAMETROS DA MUSLE ASSOCIADOS AO USO
INTEGER,PARAMETER:: FILTEX=2	!ARQUIVO DAS TEXTURAS DOS SOLOS ASSOCIADOS AO USO
CHARACTER (10) Nuso(12)			!NOMES DOS USOS DO SOLO (BLOCOS)
CHARACTER (10) CABE(300)		!CABECALHOS SEM IMPORTANCIA

real(KIND=8),allocatable :: Kusle(:),Cusle(:),Pusle(:),Rgros(:)
real(KIND=8),allocatable :: Mareia(:),Msilte(:),Margila(:),Morg(:),Mrocha(:)


end