!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ ATUALIZAÇÃO: Jun 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA ENTRADA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ ENTRADA DE DADOS PARA CÁLCULOS DOS PARÂMETROS DE MUSLE
!@
!@ *********************************************************************************************

Subroutine ENTRADA

use VAR

implicit none



!@ ***************** ESCOLHA DAS DISTÂNCIAS ENTRE PIXELS
	write(*,*) 
	write(*,*) "1.1 DEFININDO COMPRIMENTOS...:"
	write(*,*)
	write(*,*) "       Os comprimentos sao sempre calculados de pixel a pixel,"
	write(*,*) " podendo considerar:"
	write(*,*)
	write(*,*) " - distancia entre dois pixels vizinhos eh igual a resolucao"
	write(*,*) "   espacial local quando a direcao entre eles eh ortogonal,"
	write(*,*) "   ou igual a 2^(1/2) vezes a resolucao espacial quando a"
	write(*,*) "   direcao eh diagonal;"
	write(*,*)
	write(*,*) " - transformacao de distancias (Distance Transforms - DT), onde"
	write(*,*) "   sao adotados os valores 0,96194 e 1,36039 para passos nas"
	write(*,*) "   direcoes ortogonais e transversais, respectivamente, buscando"
	write(*,*) "   reduzir os erros de calculos devido a discretizacao, por pixels"
	write(*,*) "   quadrados, do trecho a ser calculado (Butt e Maragos, 1998)."
	write(*,*)
	write(*,*)
	write(*,*) "POR FAVOR, ESCOLHA...:"
	write(*,*)
	write(*,*) "(0) distancias padrao entre pixels:"
	write(*,*) "    lado=1, diagonal=1.414"
	write(*,*)
	write(*,*) "(1) distancias entre pixels com DistanceTransforms:"
	write(*,*) "    lado=0.96194, diagonal=1.36039"
	write(*,*)
	write(*,*) "Digite sua opcao: 0 ou 1"
	write(*,*)
	read(*,*) respdist
	write(*,*)

	if (respdist==0) then
		lado=1.0
		diag=1.414		!@ tamanho da diagonal em numero de lados (diagonal=1.414*lado)
	else
		lado=0.96194
		diag=1.36039	!@ tamanho da diagonal em numero de lados (diagonal=1.414*lado)
	end if



!@ ***************** ARQUIVO DE AREA ACUMULADA
	write(*,*) 
	write(*,*) "1.2 LAYER DE AREA ACUMULADA...:"
	write(*,*)
	write(*,*) "       O layer com as areas acumuladas para cada pixel eh necessário para"
	write(*,*) " o calculo do fator LS-2D e deve ser fornecido em km2. Caso jah exista este"
	write(*,*) " layer, o usuario poderah fornece-lo como entrada, caso contrario, nesta"
	write(*,*) " rotina eh possivel calcula-lo. Como o calculo deste layer pode consumir um"
	write(*,*) " tempo relativamente alto de processamento, sao fornecidas as 2 opcoes abaixo."
	write(*,*)
	write(*,*) " - opcao 0: se jah existe um layer de areas acumuladas calculado anteriormente."
	write(*,*) "   Nesta opcao serah lido o arquivo existente com o nome AREA.txt;"
	write(*,*)
	write(*,*) " - opcao 1: caso deseje calcular o layer de areas acumuladas. Nesta opcao"
	write(*,*) "   serah gerado o layer de areas acumuladas com o nome AREA.txt."
	write(*,*)
	write(*,*)
	write(*,*) "POR FAVOR, ESCOLHA...:"
	write(*,*)
	write(*,*) "(0) ler layer de areas acumuladas"
	write(*,*)
	write(*,*) "(1) calcula areas acumuladas"
	write(*,*)
	write(*,*) "Digite sua opcao: 0 ou 1"
	write(*,*)
	read(*,*) resparea
	write(*,*)



!@ ***************** ARQUIVO DO FATOR LS-2D
	write(*,*) 
	write(*,*) "1.3 ARQUIVO DO FATOR TOPOGRAFICO LS-2D...:"
	write(*,*)
	write(*,*) "       O calculo deste layer pode consumir um tempo relativamente alto"
	write(*,*) "  de processamento e armazenamento. Caso jah exista este layer, o usuario"
	write(*,*) " poderah fornece-lo como entrada, caso contrario, ele eh calculado nesta"
	write(*,*) " rotina. Por favor, escolha entre as 2 opcoes abaixo."
	write(*,*)
	write(*,*) " - opcao 0: se jah existe um layer de fator topográfico LS-2D calculado"
	write(*,*) "    anteriormente. Nesta opcao serah lido o arquivo existente com o nome"
	write(*,*) "    LS_2D_DG1996.txt;"
	write(*,*)
	write(*,*) " - opcao 1: caso deseje calcular o layer de areas acumuladas. Nesta opcao"
	write(*,*) "   serah gerado o layer do fator topografico LS-2D com o nome LS_2D_DG1996.txt."
	write(*,*)
	write(*,*)
	write(*,*) "POR FAVOR, ESCOLHA...:"
	write(*,*)
	write(*,*) "(0) ler layer do fator topografico LS-2D"
	write(*,*)
	write(*,*) "(1) calcula fator topografico LS-2D"
	write(*,*)
	write(*,*) "Digite sua opcao: 0 ou 1"
	write(*,*)
	read(*,*) respls
	write(*,*)



!@ >>>>>>>>>>>>>>>>>>>>>>> LEITURA DOS ARQUIVOS DE ENTRADA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

write(*,*)
write(*,*) "2. LENDO ARQUIVOS DE ENTRADA..."


!@ ***************** LEITURA DADOS DO MNT BRUTO - MNT.txt

	write(*,*) 
	write(*,*) "2.1 LENDO LAYER DO MNT BRUTO..."
	write(*,*)
	write(*,*) "	- ARQUIVO MNT.txt"
	write(*,*)

	open(10,FILE='.\input\MNT.txt',STATUS='UNKNOWN')

	read(10,*) TRASH,nC
	read(10,*) TRASH,nL
	read(10,*) TRASH,Xmin
	read(10,*) TRASH,Ymin
	read(10,*) TRASH,cellsize
	read(10,*) TRASH,noDataFLAG

	TOL = cellsize/1000

	write(*,*) "******************************************************************"
	write(*,*)
	write(*,*) "ATENCAO - Confira abaixo!!"
	write(*,*)
	write(*,*) "	Resolucao espacial da imagem em GRAUS = ", cellsize
	write(*,*)
	write(*,*) "	Numero de linhas  = ",nL
	write(*,*)
	write(*,*) "	Numero de Colunas = ",nC
	write(*,*)
	write(*,*) "	Coordenada Xmin (em GRAUS) = ",Xmin
	write(*,*)
	write(*,*) "	Coordenada Ymin (em GRAUS) = ",Ymin
	write(*,*)
	write(*,*) "******************************************************************"
	write(*,*)
	write(*,*) "(tecle enter)"
	read(*,*)


	!@ ALOCA VARIÁVEIS DO MNT:
	allocate(MNT(nL,nC))

	!@ Le dados:
	do i=1,nL
		read(10,*) (MNT(i,j),j=1,nC)
	enddo
	close(10)




!@ ***************** LEITURA DADOS DO ARQUIVO DE DIREÇÕES DE FLUXO - DIR.txt

	write(*,*) 
	write(*,*) "2.2 LENDO LAYER DE DIRECOES DE FLUXO..."
	write(*,*)
	write(*,*) "	- ARQUIVO DIR.txt"
	write(*,*)
	write(*,*)

	open(10,FILE='.\input\DIR.txt',STATUS='UNKNOWN')
	
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
		return
	endif

	!@ ALOCA VARIÁVEIS DO DIR:
	allocate(DIR(nL,nC))

	!@ Le dados:
	do i=1,nL
		read(10,*) (DIR(i,j),j=1,nC)
	enddo
	close(10)




!@ ***************** LEITURA DADOS DAS MINIBACIAS - MINI.txt

	write(*,*) 
	write(*,*) "2.3 LENDO LAYER DE MINIBACIAS..."
	write(*,*)
	write(*,*) "	- ARQUIVO MINI.txt"
	write(*,*)
	write(*,*)

	open(10,FILE='.\input\MINI.txt',STATUS='UNKNOWN')

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
		return
	endif
	
	!@ ALOCA VARIÁVEIS DA MINIBACIA:
	allocate(MINI(0:nL+1,0:nC+1))
	MINI=-9999

	!@ Le dados:
	do i=1,nL
		read(10,*) (MINI(i,j),j=1,nC)
	enddo
	close(10)




!@ ***************** LEITURA DADOS DOS BLOCOS - BLOCOS.txt

	write(*,*) 
	write(*,*) "2.4 LENDO LAYER DE BLOCOS..."
	write(*,*)
	write(*,*) "	- ARQUIVO BLOCO.txt"
	write(*,*)
	write(*,*)

	open(10,FILE='.\input\BLOCO.txt',STATUS='UNKNOWN')

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
		return
	endif

	!@ ALOCA VARIÁVEIS DOS BLOCOS:
	allocate(BLOCOS(nL,nC))

	!@ Le dados:
	do i=1,nL
		read(10,*) (BLOCOS(i,j),j=1,nC)
	enddo
	close(10)




!@ ***************** LEITURA DADOS DA REDE - REDE.txt

	write(*,*) 
	write(*,*) "2.4 LENDO LAYER DA REDE DE DRENAGEM..."
	write(*,*)
	write(*,*) "	- ARQUIVO REDE.txt"
	write(*,*)
	write(*,*)

	open(10,FILE='.\input\REDE.txt',STATUS='UNKNOWN')

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
		return
	endif

	!@ ALOCA VARIÁVEIS DOS BLOCOS:
	allocate(DRENAGEM(nL,nC))

	!@ Le dados:
	do i=1,nL
		read(10,*) (DRENAGEM(i,j),j=1,nC)
	enddo
	close(10)





!@ ***************** OBTENDO VARIÁVEIS

	!@ MAIOR CÓDIGO DE UMA MINIBACIA = NÚMERO DE MINIBACIAS
	numMINI=maxval(MINI)
	write(*,*)
	write(*,*) "Numero de minibacias      = ",numMINI
	write(*,*)

	!@ MAIOR CÓDIGO DE UM BLOCO = NÚMERO DE BLOCOS
	numBLOCO=maxval(BLOCOS)
	write(*,*)
	write(*,*) "Numero de blocos          = ",numBLOCO
	write(*,*)




!@ ***************** ALOCANDO VARIÁVEIS

!	!@ ALOCA VARIÁVEIS PARA A AREA ACUMULADA
!	allocate(AREAAcu(nL,nC))
!	allocate(CONEXAO(nL,nC))

!	!@ ALOCA VARIÁVEIS DE AREAS DE BLOCOS
!	allocate(Apixel(numMINI,numBLOCO))

!	!@ ALOCA VARIÁVEIS DO LS-2D:
!	allocate(LS2d(nL,nC))
!	LS2d = noDataFLAG

	!@ ALOCA VARIÁVEIS DO LS médio
	allocate(LS_acu(numMINI,numBLOCO))
	LS_acu = 0.0

	!@ ALOCA VARIÁVEIS DO CONTADOR DE PIXELS
	allocate(CONT(numMINI,numBLOCO))
	
!	!@ ALOCA VARIÁVEIS DO SDR
!	allocate(Wi(numMINI,numBLOCO),WiDi(numMINI,numBLOCO),SDR(numMINI,numBLOCO))
!	allocate(WEIGHT(nL,nC))
!	Wi = 0.
!	WiDi = 0.
!	SDR = 0.
!	WEIGHT = noDataFLAG

!	!@ ALOCA VARIÁVEIS DE SEDIMENTOS
!	allocate(SED(numMINI,numBLOCO))
!	SED = 0.0


return

end