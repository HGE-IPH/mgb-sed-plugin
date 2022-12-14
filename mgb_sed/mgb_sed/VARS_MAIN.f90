	MODULE VARS_MAIN
	!DECLARAÇÃO DE VARIÁVEIS DO PROGRAMA PRINCIPAL
	IMPLICIT NONE
	SAVE
	INTEGER,PARAMETER:: NUP=12 !NUMERO MAXIMO DE USOS !RP
	INTEGER,PARAMETER:: NBP=300 !NUMERO MAXIMO DE SUB-BACIAS
	INTEGER,PARAMETER:: NCLIP=300 !NUMERO MAXIMO DE POSTOS CLIMATOLÓGICOS !RP
	INTEGER,PARAMETER:: NUMUSK=150 !NUMERO MAXIMO DE SUBTRECHOS DE PROPAGACAO MUSKINGUM CUNGE !RP
	INTEGER,PARAMETER:: NFP=3 !NUMERO DE FUNÇÕES OBJETIVO PARA CADA POSTO
	INTEGER,PARAMETER:: NTMUP=20 !NUMERO DE PONTOS DA TABELA DO MUSKINGUN-CUNGE NAO LINEAR
	
	!NUMEROS 	 
	INTEGER NC,NU !CELULAS, USOS
	INTEGER NT !INTERVALOS DE TEMPO
	INTEGER NCLI !NUMERO DE POSTOS CLIMATOLOGICOS
	INTEGER NB !NUMERO DE BACIAS
	INTEGER NOBS,IQOBS(NBP) !NUMERO DE SERIES DE VAZAO OBSERVADA, numero das células com dados observados
	INTEGER NUMHIDG,IHIDG(NBP) !NUMERO DE HIDROGRAMAS QUE SE DESEJA GRAVAR, NUMERO DAS CÉLULAS EM QUE SE DESEJA OS HIDROGRAMAS
	INTEGER NUMSUBST,ISUBST(NBP) !NUMERO DAS CÉLULAS ONDE A Q CALCULADA DEVE SER SUBSTITUÍDA PELA LIDA DE UM ARQUIVO
	!INTEGER,ALLOCATABLE:: IOBSAUX(:)!NUMERO DE SERIES DE VAZAO OBSERVADA

	!CONTADORES
	INTEGER IT !INTERVALO DE TEMPO

	!NUMEROS DOS ARQUIVOS DE DADOS CLIMATOLÓGICOS DIÁRIOS
	INTEGER,ALLOCATABLE:: NARQ(:)

	!-------VARIAVEIS COM NUMEROS DOS ARQUIVOS------------
	!ARQUIVOS DE ENTRADA
	INTEGER,PARAMETER:: FILFIX=1 !ARQUIVO DOS PARAMETROS FIXOS
	INTEGER,PARAMETER:: FILVAR=2 !ARQUIVO DOS PARAMETROS VARIAVEIS
	INTEGER,PARAMETER:: FILUSO=3 !ARQUIVO DOS PARAMETROS ASSOCIADOS AO USO
	INTEGER,PARAMETER:: FILHIG=4 !ARQUIVO DAS CELULAS
	INTEGER,PARAMETER:: FILPLU=5 !ARQUIVO DOS DADOS DE CHUVA
	INTEGER,PARAMETER:: FILMED=6 !ARQUIVO DOS DADOS CLIMATOLÓGICOS MEDIOS MENSAIS
	INTEGER,PARAMETER:: FILOBS=7 !ARQUIVO DE VAZOES OBSERVADAS
	INTEGER,PARAMETER:: FILLIM=8 !ARQUIVO DE LIMITES DOS PARÂMETROS
	INTEGER,PARAMETER:: FILPREV=9 !ARQUIVO DE DADOS DE PREVISÃO DE CHUVA
	INTEGER,PARAMETER:: FILSUBS=10 !ARQUIVO DE HIDROGRAMAS LIDOS PARA SUBST. CALCULADOS
	INTEGER,PARAMETER:: FILOUTROS=11 !ARQUIVO DE HIDROGRAMAS OBSERVADOS NO EIXO DO SAO FRANCISCO
	INTEGER,PARAMETER:: FILCALIB=12 !ARQUIVO DE PARAMETROS DE CALIBRACAO



	!ARQUIVOS DE SAÍDA
	INTEGER,PARAMETER:: FILHID=21 !ARQUIVO DOS HIDROGRAMAS DAS SUB-BACIAS
	INTEGER,PARAMETER:: FILPRP=22 !ARQUIVO DE HIDROGRAMA COM PROPORÇÕES
	INTEGER,PARAMETER:: FILAJU=23 !ARQUIVO DE AJUSTE
	INTEGER,PARAMETER:: FILBAC=24 !ARQUIVO DE SAIDA DE DADOS MEDIOS DAS BACIAS 
	INTEGER,PARAMETER:: FILSOL=25 !ARQUIVO DO SOLO (1 BLOCO)
	INTEGER,PARAMETER:: FILSOL2=25001 !ARQUIVO DO SOLO (1 BLOCO)
	INTEGER,PARAMETER:: FILEVO=26 !ARQUIVO DE PARÂMETROS DO EVOLUTION
	INTEGER,PARAMETER:: FILBOM=27 !ARQUIVO CONTENDO OS MELHORES PARÂMETROS ATÉ O MOMENTO
	!-----------------------------------------------------

	INTEGER ICALIB
	REAL,ALLOCATABLE:: R2(:),ERRV(:),R2L(:) !COEF. NASH, ERRO NOS VOLUMES, COEF NASH LOGARITMOS
	REAL VFO(NFP) !VETOR COM VALORES DAS FUNÇÕES OBJETIVO


	!____________________________________________________
	!----TIPO E DIMENSIONAMENTO DAS VARIAVEIS------------
	!REAL,ALLOCATABLE:: CORR(:) !CORRIGE O VALOR DOS PARAMETROS
	REAL,ALLOCATABLE:: QR(:,:),QOBS(:,:) !VAZÃO NOS EXUTORIOS DAS BACIAS
	REAL,ALLOCATABLE:: QLIDO(:,:) !VAZÃO EM ALGUNS PONTOS PARA SUBSTITUIÇÃO DE VALORES CALCULADOS
	REAL,ALLOCATABLE:: QMANTEIGA(:),QSAOROM(:),QSAOFRA(:),QPMCRUZ(:),QMANGA(:),QCARI(:)
	REAL,ALLOCATABLE:: QBJLAPA(:),QGAMEL(:),QMORP(:)
	INTEGER,ALLOCATABLE:: KCB(:)!NUMERO DE CELULAS EM CADA SUB-BACIA 
	!VALORES MEDIOS POR SUB-BACIA
	REAL,ALLOCATABLE:: DSUM(:,:),DINM(:,:),DBAM(:,:)
	REAL,ALLOCATABLE:: SIM(:,:),EM(:,:),WBM(:,:),PM2(:,:)
	!VAZÃO DE ACORDO COM A ORIGEM EM CADA TEMPO PARA UMA DAS CÉLULAS
	REAL,ALLOCATABLE:: QB(:),QBI(:),QBIS(:)
	INTEGER IDIA,IMES,IANO,IDINI !DATA E DIA DO CALENDARIO JULIANO
	REAL,ALLOCATABLE:: QESP(:)!VAZÃO ESPECIFICA DE BASE (M3/S/KM2)
	!ALBEDO, INDICE DE AREA FOLIAR E ALTURA MEDIA, RESISTENCIA SUPERFICIAL
	REAL,ALLOCATABLE:: ALB(:,:),RIAF(:,:),Z(:,:),RS(:,:)
	!NUMERO DE SUBTRECHOS E INTERVALO DE TEMPO DE CALCULO MUSKINGUM-CUNGE
	INTEGER,ALLOCATABLE:: NSUBT(:)
	REAL,ALLOCATABLE:: DT(:)
	REAL,ALLOCATABLE:: CEL(:),TIND(:)!CELERIDADE NO RIO E TEMPO DE CONCENTRACAO DA CELULA
	!NOMES DOS ARQUIVOS DE DADOS CLIMATOLÓGICOS DIÁRIOS
	CHARACTER (20) ARQCLI(NCLIP)
	!NOME DO ARQUIVO DE DADOS CLIMATOLÓGICOS MÉDIOS MENSAIS
	CHARACTER (20) ACLIMED
	!NOME DO ARQUIVO DAS VAZOES OBSERVADAS
	CHARACTER (20) ARQOBS
	!NOME DO ARQUIVO DAS VAZÕES LIDAS QUE SUBSTITUIRÃO AS CALCULADAS EM ALGUMAS CÉLULAS
	CHARACTER (20) ARQSUBST
	!PARAMETROS RELACIONADOS AO USO
	REAL,ALLOCATABLE:: WM(:,:),B(:,:),KINS(:,:),KBAS(:,:)
	REAL,ALLOCATABLE:: WMOLD(:,:),BOLD(:,:),KBOLD(:,:),KIOLD(:,:)
	REAL,ALLOCATABLE:: PLAM(:,:),CAP(:,:),WC(:,:)
	REAL,ALLOCATABLE:: CI(:),CB(:),CS(:)!PARAMETROS DA PROPAGAÇÃO NA CÉLULA
	REAL,ALLOCATABLE:: CIOLD(:),CSOLD(:)!CÓPIA DOS PARAMETROS DA PROPAGAÇÃO NA CÉLULA
	CHARACTER (10) AUSO(NUP)!NOMES DOS USOS DO SOLO (BLOCOS)
	!-----------------------------------------------------
	INTEGER,ALLOCATABLE:: IEXUT(:)!INDICA CELULA DO EXUTORIO DA BACIA
	REAL,ALLOCATABLE:: BRIO(:)!LARGURA DO RIO
	REAL BRX !LARGURA DO RIO (VARIÁVEL AUXILIAR) 
	REAL,ALLOCATABLE:: PM(:) !CHUVA MEDIA E CHUVA NO INTERVALO NA CELULA
	INTEGER,ALLOCATABLE:: IBAC(:),CELJUS(:)!BACIA,CELULA DE JUSANTE
	REAL,ALLOCATABLE::  HCEL(:),LCEL(:) ! Declividade e comprimento do rio mais longo
	REAL,ALLOCATABLE:: X(:),Y(:)!COORDENADAS DO CENTRO DA CELULA
	!ÁREA DA CÉLULA, AREA DRENADA EM CELULAS, COMPRIMENTO DO RIO E DECLIVIDADE DO RIO
	REAL,ALLOCATABLE:: ACEL(:),ACUR(:),SRIO(:),DECL(:)
	!PROPORCAO DE USOS NA CELULA
	REAL,ALLOCATABLE:: PUSO(:,:)
	INTEGER,ALLOCATABLE:: ICBOM(:)!NUMERO DO POSTO CLIMATOLOGICO MAIS PROXIMO DA CELULA
	REAL,ALLOCATABLE:: QREF(:) !VAZAO DE REFERENCIA
	REAL QRX !VAZAO DE REFERENCIA (VARIAVEL AUXILIAR)
	!TEMP., UMIDADE, VENTO, INSOL. DIARIAS EM CADA POSTO
	REAL,ALLOCATABLE:: TD(:,:),UD(:,:),VD(:,:),SOLD(:,:),PAD(:,:) 
	!-------------------------------------------------
	REAL,ALLOCATABLE:: PMB(:)!CHUVA MEDIA POR BACIA
	INTEGER,ALLOCATABLE:: KPM(:) !VARIAVEL AUXILIAR
	!TEMP., UMIDADE, VENTO, INSOL., PRESSAO MEDIAS MENSAIS
	REAL,ALLOCATABLE:: TAMM(:,:),URMM(:,:),VVMM(:,:),PAMM(:,:)
	REAL,ALLOCATABLE:: SOLMM(:,:)
	REAL,ALLOCATABLE:: XYC(:,:)!COORDENADAS X E Y DOS POSTOS CLIMATOLÓGICOS
	!-----------------------------------------------------
	!CONTADORES
	INTEGER IC
	!CABECALHOS SEM IMPORTANCIA
	CHARACTER (10) CABE(300)
	!LAMINA INTERCEPTADA NA SUPERFICIE
	REAL,ALLOCATABLE:: SI(:,:)

	REAL XLADO,DH,ARX !COMPRIMENTO DO LADO DA CÉLULA, DIFERENÇA DE ALTITUDE, ÁREA (VARIÁVEL AUXILIAR)
	REAL,ALLOCATABLE:: QRG(:,:)	!ARMAZENA HIDROGRAMAS ONDE SE DESEJA GRAVAR
	REAL,ALLOCATABLE:: QRB(:,:)	!ARMAZENA HIDROGRAMAS NOS EXUTÓRIOS DAS SUB-BACIAS
	REAL,ALLOCATABLE:: QM1(:),QJ1(:),QM2(:),QJ2(:) !VAZOES A MONTANTE E A JUSANTE NA CÉLULA i NO TEMPO 1 E 2
	REAL,ALLOCATABLE:: QCEL1(:),QCEL2(:) !VAZÕES ORIGINADAS NA CELULA NOS INSTANTES t E t+1 NA CÉLULA i
	REAL,ALLOCATABLE:: PMB2(:),PMI2(:),PMS2(:),PJB2(:),PJI2(:),PJS2(:) !PROPORÇÕES DE ORIGEM DAS VAZÕES NO RIO
	REAL,ALLOCATABLE:: VRB(:),VRI(:),VRS(:) !VOLUMES DAS PROPORÇÕES NO TRECHO
	REAL,ALLOCATABLE:: ET(:,:) 	!EVAPOTRASPIRACAO TOTAL
	REAL,ALLOCATABLE:: CAF(:,:)		!FLUXO CAPILAR ASCENDENTE
	REAL,ALLOCATABLE:: W(:,:) 	!LAMINA DE AGUA NO SOLO
	REAL,ALLOCATABLE:: QBAS(:),QINT(:),QSUP(:)		!VAZAO NA CELULA
	REAL,ALLOCATABLE:: VBAS(:),VINT(:),VSUP(:)		!VOLUME NA CELULA
	REAL,ALLOCATABLE:: TONTEM(:)		!TEMPERATURA NO DIA ANTERIOR
	REAL,ALLOCATABLE:: P(:) !CHUVA NO INTERVALO NA CELULA
	REAL,ALLOCATABLE:: TA(:),UR(:),VV(:),SOL(:),PA(:)	!TEMP., UMIDADE, VENTO, INSOL., PRESSAO

	INTEGER HORAINI !HORA EEM QUE INICIA A SIMULAÇÃO NO DIA INICIAL 
	INTEGER JDIA !DIA DO CALENDÁRIO JULIANO
	CHARACTER (40) TITULO !VARIÁVEL ALFA PARA LER TÍTULOS DOS ARQUIVOS DE ENTRADA
	REAL DTP !INTERVALO DE TEMPO DE CÁLCULO PRINCIPAL (DADO EM SEGUNDOS NO ARQUIVO DE ENTRADA)

	REAL PEVAPO !PARAMETRO SUJO QUE MULTIPLICA A EVAPORAÇÃO

	!VARS DA SUB RADIAÇÃO
	REAL GLAT !LATITUDE EM GRAUS DECIMAIS ( - SUL)
	REAL SOLX,T1,TAR2 !INSOLAÇÃO (H) E TEMPERATURAS (C)
	REAL ALBX !ALBEDO
	REAL RLX !RADIAÇÃO LIQUIDA RESULTANTE (MJ/m2/dia)
	REAL RDIA !DIA JULIANO (real)
	REAL URX
	REAL,PARAMETER:: MESP=1000.0 !MASSA ESP. DA ÁGUA (KG/M3)
	REAL,PARAMETER:: PI=3.141592 !pi
	REAL,PARAMETER:: STEBOL=4.903E-9 !CONST. STEFAN BOLTZMANN
	!VARIAVEIS INTERNAS DA SUBROTINA
	REAL CLATE !calor latente de vaporização
	REAL SDECL,SANG !DECLINAÇÃO SOLAR,ANGULO AO NASCER
	REAL HIM !DURACAO MÁXIMA DA INSOLAÇÃO (HORAS)
	REAL DR !DISTANCIA REL. TERRA - SOL
	REAL GS,STO !FLUXO DE CALOR P/ SOLO, RAD. TOPO ATM.
	REAL SSUP,SN !FLUXO DE CALOR
	REAL ED,ES !PRESSÕES DE VAPOR (REAL, SATURAÇAO)
	REAL SLONG !RADIAÇÃO DE ONDAS LONGAS


	!VARS DA SUB CÉLULA
		!VARIAVEIS AUXILIARES
	REAL WX,WMX,PAX,TAX,VVX,BX,KIX,KBX,RIAFX,ZX,EX,PX,SIX
	REAL XL,T2,RSX,CAPX,ETX,WCX,VBX,VIX,VSX,PPU,WXB,SIB,ETB,DCAP


	!VARS DA SUB EVAPO

	REAL D !DEFICIT DE PRESSAO DE VAPOR
	REAL DEL !DERIVADA DA FUNÇAO ES x T
	REAL GAMA !CONSTANTE PSICROMÉTRICA
	REAL MSPA !MASSA ESPECÍFICA DO AR
	REAL RUG !RUGOSIDADE
	REAL SF !IDEM
	REAL SIL !MAXIMA LAMINA INTERCEPTADA
	REAL EIX !LAMINA EVAPORADA DA INTERCEPTAÇÃO (POTENCIAL)
	REAL,ALLOCATABLE:: EVQ(:) !EVAPORAÇÃO DIRETA DAS SUPERFÍCIES LÍQUIDAS (BLOCO ÁGUA) DA CÉLULA EM M3/S
	REAL REIX !LAMINA EVAPORADA DA INTERCEPTAÇÃO (REAL)
	REAL FDE !FRAÇÃO DE DEMANDA EVAPORATIVA DISPONIVEL
	REAL,PARAMETER:: CSPA=0.001013 	!calor esp. do ar umido MJ/kg/C
	REAL RA !RESISTÊNCIA AERODINÂMICA
	REAL WL !LIMITE DE UMIDADE ACIMA DO QUAL A UMIDADE DO SOLO NAO RESTRINGE A EVAPOTRANSPIRAÇAO (SHUTTLEWORTH)
	REAL WPM !PONTO DE MURCHA (SHUTTLEWORTH)
	REAL F4 !FATOR DE CORREÇÃO DA RESISTÊCNIA SUPERFICIAL EM FUNÇÃO DO DEFICIT DE UMIDADE DO SOLO


	INTEGER IB,IU !CONTADOR BACIA CONTADOR USO

	REAL,ALLOCATABLE:: QPREV(:,:,:)
	REAL,ALLOCATABLE:: WPREV(:,:) 	!LAMINA DE AGUA NO SOLO NO INÍCIO DO CICLO DE PREVISÃO
	REAL,ALLOCATABLE:: VBASPREV(:),VINTPREV(:),VSUPPREV(:) !VOLUME NA CELULA NO INÍCIO DO LOOP PREVISÃO
	REAL,ALLOCATABLE:: TAPREV(:) !TEMPERATURA NO DIA DE INÍCIO DO LOOP DE PREVISÃO
	REAL,ALLOCATABLE:: QM2PREV(:),QJ2PREV(:),QJ1PREV(:) !VAZOES NAS CÉLULAS NO INÍCIO DO LOOP PREVISÃO
	REAL,ALLOCATABLE:: SIPREV(:,:) !ARMAZENAMENTO NO RESERVATÓRIO DE INTERCEPTAÇÃO (PREVISÃO)
	INTEGER ITCHUVA

	REAL,ALLOCATABLE:: QCONTORM(:,:),QCONTORJ(:,:) !VAZÃO DA CONDIÇÃO DE CONTORNO DE MUSKINGUM CUNGE
	REAL,ALLOCATABLE:: QRIOINI(:,:) !VAZÃO DA CONDIC'~AO INICIAL DE MUSKINGUM CUNGE NO TRECHO

	REAL,ALLOCATABLE:: QCONTORMPREV(:,:),QCONTORJPREV(:,:) !VAZÃO DA CONDIÇÃO DE CONTORNO DE MUSKINGUM CUNGE
	REAL,ALLOCATABLE:: QRIOINIPREV(:,:) !VAZÃO DA CONDIC'~AO INICIAL DE MUSKINGUM CUNGE NO TRECHO

	REAL,ALLOCATABLE:: QCEL1PREV(:),QCEL2PREV(:) !SOMATORIO DAS VAZÕES GERADAS NA CELULA
	REAL,ALLOCATABLE:: PMB2PREV(:),PMI2PREV(:),PMS2PREV(:) !PROPORÇÕES DE ORIGEM DAS VAZÕES NO RIO
	REAL,ALLOCATABLE:: PJB2PREV(:),PJI2PREV(:),PJS2PREV(:) !PROPORÇÕES DE ORIGEM DAS VAZÕES NO RIO



	REAL,ALLOCATABLE:: BPLAN(:) !LARGURA DA PLANÍCIE DE INUNDAÇÃO (INCLUI O PROPRIO RIO)
	REAL,ALLOCATABLE:: HCALHA1(:),HCALHA2(:) !PROF EM QUE INICIA E EM QUE A PLANICIE ESTA TOT INUNDADA
	REAL,ALLOCATABLE:: QMUP(:,:),AMUP(:,:),BMUP(:,:),CMUP(:,:) !TABELA MUSKINGUN NAO LINEAR
	INTEGER,ALLOCATABLE:: ICODMUSK(:) !CODIGO QUE INDICA SE É LINEAR (0) OU NAO LINEAR (1)

	integer,allocatable:: OD(:) !RP
	real,allocatable:: CalibFLAG(:,:)
	real,allocatable:: CBOLD(:),PLAMOLD(:,:),CAPOLD(:,:),WCOLD(:,:)
	integer,allocatable::	sbtFLAG(:)	! RP



	! Variáveis para calculo de indicadores de evapotranspitacao máxima:
	real,allocatable:: E0agua(:),E0TOPO(:),E0SUP(:)
	INTEGER :: flagaclimed   !para saber se trabalha com CRU - Daniel


	! Variáveis para previsao:
	integer,allocatable:: POSTORIO(:),POSTOBAS(:),POSTOINT(:),POSTOSUP(:) !Postos utilizados para atualizar cada subbacia para cada tipo de var de estado.
	REAL,ALLOCATABLE:: QRG0(:,:)	!ARMAZENA HIDROGRAMAS ONDE SE DESEJA GRAVAR


	! Váriáveis para cálculo do armazenamento total em cada minibacia:
	real,allocatable:: TWS(:),DTWS(:),QM2in(:),QM1in(:) !RP
	real,allocatable:: TWS1(:),TWS2(:),TWS3(:),TWS4(:),TWS5(:),TWS6(:),TWS7(:) !RP
    real,allocatable:: E0media(:)
	real:: Ebacia,Pbacia
	
	
	!@ DCB_sed ###########################################################################
	INTEGER SUBini, SUBfim
	REAL,ALLOCATABLE:: Qmini(:,:)	!@ DCB set 2012 - ARMAZENA HIDROGRAMAS DAS MINI-BACIAS
	!@ DCB_sed ###########################################################################
	
	REAL,ALLOCATABLE:: RUGMAN(:) !Manning coefficient  !FMF em 06/05/2019 - leitura do novo MINI.GTP
	REAL*8,ALLOCATABLE:: HRIO(:)                  !FMF em 06/05/2019 - leitura do novo MINI.GTP
    REAL*8 HRX                                    !FMF em 06/05/2019 - leitura do novo MINI.GTP

	END MODULE
