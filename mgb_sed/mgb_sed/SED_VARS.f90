!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Mar�o de 2011
!@
!@ Atualizado: Abr 2011
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MODULO SED_VARS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@-------------------------------------------------------------------------------------
!@  VARI�VEIS
!@ Kusle	= is the USLE soil erodibility factor [(0.013*ton*m2*hr)/(m3*ton*cm)]
!@ Culse	= is the USLE cover and management factor
!@ Pusle	= is the USLE support practice factor
!@ Rgros	= is the coarse fragment factor
!@
!@ Mareia	= is the percent sand content (0.05-2.00 mm diameter particles)
!@ Msilte	= is the percent silt content (0.002-0.05 mm diameter particles)
!@ Margila	= is the percent clay content (< 0.002 mm diameter particles)
!@ Morg		= is the percent of organic mater in the soil
!@ Mrocha	= is the percent rock in the first soil layer
!@
!@ LSAcu	= somat�rio do fator LS de cada uso em cada minibacia
!@ PSC		= aporte de sedimentos da minibacia para o rio (ton)
!@ PSU		= aporte de sedimentos dos blocos da minibacia para o rio (ton)
!@ SLC		= perda de solo total da minibacia (ton)
!@ SLU		= perda de solo total do bloco da minibacia (ton)
!@ QSC		= descarga s�lida das minibacias (ton/s) no passo de tempo
!@ QSU      = descarga s�lida dos blocos das minibacias (ton/s)
!@ QSAUX	= vari�vel que armazena, para o escoamento superficial de cada minibacia, o tempo de retardo (s),
!@           o volume total (m3) e as l�minas (mm) para cada bloco
!@ NPXU		= n�mero de pixels de cada uso em cada minibacia
!@ SDR      = valores do SDR de cada bloco de cada minibacia e de cada minibacia
!@
!@ DMa		= di�metro caracter�stico para as part�culas de areia
!@ DMs		= di�metro caracter�stico para as part�culas de silte
!@ DMg		= di�metro caracter�stico para as part�culas de argila
!@
!@ APIX     = �rea m�dia dos pixels de uma minibacia (km2)
!@ FCTE     = fator constante da MUSLE para cada uso da minibacia
!@ TPIC     = taxa de pico do escoamento superficial (m3/s)
!@ SLX      = carga e sedimentos remanescente na minibacia no passo de tempo
!@ SLXU     = carga e sedimentos remanescente por bloco da minibacia no passo de tempo
!@
!@ RUGMAN(IC)   = rugosidade de Manning para MC
!@ VISC     = viscosidade cinem�tica da �gua (m2/s)
!@
!@ PSCaux   = aporte das fra��es de sedimentos que chegam � drenagem das minibacias das subbacias selecionadas
!@ DMP      = di�metro nominal das part�culas (m)
!@ WSP      = velocidade de queda das part�culas (m/s)
!@ CTS      = capacidade de transporte do trecho (ton/m3)
!@ Fform    = fator de forma das part�culas de sedimento
!@
!@ CONCTREC1= concentra��o de sedimentos no trecho no tempo 1 (ton/m3)
!@ CONCTREC2= concentra��o de sedimentos no trecho no tempo 2 (ton/m3)
!@ VolTREC1 = volume m�dio de �gua no trecho no tempo 1 (m3)
!@ VolTREC2 = volume m�dio de �gua no trecho no tempo 2 (m3)
!@
!@ DEPTREC  = carga acumulada de sedimento depositado no trecho (ton)
!@ EROSTREC = carga acumulada de sedimento erodido no trecho (ton)
!@ CARGM    = carga de sedimento saindo para o trecho de jusante(ton)
!@
!@ HTR1     = PROFUNDIDADE NA SE��O DE JUSANTE DO TRECHO NO TEMPO T
!@ HTR2     = PROFUNDIDADE NA SE��O DE JUSANTE DO TRECHO NO TEMPO T+1
!@ VFL1     = VOLUME DE �GUA NA PLAN�CIENO TEMPO T
!@ VFL2     = VOLUME DE �GUA NA PLAN�CIENO TEMPO T+1
!@ CFL1     = CONCENTRA��O NA PLAN�CIENO TEMPO T
!@ CFL2     = CONCENTRA��O NA PLAN�CIENO TEMPO T+1
!@ DFL      = DEP�SITO DE SILTE E ARGILA NA PLAN�CIE
!@-------------------------------------------------------------------------------------
!@
!@ *********************************************************************************************

MODULE SED_VARS

!@ J� LISTADAS
INTEGER,PARAMETER ::    FILPAR=300 !ARQUIVO DOS PARAMETROS DA MUSLE ASSOCIADOS AO USO
INTEGER,PARAMETER ::    FILTEX=301 !ARQUIVO DAS TEXTURAS DOS SOLOS ASSOCIADOS AO USO
INTEGER,PARAMETER ::    FILLSM=302 !ARQUIVO DOS LS ACUMULADOS DE CADA USO DE CADA MINIBACIA
INTEGER,PARAMETER ::    FILHRU=303 !ARQUIVO COM NUMERO DE PIXELS DE CADA USO DE CADA MINIBACIA
INTEGER,PARAMETER ::    FILSED=304 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA � REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILSEDa=305 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA � REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILSEDs=306 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA � REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILSEDc=307 !ARQUIVO COM CARGA DE SEDIMENTOS DAS MINIBACIAS QUE CHEGA � REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILQCE=308 !ARQUIVO COM A VAZ�O SUPERFICIAL DAS MINIBACIAS QUE CHEGA � REDE DE DRENAGEM
INTEGER,PARAMETER ::    FILCAR=319 !306 !ARQUIVO COM A CONCENTRA��O DE AREIA DOS EXUTORIOS DESEJADOS(mg/L)
INTEGER,PARAMETER ::    FILCSI=320 !307 !ARQUIVO COM A CONCENTRA��O DE SILTE DOS EXUTORIOS DESEJADOS(mg/L)
INTEGER,PARAMETER ::    FILCCL=321 !308 !ARQUIVO COM A CONCENTRA��O DE ARGILA DOS EXUTORIOS DESEJADOS(mg/L)
INTEGER,PARAMETER ::    RIOCAR=309 !ARQUIVO COM A CONCENTRA��O DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    RIOCSI=310 !ARQUIVO COM A CONCENTRA��O DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    RIOCCL=311 !ARQUIVO COM A CONCENTRA��O DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    FILSDR=312 !ARQUIVO COM A CONCENTRA��O DE ARGILA DOS EXUTORIOS DAS MINIBACIAS(mg/L)
INTEGER,PARAMETER ::    SEDSAI=313 !ARQUIVOS COM A CARGA DE SEDIMENTOS NO TRECHO DA MINIBACIA (ton/dia)
INTEGER,PARAMETER ::    SEDDEP=314 !ARQUIVO COM A CARGA DE SEDIMENTOS DEPOSITADA NO TRECHO DA MINIBACIA (ton/dia)
INTEGER,PARAMETER ::    SEDERO=315 !ARQUIVO COM A CARGA DE SEDIMENTOS ERODIDA NO TRECHO DA MINIBACIA
INTEGER,PARAMETER ::    TAXPIC=316 !ARQUIVO COM O ACUMULADO DAS TAXAS DE PICO DOS BLOCOS DAS MINIBACIAS
INTEGER,PARAMETER ::    FILVAZ=317 !ARQUIVO DOS HIDROGRAMAS DAS MINI-BACIAS
INTEGER,PARAMETER ::    FILCEL=318 !317 !ARQUIVO COM A VAZ�O DAS MINIBACIAS QUE CHEGA � REDE DE DRENAGEM 2013_02_06
INTEGER,PARAMETER ::    VAZFLP=322 !399 !@ DCB_HD_Sed
INTEGER,PARAMETER ::    VAZFLP2=323 !400 !@ DCB_HD_Sed
INTEGER,PARAMETER ::    VAZFLP3=324 !401 !@ DCB_HD_Sed
INTEGER,PARAMETER ::    BALFLP=325 !401 !@ DCB_HD_Sed
CHARACTER               (10) Nuso(12) !NOMES DOS USOS DO SOLO (BLOCOS)
REAL                    DMa, DMs, DMg
REAL                    APIX, FCTE, QPIC, SLX
REAL                    VISC !, RUGMAN !FMF em 06/05/2019 - leitura do novo MINI.GTP
REAL,ALLOCATABLE ::     Kusle(:,:),Cusle(:,:),Pusle(:,:),Rgros(:,:), Ksdr(:,:) !@ Ksdr(:) DCB 02-07-2011
REAL,ALLOCATABLE ::     Mareia(:,:),Msilte(:,:),Margila(:,:),Morg(:,:),Mrocha(:,:)
REAL,ALLOCATABLE ::     LSAcu(:,:), SDR(:,:) !@ SDR(:,:) DCB 02-07-2011
REAL,ALLOCATABLE ::     PSC(:,:),SLC(:),QSC(:),QSAUX(:,:)
REAL,ALLOCATABLE ::     PSU(:,:),QSU(:,:),SLXU(:)
REAL,ALLOCATABLE ::     SLU(:,:)
REAL,ALLOCATABLE ::     PSCaux(:,:)
REAL,ALLOCATABLE ::     DMP(:), WSP(:), CTS(:,:), Fform(:)
REAL,ALLOCATABLE ::     CSM1(:,:), CSM2(:,:), CSJ1(:,:), CSJ2(:,:)  !, CSJ2aux(:,:)
REAL,ALLOCATABLE ::     VolTREC1(:), VolTREC2(:)
REAL,ALLOCATABLE ::     DEPTREC(:,:), DEPT(:,:), EROSTREC(:,:), EROT(:,:), CARGM(:,:)
REAL,ALLOCATABLE ::     QSJ2(:,:)
INTEGER,ALLOCATABLE ::  NPXU(:,:)
REAL,ALLOCATABLE ::     CAREIA(:,:)	!ARMAZENA CONCENTRACOES DE AREIA ONDE SE DESEJA GRAVAR
REAL,ALLOCATABLE ::     CSILTE(:,:)	!ARMAZENA CONCENTRACOES DE SILTE ONDE SE DESEJA GRAVAR
REAL,ALLOCATABLE ::     CARGIL(:,:)	!ARMAZENA CONCENTRACOES DE ARGILA ONDE SE DESEJA GRAVAR

INTEGER                 iSEDaux, nSEDmini
INTEGER,ALLOCATABLE ::  CONTIC(:),CONTICSUB(:)


!@ N�O LISTADAS AINDA
REAL,ALLOCATABLE ::     FracS(:,:)

INTEGER,ALLOCATABLE ::  CellTr(:)               !@ DCB_HD_Sed - ARMAZENA C�DIGO DO TRECHO HD DE CADA MINIBACIA
REAL,ALLOCATABLE ::     HTR1(:), HTR2(:)        !@ DCB_HD_Sed
REAL,ALLOCATABLE ::     VFL1(:), VFL2(:)        !@ DCB_HD_Sed
REAL,ALLOCATABLE ::     CFL1(:,:), CFL2(:,:)    !@ DCB_HD_Sed
REAL,ALLOCATABLE ::     DFL(:,:)                !@ DCB_HD_Sed

!@ ------------------------------------------------------------------------
!@ BALAN�O DE SEDIMENTOS
REAL,ALLOCATABLE ::     iENTRA(:,:), iSAI(:,:)  !@ DCB_HD_Sed   ,iERRO(:,:)
REAL,ALLOCATABLE ::     inFL(:,:), outFL(:,:), BalQFL(:,:), VTOTAL(:,:)  !@ DCB_HD_Sed
!@ ------------------------------------------------------------------------




END MODULE