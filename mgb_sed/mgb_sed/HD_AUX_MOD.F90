	MODULE AUX_MOD
	! Declaração de variáveis globais auxiliares ao hidrodinâmico.
	! 
	!
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descrição das variáveis:  ***********************************
	!
	! hdFLAG0 = indica se MGB utiliza (Sim/Não) modelo hidrodinâmico (hdFLAG0=1/hdFLAG0=0)
	! hdFLAG(.) = indica de existe propagação hidrodinâmica em cada sub-bacia
	! 
	!
	! IDhd(.)  = indices trechos de rio com propagacao hidrodinamica (nTr X 1)
	! TrCELL(.) = subbacias de cada trecho com propagacao hidrodinamica (nTr X 1) 
	! TrL(.) = comprimento dos trechos (nTr x 1) 
	! TrMan(.) = numero de manning dos trechos (nTr x 1)
	! nTr = número de trechos com propagação hidrodinâmica
	! iTr =  indice (contador) de trechos de rio com prop. hidrodinamica
	! iX = indice secoes transversais
	!
	! 
	! nXsec(.) = numero de secoes transversais por trecho (nTr x 1)
	! TrXsec(.) = trecho das secoes transversais (nX x 1)

	! TrNST(.,.) = define secoes de montante e jusante dos trechos (nTr x 2)
	!				Obs.:
	!				TrNST(I,J) J=1 Secao de montante, J=2 Secao de jusante. 
	! bounFLAG(.) = indica tipo de condicao de contorno (nBoun x 1)
	!				1 - resultado do MGB de série de vazao em bacia de cabeceira
	!               2 - resultado do MGB de série de vazão
	!               3 - série de vazão lida em arquivo
	!               4 - série de nível d'água lida em arquivo
	!               5 - curva de descarga
	!				6 - curva de descarga tipo eq. de Manning
	! bounID(.) = identifica que série temporal ou curva de descarga deve ser utilizada (nBoun x 1).
	! bounTS(.) = armazena series temporais de condicao de contorno lidas em arquivo externo (nBoun x 1).
	! Afl(.,.) = área alagada em km2 da relacao nivel x área alagada das subbacias (max de pontos da tabela x nTr)	
	! Zfl(.,.) = niveis em m da relacao nivel x área alagada das subbacias (max de pontos da tabela x nTr)
	! nPfl(.) = numero de pontos da tabela de relacao nivel x área alagada de cada subbacia (nTr)
	! Aflaux = área alagada total de minibacia do intervalo de tempo iT
	! PUSOaux(.,.) = porcentagem de áreas dos blocos (NC x NU)
	! reacFLAG(.) = indica tipo de subtrecho (nReac x 1)
	!				reacFLAG = 1 - subtrecho normal (eq. Saint Venant)
	!				reacFLAG = 2 - subtrecho curto dx=0 (eq. Energia simplificada e continuidade)
	!				reacFLAG = 3 - condicao de contorno interna ou estrutura hidraulica ou barragem
	! Vfl(.,.) = volume na planicie de inundação em km2
	!------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	SAVE

	INTEGER,ALLOCATABLE:: hdFLAG(:)
	INTEGER hdFLAG0

	INTEGER,ALLOCATABLE:: IDhd(:),TrCELL(:)
	REAL,ALLOCATABLE:: TrL(:),TrMan(:)
	INTEGER:: nTr,iTr,iX
	INTEGER,ALLOCATABLE:: nXsec(:),TrXsec(:),TrNST(:,:),bounFLAG(:)
	REAL,ALLOCATABLE:: Afl(:,:),Zfl(:,:)
	INTEGER,ALLOCATABLE:: nPfl(:)
	
	INTEGER,ALLOCATABLE:: bounID(:)
	REAL,ALLOCATABLE:: bounTS(:,:)

	INTEGER,ALLOCATABLE:: reacFLAG(:)

	REAL:: Aflaux
	REAL,ALLOCATABLE:: PUSOaux(:,:)


	real:: dAFImax,dAFImax2,AFImax
	integer:: idAFImax,idAFImax2,iAFImax

	integer:: nHmin



	real,allocatable:: HO1(:),QO1(:)

	real,allocatable:: Tp1(:),Tp2(:)
	real:: dtHD0
	real:: Hmin
	integer,allocatable:: HminFLAG(:)
	real,allocatable:: Qmin(:)


	real:: DHmax,VcritMAX,DQmax

	integer,allocatable:: HminFLAGprev(:),Qminprev(:)

	REAL,ALLOCATABLE:: Vfl(:,:)

	END MODULE




	! Conflito entre variáveis:


	! IR - módulo HD_MAT_MOD e subrotina Ran1		! S/ problema
	! NB - módulo HD_PAR1_MOD e módulo VARS_MAIN ! Verificar se existe conflito	RESOLVIDO, MUDOU NB DO HD PARA NBO
	! DX - módulo HD_PAR1_MOD e subrotinas Funcd, Musk_NL,Newtrap,Parcunge,rede ! Verificar		RESOLVIDO. MUDOU DX DO HIDRO PARA DXT
	! TA - módulo HD_PAR1_MOD e módulo VARS_MAIN ! Verificar							RESOLVIDO. MUDOU TA DO HIDRODIN PARA BA
	! F - módulo HD_PAR1_MOD e subrotinas Funcd, Newtrap ! S/ problema
	! SF - módulo HD_PP_MOD e modulo  VARS_MAIN ! Verificar (acho que não tem problema)   OBS. NO MGB FORAM TIRADAS VARIAVEIS DECLARADAS EM SUBROTINAS (CELULA,EVAPO,RADIACAO) E PASSADAS PARA GLOBAL
			!RESOLVIDO PARCIALMENTE
				! SF DO HIDRODIN PASSA PARA SFF
				! VOLTAR A NOTACAO ORIGINAL E DECLARAR VARS DO MGB DENTRO DAS ROTINAS RADIACAO,CELULA,EVAPO

	! DT - módulo HD_TIME_MOD e modulo  VARS_MAIN ! Verificar vai dar problema	RESOLVIDO. MUDOU DT DO HIDRODIN PARA DThd
	! NT - módulo HD_TIME_MOD e modulo  VARS_MAIN ! Verificar vai dar problema  RESOLVIDO. MUDOU NT DO HIDRODIN PARA NThd
	! iT - módulo HD_TIME_MOD e modulo  VARS_MAIN ! Verificar vai dar problema  RESOLVIDO. MUDOU IT DO HIDRODIN PARA IThd

