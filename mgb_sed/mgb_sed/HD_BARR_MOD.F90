	MODULE BARR_MOD
	! Declara??o de vari?veis globais relativas as curvas de descarga das barragens.
	! 
	!
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descri??o das vari?veis:  ***********************************
	!
	!
	! ! Obs.: indices das barragens devem estar numerados em ordem crescente com indices de se??es e trechos
	! NPB(.,.) = numero de pontos referente a secao de montante (LTOT x NCMAX)
	!				 (obs.: se negativo existe efeito de jusante)
	! NPB2(.,.) = numero de pontos referente a secao de jusante (LTOT x NCMAX) 
	! NPMMAX,NPJMAX = numero m?ximo de pontos referente as sec?es de montante e jusante
	! NCMAX = m?ximo numero de curvas de descarga da barragem
	! QBR(.,.,.,.) = vaz?es da tabela das curvas de descarga das barragens (NPMMAX x NPJMAX x LTOT x NCMAX)
	!				 Ex.: QBR(I,J,K,L)=vaz?o na prof. ou nivel de montante I, prof. ou nivel 
	!									de jusante J,barragem K e curva L
	! HBR(.,.,.) = tabela de prof. ou nivel dagua da secao de montante (NPMMAX x LTOT x NCMAX)
	! H2BR(.,.,.) = tabela de prof. ou nivel dagua da secao de jusante (NPJMAX x LTOT x NCMAX)
	! JIT(.,.) = intervalo de tempo ate onde ? valida a curva de descarga (LTOT x NCMAX)
	! NIT(.) = numero de curvas de descarga de cada barragem (LTOT x 1)
	! LTOT = numero total de barragens
	! barrReac(.) = subtrecho da barragem (LTOT x 1)
	!------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	SAVE

	INTEGER NPMMAX,NPJMAX,NCMAX
	INTEGER,ALLOCATABLE:: NPB(:,:)
	REAL,ALLOCATABLE:: QBR(:,:,:,:),HBR(:,:,:),H2BR(:,:,:)
	INTEGER,ALLOCATABLE:: NPB2(:,:),JIT(:,:),NIT(:)
	INTEGER LTOT
	INTEGER,ALLOCATABLE:: barrReac(:)

	END MODULE