	MODULE PAR1_MOD
	! Declaração de variáveis globais relativas topologia e geometrica da rede simulada.
	! 
	!
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descrição das variáveis:  ***********************************
	!
	!
	! NX = numero de seções transversais 
	! NREAC = numero de subtrechos (ligação entre 2 secoes)
	! NCONF = numero de confluencias
	! NBOUN = numero de condicoes de contorno
	! NST(.,.) = define secoes de montante e jusante dos trechos (NREAC x 2)
	!				Obs.:
	!				NST(I,J) J=1 Secao de montante, J=2 Secao de jusante. 
	! NCC(.,.) = numero da secoes de cada confluencia. (NCONF x 3)
	!               Obs.:
	!				NCC(I,J) J=1 e J=2 se referem as secoes dos bracos e J=3 a secao principal. 
	!					NCC(I,3)<0 quando confluencia é convergente
	!					NCC(I,3)>0 quando confluencia é divergente.
	! NBO(.) = indice das secoes com condicoes de contorno (NBOUN x 1)
	! ZO(.) = nivel de fundo das secoes (NX x 1)
	! DXT(.) = comprimento dos subtrechos em metros (NREAC x 1)
	! DXC(.) = distancia entre secoes das confluencias em metros (NCONF*2 x 1)
	!          Obs.: DXC((I-1)*2+1) = Distancia entre seções NCC(I,1) e NCC(I,3) da confluencia I        	
	!				 DXC((I-1)*2+2) = Distancia entre seções NCC(I,2) e NCC(I,3) da confluencia I        		
	! ALFA(.) = coef. de energia das confluencias (NCONF*2 x 1)
	!          Obs.: ALFA((I-1)*2+1) = alfa entre seções NCC(I,1) e NCC(I,3) da confluencia I        	
	!				 ALFA((I-1)*2+2) = alfa entre seções NCC(I,2) e NCC(I,3) da confluencia I        		
	! NP(.) = numero de pontos da tabela cota x Area x Raio Hid. x Largura das secoes (NX x 1)
	! AR(.,.) = area molhada das secoes transversais (max(NP) x NX)
	! RR(.,.) = raio hidraulico das secoes transversais (max(NP) x NX)
	! HA(.,.) = prof. ou niveis d'agua da tabela de caracteristicas das secoes transersais (max(NP) x NX)
	! BA(.,.) = largura das secoes transversais (max(NP) x NX)
	! BETA(.) = coeficiente de Boussinesq das secoes (NX x 1)
	! ICONF = flag - considera energia cinetica e perdas de carga nas confluencias
	!					Sim - ICONF=0
	!					Nao - ICONF>0
	! LD = dados informados em termos de profundidade (LD=0) ou niveis d'água (LD>0)
	! LRO(.) = > 0 numero de secoes que nao possuem um tirante minimo. VERIFICAR SE ISSO É NECESSARIO
	! F(.,.) = coef. de rugosidade de Manning (max(N) x NX)
	! N(.) = numero de pontos da tabela de rugosidade x nivel ou prof. (NX x 1)
	!			Se N(I)=0, coef. de rugosidade da secao I é constante
	! HUX(.,.) = prof. da tabela de coef. de Manning das seções (max(N) x NX)
	! THETA(.) = coeficiente theta=dt/dx (NREAC x 1)
	!
	! CONST = variavel auxiliar
	! ITT = auxiliar
	! INIC = indica se calcula condições iniciais? (Sim - INIC>0 / Não - INIC=0)
	! HFIM = nivel d'água em t=0 s na secao de jusante para calculo de remanso (condições iniciais)
	!------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	SAVE

	INTEGER NX,NREAC,NCONF,NBOUN
	INTEGER,ALLOCATABLE:: NST(:,:),NCC(:,:),NBO(:)
	REAL,ALLOCATABLE:: ZO(:),DXT(:),DXC(:),ALFA(:)
	INTEGER,ALLOCATABLE:: NP(:)
	REAL,ALLOCATABLE:: AR(:,:),RR(:,:),HA(:,:),BA(:,:)
	REAL,ALLOCATABLE:: BETA(:)
	INTEGER ICONF,LD
	REAL,ALLOCATABLE:: F(:,:),HUX(:,:)
	INTEGER,ALLOCATABLE:: N(:)
	REAL,ALLOCATABLE:: THETA(:)
	REAL CONST
	INTEGER ITT
	INTEGER INIC
	REAL HFIM
	REAL Z01,Z02
	

	END MODULE