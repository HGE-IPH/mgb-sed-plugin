	MODULE CHEIA_MOD
	! Declaração de variáveis globais relativas as tabelas das planicies de inundação
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descrição das variáveis:  ***********************************
	!
	!
	! NPF(.) = numero de pontos da tabela de areas de inundacao (NX x 1)
	! FAF(.,.) = largura da tabela nivel ou prof x largura planicie 
	!							de inundacao. (max de pontos da tabela x NX)
	! HF(.,.) = niveis ou prof. d'agua da tabela nivel ou prof x largura planicie
	!							de inundacao. (max de pontos da tabela x NX)
	!------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	SAVE
	INTEGER,ALLOCATABLE:: NPF(:)
	REAL,ALLOCATABLE:: FAF(:,:),HF(:,:)

	END MODULE