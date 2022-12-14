	SUBROUTINE ALLOCA_CALIB(IOP)
	USE VARS_CALIB
	IMPLICIT NONE
	INTEGER IOP
    SAVE

	ALLOC_CASE: SELECT CASE (IOP) !VERIFICA SE ALLOCA OU DEALLOCA
	CASE (0) ! ALLOCA
		ALLOCATE (IRUIM(NS))
		ALLOCATE (FO(NS,NF)) !FUN??ES OBJETIVO
		ALLOCATE (IPARET(NS))
		ALLOCATE (XPAR(NPAR))
		ALLOCATE (SOMAPAR(NPAR),REFLEX(NPAR),CONTRA(NPAR)) !SOMA DAS COORD., PONTO DE REFLEX?O, PONTO DE CONTRA??O
		ALLOCATE (MEDIA(NS)) !M?DIA DE FUN??ES OBJETIVO DA POPULA??O
		ALLOCATE (PPAR(NPAR,NS)) !VALOR RELATIVO DO PARAMETRO NA FAIXA DE VALIDADE
		ALLOCATE (PAR(NPAR,NS)) !VALOR DO PARAMETRO 
!		ALLOCATE (PMIN(NPAR),PMAX(NPAR)) !FAIXA DE VALIDADE DOS PARAMETROS
		ALLOCATE (FOLD(NS)) !ANTIGO VALOR DA FUN??O EM CADA PONTO DA AMOSTRA
		ALLOCATE (PROB(NS)) !PROBAB. DE ESCOLHA PONTO 1..NS
		ALLOCATE (VMIN(NF)) !VALOR MINIMO DAS FUN??ES OBJETIVO
		ALLOCATE (PARX(NPAR)) !VALOR DO PARAMETRO 
!		
	CASE (1) ! DEALLOCA
		DEALLOCATE (IRUIM,FO,IPARET,XPAR)
		DEALLOCATE (SOMAPAR,REFLEX,CONTRA) !SOMA DAS COORD., PONTO DE REFLEX?O, PONTO DE CONTRA??O
		DEALLOCATE (MEDIA) !M?DIA DE FUN??ES OBJETIVO DA POPULA??O
		DEALLOCATE (PPAR) !VALOR RELATIVO DO PARAMETRO NA FAIXA DE VALIDADE
		DEALLOCATE (PAR) !VALOR DO PARAMETRO 
		DEALLOCATE (PMIN,PMAX) !FAIXA DE VALIDADE DOS PARAMETROS
		DEALLOCATE (FOLD) !ANTIGO VALOR DA FUN??O EM CADA PONTO DA AMOSTRA
		DEALLOCATE (PROB) !PROBAB. DE ESCOLHA PONTO 1..NS
		DEALLOCATE (VMIN) !VALOR MINIMO DAS FUN??ES OBJETIVO
	CASE DEFAULT
		STOP ' ERRO: IOP DESCONHECIDO NO ROTINA ALLOCA_CALIB!!!'
	END SELECT ALLOC_CASE

	!ESTA SUBROTINA TEM OS COMANDOS DE ALLOCA??O DE MEM?RIA DAS VARI?VEIS DA CALIBRA??O AUTOM?TICA

	RETURN
	END