    PROGRAM MGBPREV
	!******************************************************************
	!***************MODELO HIDROL�GICO DE GRANDES BACIAS***************
	!******************************************************************
	!*					Walter Collischonn                            *
	!*			Instituto de Pesquisas Hidr�ulicas                    *
	!*		Universidade Federal do Rio Grande do Sul                 *
	!*                    fevereiro de 2002                           *
	!*                modificado em agosto 2002                       *
	!******************************************************************

	USE PORTLIB !biblioteca para calcular o tempo de processamento
	USE VARS_MAIN !m�dulo que cont�m as vari�veis principais
	USE VARS_CALIB !m�dulo que cont�m as vari�veis da calibra��o
!**********************
	USE AUX_MOD
!***********************
	IMPLICIT NONE
	INTEGER JULDAY

	CALL TEMPO(0,ISEED)  ! INICIA CONTAGEM DO TEMPO	

	!____________LEITURA DE ARQUIVOS PRINCIPAIS________________
	CALL LEFIX 	!SUBROTINA DE LEITURA DE PARAMETROS FIXOS
	CALL ALLOCA_VARS(0) !ALLOCA VARI�VEIS PRINCIPAIS

	IDINI=JULDAY(IMES,IDIA,IANO) !CONVERTE O DIA DE INICIO EM CALENDARIO JULIANO
	CALL LEVAR !SUBROTINA DE LEITTURA DE PARAMETROS VARIAVEIS (ALBEDO, iaf, RS E Z MENSAIS)
	CALL LEUSO !SUBROTINA DE LEITURA DO ARQUIVO DE PARAMETROS DE USO
	CALL LECELL !SUBROTINA DE LEITURA DO ARQUIVO DAS CELULAS
	CALL LECLIMED !SUBROTINA DE LEITURA DO ARQUIVO DE DADOS CLIMATOLOGICOS MENSAIS
	CALL ARQCLISUB	  !L� DADOS DOS ARQUIVOS DE CLIMA
	CALL LEQOBS	!LEITURA DE DADOS DE VAZAO OBSERVADA 
	!CALL LEOUTROSQ !LEITURA DE DADOS DE VAZ�O OBSERVADA NO EIXO PRINCIPAL DO SAO FRANCISCO
	IF(NUMSUBST>0)THEN
		CALL LESUBST
	ENDIF

!*****************************************************************
	! Modelo hidrodin�mico:
	hdFLAG=0 ! Para desativar hidrodinamico
	!hdFLAG0=sum(hdFLAG) ! Teste

!@ DCB_sed ###########################################################################
    !@ MADEIRA   - 57 a  82
    !@ TAPAJOS   - 83 a  92
    !@ TOCANTINS - 98 a 168
    !@ TUDO      -  1 a 180
	SUBini = 1
	SUBfim = 47
!@ DCB_sed ###########################################################################

			! Teste
	IF (hdFLAG0>0) THEN
		! Subrotina de leitura e organizacao dos dados do modelo hidrodinamico:
		write(*,*) 'Le dados hidrodinamico'
		CALL INPUT2
		!Organizacao da matriz de coeficientes:
		write(*,*) 'Organiza matriz de coeficientes'

		CALL MATRIX
		call SKYLINE0
!		CALL MATRIX2
	ENDIF
!*****************************************************************


	!___________FIM DA LEITURA DOS ARQUIVOS PRINCIPAIS ______________

	!_______________PREPARACAO DE DADOS_____________________
	CALL PARCEL !CALCULA ALGUNS PAR�METROS DAS C�LULAS E DO RIO
	CALL PARCUNGE !CALCULA PARAMETROS PARA PROPAGACAO MUSKINGUM-CUNGE
	!____________FIM DA PREPARACAO DE DADOS_____________________

	!ABRE ARQUIVOS DE SAIDA  E DE ENTRADA
	OPEN(FILSOL,FILE='.\output\NOSOLO.HIG',STATUS='UNKNOWN')	!ARQUIVO DE SAIDA COM DADOS DO SOLO DE UMA CELULA
	OPEN(FILSOL2,FILE='.\output\NOSOLO2.HIG',STATUS='UNKNOWN')	!ARQUIVO DE SAIDA COM DADOS DO SOLO DE UMA CELULA

	OPEN(FILAJU,FILE='.\output\AJUSTE.HIG',STATUS='UNKNOWN') !ARQUIVO DE SA�DA COM AJUSTE DOS HIDROGRAMAS
	!OPEN(FILPLU,FILE='CHUVABIN.HIG',STATUS='old',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT') !dados de chuva interpolados
	OPEN(FILPLU,FILE='.\input\CHUVAbin.HIG',STATUS='old',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT') !dados de chuva interpolados
	!OPEN(FILPLU,FILE='DIARIA_CEL.HIG',STATUS='old',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT') !dados de chuva interpolados

!	OPEN(971,FILE='E0agua.txt',STATUS='UNKNOWN') !
!	OPEN(972,FILE='E0topo.txt',STATUS='UNKNOWN') !
!	OPEN(973,FILE='E0sup.txt',STATUS='UNKNOWN') !
!@ DCB ago/2012	OPEN(974,FILE='.\output\TWSmgb.txt',STATUS='UNKNOWN') !
!@ DCB ago/2012	OPEN(975,FILE='.\output\FluxosMedio.txt',STATUS='UNKNOWN') !
!@ DCB ago/2012	OPEN(988,FILE='.\output\ET.txt',STATUS='UNKNOWN') !
!@ DCB ago/2012	OPEN(989,FILE='.\output\P.txt',STATUS='UNKNOWN') !





	!DECIDE SE VAI SIMULAR, CALIBRAR OU FAZER PREVIS�O (ICALIB=0; ICALIB=1; ICALIB=2) 
	CALIBRA_CASE: SELECT CASE (ICALIB) !VERIFICA SE CALIBRA OU N�O
	CASE (0) ! N�O CALIBRA, APENAS SIMULA
!@ DCB set/2012		WRITE(*,*)'SIMULA'
		CALL SIMULA	
	CASE (1) ! CALIBRA
		CALL LeCalib
		CALL ALLOCA_CALIB(0) !ALLOCA VARI�VEIS DE CALIBRA��O
		CALL CALIBRA
		CALL ALLOCA_CALIB(1) !DEALLOCA VARI�VEIS DE CALIBRA��O
	CASE (2) !SIMULA PREVIS�O
		WRITE(*,*)' FAZENDO PREVISAO'
		CALL PREVISAO

	CASE (6) !SIMULA Hindcast of ESP/revESP
		WRITE(*,*)' Hindcast of ESP/revESP'
		CALL PREV_ESP

	CASE DEFAULT
		STOP ' ERRO: ICALIB DESCONHECIDO NO PROGRAMA PRINCIPAL!!!'
	END SELECT CALIBRA_CASE

!@ DCB	CALL TEMPO(1,0)  ! FINALIZA CONTAGEM DO TEMPO
	CALL TEMPO(1,ISEED)  ! FINALIZA CONTAGEM DO TEMPO
	
	CLOSE (FILPLU) !ARQUIVO DE ENTRADA CHUVA
	CLOSE (FILSOL) !ARQUIVO DE SAIDA DADOS DO SOLO
	CLOSE (FILAJU)

	CALL ALLOCA_VARS(1) !DEALLOCA VARI�VEIS PRINCIPAIS


	WRITE(*,*)'PROGRAMA TERMINOU'
	WRITE(*,*)'TECLE ENTER'
	READ(*,*)
	STOP
	END         