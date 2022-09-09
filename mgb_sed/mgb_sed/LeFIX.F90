	SUBROUTINE LEFIX
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER K
	!	SUBROUTINE LEFIX(FILFIX,NC,NU,NT,NP,NCLI,NOBS,IQOBS,NBP,NB,
	!	&IDIA,IMES,IANO,NCLIP,ARQCLI,ACLIMED,ARQOBS,ICALIB)
	!esta subrotina le o arquivo PARHIG.hig
	!onde estao os dados gerais fixos da simulaçao
	!como o numero de celulas, o numero de usos do solo
	!o numero de interalos de tempo e de postos pluviometricos

	OPEN(FILFIX,FILE='.\input\PARHIG.HIG',STATUS='OLD',ACTION='READ')		!abre arquivo de dados gerais

	READ(FILFIX,75)TITULO !LE CABEÇALHO PRINICPAL
	WRITE(*,*)TITULO
	READ(FILFIX,75)TITULO !LE LINHA EM BRANCO
	WRITE(*,*)TITULO
	READ(FILFIX,75)TITULO !LE TITULO DA APLICAÇÃO
	WRITE(*,*)TITULO
	READ(FILFIX,*)TITULO !LE OUTRA LINHA EM BRANCO
	WRITE(*,*)TITULO
!	READ(*,*)

	READ(FILFIX,72)(CABE(K),K=1,4)
	WRITE(*,72)(CABE(K),K=1,4)
	READ(FILFIX,71)IDIA,IMES,IANO,HORAINI !LÊ DATA E HORA DO INÍCIO DA SIMULAÇÃO
	WRITE(*,71)IDIA,IMES,IANO,HORAINI
!READ(*,*)
	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,72)(CABE(K),K=1,2)  
	WRITE(*,72)(CABE(K),K=1,2)
	READ(FILFIX,76)NT,DTP					!LÊ NUMERO DE INTERVALOS E TAMANHO DO INTERVALO DE TEMPO
	WRITE(*,76)NT,DTP

	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,72)(CABE(K),K=1,4) 
	WRITE(*,72)(CABE(K),K=1,4)
    !@ ************************ DCB *****************************
!Passar deste formato:
!    READ(FILFIX,71)NC,NU,NB,NCLI !LÊ NUMERO DE CÉLULAS, USOS, BACIAS E POSTOS CLIMATOLÓGICOS
!Para este formato:
 !   READ(FILFIX,'(5I10,I13)')NC,NU,NB,NCLI !LÊ NUMERO DE CÉLULAS, USOS, BACIAS E POSTOS CLIMATOLÓGICOS
	! Formato livre
    READ(FILFIX,*)NC,NU,NB,NCLI !LÊ NUMERO DE CÉLULAS, USOS, BACIAS E POSTOS CLIMATOLÓGICOS
    WRITE(*,71)NC,NU,NB,NCLI !@ DCB set/2012
    !@ ************************ DCB *****************************	 

	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,72)(CABE(K),K=1,1)
	WRITE(*,72)(CABE(K),K=1,1)
	READ(FILFIX,71)ICALIB
	WRITE(*,71)ICALIB





!READ(*,*)
	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,75)TITULO					!LE TEXTO
!@ DCB set/2012	WRITE(*,*)TITULO
!READ(*,*)
	
	
	if(NCLI<=0)then  !introduzido por Daniel para Amazonia
     !arquivo do CRU
	  NCLI=-NCLI
	  flagaclimed=1	! usar dados do CRU
	else
		flagaclimed=0
		DO K=1,NCLI
			READ(FILFIX,77)ARQCLI(K)
			WRITE(*,77)ARQCLI(K)
		ENDDO
	endif

!READ(*,*)
	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,75)TITULO					!LE TEXTO
!@ DCB set/2012	WRITE(*,*)TITULO
	READ(FILFIX,*)ACLIMED
	WRITE(*,*) 'NOME ARQUIVO CLIMA: ', ACLIMED  !@ DCB set/2012
!READ(*,*)
	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,75)TITULO					!LE TEXTO
!@ DCB set/2012	WRITE(*,*)TITULO
	READ(FILFIX,*)NOBS,ARQOBS !lê quantos postos com dados observados tem e qual o nome do arquivo dos dados
	WRITE(*,*)'NUMERO DE Qobs E NOME ARQUIVO: ', NOBS,ARQOBS
!READ(*,*)
	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,75)TITULO					!LE TEXTO
!@ DCB set/2012	WRITE(*,*)TITULO
	READ(FILFIX,*)(IQOBS(K),K=1,NOBS) !Lê número da célula em que existem dados observados
!@ DCB set/2012	WRITE(*,*)(IQOBS(K),K=1,NOBS)
!READ(*,*)
	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,75)TITULO					!LE TEXTO
!@ DCB set/2012	WRITE(*,*)TITULO
	READ(FILFIX,*)NUMHIDG !Lê em quantos pontos se deseja gravar o hidrograma
!@ DCB set/2012	WRITE(*,*)NUMHIDG

!READ(*,*)
	READ(FILFIX,*)							!LÊ LINHA EM BRANCO
	READ(FILFIX,75)TITULO					!LE TEXTO
!@ DCB set/2012	WRITE(*,*)TITULO
!READ(*,*)
	DO K=1,NUMHIDG
		READ(FILFIX,*)IHIDG(K) !lê as células correspondentes aos pontos em que se deseja gravar os hidrogramas
!		WRITE(*,*)IHIDG(K) 
	ENDDO
!@ DCB set/2012	WRITE(*,*) (IHIDG(K),K=1,NUMHIDG) 
!READ(*,*)
	READ(FILFIX,*)
	READ(FILFIX,75)TITULO					!LE TEXTO
!@ DCB set/2012	WRITE(*,*)TITULO
	READ(FILFIX,*)NUMSUBST,ARQSUBST !Lê em quantos pontos se deseja SUBSITUIR a vazão calculada pela lida de arquivo

	WRITE(*,*) 'NUMERO DE Qsubst E ARQUIVO: ', NUMSUBST,ARQSUBST
!READ(*,*)

	READ(FILFIX,*)
	READ(FILFIX,75)TITULO
!@ DCB set/2012	WRITE(*,*)TITULO

!READ(*,*)
	DO K=1,NUMSUBST
		READ(FILFIX,*)ISUBST(K) !LÊ as células cuja vazao calculada será substituída pela lida
!@ DCB set/2012		WRITE(*,*)ISUBST(K)
	ENDDO

!READ(*,*)

	CLOSE (FILFIX)

71	FORMAT(6I10)
72	FORMAT(5A10)
73	FORMAT(F10.2)
74	FORMAT(A10)
75	FORMAT(A20)
76	FORMAT(I10,F10.1)
77	FORMAT(A20)


	RETURN
	END