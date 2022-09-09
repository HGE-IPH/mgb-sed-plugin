	SUBROUTINE LECELL
	!esta subrotina le o arquivo cell.hig
	!onde estao os dados das celulas, como as coordenadas
	!e as porcentagens de uso
	USE VARS_MAIN
	use AUX_MOD !RP
	IMPLICIT NONE
	!VARIAVEIS QUE DEPENDEM DO NUMERO DE CELULAS----
	INTEGER ICELL(NC)
	INTEGER IBANT,I,K,KW
	integer aux !RP
	OPEN(FILHIG,FILE='.\input\CELL.HIG',STATUS='OLD')
	READ(FILHIG,77) !RP
	IBANT=1
	IEXUT=-999999
	DO I=1,NC
		!FMF em 06/05/2019 - leitura do novo MINI.GTP
		READ(FILHIG,*)aux,ICELL(I),X(I),Y(I),IBAC(I),ACEL(I),ACUR(I),SRIO(I),DECL(I),LCEL(I),HCEL(I),CELJUS(I),OD(I),hdFLAG(I),BRIO(I),HRIO(I),RUGMAN(I),(PUSO(I,K),K=1,NU)
		!READ(FILHIG,*)aux,ICELL(I),X(I),Y(I),IBAC(I),ACEL(I),ACUR(I),SRIO(I),DECL(I),LCEL(I),HCEL(I),CELJUS(I),OD(I),hdFLAG(I),(PUSO(I,K),K=1,NU) ! RP
		!READ(FILHIG,77)ICELL(I),X(I),Y(I),IBAC(I),ACUR(I),ACEL(I),HCEL(I),LCEL(I),SRIO(I),DECL(I),CELJUS(I),(PUSO(I,K),K=1,NU) 
		! WRITE(*,*)(PUSO(I,K),K=1,NU)
		QREF(I)=ACUR(I)*0.03 !Calculates reference initial discharge                                   !@Inlcuído por Hugo Fagundes em 18/08/17 e copiado da rotina do inercial
		! Verifica erros nas áreas dos blocos:
		if (sum(PUSO(I,:))>100.01.or.sum(PUSO(I,:))<99.99.or.minval(PUSO(I,:))<0.0.or.maxval(PUSO(I,:))>100.0) then
			write(*,*) 'Erro nas porcentagens dos HRUs da minibacia ',i
			write(*,*) 'Soma=',sum(PUSO(I,:)),'Min=',minval(PUSO(I,:)),'no HRU ',minloc(PUSO(I,:)),'Max=',maxval(PUSO(I,:)),' no HRU',maxloc(PUSO(I,:))
			!@DCBteste  read(*,*)
		endif

		!DEFINE EXUTORIOS
		!verifica qual é a celula que define o fim de uma sub-bacia
!		IF(IBAC(I).GT.IBANT)THEN
!			KW=IBAC(I)-1
!			IEXUT(KW)=I-1
!			IBANT=IBAC(I)
!		ENDIF

		! Definicao de exutorios considerando que subbacias nao estao ordenadas:
		! Considera o exutorio igual a minibacia de maior codigo.
		KW=IBAC(I)
		IF(IEXUT(KW)<ICELL(I)) IEXUT(KW)=ICELL(I)	
	

	!	*************************************************************
	!		espaço usado para modificar uso do solo

		!IF(IBAC(I)>4)THEN
		!	ACEL(I)=0.00001
		!ENDIF


	!	************************************************************
	ENDDO

!	IEXUT(KW+1)=NC !O EXUTORIO DA ULTIMA SUB-BACIA É A ULTIMA CELULA !RP verificar




	! Muda unidades: !RP
	! Passa declividade do rio para m/m
	DECL=DECL/1000.0  !RP

	! Corrige valores de comp. rio mais longo:
	do iC=1,nC
		LCEL(iC)=max(LCEL(iC),0.001)
	enddo

	CLOSE (FILHIG)
77	FORMAT(I5,2F15.3,I5,2F10.1,2I5,F10.2,F10.6,I5,6F5.1)
	RETURN
	END