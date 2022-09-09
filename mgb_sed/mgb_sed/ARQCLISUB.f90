	SUBROUTINE ARQCLISUB
	!SUBROTINA QUE LÊ DADOS DOS ARQUIVOS DE DADOS CLIMÁTICOS
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER KB,KLI,K

	!*************************************************************
	!VERIFICA O NUMERO DE CELULAS EM CADA SUB-BACIA
	KCB(1)=IEXUT(1)
	DO KB=2,NB
		KCB(KB)=IEXUT(KB)-IEXUT(KB-1)
	ENDDO
	!FIM DA VERIFICACAO DO NUMERO DE CELULAS EM CADA SUB-BACIA




	if (flagaclimed==1) then !se usa CRU - Daniel

      TD(:,:)=-9999.0
	  UD(:,:)=-9999.0
	  VD(:,:)=-9999.0
	  SOLD(:,:)=-9999.0
	  PAD(:,:)=-9999.0

	else !se nao usa CRU

		!________________________________________________________________
		!___ABRE OS ARQUIVOS DE DADOS CLIMATOLOGICOS DIARIOS______
		DO KLI=1,NCLI
			NARQ(KLI)=KLI+30	!NUMERO DO ARQUIVO
			OPEN(NARQ(KLI),FILE='.\input\'//ARQCLI(KLI),STATUS='OLD')
			READ(NARQ(KLI),733)(CABE(K),K=1,4) !LE CABECALHO
		ENDDO

		DO KLI=1,NCLI
			DO IT=1,NT
				READ(NARQ(KLI),734)TD(KLI,IT),UD(KLI,IT),VD(KLI,IT),SOLD(KLI,IT),PAD(KLI,IT)
			ENDDO
		ENDDO
		DO KLI=1,NCLI
			NARQ(KLI)=KLI+30	!NUMERO DO ARQUIVO
			CLOSE (NARQ(KLI))
		ENDDO
	endif


733	FORMAT(4A10)
734	FORMAT(5F10.2)
	!__________________________________________________________

	RETURN
	END