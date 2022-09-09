	SUBROUTINE MODELO
	!Esta subrotina comanda o loop do tempo do modelo hidrológico e chama as 
	!rotinas de balanco e propagação nas células e de propagação na rede de drenagem
	USE VARS_MAIN
	USE AUX_MOD
!	
!************ DIOGO BUARQUE ********************
	USE SED_VARS
!    USE PAR1_MOD    !@ DCB_HD_Sed
!    USE TIME_MOD    !@ DCB_HD_Sed
!    USE PP_MOD      !@ DCB_HD_Sed
!************ DIOGO BUARQUE ********************
!
	IMPLICIT NONE

	INTEGER K,KHID !CONTADORES
	INTEGER KC,JB,JC,KB,JCB,MWP 	!CONTADORES
    integer count_daniel,dt_daniel
	integer iTwrite
!	
!************ DIOGO BUARQUE ********************
	INTEGER i
	REAL INmini(3), OUTmini(3)
!************ DIOGO BUARQUE ********************
!
	

! Tudo isto passou para a rotina condinic (Paiva)
!	!ZERANDO VARIAVEIS
!	PM=0.0
!	PM2=0.0
!	SIM=0.0 !GUARDA LAMINA INTERCEPTADA
!	EM=0.0 !GUARDA LAMINA EVAPOTRANSPIRAÇÃO
!	WBM=0.0 !LAMINA ARMAZENADA MEDIA
!	DSUM=0.0 !ESCOMENTO SUPERFICIAL DA BACIA
!	DINM=0.0 !ESCOMENTO INTERMEDIÁRIO DA BACIA
!	DBAM=0.0 !ESCOMENTO SUBTERRÂNEO DA BACIA
!		
!
!	!****************************CONDIÇÕES INICIAIS*****************************************
!	!DA INTERCEPTAÇÃO
!	SI=0.0 !LAMINA INTERCEPTADA É ZERO
!	!DO SOLO
!	DO IC=1,NC
!		IB=IBAC(IC)
!		DO IU=1,NU
!			W(IC,IU)=WM(IB,IU)*0.40 !CONSIDERA QUE SOLO TEM 40% AGUA NO INÍCIO DA SIMULAÇÃO
!		ENDDO
!	ENDDO
!	!DOS RESERVATORIOS SUBTERRANEO, SUBSUPERFICIAL E SUPERFICIAL
!	DO IC=1,NC
!		IB=IBAC(IC)
!		VBAS(IC)=QESP(IB)*ACEL(IC)*CB(IB)*3600.
!		VINT(IC)=0.0
!		VSUP(IC)=0.0
!	ENDDO
!	!DE TEMPERATURA
!	TA=20.0
!   count_daniel=1
!	OPEN(31,FILE='VAZAO EM CELULAS.TXT',STATUS='UNKNOWN')
!	
!	!*****************************************************************
!		! Modelo hidrodinâmico:
!		IF (hdFLAG0>0) CALL HD_CondInic
!	!*****************************************************************
!	
!	!****************************FIM DAS CONDIÇÕES INICIAIS*********************************


	!	INICIO DO LOOP DO TEMPO

	!IT=8767 !01/01/1982
	IT=0
    Dt_daniel=INT((NT-IT)/10)

	!WRITE(*,*)' DESEJA EVITAR A SIMULACAO DE ALGUMA SUB-BACIA? QUANTAS?'
	!READ(*,*)NCONGEL
	!WRITE(*,*)' QUAIS? INDIQUE LISTA SEPARADA POR VIRGULA.'
	!READ(*,*)(IBCONGEL(KB),KB=1,NCONGEL)

	!Write(*,*)'Percentagem da simulação completada'
    !Write(*,*)'10------50--------100'
 
!@ DCB ago/2012open(316666,FILE='.\output\P_SUB.bin',STATUS='UNKNOWN',RECL=NC,FORM='UNFORMATTED',ACCESS='DIRECT')
	
    itWrite = 0 !@ DCB ago/2012 (sem isso, o valor é grande no Fortran Intel)
    WRITE(*,*)  !@ DCB set/2012
    WRITE(*,*)'SIMULANDO'   !@ DCB set/2012
    WRITE(*,*)  !@ DCB set/2012
    DO WHILE (IT<NT)
		IT=IT+1
		if (icalib==0.and.itWrite<iT) then
    	WRITE(*,*) '     Passo de tempo: ', IT   !@ DCB set/2012 ,maxval(QJ2),maxloc(QJ2),minval(QJ2),minloc(QJ2),QJ2(nC) !QJ2(6848)! RP
!@ DCB ago/2012			write(7000,*) 'iT',iT
!			write(*,*) 'iT',iT
			itWrite=itWrite+1
        endif
        
		if(it==count_daniel)then
		   WRITE(*,701)'**'
			701		FORMAT(A2,$)
			count_daniel=Count_daniel+Dt_daniel
        endif
        if(it==NT)write(*,*)

		JDIA=IDINI+INT((IT+HORAINI-1)/(86400./DTP)) !VERIFICA QUAL É O DIA DO CALENDÁRIO 
		!SUBROTINA DE LEITURA E PREPARACAO DA CHUVA
		ITCHUVA=IT
		CALL LECHUVA
	
	
	
		DO KC=1,NC
			PM(KC)=PM(KC)+P(KC) !ACUMULA CHUVA
		ENDDO
	
		!SUBROTINA DE PREPARAÇÃO DOS DADOS CLIMATOL.
		TONTEM=TA
		CALL LECLIMA
	
		!SUBROTINA DA CELULA
		CALL CELULA
		
!
!
!************ DIOGO BUARQUE ********************
		CALL SED_BACIA
		!@ Resultado é o aporte de cada fracao de sedimento para a rede
		!@ PSC(IC,1:4) = 1 = total, 2 = areia, 3 = silte, 4 = argila
!************ DIOGO BUARQUE ********************
!
!
!		
	
		!************** NOSOLO.HIG ***************************************************************
		!As linhas abaixo servem para gravar algumas variáveis detalhadas de um bloco e de
		!uma célula. Os valores de JB e JC podem ser alterados, JB indica o bloco e JC a
		!célula em que se desejam os dados. Os dados são gravados num arquivo chamado NOSOLO.HIG.
		IF(ICALIB.EQ.0)THEN !NÃO ESTÁ CALIBRANDO, PODE GRAVAR
			JB=1
			JC=1
			!JC=57
			WRITE(FILSOL,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)
			JB=1
			JC=1
			WRITE(FILSOL2,75)IT,P(JC),W(JC,JB),SI(JC,JB),ET(JC,JB),CAF(JC,JB),QBAS(JC),QINT(JC),QSUP(JC)
			

!			write(971,66) (E0agua(iC),iC=1,nC)
!			write(972,66) (E0topo(iC),iC=1,nC)
!			write(973,66) (E0sup(iC),iC=1,nC)

		ENDIF
		!Fim da saida de dados para o arquivo NOSOLO.HIG

		!******************************************************************************************
	
		!SUBROTINA DA REDE DE DRENAGEM
		!CALL REDE(NC,QBAS,QINT,QSUP,SRIO,BRIO,DECL,CEL,QREF,NSUBT,DT,CELJUS,IT,QM2,QJ2,PMB2,PMI2,PMS2,PJB2,PJI2,PJS2,VRB,VRI,VRS)

!		write(*,*) maxval(TWS),maxloc(TWS),minval(TWS), minloc(TWS)
!		read(*,*)
!		read(*,*)

		CALL REDE

!		write(*,*) maxval(TWS),maxloc(TWS),minval(TWS), minloc(TWS),DECL(maxloc(TWS))
!		read(*,*)
!		read(*,*)
!*****************************************************************
		! Modelo hidrodinâmico:
		IF (hdFLAG0>0) CALL Hidrodinamico2
!*****************************************************************
!

!DO KC=1,NC
!    IF (hdFLAG(KC)>0) THEN
!        write(*,*) 'IC, Tr     = ', KC, CellTr(KC), hdFLAG0
!        write(*,*) 'IX1, IX2   = ', TrNST(CellTr(KC),1), TrNST(CellTr(KC),2)
!        write(*,*) 'L1, Z0, H0 = ', BA(1,TrNST(CellTr(KC),1)), ZO(TrNST(CellTr(KC),1)), HO(TrNST(CellTr(KC),1))
!        write(*,*) 'L2, Z0, H0 = ', BA(1,TrNST(CellTr(KC),2)), ZO(TrNST(CellTr(KC),2)), HO(TrNST(CellTr(KC),2))
!        write(*,*) 'S, S1, S2  = ', DECL(KC), SFF(TrNST(CellTr(KC),1)), SFF(TrNST(CellTr(KC),2))
!        pause
!    ENDIF
!ENDDO

!
!************ DIOGO BUARQUE ********************
		CALL SED_REDE
!************ DIOGO BUARQUE ********************
!
!
!
		
		!ARMAZENA DADOS DE VAZÃO DAS CÉLULAS EM QUE EXISTE VAZÃO OBSERVADA
		DO K=1,NOBS
			KHID=IQOBS(K) !CÉLULA QUE CORRESPONDE AO POSTO
			QR(K,IT)=QJ2(KHID) !QR(K,IT) VAI SER COMPARADO A QOBS(K,IT) NA ROTINA FOBJ
		ENDDO
	
		IF(ICALIB.EQ.0)THEN !SÓ ARMAZENA ESTES DADOS QUANDO NÃO ESTÁ CALIBRANDO
!			DO KB=1,NB	!ARMAZENA VAZOES DAS SUB-BACIAS
!				JCB=IEXUT(KB) !CELULA DO EXUTORIO DA SUB-BACIA KB
!				QRB(KB,IT)=QJ2(JCB)
!			ENDDO
	
			DO K=1,NUMHIDG !GUARDA DADOS PARA GRAVAR HIDROGRAMAS EM LOCAIS DEFINIDOS NO ARQUIVO PARHIG 
				KHID=IHIDG(K) !IHIDG(K) É O NÚMERO DA CÉLULA EM QUE SE DESEJA O HIDROGRAMA
				QRG(K,IT)=QJ2(KHID) !QRG ARAMZENA OS HIDROGRAMAS NOS LOCAIS DESEJADOS
!
!
!************ DIOGO BUARQUE ********************
                CAREIA(K,IT)=CSJ2(KHID,1)*(10.**(6.)) !@ (mg/l)
                CSILTE(K,IT)=CSJ2(KHID,2)*(10.**(6.)) !@ (mg/l)
                CARGIL(K,IT)=CSJ2(KHID,3)*(10.**(6.)) !@ (mg/l)
!                if (isNaN(CAREIA(K,IT)) .OR. isNaN(CSILTE(K,IT)) .OR. isNaN(CARGIL(K,IT))) then
!                    write(*,*) 'IT, IC     = ', IT, IC
!                    write(*,*) 'CSJ2       = ', CSJ2(KHID,1), CSJ2(KHID,2), CSJ2(KHID,3)
!                    write(*,*) 'Ca, Cs, Cc = ', CAREIA(K,IT), CSILTE(K,IT), CARGIL(K,IT)
!                    pause
!                endif
!************ DIOGO BUARQUE ********************
!
!
			ENDDO
!
!
!************ DIOGO BUARQUE set/2012 ********************
			DO i = 1,nSEDmini
			    KHID=CONTIC(i) !CODIGO DA MINI-BACIA EM QUE SE DESEJA O HIDROGRAMA
				Qmini(i,IT)=QJ2(KHID) !QRG ARAMZENA OS HIDROGRAMAS NOS LOCAIS DESEJADOS
			ENDDO
!************ DIOGO BUARQUE set/2012 ********************			
!
!			
	
			!ARMAZENA VAZÕES SEGUNDO A ORIGEM PARA UMA CÉLULA - quando não está calibrando
			MWP=119 !CÉLULA EM QUE SE DESEJAM OS RESULTADOS DE HIDROGRAMA SEPARADO POR ORIGEM
			
			MWP=2867 !CÉLULA EM QUE SE DESEJAM OS RESULTADOS DE HIDROGRAMA SEPARADO POR ORIGEM
			
			QB(IT)=QJ2(MWP)*PJB2(MWP)
			QBI(IT)=QB(IT)+QJ2(MWP)*PJI2(MWP)
			QBIS(IT)=QBI(IT)+QJ2(MWP)*PJS2(MWP)
		ENDIF

!@ DCB ago/2012        write(316666,REC=iT)(QJ2(iC),iC=1,nC)
!@ DCB ago/2012		write(974,67) (TWS(iC),iC=1,nC)
!@ DCB ago/2012		write(988,68) (E0media(iC),iC=1,nC)
!@ DCB ago/2012		write(989,68) (P(iC),iC=1,nC)

		
		QM2in=QM2
		if (IT==1) QM1in=QM2in
		do iC=1,nC

			IF (hdFLAG(IC)==0) QM2in(iC)=QM2(iC)-(QSUP(IC)+QINT(IC)+QBAS(IC))

			if (OD(iC)==1) then
				! Minibacia de cabeceira:
				DTWS(iC)=P(iC)+0.5*(-QJ1(iC)-QJ2(iC))*DTP/(ACEL(iC)*1000.0)-E0media(iC)
			else
				DTWS(iC)=P(iC)+0.5*(QM1in(iC)+QM2in(iC)-QJ1(iC)-QJ2(iC))*DTP/(ACEL(iC)*1000.0)-E0media(iC)
			endif
			! Tira valores absurdos em minibacias muito pequenas:
			if (ACEL(iC)<=10.0) DTWS(iC)=0.0

		enddo
		QM1in=QM2in

!@ DCB ago/2012		write(975,'(2F12.4)') Pbacia,Ebacia

	
	ENDDO !FIM DO LOOP DO TEMPO





!************ DIOGO BUARQUE ********************
OPEN(SEDSAI,FILE  = '.\output\SEDRIO_SAIDA.txt',STATUS='UNKNOWN')
OPEN(SEDDEP,FILE  = '.\output\SEDRIO_DEPOSITO.txt',STATUS='UNKNOWN')
OPEN(SEDERO,FILE  = '.\output\SEDRIO_EROSAO.txt',STATUS='UNKNOWN')
OPEN(BALFLP,FILE  = '.\output\BAL_PLAN.txt',STATUS='UNKNOWN')  !@ DCB_HD_Sed - SALVA ENTRADAS E SAIDA NAS PLANICIES
    WRITE(SEDSAI,'(1A10,3A16)') 'IC', 'SAI_A', 'SAI_S', 'SAI_C'
    WRITE(SEDDEP,'(1A10,3A16)') 'IC', 'DEP_A', 'DEP_S', 'DEP_C'
    WRITE(SEDERO,'(1A10,3A16)') 'IC', 'ERO_A', 'ERO_S', 'ERO_C'
    WRITE(BALFLP,'(1A10,9A16)') 'IC', 'IN_FLPsand', 'IN_FLPsilt', 'IN_FLPclay', 'OUT_FLPsand', 'OUT_FLPsilt', 'OUT_FLPclay', 'QFL_in', 'QFL_out', 'DFL'      !@ DCB_HD_Sed
    DO i = 1,nSEDmini
        WRITE(SEDSAI,'(1I10,3ES16.7)'), CONTIC(i), iSAI(i,1), iSAI(i,2), iSAI(i,3)      !@ ton
        WRITE(SEDDEP,'(1I10,3ES16.7)'), CONTIC(i), DEPT(i,1), DEPT(i,2), DEPT(i,3)      !@ ton
        WRITE(SEDERO,'(1I10,3ES16.7)'), CONTIC(i), EROT(i,1), EROT(i,2), EROT(i,3)      !@ ton
        WRITE(BALFLP,'(1I10,9ES16.7)'), CONTIC(i), inFL(i,1), inFL(i,2), inFL(i,3), outFL(i,1), outFL(i,2), outFL(i,3), BalQFL(i,1), BalQFL(i,2), DFL(i,2)+DFL(i,3)  !@ DCB_HD_Sed (ton)
    ENDDO
CLOSE(SEDSAI)
CLOSE(SEDDEP)
CLOSE(SEDERO)
CLOSE(BALFLP)



INmini(1)  = sum(iENTRA(:,1))/8
INmini(2)  = sum(iENTRA(:,2))/8
INmini(3)  = sum(iENTRA(:,3))/8
OUTmini(1) = sum(iSAI(:,1))/8
OUTmini(2) = sum(iSAI(:,2))/8
OUTmini(3) = sum(iSAI(:,3))/8
write(*,*) 
write(*,*) 
write(*,*) 'VERFICACAO DE CONSERVACAO DE SEDIMENTOS' 
write(*,*) 
write(*,*) 'ENTRADA (ton/dia)  = ', INmini(:)
write(*,*) 'SAIDA   (ton/dia)  = ', OUTmini(:)
write(*,*) 'ERRO    (%)        = ', 100.*(INmini(1)-OUTmini(1))/INmini(1), 100.*(INmini(2)-OUTmini(2))/INmini(2), 100.*(INmini(3)-OUTmini(3))/INmini(3)
write(*,*)
write(*,*)
pause


!************ DIOGO BUARQUE ********************

75	FORMAT(I6,8F10.4)


66	FORMAT(<nC>F10.4)
67	FORMAT(<nC>F10.1)
68	FORMAT(<nC>F10.3)
69	FORMAT(<nC>F12.3)

     
	RETURN
	END