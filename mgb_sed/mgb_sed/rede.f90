	SUBROUTINE REDE
	!ESTA ROTINA FAZ A PROPAGAÇÃO NA REDE DE DRENAGEM
	USE VARS_MAIN
	USE AUX_MOD
	IMPLICIT NONE
	REAL DX		!DISCRETIZACAO NO ESPACO
	REAL XXX,AK		!PARAMETROS X E K DO MUSKINGUM CUNGE
	REAL QMX1,QMX2,QJX1,QJX2,QAUX,DTCAL,DEN,VT 	!VARIAVEIS AUXILIARES
	INTEGER JUS
	INTEGER NTRECH,NTC !NUMERO DE SUBTRECHOS, NUMERO DE SUBINTERVALOS
	REAL C1,C2,C3 	!COEFICIENTES MUSKINGUM
	!PROPORÇÕES DE ORIGEM
	REAL PMB1(NC+1),PMI1(NC+1),PMS1(NC+1)
	REAL PJB1(NC),PJI1(NC),PJS1(NC)
	INTEGER NTMUSK,LT
	INTEGER I,IS !CONTADORES

	real DECLref,QJref,QMref

	IF(DT(NC)<0.00001)STOP 'ALGO ERRADO COM INTERVALO DE TEMPO NA SUBROTINA REDE'
	NTMUSK=DTP/DT(NC) !

!Condicoes inicias passam para a rotina condinic
!********************************************************************************
!	!CONDIÇÕES INICIAIS DA VAZAO NO RIO
!	IF(IT.EQ.1)THEN
!		
!		! OBS RP Melhorar condicoes iniciais do MGB, usar QREF no lugar de vazao zero nos rios ********************
!
!		QM2=0.0
!		PMB2=0.0
!		PMI2=0.0
!		PMS2=0.0
!		QM1=0.0
!		QJ1=0.0
!		QJ2=0.0
!		QCEL1=0.0
!		QCEL2=0.0
!		QRIOINI=0.0
!		DO IC=1,NC
!*****************************************************************
			! Testa se propagação na célula é computada com modelo hidrodinâmico:
!			IF (hdFLAG(IC)==1) CYCLE
!*****************************************************************
!			
!			IF(NSUBT(IC).GT.0)THEN	!CELULA  COM RIO
!
!				!TRANSFORMA PROPORÇÕES EM VAZÕES
!				PMB2(IC)=PMB2(IC)*QM2(IC) !VAZAO SUBTERRANEA
!				PMI2(IC)=PMI2(IC)*QM2(IC) !VAZAO SUB-SUPERFICIAL
!				PMS2(IC)=PMS2(IC)*QM2(IC) !VAZÃO SUPERFICIAL
!
!				!ACUMULA VAZAO
!				QM2(IC)=QM2(IC)+QSUP(IC)+QINT(IC)+QBAS(IC)

!				!CONDIÇÃO DE CONTORNO MUSKINGUM-CUNGE
!				DO LT=1,NTMUSK+1
					!QCONTORM(IC,LT)=QM1(IC)+(LT-1)*(QM2(IC)-QM1(IC))/NTMUSK
!					QCONTORM(IC,LT)=QCONTORM(IC,LT)+QM2(IC)
!				ENDDO

				!ATUALIZA PROPORÇÕES
!				PMB2(IC)=(PMB2(IC)+QBAS(IC))/QM2(IC)
!				PMI2(IC)=(PMI2(IC)+QINT(IC))/QM2(IC)
!				PMS2(IC)=(PMS2(IC)+QSUP(IC))/QM2(IC)
!
!				JUS=CELJUS(IC)
!				QJ2(IC)=QM2(IC)	!VAZAO CONSTANTE NO TRECHO
!
!				!CONDIÇÃO DE CONTORNO MUSKINGUM-CUNGE
!				DO LT=1,NTMUSK+1
!					!QCONTORJ(IC,LT)=QJ1(IC)+(LT-1)*(QJ2(IC)-QJ1(IC))/NTMUSK
!					QCONTORJ(IC,LT)=QCONTORJ(IC,LT)+QJ2(IC)
!				ENDDO
!				
!				!A PROPORÇÃO QUE ENTRA SAI
!				PJB2(IC)=PMB2(IC)
!				PJI2(IC)=PMI2(IC)
!				PJS2(IC)=PMS2(IC)
!	            !write(*,*)jus
!				PMB2(JUS)=PMB2(JUS)*QM2(JUS) !VAZÃO SUBTERRANEA
!				PMI2(JUS)=PMI2(JUS)*QM2(JUS) !VAZÃO SUB-SUPERFICIAL
!				PMS2(JUS)=PMS2(JUS)*QM2(JUS) !VAZÃO SUPERFICIAL
!
!				!O QUE SAI DO RIO DE UMA CEL. VAI PRA OUTRA
!				QM2(JUS)=QM2(JUS)+QJ2(IC) 
!
!				!PONDERA PROPORÇÕES
!				PMB2(JUS)=(PMB2(JUS)+PJB2(IC)*QJ2(IC))/QM2(JUS)  
!				PMI2(JUS)=(PMI2(JUS)+PJI2(IC)*QJ2(IC))/QM2(JUS)
!				PMS2(JUS)=(PMS2(JUS)+PJS2(IC)*QJ2(IC))/QM2(JUS)
!
!				!TENTATIVA DE PROPORÇÕES
!				QAUX=(QM2(IC)+QJ2(IC))/2
!				VRB(IC)=0.3333*(1000.*SRIO(IC)*QAUX/CEL(IC))
!				VRI(IC)=0.3333*(1000.*SRIO(IC)*QAUX/CEL(IC))
!				VRS(IC)=0.3333*(1000.*SRIO(IC)*QAUX/CEL(IC))
!
!			ELSE !CELULA FONTE
!				JUS=CELJUS(IC)
!
!				PMB2(JUS)=PMB2(JUS)*QM2(JUS) !VAZÃO SUBTERRANEA
!				PMI2(JUS)=PMI2(JUS)*QM2(JUS) !VAZÃO SUB-SUPERFICIAL
!				PMS2(JUS)=PMS2(JUS)*QM2(JUS) !VAZÃO SUPERFICIAL
!
!				QM2(JUS)=QM2(JUS)+QSUP(IC)+QINT(IC)+QBAS(IC)
!
!				!CONDIÇÃO DE CONTORNO MUSKINGUM-CUNGE
!				DO LT=1,NTMUSK+1
!					!QCONTORM(JUS,LT)=QM1(JUS)+(LT-1)*(QM2(JUS)-QM1(JUS))/NTMUSK
!					QCONTORM(JUS,LT)=QCONTORM(JUS,LT)+QM2(JUS)
!				ENDDO
!
!				!PONDERA PROPORÇÕES
!				PMB2(JUS)=(PMB2(JUS)+QBAS(IC))/QM2(JUS)  
!				PMI2(JUS)=(PMI2(JUS)+QINT(IC))/QM2(JUS)
!				PMS2(JUS)=(PMS2(JUS)+QSUP(IC))/QM2(JUS)
!
!			ENDIF
!		ENDDO
!	ENDIF
!	!FIM DAS CONDIÇÕES INICIAIS

	!EM CADA INTERVALO DE TEMPO O QUE ERA i+1 VIRA i
	DO IC=1,NC
		QM1(IC)=QM2(IC)
		QJ1(IC)=QJ2(IC)
		QM2(IC)=0.0
		QJ2(IC)=0.0
		PMB1(IC)=PMB2(IC)
		PMI1(IC)=PMI2(IC)
		PMS1(IC)=PMS2(IC)
		PJB1(IC)=PJB2(IC)
		PJI1(IC)=PJI2(IC)
		PJS1(IC)=PJS2(IC)
		QCEL1(IC)=QCEL2(IC)
	ENDDO
	QCONTORM=0.0
	QCONTORJ=0.0

	DO IC=1,NC

		if (sbtFLAG(iC)==1) cycle ! RP nao computa celulas a montante de celula com substituicao de dados

		!write(*,*),IC
		IB=IBAC(IC)
!@ DCB_sed ###########################################################################
	    !@ Calcula apenas sub-bacia de interesse.
		!IF(IB>SUBfim.OR.IB<SUBini)CYCLE
!@ DCB_sed ###########################################################################


		QCEL2(IC)=QBAS(IC)+QINT(IC)+QSUP(IC) !SOMA VAZOES GERADAS NA CELULA


!*****************************************************************
		! Testa se propagação na célula é computada com modelo hidrodinâmico:
		IF (hdFLAG(IC)==1) CYCLE
!*****************************************************************
	


		IF(NSUBT(IC).GT.0)THEN !CELULA COM RIO			
			DTCAL=DT(IC)
			DX=SRIO(IC)*1000./NSUBT(IC)
			AK=DX/CEL(IC)
!			AK=DX/(CEL(IC)/3.0) ! Teste

			XXX=0.5*(1.0-(QREF(IC)/(BRIO(IC)*DECL(IC)*CEL(IC)*DX)))


			
!			write(6661,*) ic,xxx,ak/dtcal,decl(ic)*1000.0,cel(ic),Qref(ic),BRIO(IC),DX

			! Verificação de valores absurdos de X e K: !RP
			
			! Testa condição 0.2<=X<=0.5:
			if (XXX<0.2) then
!				write(*,*) ic,xxx,ak/dtcal
!				write(*,*) decl(ic)*1000.0,cel(ic),Qref(ic),BRIO(IC),DX
				XXX=0.2
				AK=DTCAL/(3.125*XXX**1.25) !EQUAÇÃO 4.111 DO LIVRO DO TUCCI MODELOS
			elseif(XXX<0.4) then
!				AK=DTCAL/(3.125*XXX**1.25) !EQUAÇÃO 4.111 DO LIVRO DO TUCCI MODELOS
			elseif (XXX<0.5) then
!				write(*,*) ic,xxx,ak,decl(ic),cel(ic)
!				AK=DTCAL
			else
!				write(*,*) ic,xxx,ak,decl(ic),cel(ic)
				XXX=0.5
				AK=DTCAL
			endif
			! Testa condição 2X<=dt/K<=2(1-X):
!			if (DTCAL/AK<=2*XXX) then
!				write(*,*) '1',iC,XXX,AK,DTCAL
!				read(*,*)
!				AK=DTCAL/(2*XXX)
!			elseif (DTCAL/AK>=2*(1-XXX)) then
!				write(*,*) '2',iC,XXX,AK,DTCAL
!				read(*,*)
!				AK=DTCAL/(2*(1-XXX))
!			endif		
		
		
			DEN=2.*AK*(1.-XXX)+DTCAL			
			C1=(2.*AK*XXX+DTCAL)/DEN
			C2=(DTCAL-2.*AK*XXX)/DEN
			C3=(2.*AK*(1.-XXX)-DTCAL)/DEN

			! Corrige para ter somente coeficientes positivos:
!			C1=max(C1,0.0)
!			C2=max(C2,0.0)
!			C3=max(C3,0.0)

!			if (min(C1,C2,C3)<0.0) then
!				write(*,*) ic,C1,C2,C3,xxx,ak,DECL(IC),dx
!			
!				read(*,*)
!			endif
			!AS VAZOES QUE SAO GERADAS NAS CELULAS COM RIO
			! ENTRAM A MONTANTE DA PROPRIA CELULA

			!TRANSFORMA PROPORÇÕES EM VAZÕES
			PMB2(IC)=PMB2(IC)*QM2(IC) !VAZAO SUBTERRANEA
			PMI2(IC)=PMI2(IC)*QM2(IC) !VAZAO SUB-SUPERFICIAL
			PMS2(IC)=PMS2(IC)*QM2(IC) !VAZÃO SUPERFICIAL
			
!			if (IC==6848) write(*,*) 'IC,QM2,QSUP,QINT,QBAS',IC,QM2(IC),QSUP(IC),QINT(IC),QBAS(IC)
			QM2(IC)=QM2(IC)+QSUP(IC)+QINT(IC)+QBAS(IC)

			!WRITE(*,*)'CONTORNO DE MONTANTE CELULA COM RIO'

			!CONDIÇÃO DE CONTORNO MUSKINGUM-CUNGE
			DO LT=1,NTMUSK+1
				QCONTORM(IC,LT)=QCONTORM(IC,LT)+QCEL1(IC)+(LT-1)*(QCEL2(IC)-QCEL1(IC))/NTMUSK
				!WRITE(*,*)IC,LT,QCONTORM(IC,LT)
			ENDDO

			!ATUALIZA PROPORÇÕES
			PMB2(IC)=(PMB2(IC)+QBAS(IC))/QM2(IC)
			PMI2(IC)=(PMI2(IC)+QINT(IC))/QM2(IC)
			PMS2(IC)=(PMS2(IC)+QSUP(IC))/QM2(IC)

			!ATUALIZA PROPORÇÕES NO VOLUME DO RIO
			VRB(IC)=VRB(IC)+PMB2(IC)*QM2(IC)*DTP	!VRB(IC)=VRB(IC)+PMB2(IC)*QM2(IC)*3600.*24.
			VRI(IC)=VRI(IC)+PMI2(IC)*QM2(IC)*DTP	!VRI(IC)=VRI(IC)+PMI2(IC)*QM2(IC)*3600.*24.
			VRS(IC)=VRS(IC)+PMS2(IC)*QM2(IC)*DTP	!VRS(IC)=VRS(IC)+PMS2(IC)*QM2(IC)*3600.*24.

			QMX1=QM1(IC)
			QMX2=QM2(IC)
			QJX1=QJ1(IC)
			NTRECH=NSUBT(IC)
			NTC=DTP/DT(IC)			!NTC=3600.*24./DT(IC) !NUMERO DE SUBINTERVALOS DE TEMPO

		
!			if (IC==6848) then
!				write(*,*) C1,C2,C3,ntrech,ntc,QM1(IC),QM2(IC),QJ1(IC)
!			endif
			IF(ICODMUSK(IC)==0)THEN !MUSKINGUN CUNGE LINEAR
				CALL MUSK(QMX1,QMX2,QJX1,QJX2,NTRECH,NTC,C1,C2,C3)
			ELSEIF(ICODMUSK(IC)==1)THEN !MUSKINGUN CUNGE NAO LINEAR
!				CALL MUSK_NL(QMX1,QMX2,QJX1,QJX2,NTRECH,NTC,C1,C2,C3)
			ENDIF
			QJ2(IC)=QJX2

			! Testa se nao é trecho curto:
!			if (SRIO(IC)<2000) then ! RP
!				QCONTORJ(IC,:)=QCONTORM(IC,:)
!				QJ2(IC)=QM2(IC)
!			endif	! RP
	
			!ATUALIZA PROPORÇÕES NO VOLUME DO RIO
			VT=VRB(IC)+VRI(IC)+VRS(IC)
			PJB2(IC)=VRB(IC)/VT
			PJI2(IC)=VRI(IC)/VT
			PJS2(IC)=VRS(IC)/VT
			VRB(IC)=MAX(0.0,VRB(IC)-PJB2(IC)*QJ2(IC)*DTP)	!VRB(IC)=MAX(0.0,VRB(IC)-PJB2(IC)*QJ2(IC)*3600.*24.)
			VRI(IC)=MAX(0.0,VRI(IC)-PJI2(IC)*QJ2(IC)*DTP)	!VRI(IC)=MAX(0.0,VRI(IC)-PJI2(IC)*QJ2(IC)*3600.*24.)
			VRS(IC)=MAX(0.0,VRS(IC)-PJS2(IC)*QJ2(IC)*DTP)	!VRS(IC)=MAX(0.0,VRS(IC)-PJS2(IC)*QJ2(IC)*3600.*24.)

			!AS VAZOES PROPAGADAS NA CELULA ENTRAM COMO CONTORNO
			!DE MONTANTE DA CELULA QUE ESTA A JUSANTE
			JUS=CELJUS(IC)
			
!			if (JUS==6848) write(*,*) IC,JUS,QM2(JUS),QJ2(IC)
			! Testa se mini-bacia é exutorio:
			if (JUS>0) then

				PMB2(JUS)=PMB2(JUS)*QM2(JUS) !VAZÃO SUBTERRANEA
				PMI2(JUS)=PMI2(JUS)*QM2(JUS) !VAZÃO SUB-SUPERFICIAL
				PMS2(JUS)=PMS2(JUS)*QM2(JUS) !VAZÃO SUPERFICIAL

				QM2(JUS)=QM2(JUS)+QJ2(IC)

				!WRITE(*,*)'CONTORNO DE JUSANTE CELULA COM RIO'


				!CONDIÇÃO DE CONTORNO MUSKINGUM-CUNGE
				!VERIFICA SE A VAZÃO CALCULADA NA CÉLULA DEVE SER SUBSTITUÍDA POR HIDROGRAMA LIDO
				IS=0
				DO I=1,NUMSUBST 
					IF(IT>1.AND.IC==ISUBST(I))THEN 
						IS=I
						EXIT
					ELSE
						IS=0
					ENDIF
				ENDDO

				IF(IS>0)THEN
					IF(QLIDO(IS,IT)<0.0)then !falha nos dados lidos
						!STOP 'VAZÃO LIDA PARA SUBST TEM FALHA - PROGRAMA TERMINOU'
						QM2(JUS)=QM2(JUS) !NÃO FAZ NADA, USA VALOR CALCULADO
						DO LT=1,NTMUSK+1
							QCONTORM(JUS,LT)=QCONTORM(JUS,LT)+QCONTORJ(IC,LT) !O CONTORNO DE MONTANTE DA CELULA JUS RECEBE A SAIDA DO HIDROGRAMA RECEM CALCULADO
							!WRITE(*,*)IC,LT,QCONTORJ(IC,LT),QCONTORM(JUS,LT)
						ENDDO
					ELSE
						QM2(JUS)=QLIDO(IS,IT)
						DO LT=1,NTMUSK+1 !SUBSTITUI VALORES CALCULADOS POR VALORES LIDOS
							!QCONTORM(JUS,LT)=QLIDO(IS,IT-1)+(LT-1)*(QLIDO(IS,IT)-QLIDO(IS,IT-1))/NTMUSK
							QCONTORM(JUS,LT)=QLIDO(IS,IT-1)+(LT-1)*(QLIDO(IS,IT)-QLIDO(IS,IT-1))/NTMUSK+QCONTORM(JUS,LT) ! Correção Walter
						ENDDO
					ENDIF
				ELSE
					DO LT=1,NTMUSK+1
						QCONTORM(JUS,LT)=QCONTORM(JUS,LT)+QCONTORJ(IC,LT) !O CONTORNO DE MONTANTE DA CELULA JUS RECEBE A SAIDA DO HIDROGRAMA RECEM CALCULADO
						!WRITE(*,*)IC,LT,QCONTORJ(IC,LT),QCONTORM(JUS,LT)
					ENDDO
				ENDIF

				!PONDERA PROPORÇÕES
				PMB2(JUS)=(PMB2(JUS)+PJB2(IC)*QJ2(IC))/QM2(JUS)  
				PMI2(JUS)=(PMI2(JUS)+PJI2(IC)*QJ2(IC))/QM2(JUS)
				PMS2(JUS)=(PMS2(JUS)+PJS2(IC)*QJ2(IC))/QM2(JUS)
			endif	


			! Soma volume d´água no rio considerando profundidades 
			! das extremidades de jusante e montante calculadas por Manning 
			! com declividade constante:


			! Usa uma declividade minima se parametros do MC foram corrigidos. 
			DECLref=max(0.0001,DECL(iC))
			QJref=max(0.0,QJ2(iC))
			QMref=max(0.0,QM2(iC))
			! Evita declividades muito pequenas onde Manning não funciona bem.
!			if (iC==6760) write(*,*) 'res',TWS(ic),BRIO(IC),SRIO(iC),QM2(iC),QJ2(iC),DECLref,(QJ2(iC)*0.030/BRIO(iC))**0.6*DECLref**-0.3
			TWS(iC)=TWS(iC)+0.5*BRIO(IC)*SRIO(IC)*1000.0*((QMref*0.030/BRIO(iC))**0.6*DECLref**-0.3+(QJref*0.030/BRIO(iC))**0.6*DECLref**-0.3)/(ACEL(IC)*1000.0)


			TWS5(iC)=TWS5(iC)+0.5*BRIO(IC)*SRIO(IC)*1000.0*((QMref*0.030/BRIO(iC))**0.6*DECLref**-0.3+(QJref*0.030/BRIO(iC))**0.6*DECLref**-0.3)/(ACEL(IC)*1000.0)

		ELSE  !CELULA SEM RIO
			DX=0.0
			AK=0.0
			XXX=0.0

!@ DCB 26-05-2011 **************************************************************************
!            QJ2(IC)=QBAS(IC)+QINT(IC)+QSUP(IC)	!@ Soma das vazões geradas na célula
			QJ2(IC)=max(QBAS(IC)+QINT(IC)+QSUP(IC),0.0)	!@ Soma das vazões geradas na célula
!@ DCB 26-05-2011 **************************************************************************

			!@ DETERMINA PROPORÇÕES
			PJB2(IC)=QBAS(IC)/QJ2(IC)
			PJI2(IC)=QINT(IC)/QJ2(IC)
			PJS2(IC)=QSUP(IC)/QJ2(IC)

			! AS VAZOES QUE SAO GERADAS NAS CELULAS FONTE ENTRAM
			! A MONTANTE DO RIO QUE ESTÁ A JUSANTE
			JUS=CELJUS(IC)
			! Testa se mini-bacia é exutorio:
			if (JUS>0) then

				PMB2(JUS)=PMB2(JUS)*QM2(JUS) !VAZÃO SUBTERRANEA
				PMI2(JUS)=PMI2(JUS)*QM2(JUS) !VAZÃO SUB-SUPERFICIAL
				PMS2(JUS)=PMS2(JUS)*QM2(JUS) !VAZÃO SUPERFICIAL

				QM2(JUS)=QM2(JUS)+QSUP(IC)+QINT(IC)+QBAS(IC)
				

				!CONDIÇÃO DE CONTORNO MUSKINGUM-CUNGE
				DO LT=1,NTMUSK+1
					QCONTORM(JUS,LT)=QCONTORM(JUS,LT)+(QJ1(IC)+(LT-1)*(QJ2(IC)-QJ1(IC))/NTMUSK)
					!WRITE(*,*)JUS,LT,QCONTORM(JUS,LT)
				ENDDO


				!PONDERA PROPORÇÕES
				PMB2(JUS)=(PMB2(JUS)+QBAS(IC))/QM2(JUS)  
				PMI2(JUS)=(PMI2(JUS)+QINT(IC))/QM2(JUS)
				PMS2(JUS)=(PMS2(JUS)+QSUP(IC))/QM2(JUS)
			endif
		ENDIF


	ENDDO
	RETURN
	END