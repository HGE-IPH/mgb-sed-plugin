	SUBROUTINE COEF1
	! Calcula coeficientes das equações dos trechos e condições de contorno
	!-----------------------------------------------------------------------
	! Descrição das variáveis locais:
	!
	! I,J,L,K,M = variáveis auxiliares
	! H1 = profundidade d'água na seção I
	! FINT e FINT2 = resultado de interpolação
	! DF = derivada da rugosidade em relação a prof.
	! HAUX = prof. auxiliar (H1+.01)
	! RAUX = raio hidráulico auxiliar
	! AUX = área molhada auxiliar
	! DA = derivada da área em relação a prof. dA/dy
	! DR = derivada do raio hidráulico em relação a prof. dR/dy
	! FF = coef. de rugosidade de Manning para H1 dFF/dy
	! L1,L2,L3,L4,L5 = localizacao dos coeficientes no vetor AA
	! NS = indice da seção de condição de contorno
	! DF1 = derivada da vazão em relação a prof. na relacao de curva chave dQ/dy
	! R2 = raio hidráulico para H1+0.01
	! COC = auxiiar quando curva chave é calculada por manning
	! J1 = contador
	! J2 = identifica linha da matriz de coeficientes
	! LT = indice da barragem
	! AXY = derivada da largura em relação a distancia dB/dx
	! V1 e V2 = velocidades do escoamento nas seções no tempo anterior
	! CS,CS1,CS2,CS3,CS4,CS5,CS6,CS7,CS8,CS9,CS10,CS11,CS12, CS13 = coef. auxiliares para 
	!			calculo dos coef. das eq. dinâmica e continuidade discretizadas.
	! T1,T2,T3,T4,T5,T6 =  coef. auxiliares para calculo dos coef. das eq. dinâmica 
	!			e continuidade discretizadas. 
	! H2 = profundidade na seção de jusante em trecho com CC interna 
	! JP = indice da curva de descarga da barragem
	! IPB = Numero de pontos referente a secao de montante na tabela de curva de descarga da barragem
	! DFD = derivada da vazão em relação a prof. de montante em um trecho com barragem (dQ/dH1)
	! DF2 = derivada da vazão em relação a prof. de jusante em um trecho com barragem (dQ/dH2)
	! KM = coluna do coeficiente
	! JA = indice dos pares de seções das confluencias 
	! K1, K2, K3 = números das seções das confluências 
	!----------------------------------------------------------------------
	!
	! Declaração de variáveis:
	USE PAR1_MOD
	USE MAT_MOD
	USE TIME_MOD
	USE PP_MOD
	USE BARR_MOD
	USE CHEIA_MOD
	USE AUX_MOD
	
	USE MAT_MOD ! Tirar depois
	IMPLICIT NONE

	! Variáveis locais de cálculo:

	integer:: I,J,K,L,M
	real:: H1,FINT,FINT2,DF,HAUX,RAUX,AUX,DA,DR,FF
	integer:: L1,L2,L3,L4,L5,NS
	real:: DF1,R2,COC
	integer:: J1,J2,LT
	real:: AXY,V1,V2
	real:: CS,CS1,CS2,CS3,CS4,CS5,CS6,CS7,CS8,CS9,CS10,CS11,CS12,CS13
	real:: T1,T2,T3,T4,T5,T6
	real:: H2
	real:: DFD,DF2
	integer:: JP,IPB,KM,JA,K1,K2,K3,JJ
	real:: Fr,Fr2,AFIaux(nX)

	!----------------------------------------------------------------------------

	dafimax=0.0
	dafimax2=0.0
	afimax=0.0
	idafimax=0
	idafimax2=0
	iafimax=0


	nHmin=0

	! Calcula variáveis do escoamento em cada secao

	DO I=1,NX
		
		! Corrige profundidade minima:
!		if (HO(I)<Hmin) HO(I)=Hmin 

		! Número de pontos da tabela nivel x largura x area molhada x raio hidraulico
		L=NP(I)
		! Profundidade d'água:
		H1=HO(I)
		! Inicializa largura da planicie de inundação
		AFI(I)=0.0
		! Se dados em termos de niveis d'água, corrige H1:
		
		IF(LD/=0) H1=H1+ZO(I)


		! Calcula largura, area molhada e raio hidraulico:
		
		! Testa tipo de seção:
		IF (L>0) THEN
			! Dados detalhados de seção transversal
			! Interpola largura,area molhada e raio hidraulico para H1:
			T(I)= FINT(HA(1:L,I),BA(1:L,I),L,H1)
			A(I)=FINT(HA(1:L,I),AR(1:L,I),L,H1)
			R(I)=FINT(HA(1:L,I),RR(1:L,I),L,H1)
		
			! Profundidade auxiliar para calculo de derivadas numéricas:
			HAUX =H1+.01
			! Raio hidraulico e área auxiliares:
			RAUX=FINT(HA(1:L,I),RR(1:L,I),L,HAUX)
			AUX=FINT(HA(1:L,I),AR(1:L,I),L,HAUX)
		
	
		ELSE
			! Seção tipo retangular:
			! Calcula largura,area molhada e raio hidraulico:
			!*********************
			IF (HO(I)<=0) THEN
				HO(I)=0.5*HA(1,I)
				H1=HO(I)
				IF(LD/=0) H1=H1+ZO(I)
			ENDIF
			!********************
			T(I)=BA(1,I)
			A(I)=T(I)*HO(I)
			R(I)=A(I)/(T(I)+2*min(HA(1,I),HO(I)))
		
			! Profundidade auxiliar para calculo de derivadas numéricas:
			HAUX=HO(I)+.01
			! Raio hidraulico e área auxiliares:
			AUX=T(I)*HAUX
			RAUX=AUX/(T(I)+2*min(HA(1,I),HAUX))
		ENDIF

		! Testa prof. Negativa:
!		if (HO(I)<=0) then
	
!			R(I)=0.5
!			T(I)=1.0
!			A(I)=T(I)*0.5
			! Profundidade auxiliar para calculo de derivadas numéricas:
!			HAUX=HO(I)+.01
			! Raio hidraulico e área auxiliares:
!			AUX=T(I)*0.51
!			RAUX=0.51

!		endif
!******************************************************************8
! REVER SE ISSO FICA
		! Verifica erros e corrige raio hidraulico: 
!		IF (R(I)<=0)THEN  !CORRIGE VALOR DO RAIO HIDRAULICO
!			WRITE(6,'(I10,4F10.2)')I,H1,A(I),T(I),R(I)
!			R(I)=A(I)/T(I)
!		ENDIF
		! Verifica erros e corrige A,T e R
!		IF((A(I)<=0).OR.(T(I)<=0).OR.(R(I)<=0))THEN
!			IF(HO(I)<0.)THEN
!				HO(I)=1.0
!			ENDIF		
!			T(I)=QO(I)/SQRT(9.81*HO(I)**3)
!			A(I)=T(I)*HO(I)
!			R(I)=HO(I)
!		ENDIF
!***********************************************************************8

		! Calcula coef. de rugosidade de Manning para H1:
		IF(N(I)/=0)THEN
			! Manning variavel
			FF=FINT(HUX(1:N(I),I),F(1:N(I),I),N(I),H1)
		ELSE
			! Manning constante:
			FF=F(1,I)
			DF=0.0
		ENDIF
		
		! Aumenta rugosidade para baixas profundidades:
		!if (HO(I)<=min(2.0,0.1*HA(1,I))) then
!		if (HO(I)<=1.0) then
!			FF=FF+(1.0-HO(I))**2*10000 ! Manning 
!			DF=-2*(1.0-HO(I))*10000.0
!		endif  

		
		! Testa vazao exatamente igual a zero (causa problemas numéricos): ! RP
		if (QO(I)==0.0) QO(I)=0.0000000001
		
		! Condutancia hidraulica para H1:
		CK(I)=(R(I)**(2./3.)*A(I)*CONST/FF)
		! Declividade da linha de energia:
		SFF(I)=QO(I)*ABS(QO(I))/CK(I)**2.   !@ DCB_HD_Sed (Antes o expoente era 2, inteiro)


		! Derivada da área em relação a prof. dA/dy:
		DA=(AUX-A(I))/.01
		! Derivada do coef. de rugosidade de Manning em relacao a prof. dFF/dy:
		IF(N(I)/=0)THEN
			DF=(FINT(HUX(1:N(I),I),F(1:N(I),I),N(I),HAUX)-FF)/.01
		ENDIF
		! Derivada do raio hidráulico em relação a prof. dR/dy:
		DR=(RAUX-R(I))/.01
		
		! Derivada do condutancia hidráulica em relação a prof. dK/dy:
		CKY(I)=CK(I)*(DA/A(I)+(2./3.)*DR/R(I)-DF/FF)


		! Verifica se existe planicie de inundação:
		IF(NPF(I)>1.and.sum(HF(:,I))/=0.0)THEN !RP
			! Verifica se nivel dagua está acima de minimo da planicie de inundação:
			IF(HF(1,I)<=H1)THEN
				! Interpola largura da planicie de inundação:
				AFI(I)=FINT(HF(1:NPF(I),I),FAF(1:NPF(I),I),NPF(I),H1)

				
!				! Tratamento da planicie de inundacao:
				AFI(I)=min(AFI(I),FAF(NPF(I),I))
!				if (H1<ZO(I)+HA(1,I)) then
!					AFI(I)=AFI(I)*(H1-ZO(I))/HA(1,I)
!				endif
!				write(*,*) I,AFI(I),afimax,NPF(I),FAF(1,I),FAF(NPF(I),I),H1,HF(1,I),HF(NPF(I),I)
!				read(*,*)
!				write(*,*)  I,AFI(I),afimax,dafimax,(FINT(HF(1:NPF(I),I),FAF(1:NPF(I),I),NPF(I),H1+0.01)-AFI(I))/0.01
!				read(*,*)
				
				if (AFI(I)>afimax) then
					afimax=AFI(I)
					iafimax=I
				endif
				
				if (dafimax<(FINT(HF(1:NPF(I),I),FAF(1:NPF(I),I),NPF(I),H1+0.01)-AFI(I))/0.01) then
					dafimax=(FINT(HF(1:NPF(I),I),FAF(1:NPF(I),I),NPF(I),H1+0.01)-AFI(I))/0.01
					idafimax=I
				endif
				
				if (AFI(I)/=0.0.and.dafimax2<(FINT(HF(1:NPF(I),I),FAF(1:NPF(I),I),NPF(I),H1+0.01)-AFI(I))/(0.01*AFI(I))) then
					dafimax2=(FINT(HF(1:NPF(I),I),FAF(1:NPF(I),I),NPF(I),H1+0.01)-AFI(I))/(0.01*AFI(I))
					idafimax2=I
				endif
				
				
				
			ENDIF
		ENDIF
!		AFI(i)=0.0 !Teste
		if (AFI(I)<0.0) write(*,*) 'Largura planicie negativa'

		! Se prof. menor que minima, zerar largura da planicie de inundação: 
!		if (HO(I)<=1.1*Hmin) AFI(I)=0.0 


		Fr=abs(QO(i)/A(i))/(9.81*HO(i))**2
		Fr2=Fr**2
		BETA(i)=max(1-Fr2,0.0) ! Suprime termos de inercia com aumento do numero de Froude


		if (Fr>0.5.or.abs(QO(i)/A(i))>10.0) then
!			write(*,*) 'Secao',i,'Froud=',Fr, 'Prof=',HO(i),'V=',QO(i)/A(i),'Q=',QO(i),'B=',T(i),'AFI=',AFI(i)
!@ DCB ago/2012			write(7000,*) 'Secao ',i,'Fr=',Fr,'H=',HO(i),'V=',QO(i)/A(i),'Q=',QO(i),'b=',T(i),'Bplan=',AFI(i)
!			read(*,*)
		
		endif


		! Verificacoes:
		if (T(I)<=0.0.or.A(I)<=0.0.or.R(I)<=0.0.or.isNaN(T(I)).or.isNaN(A(I)).or.isNaN(R(I))) write(*,*) 'Erro secao',I,'A,T,R',A(I),T(I),R(I)
		if (CK(I)<=0.0.or.CKY(I)<=0.0.or.isNaN(CK(I)).or.isNaN(SFF(I)).or.isNaN(CKY(I))) write(*,*) 'Erro secao',I,'CK,SFF,CKY',CK(I),SFF(I),CKY(I)


	ENDDO
	
!	do j=1,1
!		AFIaux=AFI
!		do iTr=1,nTr
!				
!			do i=TrNST(iTr,1),TrNST(iTr,2)
!				if (i==TrNST(iTr,1)) then
!					AFI(i)=(AFIaux(i+1)+AFIaux(i))/2.0
!				elseif (i==TrNST(iTr,2)) then
!					AFI(i)=(AFIaux(i)+AFIaux(i-1))/2.0
!				else
!					AFI(i)=(AFIaux(i+1)+AFIaux(i)+AFIaux(i-1))/3.0
!				endif
!			enddo
!		enddo
!	enddo





	! Calcula coeficientes das equacoes das condicoe de contorno e armazena no vetor AA
	DO J=1,NBOUN
		! Identifica linha da eq. da Condicao de Contorno J
		K=JLIN(J)

		! Verifica tipo de CC:
		IF(bounFLAG(J)<=4)THEN
			! Série temporal de nivel d'água, prof. ou vazão
			! Coef. independente BB:	
			BB(K)=HQB(iThd,J)
			! Localização no vetor AA:
			K=ICOL(K,2)

			! Coef. dependente AA:
			AA(K)=1.0
		ELSE
			! Curva de descarga:
			! Localização dos coeficientes no vetor AA:
			L1=ICOL(K,2)
			L2=ICOL(K,3)
			! Seção da condição de contorno:
			NS=NBO(J)
			! Coef. independente AA relacionado a QO:
			AA(L2)=1.0
			! Prof. na secao:
			H1=HO(NS)
			! Se dados em nivel d'água, corrige H1:
			IF(LD/=0) H1=H1+ZO(NS)
			
			! Verifica tipo de curva de descarga: 
			IF(bounFLAG(J)==5)THEN
				! Tabela:
				JJ=bounID(J)
				! Vazão em funcao de H1:
				DF=FINT(HT(1:NPX(JJ),JJ),QT(1:NPX(JJ),JJ),NPX(JJ),H1)
				! Derivada do vazão em relação a prof. dQ/dy:
				HAUX=H1+.01
				DF1=FINT(HT(1:NPX(JJ),JJ),QT(1:NPX(JJ),JJ),NPX(JJ),HAUX)
				DF1=(DF1-DF)/.01
				! Coef. independente AA relacionado a HO:
				AA(L1)=-DF1
				! Coef. independente BB:
				BB(K)=DF-DF1*HO(NS)
			ELSE
				! Eq. de Manning:
				! Numero de pontos na tabela das secoes:
				L=NP(NS)
				
				! Testa tipo de seção:
				IF (L>0) THEN
					! Dados detalhados de seção transversal
					! Prof. auxiliar:
					H1=H1+.01
					! Raio hidraulico para H1+0.01:
					R2=FINT(HA(1:L,NS),RR(1:L,NS),L,H1)
				ELSE
					! Seção tipo retangular:
					! Profundidade auxiliar:
					H1 =HO(NS)+.01
					! Raio hidraulico para H1+0.01:
					R2=T(NS)*H1/(T(NS)+2*min(HA(1,NS),H1))
				ENDIF
	
				! Derivada do raio hidráulico em relação a prof. dR/dy:
				DR=(R2-R(NS))/.01
				! Derivada em relacao a condutancia hidraulica: (Obs.: nao considera coef. Manning variável)
				COC=2/3.*DR/R(NS)+T(NS)/A(NS)
				! Coef. independente AA relacionado a HO:
				AA(L1)=-QO(NS)*COC
				! Coef. independente BB:
				BB(K)=QO(NS)*(1.-HO(NS)*COC)
			ENDIF
		ENDIF
	ENDDO

	!!$OMP END DO
	!!$OMP END PARALLEL

	J1=NBOUN+1

	!!$OMP PARALLEL NUM_THREADS(4)
	!!$OMP DO PRIVATE(J2,M,J,CS13,CS1,K,AXY,L2,L3,L4,L5,V1,V2,CS2,CS3,CS4,CS5,CS6,CS7,CS8,CS9,CS10,CS11,CS12,T1,T2,T3,T4,T5,T6,H1,H2,JP,DF,HAUX,DF1,IPB,DFD,DF2)
	! Obs.: Se existir mais de uma barragem o programa em paralelo não funciona.

	! Calcula coeficientes das equações dos trechos e armazena no vetor AA:
	DO I=1,NREAC
		! Linha da eq. continuidade na ordem original:
		J2=J1+(I-1)*2 +1
		! Seção de montante:
		M=NST(I,1)
		! Seção de jusante:
		J=NST(I,2)
		! Linha da eq. da continuidade no sistema ordenado:
		K=JLIN(J2)
		! Localizacao dos coeficientes dependentes da eq. continuidade no vetor AA
		L2=ICOL(K,2)
		L3=ICOL(K,3)
		L4=ICOL(K,4)
		L5=ICOL(K,5)

		! Verifica tipo de trecho:
		SELECT CASE (reacFLAG(I))
		CASE(1)
			! Equações de Saint Venant (continuidade e dinâmica)
			
			! Eq. Continuidade:

			! Aux.:
			CS13=HO(M)+HO(J)
			!CS1=T(M)+T(J)+2*AFI(M) ! So entra largura da planicie de montante? !RP
			!CS1=T(M)+T(J)+AFI(M)+AFI(J)
			CS1=T(M)+T(J)+2*AFI(J) ! Usa largura de jusante

			! Desconsidera armazenamento na planicie de inundacao em profundidade muito baixa:
!			if (HO(M)<Hmin*1.05.or.HO(J)<Hmin*1.05) CS1=T(M)+T(J)
			
			! Derivada da largura em relacao a distancia dT/dx
			AXY=(T(J)-T(M))/DXT(I)
			
			! Coef. dependentes AA:
			AA(L3)=-4.*THETA(I)
			AA(L5)=-AA(L3)
			AA(L2)=CS1
			AA(L4)=CS1

			if (isNaN(AA(L3)).or.isNaN(AA(L4)).or.isNaN(AA(L2)).or.isNaN(AA(L5))) then
				write(*,*) 'Eq Continuidade: digite ENTER',I,AA(L3),AA(L5),AA(L2),AA(L4),theta(I),T(M),T(J),AFI(M),AFI(J)
				!@DCBteste  read(*,*)
			endif

			! Coef. independente BB:
!			BB(K)=DThd*(QL2(I,1)+QL2(I,2))*2+CS1*CS13
			BB(K)=DThd*(QL2(I,1)+QL2(I,2))*2/DXT(I)+CS1*CS13	!RP		
			! Eq. Dinâmica:

			! Velocidades médias no tempo anterior:
			V1=QO(M)/A(M)
			V2=QO(J)/A(J)
			
			! Auxiliares para calculo dos coeficientes:
			CS1=THETA(I)*(V1+V2)
			CS2=THETA(I)*(BETA(M)*V1+BETA(J)*V2)
			CS3=THETA(I)*(BETA(M)*V1*V1*T(M)+BETA(J)*V2*V2*T(J))
			CS4=THETA(I)*G*(A(M)+A(J))
			CS5=2.*CSO*A(M)*SFF(M)/CK(M)
			CS6=2.*CSO*A(J)*SFF(J)/CK(J)
			CS7=0.0
			CS8=THETA(I)*G*(A(M)+A(J))
			CS9=SFF(M)/V1*CSO
			CS10=SFF(J)/V2*CSO
			CS11=SFF(M)*T(M)*CSO
			CS12=SFF(J)*T(J)*CSO
			T1=1.+2.*CS9
			! Linha da eq. dinamica na ordem original:
			J2=J2-1
			! Linha da eq. dinamica no sistema ordenado:
			K=JLIN(J2)
			! Localizacao dos coeficientes dependentes da eq. dinamica no vetor AA
			L2=ICOL(K,2)
			L3=ICOL(K,3)
			L4=ICOL(K,4)
			L5=ICOL(K,5)

			! Coef. dependentes AA:
			AA(L3)=T1-CS2-CS1*BETA(M)
			T1=CS4+CS5*CKY(M)-CS11
			AA(L2)=CS3-T1+CS7
			T1=1.+2.*CS10
			AA(L5)=T1+CS2+CS1*BETA(J)

!			if (isNaN(AA(L5))) write(*,*) T1,CS2,CS1,BETA(J),BETA(M),CS10,SFF(J),SFF(M),V2,V1,QO(M),QO(J),A(M),A(J)

			T1=CS4-CS6*CKY(J)+CS12
			AA(L4 )=-CS3+T1+CS7
			T1=QO(M)*(1.+CS9)
			T2=QO(J)*(1.+CS10)
			T3=HO(M)*(CS11-CS5*CKY(M))
			T4=HO(J)*(CS12-CS6*CKY(J))
			T5=CS7*CS13
			T6=DThd*(BETA(M)*V1*V1+BETA(J)*V2*V2)*AXY
			! Coef. independente BB:
			BB(K)=T1+T2+T3+T4+T5+T6+CS8*(ZO(M)-ZO(J))

			if (isNaN(AA(L3)).or.isNaN(AA(L4)).or.isNaN(AA(L2)).or.isNaN(AA(L5))) then
				write(*,*) 'Eq Dinamica: digite ENTER',I,AA(L3),AA(L5),AA(L2),AA(L4)
				!@DCBteste  read(*,*)
			endif


			! Condição de contorno interna caso profundidade d'água menor que a minima:
			if (HminFLAG(M)==1) then
!			if (1==2) then
!				HO(M)=Hmin
				nHmin=nHmin+1
!				write(*,*) 'Produndidade minima',M,J
				
				! Modificação na ordem das equações para evitar termos nulos na diagonal principal

				! Linha da eq. continuidade na ordem original:
				J2=J2
				! Linha da eq. da continuidade no sistema ordenado:
				K=JLIN(J2)

!				write(*,*) 'linha eq:',k,'coluna coef.:',(M-1)*2+1,(M-1)*2+2,(J-1)*2+1,(J-1)*2+2	
!				read(*,*)		

				! Localizacao dos coeficientes dependentes da eq. continuidade no vetor AA
				L2=ICOL(K,2)
				L3=ICOL(K,3)
				L4=ICOL(K,4)
				L5=ICOL(K,5)
				! Eq. Continuidade QO(M)-QO(J)=-QL2(I,2)
				! Aux.:
				CS13=HO(M)+HO(J)
				CS1=T(M)+T(J)
!				CS1=T(M)+T(J)+AFI(M)+AFI(J)

				! Derivada da largura em relacao a distancia dT/dx
				AXY=(T(J)-T(M))/DXT(I)
			
				! Coef. dependentes AA:
				AA(L3)=-4.*THETA(I)
				AA(L5)=-AA(L3)
				AA(L2)=CS1
				AA(L4)=CS1
				! Coef. independente BB:
				BB(K)=DThd*(QL2(I,1)+QL2(I,2))*2/DXT(I)+CS1*CS13	!RP		

!				AA(L3)=1.0
!				AA(L5)=-1.0
!				AA(L2)=0.0
!				AA(L4)=0.0

				! Coef. independente BB:
!				BB(K)=-QL2(I,2)

		
				! Linha da eq. dinamica na ordem original:
				J2=J2+1
				! Linha da eq. dinamica no sistema ordenado:
				K=JLIN(J2)

!				write(*,*) 'linha eq:',k,'coluna coef.:',(M-1)*2+1,(M-1)*2+2,(J-1)*2+1,(J-1)*2+2	
!				read(*,*)		

				! Localizacao dos coeficientes dependentes da eq. dinamica no vetor AA
				L2=ICOL(K,2)
				L3=ICOL(K,3)
				L4=ICOL(K,4)
				L5=ICOL(K,5)

				! Coef. dependentes AA:
				AA(L3)=0.0000009
				AA(L5)=0.0000011
				AA(L2)=1.0
				AA(L4)=0.0000001 !/HO(J) !-1.0 !0.0
				! Coef. independente BB:
				BB(K)=Hmin+0.0000001*HO(J)+0.0000009*QO(M)+0.0000011*QO(J) !min(ZO(J)-ZO(M),0.0)  

				! Usar curva de descarga considerando equação de manning na 
				!seção de montante e equação de manning ponderada na seção de jusante:
				! Coef. dependentes AA:
!				DF1=QO(M)
!				DFD=CKY(I)*SFF(M)**0.5
!				DF2=DFD*-0.001
!				AA(L5)=0.0
!				AA(L3)=1.0
!				AA(L2)=-DFD
!				AA(L4)=-DF2
				! Coef. independente BB:
!				BB(K)=DF1-DFD*HO(M)-DF2*HO(J)




			endif

		CASE(2)
			! Trecho curto:
			! RP talvez teremos problema neste tipo de abordagem a depender da organizacao das equacoes na rotina matrix
			! RP pode haver elemento nulo na diagonal principal
			! RP Muito dificil de acontecer (somente quando trecho possui condicao de contorno em uma das secoes)
			! Eq. Continuidade:
			! Coef. dependentes AA:
			
			! A primeira equação é a da continuidade:
			J2=J2-1 ! Modificação para evitar termos nulos na diagonal principal:
		
			! Linha da eq. da continuidade no sistema ordenado:
			K=JLIN(J2)
			! Localizacao dos coeficientes dependentes da eq. continuidade no vetor AA
			L2=ICOL(K,2)
			L3=ICOL(K,3)
			L4=ICOL(K,4)
			L5=ICOL(K,5)
			
			AA(L3)=1.0
			AA(L5)=-1.0
			AA(L2)=0.0
			AA(L4)=0.0

			! Coef. independente BB:
			BB(K)=-QL2(I,2)

			! Eq. Energia simplificada (esc. permanente e s/ termos de energia cinética e perda de carga):
		
			! Linha da eq. dinamica na ordem original:
			!J2=J2-1
			
			! Segunda equacão é a dinamica
			J2=J2+1 ! Modificação para evitar termos nulos na diagonal principal:
			
			! Linha da eq. dinamica no sistema ordenado:
			K=JLIN(J2)
			! Localizacao dos coeficientes dependentes da eq. dinamica no vetor AA
			L2=ICOL(K,2)
			L3=ICOL(K,3)
			L4=ICOL(K,4)
			L5=ICOL(K,5)
			

			! Coef. dependentes AA:
			AA(L3)=0.0
			AA(L5)=0.0
			AA(L2)=1.0
			AA(L4)=-1.0
			! Coef. independente BB:
			BB(K)=(ZO(J)-ZO(M))


			! Correcao caso profundidade minima:
			!if (HminFLAG(M)==1) then
			!	AA(L3)=0.0
			!	AA(L5)=0.0
			!	AA(L2)=1.0
			!	AA(L4)=0.0000001 !/HO(J) !-1.0 !0.0
				! Coef. independente BB:
			!	BB(K)=Hmin+0.0000001*HO(J)
			!endif		


		CASE(3)
		! Trecho com barragem ou condicao de contorno interna
		! Equação da continuidade + curva de descarga da barragem
			
			! Eq. Continuidade:
			! Coef. dependentes AA:
			AA(L2)=0.0
			AA(L3)=1.0
			AA(L4)=0.0
			AA(L5)=-1.0
			BB(K)=0.0

			! Eq. Dinâmica:
			! Linha da eq. dinamica na ordem original:		
			J2=J2-1
			! Linha da eq. dinamica no sistema ordenado: 
			K=JLIN(J2)
			
			! Indice da barragem
			LT=0
			DO
				! Percorre barragens
				LT=LT+1
				! Testa se indice trecho é igual a indice trecho da barragem:
				IF (barrReac(LT)==I) EXIT
			ENDDO
			

			! Localizacao dos coeficientes dependentes da eq. dinamica no vetor AA
			L2=ICOL(K,2)
			L3=ICOL(K,3)
			L4=ICOL(K,4)
			L5=ICOL(K,5)
			! Prof. dagua
			H1=HO(M)
			H2=HO(J)
			
			! Se dados em nivel d'água, corrige H1 e H2: 
			IF(LD/=0)THEN
				H1=H1+ZO(M)
				H2=H2+ZO(J)
			ENDIF
			! Verifica curva de descarga utilizada dependento do intervalo de tempo:
			DO JP=1,NIT(LT)
				IF (JIT(LT,JP)>=iThd)EXIT
			ENDDO
			

			! RP Verificar coeficientes: pode dar problema - zeros em elementos da diagonal principal

			! Verifica se existe efeito de jusante:
			IF(NPB(LT,JP)>=0)THEN
				! S/ efeito de jusante:

				! Interpola vazão em funcao de H1:			
				DF=FINT(HBR(1:NPB(LT,JP),LT,JP),QBR(1:NPB(LT,JP),1,LT,JP),NPB(LT,JP),H1)
				! Prof. auxiliar:
				HAUX= H1+.01
				! Interpola vazão em funcao de HAUX:
				DF1=FINT(HBR(1:NPB(LT,JP),LT,JP),QBR(1:NPB(LT,JP),1,LT,JP),NPB(LT,JP),HAUX)
				! Derivada da vazão em função de H1:
				DF1=(DF1-DF)/.01
				! Coef. dependentes AA:
				AA(L2)=-DF1
				AA(L3)=1.0
				AA(L4)=0.0
				AA(L5)=0.0
				! Coef. independente BB:
				BB(K)=DF-DF1*HO(M)
			ELSE
				! C/ efeito de jusante:
				! Numero de pontos da tabela de montante:
				IPB=IABS(NPB(LT,JP))
				! Interpola vazão em funcao de H1 e H2:
				DF1=FINT2(HBR(1:IPB,LT,JP),H2BR(1:NPB2(LT,JP),LT,JP),QBR(1:IPB,1:NPB2(LT,JP),JP,LT),IPB,NPB2(LT,JP),H1,H2)
				
				! Interpola vazão em funcao de H1+0.01 e H2:
				! Calcula derivada da vazão em relacao a prof. de montante dQ/dH1:
				HAUX= H1+0.01
				DFD=(FINT2(HBR(1:IPB,LT,JP),H2BR(1:NPB2(LT,JP),LT,JP),QBR(1:IPB,1:NPB2(LT,JP),JP,LT),IPB,NPB2(LT,JP),HAUX,H2)-DF1)/0.01
				! Interpola vazão em funcao de H1 e H2+0.01:
				! Calcula derivada da vazão em relacao a prof. de jusante dQ/dH2:
				HAUX= H2+0.01
				DF2=(FINT2(HBR(1:IPB,LT,JP),H2BR(1:NPB2(LT,JP),LT,JP),QBR(1:IPB,1:NPB2(LT,JP),JP,LT),IPB,NPB2(LT,JP),H1,HAUX)-DF1)/0.01
				! Coef. dependentes AA:
				AA(L5)=1.0
				AA(L3)=0.0
				AA(L2)=-DFD
				AA(L4)=-DF2
				! Coef. independente BB:
				BB(K)=DF1-DFD*HO(M)-DF2*HO(J)
			ENDIF
		ENDSELECT
		
	ENDDO

	!!$OMP END DO
	!!$OMP END PARALLEL

	! Calcula coeficientes das equações das confluencias e armazena no vetor AA:
	IF(NCONF/=0)THEN
    
		J1=NBOUN+NREAC*2+1
			
		!!$OMP PARALLEL NUM_THREADS(4)
		!!$OMP DO PRIVATE(J2,K,KM,JA,K1,K2,K3,ITT)
		DO I=1,NCONF
			! Eq. Continuidade:

			! Linha da eq. continuidade na ordem original:
			J2=J1+(I-1)*3
			! Linha da eq. da continuidade no sistema ordenado:
			K=JLIN(J2)
			! Localizacao dos coeficientes dependentes da eq. continuidade no vetor AA
			! E coef. dependentes:
			KM=ICOL(K,2)
			AA(KM)=1.
			KM=ICOL(K,3)
			AA(KM)=1.0
			KM=ICOL(K,4)
			AA(KM)=-1.0
			JA=(I-1)*2+1
			! Indices das secoes da  confluencia:		
			K1=NCC(I,1)
			K2=NCC(I,2)
			K3=IABS(NCC(I,3))
			ITT=1 ! Verificar esse lance de IT

			! Verifica tipo de confluencia:
			IF(NCC(I,3)>=0)THEN
				! Confluencia divergente:
				! Coef. dependente da eq. continuidade: (Obs.: considera contribuição lateral)
				!BB(K)=(DXC(JA)+DXC(JA+1))*(QL2(K3,1)+QL2(K3,2))/2			! Tirar contribuicao lateral daqui RP
				BB(K)=0.0
				! Coeficientes das equações de energia:			
				J2=J2+1
				CALL COEF2(K3,K1,DXC(JA),ALFA(JA),J2)
				J2=J2+1
				CALL COEF2(K3,K2,DXC(JA+1),ALFA(JA+1),J2)
			ELSE
				! Confluencia convergente:
				! Coef. dependente da eq. continuidade: (Obs.: considera contribuição lateral)
				!BB(K)=-DXC(JA)*(QL2(K1,1)+QL2(K1,2))/2-DXC(JA+1)*(QL2(K2,1)+QL2(K2,2))/2	! Tirar contribuicao lateral daqui RP
				BB(K)=0.0
				! Coeficientes das equações de energia:
				J2=J2+1
				CALL COEF2(K1,K3,DXC(JA),ALFA(JA),J2)
				ITT=0
				J2=J2+1
				CALL COEF2(K2,K3,DXC(JA+1),ALFA(JA+1),J2)
			ENDIF
		ENDDO
		!!$OMP END DO
		!!$OMP END PARALLEL
	ENDIF






	! *************
	! Verificacoes
!	do j=1,num
!		if (isNaN(XI(j))) then
!			write(*,*) 'NaN nas vars de estado',j
!			write(11222,*) 'NaN nas vars de estado',j

!		endif
!	enddo
	
	do j=1,ndimAA
		if (isNaN(AA(j))) then
			
!@ DCB ago/2012			write(7000,*) 'NaN na matriz de coeficientes AA',j
!			do i=1,num
!				do k=2,ICOL(i,1)
!					if (ICOL(i,k)==j) then
!						write(*,*) i,ICOL(i,1),k-1
!						! Procura equacao:
!						do L=1,num
!							if (JLIN(L)==i) then
!								if (L-nboun-2*nreac>0) then
!									write(*,*) 'Equacao das confluencias:',L-nboun-2*nreac
!								elseif (L-nboun>0) then
!									write(*,*) 'Equacao dos subtrechos:',L-nboun
!								else
!									write(*,*) 'Cond.Contorno:',L
!								endif
!							endif
!						enddo
!					endif
!				enddo
!			enddo
!			read(*,*)
		endif
		
	enddo
	
	do j=1,num
		if (isNaN(BB(j))) then
!			write(*,*) 'NaN na matriz de coeficientes BB',j
!@ DCB ago/2012			write(7000,*) 'NaN na matriz de coeficientes BB',j
!			do L=1,num
!				if (JLIN(L)==j) then
!					if (L-nboun-2*nreac>0) then
!						write(*,*) 'Equacao das confluencias:',L-nboun-2*nreac
!					elseif (L-nboun>0) then
!						write(*,*) 'Equacao dos subtrechos:',L-nboun
!					else
!						write(*,*) 'Cond.Contorno:',L
!					endif
!				endif
!
!			enddo
		endif
	enddo
	
	
	!****************	



	RETURN
	END
