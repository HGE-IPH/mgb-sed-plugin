	SUBROUTINE Hidrodinamico2
	!********************************************************************************************
	!										    IPH4 
	!										   ------
	! 
	!						Modelo Hidrodinamico IPH4 Versão Fortran 90
	! 
	!							INSTITUTO DE PESQUISAS HIDRAULICAS
	!					da Universidade Federal do Rio Grande do Sul- UFRGS 
	!
	!  
	!  - Versão Original em Fortran 77 (mar/1996):
	!
	!		Carlos E.M. Tucci
	!
	!  - Versão Atual em Fortran 90:
	!
	!		Rodrigo Paiva e Juan Martin Bravo
	!  
	!  - Detalhes sobre aspectos teóricos:
	!       
	!    Tucci, C.E.M. 1978. Hydraulic and Water Quality Model for a River Network. 
	!			PhD dissertation, Colorado State University, Fort Collins, USA.
	!
	!    Tucci, C.E.M. 2005. Modelos Hidrológicos. Ed. UFRGS/ABRH, Porto Alegre.
	!
	!********************************************************************************************
	!-----------------------------------------------------------------------
	! Descrição das variáveis locais:
	!
	! I,J,L = variáveis auxiliares
	! TIME1,TIME2 = auxiliares para calculo do tempo de processamento
	!----------------------------------------------------------------------
	! Declaração de variáveis:
	USE PORTLIB !biblioteca para calcular o tempo de processamento
	USE PAR1_MOD
	USE MAT_MOD
	USE TIME_MOD
	USE BARR_MOD
	USE PP_MOD
	USE AUX_MOD
	USE VARS_MAIN
	IMPLICIT NONE

	! Variáveis locais de cálculo:
	integer I,J,L
	integer IV
	integer iWrite,iThd1,iThd2
	real aux,aux2
	real courant
	real velcrit
	integer ijus,imon,ijus2,imon2,i1,i2,i3,iX2
	real H1,FINT
	real Thd
	real DHreal,DQreal
	logical NaNflag
	integer erroFLAG,nItera
	real dtHD0_1
	real,allocatable:: HO2(:),QO2(:),QL2_1(:,:)
	!----------------------------------------------------------------------------


	! Alteração para recomeçar o loop quando o intervalo de tempo é muito pequeno.
	!********************************************************
	allocate(HO2(nX),QO2(nX),QL2_1(NREAC,2))
	HO2=HO
	QO2=QO
	QL2_1=QL2
	dtHD0_1=dtHD0
	!****************


	aux=DT(NC)
	Thd=dtHD
	J=1

	G=9.81
	CONST=1.


	
	HO1=HO
	QO1=QO

	! Intervalos de tempo inicial e final na condicao de contorno lida de arquivo:
	iThd1=(nThd-1)*(iT-1)+1
	iThd2=iThd1+nThd-1


	! Valores padrao de máxima variação dos niveis d'água por intervalo de tempo e maximo criterio da velocidade:
	if (INIC==0) then
		DHmax=0.50
		VcritMAX=2.5
		DQmax=25000.0
		Hmin=0.20
	endif


	iThd=1
	nItera=0

!	DO iThd=2,nTHD
	do while (Thd<=DTP)
		! Zera vetor dos coeficientes:
		AA=0.0
		BB=0.0

		! Inicializa HO e QO:
!		HO1=max(HO1,Hmin)
		HO=HO1
		QO=QO1

		CSO=G*DThd

		! Calcula coef. theta:
		DO I=1,NREAC
			IF(DXT(I)/=0)THEN
				THETA(I)=DThd/DXT(I)
			ENDIF
		ENDDO


		nItera=nItera+1	
		write(*,*) nItera, nhmin,100*(Thd-DThd)/DTP,'%',DThd,minval(HO)
!@ DCB ago/2012		write(7000,*) nItera, nhmin,100*(Thd-DThd)/DTP,'%',DThd


		! Critério para selecao do dt:
		! Máxima variação na profundidade:

		! Critério da velocidade

		! Condições de contorno:

		! Intervalos de tempo inicial e final:
		iThd1=(nThd-1)*(iT-1)+1
		iThd2=iThd1+nThd-1

		DO I=1,NBOUN
			SELECTCASE (bounFLAG(I))
			CASE(1)
				! Resultado do MGB de série de vazao em bacia de cabeceira:
				! Secao:
				iX=NBO(I)
				! Trecho:
				iTr=TrXSec(iX)
				! Subbacia:
				iC=TrCELL(iTr)
				! Interpola contribuição da subbacia:
				HQB(1,I)=QCEL1(IC)+(QCEL2(IC)-QCEL1(IC))*Thd/DTP
				i1=TrNST(iTr,1)
				i2=TrNST(iTr,2)
				if ((HO(i1)+ZO(i1)-HO(i2)+ZO(i2))/TrL(iTr)>0.05) HQB(1,I)=max(HQB(1,I),0.0) ! Se declividade da linha dágua é muito alta não considera vazão negativa					


			CASE(2)
				! Resultado do MGB de série de vazão (confluencia entre bacias s/ hidrodinamico)
				! Secao:
				iX=NBO(I)
				! Trecho:
				iTr=TrXSec(iX)
				! Subbacia:
				iC=TrCELL(iTr)
				! Interpola dados para intervalo de tempo do hidrodinâmico:
				HQB(1,I)=FINT(TP1,QCONTORM(IC,:),nThd,Thd)
				HQB(1,I)=FINT(TP1,QCONTORM(IC,:),25,Thd) ! Correção para considerar que dt musk = 1 h ! RP fev 2011.
			CASE(3,4)
				! Séries temporais lidas em arquivo externo:
				J=bounID(I)
				HQB(1,I)=FINT(TP2(iThd1:iThd2),bounTS(iThd1:iThd2,J),iThd2-iThd1+1,Thd+(iT-1)*DTP) ! RP Verificar bem esta parte.

			ENDSELECT
			
			! Condicao de contorno no aquecimento nas condicoes iniciais.
			if (INIC==1.and.bounFLAG(I)==4) then
				iX=NBO(I)
				HQB(1,I)=Z01+(Z02-Z01)*Thd/DTP-ZO(iX)
				HQB(1,I)=max(bounTS(1,J),HQB(1,I))
			endif
		ENDDO


		! Contribuição lateral:

		DO I=1,NREAC
			! Testa se é contribuição lateral concentrada ou distribuida:
			IF (TrXSec(NST(I,2))/=TrXSec(NST(I,1))) THEN
				! Secao:
				iX=NST(I,2)
				! Trecho:
				iTr=TrXSec(iX)
				! Subbacia:
				iC=TrCELL(iTr)		
	
				QL2(I,2)=FINT(TP1,QCONTORM(IC,:),nthd,Thd)
				QL2(I,2)=FINT(TP1,QCONTORM(IC,:),25,Thd) ! Correção para considerar que dt musk = 1 h ! RP fev 2011.
			ELSE
				! Secao:
				iX=NST(I,1)
				! Trecho:
				iTr=TrXSec(iX)
				! Subbacia:
				iC=TrCELL(iTr)		
				
				IF (OD(iC)==1) CYCLE ! Cuidar para não repetir vazao em cond.contorno e qlat em bacias de cabeceira
				! IMPORTANTE********************************************************************************************

				! Interpola e pondera pela relacao entre dx subtrecho e comprimento total do trecho:

				QL2(I,2)=((QCEL1(IC)+(QCEL2(IC)-QCEL1(IC))*Thd/DTP))*DXT(I)/(TrL(iTr)*1000)
			ENDIF
		ENDDO

		

		!Calculo dos coeficientes das matrizes A e B:
		CALL COEF1
		
		!Solucao do sistema de equações A*XI=B:
!		call GradBiConj
!		CALL SKYLIN
		call skyotm

		! Armazena variáveis de estado XI em HO e QO:
		DO I=1,NX
			! Identifica posição da seção no vetor XI
			IV=(I-1)*2+1
			HO(I)=XI(IV)
			QO(I)=XI(IV+1)

		ENDDO

		
		! Máxima variação de nível d'água:
!		DHreal=maxval(abs(HO-HO1))
		i2=-1
		DHreal=-1.0
		do iX=1,nX
			if (abs(HO(iX)-HO1(iX))>DHreal) then
				! Se passa para prof. abaixo da minima e variacao nao é muito brusca, nao considera:
				if (HO(iX)<Hmin.and.abs(HO(iX)-HO1(iX))<2*DHmax) cycle
				i2=iX
				DHreal=abs(HO(iX)-HO1(iX))
				
			endif
		enddo	
!		i=maxloc(abs(HO-HO1))
		! Máxima variação de vazão:
!		DQreal=maxval(abs(QO-QO1))	
		DQreal=-1.0
		do iX=1,nX
			if (abs(QO(iX)-QO1(iX))>DQreal) then
				i1=iX
				DQreal=abs(QO(iX)-QO1(iX))
			endif
		enddo	
		! Velocidade:
		velcrit=-1000.0
		do i=1,nReac
			if (DXT(i)==0.0) cycle ! Pula trechos curtos
			ijus=NST(i,2)
			imon=NST(i,1)
			aux=max(abs(QO(ijus)/A(ijus)),abs(QO(imon)/A(imon)))*dtHD/DXT(i)
			if (velcrit<aux) then
				ijus2=ijus
				imon2=imon
				velcrit=aux
			endif
		enddo

		! Verifica se existe NaN na solução:
		NaNflag=.false.
		do iX=1,nX
			NaNflag=NaNflag.or.isNaN(QO(iX)).or.isNaN(HO(iX))
		enddo
		

if (1==1) then		
		! Tratamento a depender do tipo de erro:	
		if (NaNflag) then
			! NaN nas variaveis:
			write(*,*) 'Erro. NaN nas vars de estado do modelo hidrodinamico'
			write(*,*) 'Diminui intervalo de tempo'
!@ DCB ago/2012			write(7000,*) 'Erro. NaN nas vars de estado do modelo hidrodinamico'

			Thd=Thd-dtHD
			dtHD=dtHD*0.1
			Thd=Thd+dtHD
			erroFLAG=1
		elseif (DQreal>DQmax) then
			! Variação brusca na vazao:
			! Melhorar criterio da vazao:	
			write(*,*) 'Alta variacao na vazao na secao',i1
			write(*,*) 'dQmax=',DQreal

!@ DCB ago/2012			write(7000,*) 'Alta variacao na vazao na secao',i1
!@ DCB ago/2012			write(7000,*) 'dQmax=',DQreal
			
			write(*,*) 'Diminui intervalo de tempo'
			
			Thd=Thd-dtHD
			dtHD=dtHD*0.1
			Thd=Thd+dtHD
			erroFLAG=2
				
		elseif (DHreal>DHmax) then
			! Variação brusca no nível d'água:

			write(*,*) 'Alta variacao nos niveis dagua na secao',i2
			write(*,*) 'dHmax=',DHreal,'Q=',QO(i2)

!@ DCB ago/2012			write(7000,*) 'Alta variacao nos niveis dagua na secao',i2
!@ DCB ago/2012			write(7000,*) 'dHmax=',DHreal,'Q=',QO(i2),'H=',HO(i2),'H1=',HO1(i2)
			
			write(*,*) 'Diminui intervalo de tempo'
			Thd=Thd-dtHD
			dtHD=0.8*dtHD*DHmax/DHreal ! Escolhe dt para que DHreal=DHmax
			Thd=Thd+dtHD
			erroFLAG=3	
		elseif (velcrit>VcritMAX) then
			! Velocidade muito alta:
			write(*,*) 'Velocidade muito alta'
			write(*,*) 'v*dt/dx=',velcrit,'Secoes=',imon2,ijus2

!@ DCB ago/2012			write(7000,*) 'Velocidade muito alta'
!@ DCB ago/2012			write(7000,*) 'v*dt/dx=',velcrit,'Secoes=',imon2,ijus2
			
			write(*,*) 'Diminui intervalo de tempo'
			Thd=Thd-dtHD
			dtHD=0.8*dtHD*VcritMAX/velcrit ! Escolhe dt para que v*dt/dx=VcritMAX
			Thd=Thd+dtHD
			erroFLAG=4	
		else

			! Não houveram erros:

			

			! Verifica profundidade minima:
			do i=1,nReac
				iX=NST(i,1)
				iX2=NST(i,2)
				if (HO(iX)<0.999*Hmin) then
					if (HO(iX)+ZO(iX)<HO(iX2)+ZO(iX2)-0.001) cycle ! Nao entra na condicao de contorno se profundidade da secao de jusante é maior:
					!TrXsec(iX)==TrXsec(iX+1).and.	
					HminFLAG(iX)=1
					HO(iX)=Hmin 
					if (Qmin(iX)==-1.0) Qmin(iX)=QO(iX)
				elseif (HminFLAG(iX)==1.and.(QO(iX)>1.1*Qmin(iX).or.(HO(iX)+ZO(iX))<(HO(iX2)+ZO(iX2))-0.001)) then
					! Somente para condicao de contorno interna se vazao for crescente, maior que x vezes 
					! a vazao quando comecou CC interna ou nivel do fundo de jusante for maior.
					HminFLAG(iX)=0
					Qmin(iX)=-1.0
				endif
			enddo

			! Corrige niveis em confluencias:
!			do i=1,nconf
!				i1=abs(NCC(I,1))
!				i2=abs(NCC(I,2))
!				i3=abs(NCC(I,3))
				! Conflu, secoes 1, 2 e 3, trechos,sentido do escoamento:
!				if (HminFLAG(i3)==1.or.HminFLAG(i2)==1.or.HminFLAG(i1)==1) then
!					HO(i1)=HO(i3)+ZO(i3)-ZO(i1)
!					HO(i2)=HO(i3)+ZO(i3)-ZO(i2)
!					write(*,*) 'conflu', i1,i2,i3,HminFLAG(i3),HminFLAG(i2),HminFLAG(i1)
!				endif
!			enddo
			
			do iX=1,nX
				if (HO(iX)<Hmin) HO(iX)=Hmin ! Para corrigir niveis quando prof. é negativa.
			enddo	
			HO1=HO
			QO1=QO

			! Passa contribuição lateral tempo IT para IT-1:
			QL2(:,1)=QL2(:,2)
			QL2(:,2)=0.0

			! Verifica se intervalo de tempo igual ao máximo:
			if (dtHD/=dtHD0_1) then
				! Aumenta intervalo de tempo:
				selectcase(erroFLAG)
				case(1,2)				
					aux=dtHD*1.5	! Alterei proporcoes 2.0 antigo 2/2011
				case(3)
					aux=0.7*dtHD*DHmax/DHreal ! Aumenta gradativamente ! Alterei proporcoes 0.8 antigo 2/2011
				case(4)	
					aux=0.7*dtHD*VcritMAX/velcrit ! Aumenta gradativamente  ! Alterei proporcoes 0.8 antigo 2/2011
				endselect				
	
				dtHD=min(aux,dtHD*VcritMAX/velcrit,dtHD*DHmax/DHreal)	
				dtHD=min(dtHD,dtHD0_1)
			endif

			! Tempo:

			! Verifica se chegou ao final do intervalo de tempo	
			if (Thd<DTP*0.999) then
				Thd=Thd+dtHD
				if (Thd>DTP) then
					dtHD=dtHD-(Thd-DTP)
					Thd=DTP
				endif
			else
				Thd=Thd+dtHD
			endif			
		endif
		! Para se dt é muito pequeno
		if (dThd<0.01) then
			HO1=HO2
			QO1=QO2
			QL2=QL2_1
			dtHD0_1=dtHD0_1*0.5
			dtHD=dtHD0_1
			Thd=dtHD
			write(*,*) 'Intervalo de tempo muito pequeno. Reinicia loop do dia com dt=',dtHD0_1 
			if (dtHD0_1<0.01*dtHD0) then
				dtHD0_1=dtHD0*1.5
				dtHD=dtHD0_1	
				Thd=dtHD
				!@DCBteste  read(*,*)
				!@DCBteste  read(*,*)
				!@DCBteste  read(*,*)
			endif
		endif
else
	write(*,*) velcrit,NaNflag,DHreal
!	Thd=Thd-dtHD
!	dtHD=dtHD/2.0
	Thd=Thd+dtHD
	! Não houveram erros:
	HO1=HO
	QO1=QO
	! Passa contribuição lateral tempo IT para IT-1:
	QL2(:,1)=QL2(:,2)
	QL2(:,2)=0.0
endif

	ENDDO
	
	! Verifica instabilidade:
	do iX=1,nX
		if (isNaN(HO(iX)).or.HO(iX)<0.or.isNaN(QO(iX)).or.QO(iX)<-1000000.or.QO(iX)>1000000) then
			write(*,*) 'Instabilidade secao ',iX,',H=',HO(iX),',Q=',QO(iX)
			write(11222,*) 'Instabilidade secao ',iX,',H=',HO(iX),',Q=',QO(iX)
!			read(*,*)
		endif
	enddo
	courant=0.0
	velcrit=0.0
	do i=1,nreac
			if (dxt(i)==0.0) cycle
			if (dthd/dxt(i)*(9.81*(HO(NST(i,1))+HO(NST(i,2)))/2)**0.5>courant) then
				courant=dthd/dxt(i)*(9.81*(HO(NST(i,1))+HO(NST(i,2)))/2)**0.5
				ijus=NST(i,2)
				imon=NST(i,1)
				i1=i
			endif
			if (dthd/dxt(i)*max(abs(QO(NST(i,1))/(HO(NST(i,1))*BA(1,NST(i,1)))),abs(QO(NST(i,2))/(HO(NST(i,2))*BA(1,NST(i,2)))))>velcrit) then
				velcrit=dthd/dxt(i)*max(abs(QO(NST(i,1))/(HO(NST(i,1))*BA(1,NST(i,1)))),abs(QO(NST(i,2))/(HO(NST(i,2))*BA(1,NST(i,2)))))
				ijus2=NST(i,2)
				imon2=NST(i,1)
				i2=1
			endif
	enddo
!@ DCB ago/2012	write(7000,*) 'Courant max=',courant,'Secoes',imon,ijus,'dx=',dxt(i1)
!@ DCB ago/2012	write(7000,*) 'Maximo V*dt/dx = ',velcrit,'Secoes',imon2,ijus2,'dx=',dxt(i1),'V=',QO(imon2)/(HO(imon2)*BA(1,imon2)),QO(ijus2)/(HO(ijus2)*BA(1,ijus2))

!@ DCB ago/2012	write(7000,*) 'Maximo V= ',QO(maxloc(abs(QO/(HO*BA(1,:)))))/(HO(maxloc(abs(QO/(HO*BA(1,:)))))*BA(1,maxloc(abs(QO/(HO*BA(1,:)))))),'I=',maxloc(abs(QO/(HO*BA(1,:)))),'Minimo H',minval(HO),'I=',minloc(HO)

		!Escreve nos arquivos Qout e Yout
!@ DCB ago/2012if (ICALIB/=6)	write(700,700) (QO(iWrite),iWrite=1,NX)
!@ DCB ago/2012if (ICALIB/=6)	write(701,700) (HO(iWrite),iWrite=1,NX)
!@ DCB ago/2012if (ICALIB/=6)	write(702,700) ((HO(iWrite)+ZO(iWrite)),iWrite=1,NX)
!		write(703,701) ((QL2(iWrite,1)),iWrite=1,NREAC)

!@ DCB ago/2012if (ICALIB/=6)	write(8000,700) ((AFI(iWrite)),iWrite=1,NX)
!		write(8001,700) ((T(iWrite)),iWrite=1,NX)
!@ DCB ago/2012if (ICALIB/=6)	write(9000,*) iT,iafimax,afimax,idafimax,dafimax,idafimax2,dafimax2
		!OOOOO! Armazena vazoes e niveis nas extremidades de montante e jusante dos trechos:
		! Armazenar em intervalo de tempo do MGB ou intervalo de tempo do hidrodinamico?????
	! Armazena vazões nos exutórios das subbacias:
	do iTr=1,nTr
		! Subbacia:
		iC=TrCELL(iTr)		
		! Montante:
		iX=TrNST(iTr,1)
		! Vazão:
		QM2(iC)=QO(iX)
		! Jusante:
		iX=TrNST(iTr,2)
		! Vazão:
		QJ2(iC)=QO(iX)
	enddo


	! Calcula área alagada por bacia com trecho hidrodinâmico:
if (1==2) then	
	
	do iTr=1,nTr
		! Inicializa vars.::
		Aflaux=0.0

		! Subbacia:
		iC=TrCELL(iTr)
		
		! Soma áreas alagadas de cada secao do trecho:
		do iX=TrNST(iTr,1),TrNST(iTr,2)
			! Nivel dagua:
			H1=HO(iX)
			if (LD/=0) H1=H1+ZO(iX)
			
			! Interpola área alagada:	
			if (Zfl(1,iX)<=H1.and.nPfl(iX)>0) then
				! Interpola area alagada:
				! Obs.: Considera minimo entre área interpolada e área no nivel máximo
				Aflaux=Aflaux+min(FINT(Zfl(1:nPfl(iX),iX),Afl(1:nPfl(iX),iX),nPfl(iX),H1),Afl(nPfl(iX),iX))
			endif
		enddo	

		! Considera máximo entre área alagada das relacoes Z x Afl e área calculada por Afl=B*L
		Aflaux=max(Aflaux,0.5*0.001*(BA(1,TrNST(iTr,1))+BA(1,TrNST(iTr,2)))*TrL(iTr))		
		
		! Passa para percentagem:
		Aflaux=100.0*Aflaux/ACEL(IC)

		if (Aflaux>=100.0) write(*,*) 'Aflaux=',Aflaux,'%' 
		
		Aflaux=min(Aflaux,100.0)
		PUSO(IC,NU)=Aflaux
		! Preenche demais usos:
		Aflaux=(Aflaux-PUSOaux(IC,NU))
		
		if (Aflaux<0.0) then
			! Passa resto do bloco água para o primeiro bloco
			
			!***
			! Corrigir, pode dar puso>100
			PUSO(IC,1)=PUSOaux(IC,1)-Aflaux			
		else
			i=1	
			do
				if (Aflaux==0.0.or.i==NU) exit
				
				if (Aflaux>PUSOaux(IC,i)) then
					PUSO(IC,i)=0.0
					Aflaux=(Aflaux-PUSOaux(IC,i))
				else			
					PUSO(IC,i)=PUSOaux(IC,i)-Aflaux
					Aflaux=0.0
				endif
				i=i+1
			enddo
		endif
		if (sum(PUSOaux(iC,:))>100.01.or.sum(PUSOaux(iC,:))<99.99.or.minval(PUSOaux(iC,:))<0.0.or.minval(PUSOaux(iC,:))>100.0) write(*,*) 'PUSOaux',iC,iTr,sum(PUSOaux(iC,:)),minval(PUSOaux(iC,:)),minloc(PUSOaux(iC,:)),maxval(PUSOaux(iC,:)),maxloc(PUSOaux(iC,:))
		if (sum(PUSO(iC,:))>100.01.or.sum(PUSO(iC,:))<99.99.or.minval(PUSO(iC,:))<0.0.or.minval(PUSO(iC,:))>100.0) then
	!		write(*,*) 'PUSO',iC,iTr,sum(PUSO(iC,:)),minval(PUSO(iC,:)),minloc(PUSO(iC,:)),maxval(PUSO(iC,:)),maxloc(PUSO(iC,:))
			write(11001,*) 'PUSO',iT,iC,iTr,sum(PUSO(iC,:)),minval(PUSO(iC,:)),minloc(PUSO(iC,:)),maxval(PUSO(iC,:)),maxloc(PUSO(iC,:))
		endif
	enddo		

	write(11000,804)  (PUSO(TrCell(iTr),NU),iTr=1,nTr)
endif
	


	! Calcula volume na planicie de inundação e no rio:
	do i=1,nReac
		! Inicializa vars.::
		Aflaux=0.0
		iX=NST(i,2)
		iTr=TrXSec(iX)
		if (DXT(i)==0.0) cycle ! Pula trecho curto

		! Minibacia:
		iC=TrCELL(iTr)

		if (ACEL(iC)<=5.0) cycle ! Não calcula para minibacias muito pequenas.
		
		! Nivel dagua:
		H1=0.5*(HO(iX)+HO(iX-1))
		if (LD/=0) H1=H1+ZO(iX)

		! Interpola volume:	
		if (Zfl(1,iX)<=H1.and.nPfl(iX)>0) then
			! Interpola area alagada:
			! Obs.: Considera minimo entre volume interpolado e volume no nivel máximo para evitar erro na extrapolacao.
			Aflaux=Aflaux+min(FINT(Zfl(1:nPfl(iX),iX),Vfl(1:nPfl(iX),iX),nPfl(iX),H1),Vfl(nPfl(iX),iX))
		endif

		TWS6(iC)=TWS6(iC)+Aflaux/(ACEL(IC)*1000.0)

		TWS7(iC)=TWS7(iC)+0.5*(BA(1,iX)*HO(iX)+BA(1,iX-1)*HO(iX-1))*DXT(i)/(ACEL(IC)*1000.0)

		
		! Soma volume no rio tambem:
		Aflaux=Aflaux+0.5*(BA(1,iX)*HO(iX)+BA(1,iX-1)*HO(iX-1))*DXT(i)

		TWS(iC)=TWS(iC)+Aflaux/(ACEL(IC)*1000.0)

	enddo		







deallocate(HO2,QO2,QL2_1)


	700 FORMAT(<nX>F15.3)
	701 FORMAT(<nX>F15.6)
	804 FORMAT(<nTr>F10.3)

	END