	SUBROUTINE Hidrodinamico
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
	integer ijus,imon,ijus2,imon2,i1,i2
	real H1,FINT
	!----------------------------------------------------------------------------




	

	! OBS.: VAMOS CONSIDERAR TRECHOS NUMERADOS DE 1 A NTR

	! VERIFICAR TODAS ESTAS VARIAVEIS, ISTO É SO UM RASCUNHO
	

!0 1 2 3 4 5 6 7 8 9  10 11 12 13 14 15 16 17 18 19 20 21 22 23  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23  0  
!1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49

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
			DO iThd=1,nThd
				!HQB(iThd,I)=(QJ1(IC)+(iThd-1)*(QJ2(IC)-QJ1(IC))/nThd)	! Verificar se nao pega qcel
				HQB(iThd,I)=(QCEL1(IC)+(iThd-1)*(QCEL2(IC)-QCEL1(IC))/nThd)
			ENDDO
		CASE(2)
			! Resultado do MGB de série de vazão (confluencia entre bacias s/ hidrodinamico)
			! Secao:
			iX=NBO(I)
			! Trecho:
			iTr=TrXSec(iX)
			! Subbacia:
			iC=TrCELL(iTr)
			! Interpola dados para intervalo de tempo do hidrodinâmico:
			! Ou considera intervalo de tempo igual ao muskingun cunge:
			if (dtHD==DT(IC)) then
				HQB(:,I)=QCONTORM(IC,:)
			else
				HQB(1,I)=QCONTORM(IC,1)
				aux=DT(IC)
				aux2=dtHD
				J=1
				do iThd=2,nTHD
					if (aux2>aux) then
						aux=aux+DT(IC)
						J=J+1
					endif
					HQB(iThd,I)=QCONTORM(IC,J)+(QCONTORM(IC,J+1)-QCONTORM(IC,J))*(aux2-(aux-DT(IC)))/DT(IC)
					aux2=aux2+dtHD
				enddo
			endif
		CASE(3,4)
			! Séries temporais lidas em arquivo externo:
			J=bounID(I)
			HQB(:,I)=bounTS(iThd1:iThd2,J)
		ENDSELECT
	ENDDO
	!********************************************************

	! Contribuicao no intervalo de tempo anterior:
!	QL2(:,1)=QL2(:,1)

	aux=DT(NC)
	aux2=dtHD
	J=1

	write(3000,701) (HQB(ntHD,i),i=1,nboun)

	DO iThd=2,nTHD
		! Zera vetor dos coeficientes:
		AA=0.0
		BB=0.0

		write(*,*) ithd, nhmin,maxval(QO)
!@ DCB ago/2012		write(7000,*) ithd, nhmin
		!Contribuicao lateral:
		if (aux2>aux) then
			aux=aux+DT(NC)
			J=J+1
		endif

		DO I=1,NREAC
			! Testa se é contribuição lateral concentrada ou distribuida:
			IF (TrXSec(NST(I,2))/=TrXSec(NST(I,1))) THEN
				! Secao:
				iX=NST(I,2)
				! Trecho:
				iTr=TrXSec(iX)
				! Subbacia:
				iC=TrCELL(iTr)		
				! OBS.: Se dt hidrodin=dt muskingum cunge:
			!	QL2(I,2)=QCONTORM(IC,iThd)

				QL2(I,2)=QCONTORM(IC,J)+(QCONTORM(IC,J+1)-QCONTORM(IC,J))*(aux2-(aux-DT(IC)))/DT(IC)

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
				!QL2(I,2)=((QJ1(IC)+(iThd-1)*(QJ2(IC)-QJ1(IC))/nThd))*DXT(I)/(TrL(iTr)*1000)
				QL2(I,2)=((QCEL1(IC)+(iThd-1)*(QCEL2(IC)-QCEL1(IC))/nThd))*DXT(I)/(TrL(iTr)*1000)
			ENDIF
		ENDDO
		aux2=aux2+dtHD

		!Calculo dos coeficientes das matrizes A e B:
		CALL COEF1
		
		!Solucao do sistema de equações A*XI=B:
!		call GradBiConj
		CALL SKYLIN

		! Armazena variáveis de estado XI em HO e QO:
		DO I=1,NX
			! Identifica posição da seção no vetor XI
			IV=(I-1)*2+1
			HO(I)=XI(IV)
			QO(I)=XI(IV+1)

		ENDDO

		! Verifica instabilidade:
		do iX=1,nX
			if (isNaN(HO(iX)).or.HO(iX)<0.or.isNaN(QO(iX)).or.QO(iX)<-1000000.or.QO(iX)>1000000) then
				write(*,*) 'Instabilidade secao ',iX,',H=',HO(iX),',Q=',QO(iX), 'iThd=',iThd,'iT=',iT
				write(11222,*) 'Instabilidade secao ',iX,',H=',HO(iX),',Q=',QO(iX), 'iThd=',iThd,'iT=',iT
!				read(*,*)
			endif
		enddo

		! Passa contribuição lateral tempo IT para IT-1:
		QL2(:,1)=QL2(:,2)
		QL2(:,2)=0.0
	
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
!@ DCB ago/2012		write(700,700) (QO(iWrite),iWrite=1,NX)
!@ DCB ago/2012		write(701,700) (HO(iWrite),iWrite=1,NX)
!@ DCB ago/2012		write(702,700) ((HO(iWrite)+ZO(iWrite)),iWrite=1,NX)
!		write(703,701) ((QL2(iWrite,1)),iWrite=1,NREAC)

!@ DCB ago/2012		write(8000,700) ((AFI(iWrite)),iWrite=1,NX)
!		write(8001,700) ((T(iWrite)),iWrite=1,NX)
!@ DCB ago/2012		write(9000,*) iT,iafimax,afimax,idafimax,dafimax,idafimax2,dafimax2
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
	


	700 FORMAT(<nX>F15.3)
	701 FORMAT(<nX>F15.6)
	804 FORMAT(<nTr>F10.3)

	END