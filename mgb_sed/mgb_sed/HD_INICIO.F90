	SUBROUTINE INICIO
	! Calcula condicoes iniciais em regime permanente
	!-------------------------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	!
	! ARK(.,.) = condutancia A*Rh^0.6667/Manning em cada seção (max(NP) x NX)
	! FLAG = auxiliar
	! FLAG1 = auxiliar
	! FX = coef. de manning
	! J1,J2 = secoes de montante e jusante
	! J,I,JJ,LL,NN,M,MT,MM,LTT,JI,K,I1,I2 = auxiliares
	! N1,N2,N3 = indices das secoes 1,2 e 3 da confluencia
	! HAU1,HAU2 = nivel (ou prof.) d´'agua auxiliar
	! AUX1,AUX2,AUX3 = auxiliares para calculo em barragem c/ efeito de jusante
	! XK1,XK2 = condutancia hidraulica nas secoes de confluencia divergente
	! SS = declividade da linha de energia (m/m)
	! ARM = condutancia hidraulica na secao de jusante em confluencia convergente
	! SO = declividade de fundo (m/m)
	! RK = condutancia estimada pela eq. Manning na estimativa da cond. de contorno jusante
	! FINT = auxiliar rotina de interpolacao
	!--------------------------------------------------------------------------------------

	! Declaração de variáveis:
	USE PAR1_MOD
	USE TIME_MOD
	USE BARR_MOD
	USE AUX_MOD

	IMPLICIT NONE

	! Variáveis locais de cálculo:

	real,allocatable:: ARK(:,:)
	integer:: FLAG,FLAG1
	real:: FX
	integer:: J1,J2,J,I,JJ,LL,NN,M,MT,MM,LTT,JI,K,I1,I2,N1,N2,N3
	real:: HAU1,HAU2,AUX1,AUX2,AUX3,XK1,XK2,SS,ARM,SO,RK,FINT
	!-------------------------------------------------------------------------------------

	! Inicializa variáveis:
	J=maxval(NP)
!@ DCB set/2012	write(*,*) 'maxvalNP',J
	allocate(ARK(J,NX))
	HO=0.0
	QO=0.0


	! Calcula relação nivel d'água (ou prof.) versus condutancia nas secoes com dados detalhados:
	do J=1,NX
		! Testa tipo de secao:
		if (NP(J)==0) cycle
		do I=1,NP(J)   
			! Verifica se Manning é constante:
			if (N(J)==0) then
				FX=F(1,J)
			else
				FX=FINT(HUX(1:N(J),J),F(1:N(J),J),N(J),HA(I,J))
			endif
			ARK(I,J)=AR(I,J)*RR(I,J)**.6667/FX
		enddo
	enddo



	!Computa condicao inicial de vazao nas secoes de condicao de contorno:
	do J=1,NBOUN
		! Verifica se condição de contorno é do tipo vazão:
		if (bounFLAG(J)<=3) then
			LL=NBO(J)
			! Armazena vazão no primeiro intervalo de tempo:
			QO(LL)=HQB(1,J)
		endif
	enddo




	!Computa condição inicial de vazão QO no restante das secoes (montante para jusante)

	! Loop dos trechos:
do iTr=1,nTr	! Para considerar subtrechos numerados nao respeitando relacao de montante para jusante

	do I=1,NREAC

		! Indice da secao de montante
		J1=NST(I,1)
if (TrXsec(J1)/=iTr) cycle	! Para considerar subtrechos numerados nao respeitando relacao de montante para jusante

		! Testa se trecho inicia em confluencia
		if (QO(J1)==0) then 
			! Se trecho inicia em confluencia a vazão em J1 ainda não foi definida 
		
			FLAG=0
			! Identifica confluencia da secao J1:
			do J=1,NCONF
				do JJ=1,3
					if (J1==iabs(NCC(J,JJ))) then
						FLAG=1
					endif
					if (FLAG==1) exit
				enddo
				if (FLAG==1) exit
			enddo

			! Seções da confluencia a qual J1 faz parte:
			N1=iabs(NCC(J,1))
			N2=iabs(NCC(J,2))
			N3=iabs(NCC(J,3))
			
			! Calcula vazão na secao J1:

			! Verifica se confluencia é convergente:
			if (NCC(J,3)<0) then
				! Confluencia convergente:
				QO(J1)=QO(N1)+QO(N2)
			else 
				! Confluencia divergente:
				
				! Adota niveis nas secoes de jusante: 
				HAU1=ZO(N1)+3
				HAU2=ZO(N2)+3
				
				! Calcula condutancia nas secoes de jusante:

				! Testa tipo de secao:
				if (NP(N1)==0) then
					HAU1=HAU1-ZO(N1)
					XK1=BA(1,N1)*(HAU1**(5./3.))/F(1,N1)
				else
					! Corrige se dados estão em termos de prof. d'água:
					if (LD==0) then
						HAU1=HAU1-ZO(N1)
					endif
					XK1=FINT(HA(1:NP(N1),N1),ARK(1:NP(N1),N1),NP(N1),HAU1)
				endif
	
				! Testa tipo de secao:
				if (NP(N2)==0) then
					HAU2=HAU2-ZO(N2)
					XK2=BA(1,N2)*(HAU2**(5./3.))/F(1,N2)
				else
					! Corrige se dados estão em termos de prof. d'água:
					if (LD==0) then
						HAU2=HAU2-ZO(N2)
					endif
					XK2=FINT(HA(1:NP(N2),N2),ARK(1:NP(N2),N2),NP(N2),HAU2)
				endif
				

				! Calcula vazões nas secoes de jusante, 
				! baseado nos niveis d'água adotados e condutancia das secoes:
				! Vazão é ponderada pela condutancia
				QO(N1)=QO(N3)*XK1/(XK1+XK2)
				QO(N2)=QO(N3)*XK2/(XK1+XK2)
				
				! Obs.: Este critério poderia ser melhorado
			endif
		endif

		! Indice da secao de jusante:
		J2=NST(I,2)

		! Contribuicao lateral ao trecho:
		
		! Vazão na secao de jusante:
		QO(J2)=QL2(I,2)+QO(J1)
	enddo
enddo


	
	
	! Nivel d'água da seção de jusante:

	! Secao:
	iX=NBO(NBOUN)
	! Considera que secao de jusante esta na ultima condicao de contorno
	! Considera que subtrechos estao numerados de montante para jusante

	if (bounFLAG(NBOUN)==4) then
		! Cond. de contorno série de nivel:
		HFIM=HQB(1,NBOUN)
		! Passa para nivel d'água:
		HFIM=HFIM+ZO(iX)
	elseif (bounFLAG(NBOUN)==5) then	
		! Cond. de contorno curva de descarga:
		J=bounID(NBOUN)
		HFIM=FINT(QT(1:NPX(J),J),HT(1:NPX(J),J),NPX(J),QO(iX))
		! Corrige se dados estão em termos de prof. d'água:
		if (LD==0) HFIM=HFIM+ZO(NX)
	elseif (bounFLAG(NBOUN)==6) then
		! Cond. de contorno eq. Manning:
		
		! Secao a montante:
		J=0
		do I=1,NREAC
			J=J+1
			if (NST(I,2)==iX) exit
		enddo
		I=NST(J,2)
		! Declividade fundo:
		SO=(ZO(I)-ZO(iX))/DXT(J) 
		

		! Verifica tipo de secao:
		if (NP(iX)==0) then
			! Seção tipo retangular:
			HFIM=(QO(iX)*F(1,iX)/(BA(1,iX)*SO**0.5))**(3./5.)										
			! Passa para nivel d'água:
			HFIM=HFIM+ZO(iX)
		else
			! Condutancia:
			RK=QO(iX)/SO**0.5

			! Nivel d'água para condutancia RK:
			HFIM=FINT(ARK(1:NP(iX),iX),HA(1:NP(iX),iX),NP(iX),RK)
			! Corrige se dados estão em termos de prof. d'água:
			if (LD==0) HFIM=HFIM+ZO(NX)
		endif
	endif

	
	
	HO(iX)=HFIM
!@ DCB set/2012	write(*,*) 'HFIM',HO(iX)
!	read(*,*)
	!Computa condição inicial de nivel HO no restante das secoes (jusante para montante):
	MT=0


do iTr=nTr,1,-1	! Para considerar subtrechos numerados nao respeitando relacao de montante para jusante


	do I=1,NREAC
		J=NREAC-I+1

if (TrXsec(NST(J,1))/=iTr) cycle	! Para considerar subtrechos numerados nao respeitando relacao de montante para jusante


		! Seção de montante
		J1=NST(J,1)
		! Secao de jusante
		J2=NST(J,2)
		! Obs.: Valores nao computados são HO=0.0 (pode causar problemas em alguns casos)
!		write(*,*) 'j1 j2',j1,j2
		! Verifica de secao de jusante faz parte de confluencia:
		if (HO(J2)<=0) then
			FLAG1=0
			! Identifica confluencia da secao J2:
			do K=1,NCONF
				do JI=1,3
					if (J2==NCC(K,JI)) then
						FLAG1=1
					endif
					if (FLAG1==1) exit
				enddo
				if (FLAG1==1) exit
			enddo
			
			
			! Seções da confluencia a qual J2 faz parte:
!			write(*,*) k,HO(J2),J1,J2
			N1=NCC(K,1)
			N2=NCC(K,2)
			N3=iabs(NCC(K,3))
		
			! Verifica se confluencia é convergente:	
			if (NCC(K,3)<=0) then
				! Convergente:
		
				! Estima niveis nas secoes de montante da 
				! confluencia pela equacao da energia s/ termos de energia cinética:
				
				HAU1=HO(N3)
				
				! Condutancia da secao de jusante:
				! Testa tipo de secao:
				if (NP(N3)==0) then
					HAU1=HAU1-ZO(N3)
					ARM=BA(1,N3)*(HAU1**(5./3.))/F(1,N3)
				else
					! Corrige se dados estão em termos de prof. d'água:
					if (LD==0) HAU1=HAU1-ZO(N3)
					
					ARM=FINT(HA(1:NP(N3),N3),ARK(1:NP(N3),N3),NP(N3),HAU1)
				endif
				! Declividade da linha de energia:
				! Obs.: Estimada na seção de jusante
				SS=(QO(N3)/ARM)**2
!				write(*,*) '**********************'
!				write(*,*) 'hau1 arm ss',hau1,arm,ss 
				! Equação da energias s/ termos de energia cinética:
				HO(N1)=HO(N3)+DXC(2*(K-1)+1)*SS
				HO(N2)=HO(N3)+DXC(2*(K-1)+2)*SS

				call ITERA2(N2,N3,ARK(1:NP(N2),N2),ARK(1:NP(N3),N3),DXC(2*(K-1)+2)) ! RP
				call ITERA2(N1,N3,ARK(1:NP(N1),N1),ARK(1:NP(N3),N3),DXC(2*(K-1)+1)) ! RP

!				write(*,*) 'n1,n2',n1,n2
!				write(*,*) 'zo(n1),zo(n2),ho(n1),ho(n2)',zo(n1),zo(n2),ho(n1),ho(n2)
!				write(*,*) '**********************'
			else
				! Confluencia divergente:
!				call ITERA(N3,N1,ARK(1:NP(N3),N3),ARK(1:NP(N1),N1),DXC(2*(K-1)+1))
				call ITERA2(N3,N1,ARK(1:NP(N3),N3),ARK(1:NP(N1),N1),DXC(2*(K-1)+1))
				! Obs.: Como garantir que HO em N1 ou N2 ja foram computados? (Paiva)
			endif
		endif

		! Verifica se existe barragem no trecho
		if (reacFLAG(J)==3) then
			! Trecho com barragem:
			! Indice da barragem
			LTT=LTOT-MT
			MT=MT+1

			! Verifica se existe efeito de jusante:
			if (NPB(LTT,1)>=0) then
				! s/ efeito jusante:
				HO(J1)= FINT(QBR(1:NPB(LTT,1),1,LTT,1),HBR(1:NPB(LTT,1),LTT,1),NPB(LTT,1),QO(J1))

				! Corrige se dados estão em termos de prof. d'água:
				if (LD==0) HO(J1)=HO(J1)+ZO(J1)

			else
				! c/ efeito jusante:
				! Procura localização na tabela do nivel na secao de jusante:
				I1=1
				I2=1
				HAU1=HO(J2)
				
				! Corrige se dados estão em termos de prof. d'água:
				if (LD==0) HAU1=HAU1-ZO(J2)
				
				
				do I2=2,NPB2(LTT,1)
					if (HAU1<=H2BR(I2,LTT,1).and.HAU1>H2BR(I2-1,LTT,1)) then 
						I1=I2-I1
						exit
					endif	
				enddo
				if (I1==1) I1=I2 ! HO maior ou menor que valores maximo e minimo de H na secao de jusante

				! Calcula para altura de jusante I1 e I2 e usa media ponderada:
				HAU1=FINT(QBR(1:NPB(LTT,1),I1,LTT,1),HBR(1:NPB(LTT,1),LTT,1),NPB(LTT,1),QO(J1))
				HAU2=FINT(QBR(1:NPB(LTT,1),I2,LTT,1),HBR(1:NPB(LTT,1),LTT,1),NPB(LTT,1),QO(J1))
				
				AUX1=ABS(HO(J2)-H2BR(I2,LTT,1))
				AUX2=ABS(HO(J2)-H2BR(I1,LTT,1))
				AUX3=AUX1+AUX2
				AUX1=AUX1/AUX3
				AUX2=AUX2/AUX3

				HO(J1)=HAU1*AUX1+HAU2*AUX2

				! Corrige se dados estão em termos de prof. d'água:
				if (LD==0) HO(J1)=HO(J1)+ZO(J1)
			endif
		else
			! Trecho sem barragem:
			! Eq. energia:
!			write(*,*) 'Trecho:', iTr
!			write(*,*) 'IteraAntes',J1,J2,HO(J1),HO(J2),QO(J1),QO(J2) 
!			call ITERA(J1,J2,ARK(1:NP(J1),J1),ARK(1:NP(J2),J2),DXT(J))
			call ITERA2(J1,J2,ARK(1:NP(J1),J1),ARK(1:NP(J2),J2),DXT(J))
!			write(*,*) 'IteraDepois',J1,J2,HO(J1),HO(J2),HO(J1)-ZO(J1),HO(J2)-ZO(J2)
			
!			read(*,*)
		endif
	enddo

enddo 	! Para considerar subtrechos numerados nao respeitando relacao de montante para jusante

	RETURN
	END