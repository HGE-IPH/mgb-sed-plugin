	SUBROUTINE ITERA(J1,J2,ARK1,ARK2,DXaux)
	! Resolve equacao de energia. Calcula HO da secao de montante em trecho ou confluencia divergente
	!------------------------------------------------------------------------------------------------
	! Descrição das variáveis da entrada:
	!  
	! J1,J2 = secoes de montante e jusante
	! ARK1,ARK2 = relacoes nivel x condutancia das secoes de montante e jusante (NP x 1)
	! DXaux = distancia entre secoes
	!
	! Descrição das variáveis locais:
	!
	! AM = area molhada da secao de montante
	! RM = raio hidraulico da secao de montante
	! FX = coef. de manning da secao de montante
	! HAUX = auxiliar
	! XKJ = condutancia da secao de jusante para estimativa de nivel d'água
	! R1,R2 = ponderadores
	! DD = diferenca entre niveis de fundo das secoes
	! HHH = estimativa do nivel d'água da secao de montante
	! ICON = contador de numero de iteracoes
	! XKM = condutancia da secao de montante
	! XKM2 = ((XKM+XKJ)/2)**2
	! VZ = estimativa da condutancia pela vazão e declividade do fundo (usada quando solucao não converge)
	! FINT = auxiliar interpolação
	! FLAG = auxiliar
	!
	!
	!-----------------------------------------------------------------------------------------------

	! Declaração de variáveis:
	USE PAR1_MOD
	USE TIME_MOD

	IMPLICIT NONE

	! Variáveis de entrada:
	integer,intent(in):: J1,J2
	real,intent(in):: ARK1(NP(J1)),ARK2(NP(J2))
	real,intent(in):: DXaux
	! Variáveis locais de cálculo:
	!
	real:: AM,RM,FX,XKJ,R1,R2,DD,HHH,HAUX
	integer::ICON 
	real:: XKM,XKM2,VZ,FINT
	integer:: FLAG
	!---------------------------------------------------------------------------------------------




	! Calcula área molhada, raio hidraulico, e coef. de rugosidade de manning na seção de jusante:

	HAUX=HO(J2)
	! Testa tipo de secao:
	if (NP(J2)==0) then
		HAUX=HAUX-ZO(J2)
		AM=BA(1,J2)*HAUX
		RM=HAUX
		FX=F(1,J2)
	else
		! Corrige se dados estão em termos de prof. d'água:
		if (LD==0) HAUX=HAUX-ZO(J2)

		AM=FINT(HA(1:NP(J2),J2),AR(1:NP(J2),J2),NP(J2),HAUX)
		RM=FINT(HA(1:NP(J2),J2),RR(1:NP(J2),J2),NP(J2),HAUX)

		! Verifica se Manning é constante:
		if (N(J2)==0) then
			FX=F(1,J2)
		else
			FX=FINT(HUX(1:N(J2),J2),F(1:N(J2),J2),N(J2),HAUX)
		endif
	endif

	
	! Calcula condutancia na secao de jusante
	XKJ=AM*RM**.6667/FX
!	if ((J2)<=30) write(*,*) 'am,rm,fx,XKJ',am,rm,fx,XKJ
	!Define ponderadores:
	r1=0.7
	r2=0.3
	!WRITE(6,*)XKJ,HO(J2)

	! Calcula diferenca entre niveis de fundo:
	if (ZO(J1)>ZO(J2)) then
		DD=ZO(J1)-ZO(J2)
	else
		DD=0.0
	endif

	!Estimativa inicial do nivel de montante: 
	HHH=HO(J2)+DD
	!if ((J2)<=30) write(*,*) 'DD','HHH',DD,HHH
	! Inicializa vars.:
	ICON=0
	flag=0

	! Resolve equacao da energia (iterativo):
	do while (flag==0)
		
		! Estimativa da área, raio hidraulico e coef. manning:
		HAUX=HHH
		
		! Testa tipo de secao:
		if (NP(J1)==0) then
			HAUX=HAUX-ZO(J1)
			AM=BA(1,J1)*HAUX
			RM=HAUX
			FX=F(1,J2)
		else
			! Corrige se dados estão em termos de prof. d'água:
			if (LD==0) HAUX=HAUX-ZO(J1)
				AM=FINT(HA(1:NP(J1),J1),AR(1:NP(J1),J1),NP(J1),HAUX)
				RM=FINT(HA(1:NP(J1),J1),RR(1:NP(J1),J1),NP(J1),HAUX)
    	
			! Verifica se Manning é constante:
			if (N(J1)==0) then
				FX=F(1,J2)
			else
				FX=FINT(HUX(1:N(J1),J1),F(1:N(J1),J1),N(J1),HAUX)
			endif
		endif

!		if ((J2)<=30) write(*,*) 'am rm fx',am,rm,fx
		! Condutancia da secao de montante
		XKM=AM*RM**0.6667/FX
		XKM2=((XKM+XKJ)/2)**2

		! Estimativa do nivel na secao de montante (eq. da energia s/ termos de energia cinética):
		HO(J1)=((QO(J1)+QO(J2))/2.)**2*DXaux/XKM2+HO(J2)
		
		if (HO(J1)<=ZO(J1)) then
			HO(J1)=0.1*ICON+0.1+max(HO(J2)+((QO(J1)+QO(J2))/2.)**2/XKJ**2,HO(J2)+ZO(J1)-ZO(J2))
			write(*,*) 'Corrige prof. negativa no calculo de remanso'
		endif
		
!		if ((J2)<=30) write(*,*)'XKM XKJ',' HO ZO',XKM,XKJ,HO(J1),ZO(J1)
!		if ((J2)<=30) write(*,*)'QO(J1) QO(J2) DXaux ',QO(J1), QO(J2), DXaux 
		! Conta numero de iterações:
		ICON=ICON+1
		!WRITE(6,1)ICON,HO(J1),HO(J2),DD
		
		if(ICON<=10) then
			! Testa convergencia:
			if (abs(HO(J1)-HHH)>0.05) then
				! Não convergiu.
				! Novo valor de HHH:
				HHH=r1*HO(J1)+r2*HHH
				!WRITE(6,*)ICON,HHH,HO(J1)
				
				
			else
				! Convergiu:
				flag=1
			endif
		else
			! A cada 10 iteracoes muda pesos r1 e r2:
			r1=r1-0.1
			r2=r2+0.1
			if(r1>=0.2) then
				ICON=0
			else	
				! Solucao nao converge.
				! Estima condutancia pela vazão e declividade do fundo:
				VZ=QO(J1)/(DD/DXaux)**0.5
				
				! Testa tipo de secao
				if (NP(J1)==0) then
					HO(J1)=(VZ*F(1,J1)/BA(1,J1))**(3./5.)
					write(*,*) 'Prof nao converge',HO(J1)
					HO(J1)=HO(J1)+ZO(J1)
				else		
					! Interpola nivel na tabela nivel x condutancia:
					HO(J1)=FINT(ARK1,HA(1:NP(J1),J1),NP(J1),VZ)
					! Corrige se dados estão em termos de prof. d'água:
					if (LD==0) HO(J1)=HO(J1)+ZO(J1)
				endif
				
				! Testa tipo de secao
				if (NP(J2)==0) then
					HO(J2)=(VZ*F(1,J2)/BA(1,J2))**(3./5.)
					write(*,*) 'Prof nao converge',HO(J2)
					HO(J2)=HO(J2)+ZO(J2)
				else		
					! Interpola nivel na tabela nivel x condutancia:
					HO(J2)=FINT(ARK2,HA(1:NP(J2),J2),NP(J2),VZ)
					! Corrige se dados estão em termos de prof. d'água:
					if (LD==0) HO(J2)=HO(J2)+ZO(J2)
				endif

				flag=1
			endif
		endif	
	enddo

	RETURN
	END
