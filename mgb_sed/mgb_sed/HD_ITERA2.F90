	SUBROUTINE ITERA2(J1,J2,ARK1,ARK2,DXaux)
	! Resolve equacao de energia. Calcula HO da secao de montante em trecho ou confluencia divergente.
	! Utiliza método numérico para encontrar raiz de equação da Falsa Posição 
	! (ver:Chapra e Canale (2002) - Numerical Methods for Engineers) 
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
	! nItMax = numero máximo de iteracoes
	! xl e xu =  limites inferior e superior de busca
	! errorMax = erro relativo admitido em xx 
	! xx =  valor atual de xx
	! xold = valor de x na iteração anterior
	! fxl e fxu  =  valor da função nos limites inferior e superior de busca
	! fxx = valor da função em x
	! fxxold = valor da funcao em x
	! nIt = número de iterações  
	! error =   	! 
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
		
	integer:: nItMax,nIt
	real:: xl,xu,errorMax,xx,xold,fxl,fxu,fxx,error,fxxold
	
	!---------------------------------------------------------------------------------------------


	if (DXaux==0.0) then
		! Estimativa inicial:
		! Calcula diferenca entre niveis de fundo a montante e nivel dagua a jusante:
		if (ZO(J1)>HO(J2)) then
			DD=ZO(J1)-ZO(J2)
		else
			DD=0.0
		endif
		!Estimativa inicial do nivel de montante: 
		HO(J1)=HO(J2)+DD
	else


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

	! Estimativa inicial:
	! Calcula diferenca entre niveis de fundo:
	if (ZO(J1)>ZO(J2)) then
		DD=ZO(J1)-ZO(J2)
	else
		DD=0.0
	endif

	!Estimativa inicial do nivel de montante: 
	HHH=HO(J2)+DD
	xx=HHH


	! Estima limites inferior e superior de busca:
	! xl<xx<xu e fxx(xu)*fxx(xl)<0

	call FUNC(fxx,xx,XKJ,J1,J2,DXaux)
	if (fxx<0) then
		flag=1
		xl=xx
		fxl=fxx
	else
		flag=-1
		xu=xx
		fxu=fxx
	endif
	fxxold=fxx
	! Testa valoes de xx até encontrar outro limite de busca:
	nIt=0
	do
		xx=xx+0.10*flag
	
		if (xx<ZO(J1)) then
			xx=ZO(J1)+0.0000001
			call FUNC(fxx,xx,XKJ,J1,J2,DXaux)
			if (fxx*fxxold<0) exit	
			
			flag=1
			xl=xu
			fxl=fxu
			fxxold=fxl
			xx=xl
			xx=xx+0.10*flag
		endif

		call FUNC(fxx,xx,XKJ,J1,J2,DXaux)

!		write(*,*) xx,fxx,nIt
!		read(*,*)

		if (fxx*fxxold<0) exit
		nIt=nIt+1
		if (nIt==2000) then
			HO(J1)=HHH
			return
		endif

	enddo

	if (flag==1) then
		xu=xx
		fxu=fxx
	else
		xl=xx
		fxl=fxx
	endif
!	write(*,*) xl,xu,fxl,fxu

	! Inicializa:
	xx=-99999999.0
	fxx=-99999999.0
	errorMax=0.001
	nIt=0
	nItMax=200

	! Busca solução:
	do
	    xold=xx
		fxxold=fxx
	    xx=xu-fxu*(xl-xu)/(fxl-fxu)
	    
	    
	  
	    call FUNC(fxx,xx,XKJ,J1,J2,DXaux)
	    
		error=abs(fxx)

	    if (fx*fxl<0) then
	        xu=xx
	        fxu=fxx
	    elseif (fx*fxl>0) then
	        xl=xx
	        fxl=fxx
	    else
	        error=0
	    endif

		if (error<errorMax.or.nIt==nITMax) exit
		nIt=nIt+1
!		write(*,*) 'iteracao:',nIt,error,fxx,xx
!		read(*,*)
	enddo

	HO(J1)=xx

	endif
	RETURN
	END
