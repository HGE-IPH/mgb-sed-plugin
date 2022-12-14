	SUBROUTINE FUNC(fxx,xx,XKJ,J1,J2,DXaux)
	! Fun??o auxiliar para calculo de remanso com m?todo da Falsa Posicao
	!-----------------------------------------------------------------------
	! Descri??o das vari?veis de entrada:
	!
	! XKJ = condutancia da secao de jusante 
	! xx = estimativa do nivel na secao de montante
	! J1,J2 = secoes de montante e jusante
	! DXaux = distancia entre secoes
	!
	! Descri??o das vari?veis de sa?da:
	!
	! FINT = valor interpolado de y
	!
	! Descri??o das vari?veis locais:
	!
	! AM = area molhada da secao de montante
	! RM = raio hidraulico da secao de montante
	! FX = coef. de manning da secao de montante
	! HAUX = auxiliar
	! XKM = condutancia da secao de montante
	! XKM2 = ((XKM+XKJ)/2)**2
	!
	!----------------------------------------------------------------------
	!
	! Declara??o de vari?veis:
	USE PAR1_MOD
	USE TIME_MOD
	IMPLICIT NONE
	!
	! Vari?veis de entrada:
	real, intent(in):: xx,XKJ,DXaux
	integer,intent(in):: J1,J2
	! Vari?veis de sa?da:
	real,intent(out) :: fxx
	! Vari?veis locais de c?lculo:
	real:: HAUX,AM,RM,FX
	real:: XKM,XKM2,FINT
	!----------------------------------------------------------------------------
	


	! Estimativa da ?rea, raio hidraulico e coef. manning:
	HAUX=xx
		
	! Testa tipo de secao:
	if (NP(J1)==0) then
		HAUX=HAUX-ZO(J1)
		AM=BA(1,J1)*HAUX
		RM=HAUX
		FX=F(1,J2)
	else
		! Corrige se dados est?o em termos de prof. d'?gua:
		if (LD==0) HAUX=HAUX-ZO(J1)
		AM=FINT(HA(1:NP(J1),J1),AR(1:NP(J1),J1),NP(J1),HAUX)
		RM=FINT(HA(1:NP(J1),J1),RR(1:NP(J1),J1),NP(J1),HAUX)
    	! Verifica se Manning ? constante:
		if (N(J1)==0) then
			FX=F(1,J2)
		else
			FX=FINT(HUX(1:N(J1),J1),F(1:N(J1),J1),N(J1),HAUX)
		endif
	endif
	
	! Conducancia da secao de montante:
	XKM=AM*RM**0.6667/FX
	XKM2=((XKM+XKJ)/2)**2
	
	! Eq. energia s/ termos de energia cin?tica:
	fxx=xx-(((QO(J1)+QO(J2))/2.)**2*DXaux/XKM2+HO(J2))


	ENDSUBROUTINE FUNC
