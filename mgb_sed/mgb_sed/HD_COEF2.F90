	SUBROUTINE COEF2(K1,K3,DX1,ALF,J)
	! Calcula coeficientes das equacoes de energia das confluencias
	!-------------------------------------------------------------------------------
	! Descrição das variáveis da entrada:
	!  
	! K1,K3 = indices das secoes da confluencia
	! DX1 = distancia entre secoes da confluencia 
	! ALF = coef. de energia da confluencia entre secoes K1 e K3
	! J = Linha da eq. energia na ordem original
	!
	! Descrição das variáveis locais:
	!
	! V1,V3 = velocidade nas secoes K1 e K3
	! P1P,P3,P5,P6,P8,P9,P11,P12,P13,P14,P15,P17,P18,P20 = variáveis auxiliares para calculo dos coeficientes.
	! L2,L3,L4,L5 = colunas dos coeficientes da equação.
	! K = indice que identifica linha da equação nas matrizes AA e BB
	!-------------------------------------------------------------------------------
	! Declaração de variáveis:
	USE PAR1_MOD
	USE TIME_MOD
	USE MAT_MOD
	USE PP_MOD

	IMPLICIT NONE

	! Variáveis de entrada:
	integer,intent(in):: K1,K3,J
	real, intent(in):: DX1,ALF

	! Variáveis locais de cálculo:

	real:: V1,V3 
	real:: P1P,P3,P5,P6,P8,P9,P11,P12,P13,P14,P15,P17,P18,P20
	integer:: L2,L3,L4,L5
	integer:: K
	!----------------------------------------------------------------------------

	! Linha da eq. da energia no sistema ordenado:
	K=JLIN(J)

	! Calcula coeficientes:

	! Verifica tipo de eq. da energia:
	IF (ICONF/=0) THEN
		! Eq. da energia sem termos de energia cinética e perda de carga:
		P14=0.0
		P15=1.0
		P17=0.0
		P18=-1.0
		P20=ZO(K3)-ZO(K1)
	ELSE
		! Eq. da energia com termos de energia cinética e perda de carga:
		IF(ITT/=0)THEN
			V3=QO(K3)/A(K3)
			P5=V3*V3*T(K3)/A(K3)/G
			P8=SFF(K3)/CK(K3)
			P11=HO(K3)+ZO(K3)
			P12=V3/A(K3)/G
			P13=SFF(K3)/QO(K3)
			P1P=V3*V3/2./G
		ENDIF
		V1=QO(K1)/A(K1)
		P3=V1*V1*T(K1)/A(K1)/G
		P6=DX1*SFF(K1)/CK(K1)
		P9=HO(K1)+ZO(K1)+V1*V1/2./G
		P14=V1/A(K1)/G-DX1*SFF(K1)/QO(K1)
		P15=1.-P3+P6*CKY(K1)
		P17=-ALF*P12-DX1*P13
		P18=ALF*P5-1.+DX1*P8*CKY(K3)
		P20=-P9+P11+ALF*P1P+DX1/2.*(SFF(K1)+SFF(K3))
		P20=P14*QO(K1)+P15*HO(K1)+P17*QO(K3)+P18*HO(K3)+P20
	ENDIF

	! Localizacao dos coeficientes dependentes da eq. da energia no vetor AA:
	L2=ICOL(K,2)
	L3=ICOL(K,3)
	L4=ICOL(K,4)
	L5=ICOL(K,5)


	! Armazena coef. dependentes AA:
	AA(L3)=P14
	AA(L2)=P15
	AA(L5)=P17
	AA(L4)=P18
	! Armazena coef. independente BB:
	BB(K)=P20


	RETURN
	END
