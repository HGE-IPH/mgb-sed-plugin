	MODULE TIME_MOD
	! Declaração de variáveis globais relativas vazões, niveis nas seções, contribuição lateral, etc.
	! 
	!
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descrição das variáveis:  ***********************************
	!
	!
	! HO(.) = nível ou profundidades d'agua nas secoes (NX x 1) 
	! QO(.) = vazão nas secoes (NX x 1)
	!rever HQB(.,.) = series temporais de niveis ou vazões das condições de contorno (IFIN x NBOUN)
	! QL2(.,.) = contribuição lateral (NREAC x 2)
	!			QL2(.,1) = intervalo de tempo i
	!			QL2(.,2) = intervalo de tempo i+1
	! DThd= intervalo de tempo dos dados (s)
	! DT1 = intervalo de tempo de calculo
	! NThd = numero de intervalos do modelo hidrodinâmico para 1 intervalo de tempo do MGB
	! NThd2 = numero de intervalos do modelo hidrodinâmico total
	! QT(.,.) = vazões da tabela da curva de descarga (max pontos na tabela x NBOUNRC)
	! HT(.,.) = niveis ou prof. da tabela da curva de descarga (max pontos na tabela x NBOUNRC)
	! NPX(.) = numero de pontos da tabela de curva de descarga (NBOUNRC)


	! iThd = intervalo de tempo
	! G = aceleração da gravidade (m/s2)
	! CSO = DT*G
	! AFI(.) = largura da planicie de inundação em cada secao em metros (NX x 1)
	!------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	SAVE

	REAL,ALLOCATABLE:: HO(:),QO(:),HQB(:,:),QL2(:,:)
	REAL DThd,DT1
	INTEGER NThd,NThd2
	REAL,ALLOCATABLE:: QT(:,:),HT(:,:)
	INTEGER,ALLOCATABLE:: NPX(:)
	INTEGER iThd
	REAL G
	REAL CSO
	REAL,ALLOCATABLE:: AFI(:)


	REAL,ALLOCATABLE:: HOprev(:),QOprev(:)


	END MODULE