	MODULE TIME_MOD
	! Declara??o de vari?veis globais relativas vaz?es, niveis nas se??es, contribui??o lateral, etc.
	! 
	!
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descri??o das vari?veis:  ***********************************
	!
	!
	! HO(.) = n?vel ou profundidades d'agua nas secoes (NX x 1) 
	! QO(.) = vaz?o nas secoes (NX x 1)
	!rever HQB(.,.) = series temporais de niveis ou vaz?es das condi??es de contorno (IFIN x NBOUN)
	! QL2(.,.) = contribui??o lateral (NREAC x 2)
	!			QL2(.,1) = intervalo de tempo i
	!			QL2(.,2) = intervalo de tempo i+1
	! DThd= intervalo de tempo dos dados (s)
	! DT1 = intervalo de tempo de calculo
	! NThd = numero de intervalos do modelo hidrodin?mico para 1 intervalo de tempo do MGB
	! NThd2 = numero de intervalos do modelo hidrodin?mico total
	! QT(.,.) = vaz?es da tabela da curva de descarga (max pontos na tabela x NBOUNRC)
	! HT(.,.) = niveis ou prof. da tabela da curva de descarga (max pontos na tabela x NBOUNRC)
	! NPX(.) = numero de pontos da tabela de curva de descarga (NBOUNRC)


	! iThd = intervalo de tempo
	! G = acelera??o da gravidade (m/s2)
	! CSO = DT*G
	! AFI(.) = largura da planicie de inunda??o em cada secao em metros (NX x 1)
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