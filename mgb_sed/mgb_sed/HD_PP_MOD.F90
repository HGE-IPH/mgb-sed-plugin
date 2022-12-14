	MODULE PP_MOD
	! Declara??o de vari?veis globais relativas as caracteristicas hidraulicas das se??es transversais
	! em determinado intervalo de tempo.
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descri??o das vari?veis:  ***********************************
	!
	!
	! T(.) = largura em cada se??o transversal para H1 (NX x 1)
	! A(.) = area molhada da se??o transversal para H1 (NX x 1)
	! R(.) = raio hidraulico da se??o transversal para H1 (NX x 1)
	! CK(.) = condutancia hidraulica da se??o transversal para H1 (NX x 1)
	! CKY(.) = condutancia hidraulica da se??o transversal para HAUX (NX x 1)
	!			Obs.: HAUX=H1+0.01
	! SFF(.) = declividade da linha de energia na se??o transversal (NX x 1)

	!------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	SAVE


	REAL,ALLOCATABLE:: T(:),A(:),R(:),CK(:),CKY(:),SFF(:)

	END MODULE