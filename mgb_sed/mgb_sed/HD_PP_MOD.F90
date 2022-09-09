	MODULE PP_MOD
	! Declaração de variáveis globais relativas as caracteristicas hidraulicas das seções transversais
	! em determinado intervalo de tempo.
	!------------------------------------------------------------------------------------------------
	!
	! ********************************  Descrição das variáveis:  ***********************************
	!
	!
	! T(.) = largura em cada seção transversal para H1 (NX x 1)
	! A(.) = area molhada da seção transversal para H1 (NX x 1)
	! R(.) = raio hidraulico da seção transversal para H1 (NX x 1)
	! CK(.) = condutancia hidraulica da seção transversal para H1 (NX x 1)
	! CKY(.) = condutancia hidraulica da seção transversal para HAUX (NX x 1)
	!			Obs.: HAUX=H1+0.01
	! SFF(.) = declividade da linha de energia na seção transversal (NX x 1)

	!------------------------------------------------------------------------------------------------
	IMPLICIT NONE
	SAVE


	REAL,ALLOCATABLE:: T(:),A(:),R(:),CK(:),CKY(:),SFF(:)

	END MODULE