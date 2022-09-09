!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA DEVIVA ArcGIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ CALCULA AS DERIVADAS PARCIAIS POR DIFERENÇAS FINITAS DE ACORDO COM:
!@ 
!@		ARCINFO do ArcGIS 9.2
!@
!@
!@ *********************************************************************************************

Subroutine DERIVA_ArcGIS(Zxaux,Zyaux)

use VAR

implicit none

real(KIND=8) Zxaux, Zyaux

integer(2) MNTaux(3,3), ii, jj, iaux, jaux


MNTaux = 0
iaux = 0
do ii = lin-1,lin+1
	iaux = iaux + 1
	jaux = 0
	do jj = col-1,col+1
		jaux = jaux + 1

		if (ii<1 .or. ii>nL .or. jj<1 .or. jj>nC .or. MNT(ii,jj)==noDataFLAG) then
		!@ Atribui a cota do pixel em analise aos pixels vizinhos fora da área de trabalho ou com valores noData
			MNTaux(iaux,jaux) = MNT(lin,col)
		else
			MNTaux(iaux,jaux) = MNT(ii,jj)
		endif

	enddo
enddo


!@ DERIVADAS PARCIAIS NA DIREÇÃO X
Zxaux =      (MNTaux(3,3) + 2*MNTaux(2,3) + MNTaux(1,3))
Zxaux = Zxaux - (MNTaux(3,1) + 2*MNTaux(2,1) + MNTaux(1,1))
Zxaux = Zxaux / (8.0*DeltaX*1000.0)

!@ DERIVADAS PARCIAIS NA DIREÇÃO Y
Zyaux =      (MNTaux(1,1) + 2*MNTaux(1,2) + MNTaux(1,3))
Zyaux = Zyaux - (MNTaux(3,1) + 2*MNTaux(3,2) + MNTaux(3,3))
Zyaux = Zyaux / (8.0*DeltaY*1000.0)


return

end