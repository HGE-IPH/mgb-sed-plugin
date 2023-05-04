!@ *********************************************************************************************
!@ Criado por Diogo Buarque
!@ Data: Setembro de 2010
!@
!@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SUBROTINA DEVIVA WG2000 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!@
!@ CALCULA AS DERIVADAS PARCIAIS POR DIFERENÇAS FINITAS DE ACORDO COM:
!@ 
!@		Wilson e Gallant 2000 - Terrain Analysis: Principles and Applications.
!@
!@
!@ *********************************************************************************************

Subroutine DERIVA_WG2000(Zxaux,Zyaux)

use VAR

implicit none

real(KIND=8) Zxaux, Zyaux


							!@ DERIVADAS PARCIAIS NA DIREÇÃO X

		!@ verificando derivada em X fora dos extremos e com vizinho esquerdo e direito com valor NoData
			if(col>1 .and. col<nC .and. MNT(lin,col-1)==noDataFLAG .and. MNT(lin,col+1)==noDataFLAG) then
				Zxaux = 0.0

		!@ verificando derivada em X no extremo esquerdo ou com vizinho esquerdo com valor NoData
			elseif(col == 1 .or. MNT(lin,col-1) == noDataFLAG) then
				Zxaux = (MNT(lin,col+1) - MNT(lin,col))/(DeltaX*1000.0)

		!@ verificando derivada em X no extremo direito ou com vizinho direito com valor NoData
			elseif(col == nC .or. MNT(lin,col+1) == noDataFLAG) then
				Zxaux = (MNT(lin,col) - MNT(lin,col-1))/(DeltaX*1000.0)

		!@ obtendo derivada em X centrada
			else
				Zxaux = (MNT(lin,col+1) - MNT(lin,col-1))/(2.0*DeltaX*1000.0)
			endif
		

							!@ DERIVADAS PARCIAIS NA DIREÇÃO Y

		!@ verificando derivada em Y fora dos extremos e com vizinho superior e inferior com valor NoData
			if(lin>1 .and. lin<nL .and. MNT(lin-1,col)==noDataFLAG .and. MNT(lin+1,col)==noDataFLAG) then
				Zyaux = 0.0

		!@ verificando derivada em Y no extremo inferior ou com vizinho inferior com valor NoData
			elseif(lin == 1 .or. MNT(lin-1,col) == noDataFLAG) then
				Zyaux = (MNT(lin+1,col) - MNT(lin,col))/(DeltaY*1000.0)

		!@ verificando derivada em Y no extremo superior ou com vizinho superior com valor NoData
			elseif(lin == nL .or. MNT(lin+1,col) == noDataFLAG) then
				Zyaux = (MNT(lin,col) - MNT(lin-1,col))/(DeltaY*1000.0)

		!@ obtendo derivada em Y centrada
			else
				Zyaux = (MNT(lin+1,col) - MNT(lin-1,col))/(2.0*DeltaY*1000.0)
			endif


return

end