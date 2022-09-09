	subroutine CalibParam
	! Subrotina para passar parametros do MOCOM-UA para MGB-IPH
	!-----------------------------------------------------------------------
	!
	!
	use VARS_MAIN
	use VARS_CALIB
	IMPLICIT NONE
	! Declaração de variáveis:
	! Variáveis locais:
	integer:: i1,i2,i3,L,i,j,k
	!----------------------------------------------------------------------------

		! Linhas com novo codigo:
		WM=WMOLD
		B=BOLD
		KINS=KIOLD
		KBAS=KBOLD
		PLAM=PLAMOLD
		CAP=CAPOLD
		WC=WCOLD
		CS=CSOLD
		CI=CIOLD
		CB=CBOLD
		L=0
		do i1=1,NU
			do i2=1,7
				if (p_Calib(i1,i2)==0) cycle

				if (p_Calib(i1,i2)==1) then
					L=L+1
					selectcase (i2)
					case(1)
						WM(:,i1)=WMOLD(:,i1)*PARX(L)
					case(2)
						B(:,i1)=BOLD(:,i1)*PARX(L)
					case(3)
						KBAS(:,i1)=KBOLD(:,i1)*PARX(L)
						!KINS(:,i1)=KIOLD(:,i1)*PARX(L)
					case(4)
						KINS(:,i1)=KIOLD(:,i1)*PARX(L)
						!KBAS(:,i1)=KBOLD(:,i1)*PARX(L)
					case(5)
						PLAM(:,i1)=PLAMOLD(:,i1)*PARX(L)
					case(6)
						CAP(:,i1)=CAPOLD(:,i1)*PARX(L)
					case(7)
						WC(:,i1)=WCOLD(:,i1)*PARX(L)
					endselect
				elseif (p_Calib(i1,i2)<0) then
					i3=-p_Calib(i1,i2)
					selectcase (i2)
					case(1)
						WM(:,i1)=WMOLD(:,i1)*WM(:,i3)/WMOLD(:,i3)
					case(2)
						B(:,i1)=BOLD(:,i1)*B(:,i3)/BOLD(:,i3)
					case(3)
						KBAS(:,i1)=KBOLD(:,i1)*KBAS(:,i3)/KBOLD(:,i3)
						!KINS(:,i1)=KIOLD(:,i1)*KINS(:,i3)/KIOLD(:,i3)
					case(4)
						KINS(:,i1)=KIOLD(:,i1)*KINS(:,i3)/KIOLD(:,i3)
						!KBAS(:,i1)=KBOLD(:,i1)*KBAS(:,i3)/KBOLD(:,i3)

					case(5)
						PLAM(:,i1)=PLAMOLD(:,i1)*PLAM(:,i3)/PLAMOLD(:,i3)
					case(6)
						CAP(:,i1)=CAPOLD(:,i1)*CAP(:,i3)/CAPOLD(:,i3)
					case(7)
						WC(:,i1)=WCOLD(:,i1)*WC(:,i3)/WCOLD(:,i3)
					endselect
				endif
			enddo
		enddo
		if (p_Calib(NU+1,1)==1) then
			L=L+1
			CS=CSOLD*PARX(L)
		endif
		if (p_Calib(NU+2,1)==1) then
			L=L+1
			CI=CIOLD*PARX(L)
		endif
		if (p_Calib(NU+3,1)==1) then
			L=L+1
			CB=CBOLD*PARX(L)
		endif

		! Linhas para fixar parametros da algumas subbacias:
		if (NCONGEL>0) then
			do i1=1,NCONGEL
				IB=IBCONGEL(i1)
				WM(IB,:)=WMOLD(IB,:)
				B(IB,:)=BOLD(IB,:)
				KINS(IB,:)=KIOLD(IB,:)
				KBAS(IB,:)=KBOLD(IB,:)
				PLAM(IB,:)=PLAMOLD(IB,:)
				CAP(IB,:)=CAPOLD(IB,:)
				WC(IB,:)=WCOLD(IB,:)
				CS(IB)=CSOLD(IB)
				CI(IB)=CIOLD(IB)
				CB(IB)=CBOLD(IB)
			enddo
		endif
		



!		write(*,*) 'Teste Passagem parametros'
!		write(*,*) PARX
!		read(*,*)
!		IB=1
!		write(*,*) 'subbacia',IB
!		do i1=1,NU
!				write(*,*) 'Bloco',i1
!				write(*,*) WM(IB,i1)/WmOLD(IB,i1),B(IB,i1)/BOLD(IB,i1),KINS(IB,i1)/KIOLD(IB,i1),KBAS(IB,i1)/KBOLD(IB,i1),PLAM(IB,i1)/PLAMOLD(IB,i1),CAP(IB,i1)/CAPOLD(IB,i1),WC(IB,i1)/WCOLD(IB,i1)  
!
!		enddo
!			write(*,*) CS(IB)/CSOLD(IB)
!			write(*,*)	CI(IB)/CIOLD(IB)
!			write(*,*)	CB(IB)/CBOLD(IB)
!
!		read(*,*) 
!
!		IB=2
!		write(*,*) 'subbacia',IB
!		do i1=1,NU
!			
!				write(*,*) 'Bloco',i1
!				write(*,*) WM(IB,i1)/WmOLD(IB,i1),B(IB,i1)/BOLD(IB,i1),KINS(IB,i1)/KIOLD(IB,i1),KBAS(IB,i1)/KBOLD(IB,i1),PLAM(IB,i1)/PLAMOLD(IB,i1),CAP(IB,i1)/CAPOLD(IB,i1),WC(IB,i1)/WCOLD(IB,i1)  
!			
!		enddo
!		write(*,*) CS(IB)/CSOLD(IB)
!		write(*,*)	CI(IB)/CIOLD(IB)
!		write(*,*)	CB(IB)/CBOLD(IB)

!		read(*,*) 

			

	RETURN
	END