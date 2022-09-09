	subroutine LeCalib
	! Subrotina para leitura do arquivo ParCalib.mgb
	!-----------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	!
	! I,J = variáveis auxiliares
	!
	!----------------------------------------------------------------------
	!
	use VARS_MAIN
	use VARS_CALIB
	IMPLICIT NONE
	! Declaração de variáveis:
	! Variáveis locais:
	integer:: i,j,k
	character(60):: trashTex
	character(10):: trashTex2
	!----------------------------------------------------------------------------


	WRITE(*,*)'Le arquivo de parametros de calibracao:'

	OPEN(FILCALIB,FILE='.\input\ParCalib.mgb',STATUS='OLD',ACTION='READ')

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) NS !Numero de individuos na populacao do MOCOM-UA
	WRITE(*,*)NS
read(*,*)



	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) NF !Numero de funcoes objetivo
	WRITE(*,*) NF
	read(*,*)
	
	allocate(iFO(nF))

	do i=1,4
		READ(FILCALIB,75) trashTex
		WRITE(*,*)trashTex
	enddo	
	READ(FILCALIB,*) (iFO(i),i=1,NF) !Funcoes objetivo utilizadas
	WRITE(*,*) (iFO(i),i=1,NF)
	read(*,*)

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) iMaxGen
	WRITE(*,*)iMaxGen

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) NCONGEL !Numero de subbacias com parametros fixos
	WRITE(*,*) NCONGEL
	
	allocate(IBCONGEL(NCONGEL))

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,*) (IBCONGEL(i),i=1,NCONGEL) !Numero das subbacias com parametros fixos
	WRITE(*,*) (IBCONGEL(i),i=1,NCONGEL)

	read(*,*)



	allocate(p_Calib(NU+3,7))
	read(*,*)

	p_Calib=0

	do i=1,5
		READ(FILCALIB,75) trashTex
		WRITE(*,*)trashTex
	enddo

	do i=1,NU
		READ(FILCALIB,73) trashTex2,(p_Calib(i,j),j=1,7)
		WRITE(*,73) trashTex2,(p_Calib(i,j),j=1,7)
	enddo

	do i=NU+1,NU+3
		READ(FILCALIB,73) trashTex2,p_Calib(i,1)
		WRITE(*,73) trashTex2,p_Calib(i,1)
	enddo

	NPAR=0
	do i=1,NU+3
		do j=1,7
			if (p_Calib(i,j)==1) NPAR=NPAR+1
		enddo
	enddo
	WRITE(*,*)'Numero de parametros:'
	WRITE(*,*) NPAR
	read(*,*)
	

	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	
	ALLOCATE (PMIN(NPAR),PMAX(NPAR)) ! Limites dos parametros.
	do i=1,NPAR
		READ(FILCALIB,*) PMIN(i),PMAX(i),trashTex2
		WRITE(*,*) PMIN(i),PMAX(i),trashTex2
	enddo

	allocate (CalibFLAG(NOBS,nf)) ! RP
	
	READ(FILCALIB,75) trashTex
	WRITE(*,*)trashTex
	do i=1,NOBS
		READ(FILCALIB,*) (calibFLAG(i,j),j=1,NF)
	enddo










	CLOSE (FILCALIB)

71	FORMAT(6I10)
72	FORMAT(5A10)
73	FORMAT(A10,7I6)
74	FORMAT(A20)
75	FORMAT(A60)
76	FORMAT(I10,F10.1)
77	FORMAT(A20)
78	FORMAT(A10,1I6)
!79	FORMAT(A10,2F6)

	RETURN
	END