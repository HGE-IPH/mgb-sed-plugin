	SUBROUTINE IPH4_F90
	!PROGRAM IPH4_F90
	!********************************************************************************************
	!										    IPH4 
	!										   ------
	! 
	!						Modelo Hidrodinamico IPH4 Vers?o Fortran 90
	! 
	!							INSTITUTO DE PESQUISAS HIDRAULICAS
	!					da Universidade Federal do Rio Grande do Sul- UFRGS 
	!
	!  
	!  - Vers?o Original em Fortran 77 (mar/1996):
	!
	!		Carlos E.M. Tucci
	!
	!  - Vers?o Atual em Fortran 90:
	!
	!		Rodrigo Paiva e Juan Martin Bravo
	!  
	!  - Detalhes sobre aspectos te?ricos:
	!       
	!    Tucci, C.E.M. 1978. Hydraulic and Water Quality Model for a River Network. 
	!			PhD dissertation, Colorado State University, Fort Collins, USA.
	!
	!    Tucci, C.E.M. 2005. Modelos Hidrol?gicos. Ed. UFRGS/ABRH, Porto Alegre.
	!
	!********************************************************************************************
	!-----------------------------------------------------------------------
	! Descri??o das vari?veis locais:
	!
	! I,J,L,IVK,M = vari?veis auxiliares
	! IS = numero de linhas vetor AA
	! TIME1,TIME2 = auxiliares para calculo do tempo de processamento
	!----------------------------------------------------------------------
	! Declara??o de vari?veis:
	USE PORTLIB !biblioteca para calcular o tempo de processamento
	USE PAR1_MOD
	USE MAT_MOD
	USE TIME_MOD
	USE BARR_MOD
	USE PP_MOD
	IMPLICIT NONE

	! Vari?veis locais de c?lculo:
	integer I,J,L
	integer IV
	integer IS
	integer iWrite
	INTEGER TIME1,TIME2
	REAL TIME3
	!----------------------------------------------------------------------------


	! Contador de tempo:
	TIME1=TIME()

	!CALL CPU_TIME (TIME1)
	write(*,*) time1



	! Inicializa vari?veis:
	G=9.81
	CONST=1.

	! Leitura e organiza??o de dados de entrada:
	CALL INPUT2

	!Organizacao da matriz de coeficientes:
	CALL MATRIX

	! Inicializa vari?veis:
	CSO=G*DThd
	QL2=0.0  ! Zera contribui??o lateral


	! Entrada da contribuicao lateral do primeiro intervalo de tempo:
!	IF (NQS>0) THEN 
!		DO I=1,NX
			! Indice da contribuicao lateral:
			!L=LQ(I)
			!IF (L>0) QL2(I,1)=QWL(1,L)+QL2(I,1)
!		ENDDO
!	ENDIF

	! Armazena condi??es iniciais QO e HO no vetor de vari?veis de estado XI:
	DO I=1,NX
		IV=(I-1)*2+1
		XI(IV)=HO(I)
		XI(IV+1)=QO(I) 
	ENDDO
	DO I=1,NREAC
		! Calcula coef. theta:
		IF(DXT(I)/=0)THEN
			THETA(I)=DThd/DXT(I)
		ENDIF
	ENDDO


	! Calcula numero de elementos do vetor AA:
	IS=0
	DO J=1,NUM
		! Soma elementos a esquerda e acima da diagonal princial:
		IS=IS+IR(J)+IHIGH(J)
	ENDDO


	!Escreve arquivos Qout.txt e Yout.txt
!@ DCB ago/2012	open(700,FILE='.\output\Qout.txt')
!@ DCB ago/2012	open(701,FILE='.\output\Yout.txt')
!@ DCB ago/2012	write(700,700) (QO(iWrite),iWrite=1,NX)
!@ DCB ago/2012	write(701,700) (HO(iWrite),iWrite=1,NX)

	DO iThd=2,nThd
		! Zera vetor dos coeficientes:
		AA=0.0
		BB=0.0
		!Contribuicao lateral
!		IF(NQS>0)THEN
!			DO I=1,NX
				! Indice da contribuicao lateral:
!				J=LQ(I)
				!IF (J>0) QL2(I,2)=QWL(iThd,J)+QL2(I,2)
!			ENDDO
!		ENDIF

		!Calculo dos coeficientes das matrizes A e B
		CALL COEF1

		!Solucao do sistema de equa??es A*XI=B:
		CALL SKYLIN

		! Armazena vari?veis de estado XI em HO e QO:
		DO I=1,NX
			! Identifica posi??o da se??o no vetor XI
			IV=(I-1)*2+1
			HO(I)=XI(IV)
			QO(I)=XI(IV+1)
		ENDDO

		! Passa contribui??o lateral tempo IT para IT-1:
!		IF(NQS>0) THEN
!			QL2(:,1)=QL2(:,2)
!			QL2(:,2)=0.0
!		ENDIF	

		!Escreve nos arquivos Qout e Yout
!@ DCB ago/2012		write(700,700) (QO(iWrite),iWrite=1,NX)
!@ DCB ago/2012		write(701,700) (HO(iWrite),iWrite=1,NX)

	ENDDO

	! Verifica??o tempo processamento:
	TIME2=TIME()
	TIME3=REAL(TIME2-TIME1)

	WRITE(*,*) 'Tempo total: ', TIME3,' s'
	read(*,*)

	700 FORMAT(200F10.3)


	WRITE(*,*)'Fim da simulacao'
	read(*,*)

	!END PROGRAM IPH4_F90
	END