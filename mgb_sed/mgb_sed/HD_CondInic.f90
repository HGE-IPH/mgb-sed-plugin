	SUBROUTINE HD_CondInic
	! Inicializa variáveis do modelo hidrodinamico e calcula condições iniciais
	!-------------------------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	! I,J,K,IV = contadores
	! deltaaux
	!--------------------------------------------------------------------------------------
	! Declaração de variáveis:
	USE PAR1_MOD
	USE MAT_MOD
	USE TIME_MOD
	USE BARR_MOD
	USE PP_MOD
	USE AUX_MOD
	USE VARS_MAIN
	USE CHEIA_MOD

	IMPLICIT NONE
	! Variáveis locais de cálculo:
	integer I,J,K,IV,nitera,nitera2,nitera0
	real deltaaux,deltaaux2
	real courant,Z0inicio,Z0fim
	real,allocatable:: FAFaux(:,:)
	!-------------------------------------------------------------------------------------


	Z0inicio=300.0	
	nitera=500
	nitera2=50
	nitera0=10

!***********************
	! Inicializa variáveis:


	G=9.81
	CONST=1.
	CSO=G*DThd

	! Calcula coef. theta:
	DO I=1,NREAC
		IF(DXT(I)/=0)THEN
			THETA(I)=DThd/DXT(I)
		ENDIF
	ENDDO

	ALLOCATE(PUSOaux(NC,NU))
	PUSOaux=PUSO	
	!PENSAR ONDE FICA ESTA INICIALIZACAO DE VARIAVEIS:
!*************



	write(*,*) 'Condicoes Contorno:'
	! Condições de contorno:
	DO I=1,NBOUN
		SELECTCASE (bounFLAG(I))
		CASE(1)
			! Resultado do MGB de série de vazao em bacia de cabeceira:
			! Secao:
			iX=NBO(I)
			! Trecho:
			iTr=TrXSec(iX)
			! Subbacia:
			iC=TrCELL(iTr)
			! Contribuição da subbacia:
			HQB(:,I)=QCEL2(IC)
		CASE(2)
			! Resultado do MGB de série de vazão (confluencia entre bacias s/ hidrodinamico)
			! Secao:
			iX=NBO(I)
			! Trecho:
			iTr=TrXSec(iX)
			! Subbacia:
			iC=TrCELL(iTr)
			! Interpola dados para intervalo de tempo do hidrodinâmico:
			! Ou considera intervalo de tempo igual ao muskingun cunge:
			HQB(:,I)=QM2(IC)

		CASE(3,4)
			! Séries temporais lidas em arquivo externo:
			J=bounID(I)
			HQB(:,I)=bounTS(1,J)
		ENDSELECT
	ENDDO
	
	! Contribuição lateral:
	QL2=0.0

	DO I=1,NREAC
		! Testa se é contribuição lateral concentrada ou distribuida:
		IF (TrXSec(NST(I,2))/=TrXSec(NST(I,1))) THEN
			! Secao:
			iX=NST(I,2)
			! Trecho:
			iTr=TrXSec(iX)
			! Subbacia:
			iC=TrCELL(iTr)		
			! OBS.: Se dt hidrodin=dt muskingum cunge:
			QL2(I,2)=QM2(IC)	! OBS.: Verificar esta parte, acho que vai ficar errado RP 
		ELSE
			! Secao:
			iX=NST(I,1)
			! Trecho:
			iTr=TrXSec(iX)
			! Subbacia:
			iC=TrCELL(iTr)		
			IF (OD(iC)==1) CYCLE ! Cuidar para não repetir vazao em cond.contorno e qlat em bacias de cabeceira
			! Pondera pela relacao entre dx subtrecho e comprimento total do trecho:
			QL2(I,2)=(QCEL2(IC))*DXT(I)/(TrL(iTr)*1000)	! Verificar o que pegar aqui. Acho que Qcell é suficiente
		ENDIF
		
	ENDDO
	QL2(:,1)=QL2(:,2)
	IF (INIC==0) RETURN


	! Desativa planicie de inundacao:
	I=MAXVAL(NPF)
	allocate(FAFaux(I,nX))
	FAFaux=FAF
	FAF=0.0






	! Considera inicialmente rios cheios:
	if (bounFLAG(NBOUN)==4) then
		! Cond. de contorno série de nivel:
		iX=NBO(NBOUN)
		HQB(:,NBOUN)=Z0inicio-ZO(iX)
	endif	

	! Calcula condicoes iniciais baseado na eq. da continuidade e energia no regime permanente:
	IF (INIC/=0) THEN
		CALL INICIO ! RP REVER ESTA ROTINA INICIO e ITERA (MUITA CHANCE DE DAR PAU)
		! Passa dados de níveis d'água para profundidade
		HO=HO-ZO
		
	ENDIF



!	HO=Z0inicio-ZO ! Considera nivel d'água horizontal.	

	do iX=1,nX
		! Condicao inicial de nivel ou profundidade e vazão em cada secao:
		write(3001,'(I8,3F15.3)') iX,HO(iX),HO(iX)+ZO(iX),QO(iX) 
	enddo


	! Aquece modelo a partir de condição inicial calculada pela equação da energia:
	
	! Armazena condições iniciais QO e HO no vetor de variáveis de estado XI:
	DO I=1,NX
		IV=(I-1)*2+1
		XI(IV)=HO(I)
		XI(IV+1)=QO(I) 
	ENDDO
	
	Z0fim=Z0inicio
	! Define Z0fim
	do I=1,NBOUN
		if (bounFLAG(I)/=4) cycle
		! Secao:
		iX=NBO(I)
		! Séries temporais lidas em arquivo externo:
		J=bounID(I)
		Z0fim=min(Z0fim,bounTS(1,J)+ZO(iX))
	enddo
	if (Z0fim==Z0inicio) then
		write(*,*) 'Para calcular condicoes iniciais configure no '
		write(*,*) 'minimo uma condicao de contorno de prodf. ou nivel dagua.'
		stop ''
	endif


	Z02=Z0inicio
	do iT=1,nitera0+nitera+nitera2
		
		Z01=Z02
		if (iT<=nitera0+nitera.and.iT>nitera0) then
			Z02=Z02-(Z0inicio-Z0fim)/nitera
		endif

		write(*,*) 'iT,minH,secao minH,Q min secao,H(nx),Hboun:',it,minval(HO),minloc(ho),qo(minloc(ho)),HO(nX),HQB(1,NBOUN)
!@ DCB ago/2012		write(7000,*) 'iT,minH,secao minH,Q min secao,H(nx),Hboun:',it,minval(HO),minloc(ho),qo(minloc(ho)),HO(nX),HQB(1,NBOUN)
		write(*,*) minval(QO),minloc(Qo),maxval(QO),maxloc(Qo)


        DHmax=2.00+10.0*(nitera0+nitera+nitera2-iT)/(nitera0+nitera+nitera2)
        DHmax=max(2.0,DHmax)
        VcritMAX=2.0+7.0*(nitera0+nitera+nitera2-iT)/(nitera0+nitera+nitera2)
        VcritMAX=max(2.0,VcritMAX)
        DQmax=30000.0+100000.0*(nitera0+nitera+nitera2-iT)/(nitera0+nitera+nitera2)
        DQmax=max(DQmax,30000.0)
        Hmin=0.30
		
		call hidrodinamico2
	enddo

	 
	! Reativa planicie de inundacao:
	FAF=FAFaux
	!***
	
	write(*,*) 'Fim do calculo de condicoes iniciais do mod. Hidrodinamico'

	! Escreve arquivo de condicoes iniciais:
	OPEN(301,FILE='.\output\HD_CondInic.txt')
	WRITE(301,*) 'Secao     H      Q'
	do iX=1,nX
		! Condicao inicial de nivel ou profundidade e vazão em cada secao:
		write(301,'(I8,2F15.3)') iX,HO(iX),QO(iX) 
	enddo
	close(301)
		


	700 FORMAT(200F15.3)
	701 FORMAT(200F15.6)

	
	RETURN
	END