	SUBROUTINE SOLO(PX,EX,WX,WMX,BX,KIX,KBX,XL,DSUP,DINT,DBAS,CAPX,WCX,DTP)
	!esta subrotina faz o balanço de água
	!na camada de solo, gerando
	!1 drenagem superficial ou direta
	!2 drenagem sub-superficial ou intermediária
	!3 drenagem subterrânea ou de base
	IMPLICIT NONE

	!VARIAVEIS DE ESTADO
	REAL WX		! ESTADO DE ARMAZENAMENTO DO SOLO
	!PARÂMETROS
	REAL BX,WMX ! ARMAZENAMENTO DO SOLO
	REAL KIX,KBX !CONDUTIVIDADE HIDRÁULICA
	REAL WZ,WB !LIMITES PARA DRENAGEM INT E BAS
	REAL XL !PARÂMETRO DE FORMA DA DRENAGEM INT
	REAL CAPX,WCX,WC !FLUXO CAPILAR ASCENDENTE
	!PRECIPITAÇÃO, EVAPORAÇÃO
	REAL PX,EX
	!DRENAGEM SUPERFICIAL, INTERMEDIARIA, BASICA
	REAL DSUP,DINT,DBAS
	!VARIÁVEIS AUXILIARES
	REAL VTEMP
	REAL WXOLD
	REAL DTP !INTERVALO DE TEMPO PRINCIPAL DE CÁLCULO (EM SEGUNDOS)


	!			PARAMETROS 

	!XL É O PARAMETRO DA RELAÇÀO ENTRE A DRENAGEM SUB-SUPERCIAL E A UMIDADE DO SOLO 
	!DE BROOKS E COREY (VER LIVRO DO MAIDMENT)
	!XL É O pore size index E TEM VALORES ENTRE 0.165 PARA ARGILAS E 0.694 PARA AREIA, CONFORME RAWLS
	!XL FOI DEFINIDO NO ARQUIVO PARUSO.HIG

	!WB É O LIMITE PARA OCORRER ESCOAMENTO SUBTERRÂNEO 
	WB=WCX*WMX !WCX É A VARIÁVEL TEMPORARIA DE WC(IB,IU) QUE É UM PARAMETRO DO ARQUIVO PARUSO.HIG

	!WZ É O LIMITE PARA OCORRER ESCOAMENTO SUB-SUPERFICIAL
	WZ=WB !POR SIMPLICIDADE É CONSIDERADO IGUAL A WB

	!WC É O LIMITE PARA INICIAR O FLUXO CAPILAR ASCENDENTE
	WC=WCX*WMX

	!******CALCULA O ESCOAMENTO DIRETO (SUPERFICIAL)***
	IF(PX.GT.0.0) THEN
		VTEMP=(1.-WX/WMX)**(1./(BX+1.))-(PX/((BX+1.)*WMX))
		IF(VTEMP.LE.0.0) THEN ! A CHUVA SATURA TUDO
!@ DCB 26-05-2011 **************************************************************************
!				DSUP=PX-(WMX-WX)
				DSUP=max(PX-(WMX-WX),0.0)
		ELSE                  ! A CHUVA NAO SATURA TUDO
!				DSUP=PX-(WMX-WX)+WMX*(VTEMP)**(BX+1.)
				DSUP=max(PX-(WMX-WX)+WMX*(VTEMP)**(BX+1.),0.0)
!@ DCB 26-05-2011 **************************************************************************
		ENDIF
	ELSE
		DSUP=0.0 !SE NÃO CHOVE NAO TEM ESCOAMENTO DIRETO
	ENDIF
	!**********FIM DO ESCOAMENTO DIRETO***************

	!**CALCULA O ESCOAMENTO MUITO LENTO (SUBTERRANEO)**        
	IF(WX.GE.WB) THEN
		DBAS=KBX*(WX-WB)/(WMX-WB)!RELACAO PROVISORIA
	ELSE
		DBAS=0.0
	ENDIF
	!******FIM DO ESCOAMENTO MUITO LENTO (SUBTERRANEO)*

	!**CALCULA ASCENÇÃO POR CAPILARIDADE**************
	CAPX=MAX(CAPX*(WC-WX)/WC,0.0) !BREMICKER
	!***FIM DO CALCULO DA ASCENÇÃO POR CAPILARIDADE***

	!****CALCULA O ESCOAMENTO INT (SUB-SUPERFICIAL)**
	!UTILIZA UMA VERSAO DA RELAÇAO DE BROOKS E COREY (VER PÁGINA 5.6 DO MAIDMENT Handbook of Hydrology)
	IF(WX.LE.WZ) THEN
		DINT=0.0 !DRENAGEM INT É NULA SE UMIDADE BAIXA
	ELSE
		DINT=KIX*((WX-WZ)/(WMX-WZ))**(3.+2./XL)	
	ENDIF

	!ESTA PARTE LIMITA O VALOR DE DINT PARA NAO OCORRER W NEGATIVO 
	DINT=MIN(DINT,WX/10.0) !este ajuste se justifica porque 
	DINT=MAX(DINT,0.0)     !a cond. hidr. se reduz muito rapido
	!**********FIM DO ESCOAMENTO INT (SUB-SUPERF.)***

	WXOLD=WX
	CAPX=CAPX*DTP/86400. !CONVERTE PARA MM/DTP
	DSUP=DSUP !CONVERTE PARA MM/DTP
	DINT=DINT*DTP/86400. !CONVERTE PARA MM/DTP
	DBAS=DBAS*DTP/86400.	!CONVERTE PARA MM/DTP

	DO !SÓ SAI DESTE LOOP QUANDO A EQUAÇÃO NÃO RESULTA EM ARMAZENAMENTO NEGATIVO
		!WC DTP WX=WX+PX-EX-DSUP-DINT-DBAS+CAPX  !FAZ O BALANÇO HIDRICO
		WX=WX+PX+(CAPX-EX-DSUP-DINT-DBAS)  !FAZ O BALANÇO HIDRICO


!		! Ajuste de balanço antigo (Walter)
!		WX=MIN(WX,WMX) !IMPEDE ERRO DE ARREDONDAMENTO
!		IF(WX.LT.0.0)THEN
!!RP			WRITE(*,*) 'SOLO SECOU'
!			!REDUZ A EVAPORAÇÃO PORQUE SOLO ESTÁ SECANDO
!			EX=WXOLD+PX-DSUP-DINT-DBAS 
!			WX=WXOLD
!			CYCLE
!		ELSE
!			EXIT
!		ENDIF

		! Ajuste de balanço novo (Paiva)
		! Considera correcao na Evapo e escoamentos sup, int e sub.
		if (WX<WMX+0.01.and.WX>WMX) WX=WMX
		if (WX>-0.01.and.WX<0.0) WX=0.0
		if(WX>WMX)then
			! Se armazenamento no solo é maior que o máximo
			! aumenta escoamento superficial:
			DSUP=DSUP+WX-WMX
			WX=WXOLD
!write(*,*) '4'

		elseif (WX<0.0) then
!write(*,*) '4.5'
			! Diminui evapotranspiração:

			EX=EX-(-WX)
			if (EX>=0.0) then
				WX=WXOLD

				cycle
				
			endif
			! Evapo zerou e balanco ainda nao fecha:
			WX=EX ! Valor de WX é igual ao residuo
			EX=0.0
!write(*,*) '5'

			! Diminui percolação:
			DBAS=DBAS-(-WX)
			if (DBAS>=0.0) then
				WX=WXOLD
				cycle
			endif
			! Percolacao zerou e balanco ainda nao fecha:
			WX=DBAS ! Valor de WX é igual ao residuo					
			DBAS=0.0
!write(*,*) '6'

			! Diminui esc. subsuperficial:
			DINT=DINT-(-WX)
			if (DINT>=0.0) then
				WX=WXOLD
				cycle
			endif
			! Esc. subsup. zerou e balanco ainda nao fecha:
			WX=DINT ! Valor de WX é igual ao residuo					
			DINT=0.0
!write(*,*) '7'

			! Diminui esc. superficial:
			DSUP=DSUP-(-WX)
			if (DSUP>=0.0) then
				WX=WXOLD
				cycle
			endif
			DSUP=0.0
!write(*,*) '8'

			! Esc. sup. zerou e balanco ainda nao fecha:
			write(*,*) 'Solo secou'
			WX=0.001
			exit
		else
			exit
		endif
	
!		write(*,*) '9'

	ENDDO
!		write(*,*) 'Saiu'

	RETURN
	END