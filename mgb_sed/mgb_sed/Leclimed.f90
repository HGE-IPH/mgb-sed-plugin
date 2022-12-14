	SUBROUTINE LECLIMED
	!SUBROTINA QUE L? OS DADOS MEDIOS MENSAIS DE VARIAVEIS CLIMATICAS
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER IPCL,K

	!VARIAVEIS ALFANUMERICAS
	CHARACTER (20) ANOME
	CHARACTER (10) AMES(12)

	OPEN(FILMED,FILE='.\input\'//ACLIMED,STATUS='OLD')

	!LE COORDENADAS DOS POSTOS
	DO IPCL=1,NCLI
!		READ(FILMED,710)XYC(IPCL,1),XYC(IPCL,2)
		READ(FILMED,*)XYC(IPCL,1),XYC(IPCL,2) ! RP formato livre
	ENDDO

	!temperatura media (graus C)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
	DO IPCL=1,NCLI        
		READ(FILMED,713)ANOME,(TAMM(IPCL,K),K=1,12)        
	ENDDO
	!umidade relativa (%)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(URMM(IPCL,K),K=1,12)
	ENDDO
	!insola??o (horas/dia)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
 	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(SOLMM(IPCL,K),K=1,12)
        !write(*,*)(SOLMM(IPCL,K),K=1,12)
	ENDDO
	!velocidade do vento (m/s)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
  	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(VVMM(IPCL,K),K=1,12)
        
	ENDDO
	!pressao atmosferica (kPa)
	READ(FILMED,711)ANOME
	READ(FILMED,712)ANOME,(AMES(K),K=1,12)
  	DO IPCL=1,NCLI
		READ(FILMED,713)ANOME,(PAMM(IPCL,K),K=1,12)        
	ENDDO

	CLOSE(FILMED)

710	FORMAT(2F10.1)
711	FORMAT(A20)
712	FORMAT(A20,12A10)
713	FORMAT(A20,12F10.2)

	RETURN
	END