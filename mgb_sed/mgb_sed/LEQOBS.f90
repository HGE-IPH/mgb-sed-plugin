	SUBROUTINE LEQOBS
	!SUBROTINA QUE L� DADOS OBSERVADOS DE VAZ�O
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER I,J,K,L

	OPEN(FILOBS,FILE='.\input\'//ARQOBS,STATUS='OLD')
	READ(FILOBS,701)(CABE(K),K=1,NOBS)
	DO IT=1,NT
!		READ(FILOBS,702)CABE(1),(QOBS(K,IT),K=1,NOBS)        
		READ(FILOBS,*)I,J,L,(QOBS(K,IT),K=1,NOBS)
	ENDDO
	CLOSE (FILOBS)
701	FORMAT(<NOBS>A10)
702	FORMAT(A20,43F10.2)
	!FIM DA LEITURA DE DADOS DE VAZAO OBSERVADA --------------------!


	RETURN
	END