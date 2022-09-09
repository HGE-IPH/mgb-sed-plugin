	subroutine GradBiConj
	! Resolve sistema de equações lineares pelo método iterativo para matrizes esparsas dos
	! gradientes biconjugados. Utiliza rotinas do Numerical Recipes in Fortran
	!
	! 
	!-------------------------------------------------------------------------------
	!
	! Descrição das variáveis locais:
	!
	! itol =  criterio de parada utilizado
	! tol = erro maximo tolerado
	! iter = numero de iteracoes
	! err = erro 
	! itmax = numero máximo de iteracoes
	!--------------------------------------------------------------------------------
	! Declaração de variáveis:
	use SPAMAT
	use MAT_MOD

	implicit none

	! Variáveis locais de cálculo:

	integer:: itol,iter,itmax,i
	real(8):: tol,err

	!----------------------------------------------------------------------------

	! Inicializa variáveis:
	tol=0.000005
	itol=2
	itmax=NUM
	sa=dble(AA)
	bbb=dble(BB)
	xxi=dble(XI)

!	write(*,*) ija(1),NUM
!	call dsprsax(sa,ija,XI,bbb,NUM,dimSA)
!	do i=1,NUM
!		write(*,*) ICOL(i,1),bb(i),bbb(i)
!		read(*,*)
!	enddo
!	write(*,*) minval(sa(1:dimSA)),minval(abs(sa(1:NUM))),minloc(abs(sa(1:NUM))),dimSA
!	read(*,*) 
!	write(*,*) XI
!	read(*,*)
	!call linbcg(NUM,BB,XI,itol,tol,itmax,iter,err,sa,ija,dimSA)
	call linbcg(NUM,bbb,xxi,itol,tol,itmax,iter,err,dimSA)
	write (*,*) ' iter=',iter,' err=',err
!	read(*,*) 
!	write(*,*) XI
!	read(*,*)

!	call dsprsax(sa,ija,xxi,bbb,NUM,dimSA)
!	do i=1,NUM
!		write(*,*) ICOL(i,1),bb(i),bbb(i)
!		read(*,*)
!	enddo

	XI=real(xxi)

	endsubroutine