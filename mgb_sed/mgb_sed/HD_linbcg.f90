      SUBROUTINE linbcg(N,b,x,itol,tol,itmax,iter,err,NMAX)

	  IMPLICIT NONE
      INTEGER iter,itmax,itol,n,NMAX
	  !DOUBLE PRECISION err,tol,b(N),x(N),EPS,SA(NMAX)
      REAL(8) err,tol,b(N),x(N),EPS

      PARAMETER (EPS=1.d-14)
!CU    USES atimes,asolve,snrm
      INTEGER j,i
	  !DOUBLE PRECISION ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm
	  !DOUBLE PRECISION znrm,p(NMAX),pp(NMAX),r(NMAX),rr(NMAX),z(NMAX),zz(NMAX),snrm
	  REAL(8) ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm
	  REAL(8) znrm,p(N),pp(N),r(N),rr(N),z(N),zz(N),snrm
      iter=0
 
 	  write(*,*) 'inicio' !     
!	  write(*,*) b !
!	  read(*,*)  !
 ! write(*,*) x !
!	  read(*,*) !

     call atimes(n,x,r,0,NMAX)
!     write(*,*) 'r:',r !
!	  read(*,*) !
	  do 11 j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
11    continue
	  call atimes(n,r,rr,0,NMAX)
!      write(*,*) 'rr:',rr !
!	  read(*,*) !
	  if(itol.eq.1) then
        bnrm=snrm(n,b,itol,NMAX)
        call asolve(n,r,z,0,NMAX)  !PRIMEIRA

      else if (itol.eq.2) then
        call asolve(n,b,z,0,NMAX)
        bnrm=snrm(n,z,itol,NMAX)
        call asolve(n,r,z,0,NMAX)
      else if (itol.eq.3.or.itol.eq.4) then
        call asolve(n,b,z,0,NMAX)
        bnrm=snrm(n,z,itol,NMAX)
        call asolve(n,r,z,0,NMAX)
        znrm=snrm(n,z,itol,NMAX)
      elseif (itol.eq.5) then ! Criterio de convergencia
		call asolve(n,r,z,0,NMAX)
	  else
        pause 'illegal itol in linbcg'
      endif
!	  write(*,*) 'z:',z !
!	   write(*,*) 'znrm',znrm !
!	  read(*,*) !
		write(*,*) 'inicio2' 
100   if (iter.le.itmax) then
!      write(*,*) x
!		read(*,*)
		iter=iter+1
        call asolve(n,rr,zz,1,NMAX)	  !SEGUNDA
        bknum=0.d0
        do 12 j=1,n
          bknum=bknum+z(j)*rr(j)
12      continue
        if(iter.eq.1) then
          do 13 j=1,n
            p(j)=z(j)

            pp(j)=zz(j)
13        continue
        else
          bk=bknum/bkden
          do 14 j=1,n
            p(j)=bk*p(j)+z(j)
            pp(j)=bk*pp(j)+zz(j)
14        continue
        endif
        bkden=bknum
        call atimes(n,p,z,0,NMAX)
        akden=0.d0
        do 15 j=1,n
          akden=akden+z(j)*pp(j)
15      continue
        ak=bknum/akden
        call atimes(n,pp,zz,1,NMAX)
 !      write(*,*) 'akden,bknum,ak:',akden,bknum,ak
!		read(*,*)
!		write(*,*) 'p:',p
!		read(*,*) 
		do 16 j=1,n
          x(j)=x(j)+ak*p(j)
		  
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
16      continue
!       write(*,*) 'iteracao e x depois:',iter,x
!		read(*,*)
		
		call asolve(n,r,z,0,NMAX)
        if(itol.eq.1)then

          err=snrm(n,r,itol,NMAX)/bnrm
        else if(itol.eq.2)then
          err=snrm(n,z,itol,NMAX)/bnrm
        else if(itol.eq.3.or.itol.eq.4)then
          zm1nrm=znrm
          znrm=snrm(n,z,itol,NMAX)
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then
            dxnrm=abs(ak)*snrm(n,p,itol,NMAX)
            err=znrm/abs(zm1nrm-znrm)*dxnrm
		  else
            err=znrm/bnrm
            goto 100
          endif
          xnrm=snrm(n,x,itol,NMAX)
          if(err.le.0.5d0*xnrm) then
            err=err/xnrm

          else
            err=znrm/bnrm
            goto 100
          endif
        elseif(itol.eq.5) then
			if (iter==1) then
				err=tol+1.0
			else
				err=0.0
				do i=1,n
					if (abs(ak*p(i))>err) err=abs(ak*p(i))
				enddo
			endif
		endif
 !       write (*,*) ' iter=',iter,' err=',err
      if(err.gt.tol) goto 100
      endif
      return
      END