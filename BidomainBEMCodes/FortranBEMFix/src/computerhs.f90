program tester
  use globalconstants
  use inputoutput
  use solver
type(coilarr) :: coil
type(trianglemesh) :: head
real(kind=dp),dimension(:),allocatable :: rhs,xval
real(kind=dp) :: del,relres,st,en
integer(kind=spint) :: iter
complex(kind=dp) :: ima
data ima/(0.0d0,1.0d0)/
!!!!!!!!!!!!!!!load mesh


coil%npts=193536
!head%nt=1202
!head%np=603
head%nt=62706
head%np=31355
call CPU_TIME(st)
 call initmeshfromfile('/home/ljg/Downloads/sphere.txt',head,head%nt*4+3*head%np)
del=sqrt(2.0_dp*maxval(head%vol))/2.0_dp
 call initcoilfromfile('/home/ljg/Downloads/fig8.txt',coil,coil%npts)
call CPU_TIME(en)
write(*,*) 'read mesh time=',en-st
call CPU_TIME(st)
 allocate(rhs(head%nt));
 call computerhs(coil,head,rhs)
call CPU_TIME(en)
write(*,*) 'compute rhs time=',en-st
call CPU_TIME(st)
 call initializematvec(head,del)
call CPU_TIME(en)
write(*,*) 'setup mattime=',en-st
call CPU_TIME(st)
 allocate(xval(head%nt))
xval=rhs
relres=1.0d-7
iter=100

 call ztfqmr(head%nt,rhs,xval,relres,iter,head)
call CPU_TIME(en)
write(*,*) 'solve time=',en-st
 call rwritetofile('/home/ljg/Downloads/FortranBEM/Ex',dble(rhs),head%nt)
 call rwritetofile('/home/ljg/Downloads/FortranBEM/Ey',dble(xval),head%nt)









end program tester
