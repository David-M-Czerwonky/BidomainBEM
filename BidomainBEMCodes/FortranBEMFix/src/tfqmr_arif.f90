!*********************************************************************
!  this subroutine solves the system ax=b using transpose free quasi-
!  minimal residual (tfqmr) algorithm. reference:
!    siam j. sci. compt. vol.14, no.2, pp. 470-482, march 93
!       by roland. w. freund
!
!  the program terminates when the required precision is reached.
!  if the required precision is not established only n iterations
!  are carried out.     a.a.ergin may 1995
!
! Modified by a.e. yilmaz 2003
! needs mat_vec_mult and initialize_r0 routine
! that computes an initial random vector r0_initial
!*********************************************************************
subroutine ztfqmr(ntotal,b,x,errf,iter,triamesh)
  use solver
  implicit none
  type(trianglemesh) :: triamesh
  integer(kind=spint),intent(in)::ntotal
  real(kind=dp),dimension(1:ntotal)::x,bb,b
  real(kind=dp)::errf,rerr
  integer(kind=spint) ::iter,itmax,it
  real(kind=dp),dimension(1:ntotal)::w,yo,ayo,ye,aye,r,d,v
  real(kind=dp)::ta,we,cm
  real(kind=dp)::etha,rho,amgis,ahpla,dum,beta
  real(kind=dp)::bmag
  integer(kind=spint) :: ndim
  integer(kind=spint) :: mpierr

  ndim=ntotal
  itmax=iter
  if (iter.eq.0) itmax=ntotal
  if(diag_precon) call precon(bb(1:ndim),ndim)

bb(1:ndim) = b(1:ndim)
  !
  !  set initial values
  !
  r(1:ndim)=0.0_dp
  d(1:ndim)=0.0_dp
  call matvec(x(1:ndim), r(1:ndim),triamesh)
  if(diag_precon) call precon(r(1:ndim),ndim)
  r(1:ndim)=bb(1:ndim)-r(1:ndim) !residual from the initial guess
  w(1:ndim)=r(1:ndim)
  yo(1:ndim)=r(1:ndim)
  call matvec(yo(1:ndim), ayo(1:ndim),triamesh)
  if(diag_precon) call precon(ayo(1:ndim),ndim)
  v(1:ndim)=ayo(1:ndim)
  we =0.0_dp
  etha=0.0_dp

  ta=sqrt(dot_product(r(1:ndim),r(1:ndim)))
  rho=dot_product(r0_initial(1:ndim),r(1:ndim))
  bmag=sqrt(dot_product(bb(1:ndim),bb(1:ndim)))
  rerr=ta/bmag

  iters: do it=1,itmax
     amgis=dot_product(r0_initial(1:ndim),v(1:ndim))
     ahpla=rho/amgis
     ye(1:ndim)=yo(1:ndim)-ahpla*v(1:ndim)
     call matvec(ye(1:ndim), aye(1:ndim),triamesh)
     if(diag_precon) call precon(aye(1:ndim),ndim)
   !  start odd (2n-1) m loop
     d(1:ndim)=yo(1:ndim)+(we*we*etha/ahpla)*d(1:ndim)
     w(1:ndim)=w(1:ndim)-ahpla*ayo(1:ndim)
     we=sqrt(abs(dot_product(w(1:ndim),w(1:ndim))))/ta
     cm=1.0d0/sqrt(1.0d0+we*we)
     ta=ta*we*cm
     etha=ahpla*cm*cm
     x(1:ndim)=x(1:ndim)+etha*d(1:ndim)
!  check if the result has converged.
!a        if (err*bmag .gt. ta*sqrt(2.*it)) then
!
!  start even (2n)  m loop
     d(1:ndim)=ye(1:ndim)+(we*we*etha/ahpla)*d(1:ndim)
     w(1:ndim)=w(1:ndim)-ahpla*aye(1:ndim)
     we=sqrt(abs(dot_product(w(1:ndim),w(1:ndim))))/ta
     cm=1.0d0/sqrt(1.0d0+we*we)
     ta=ta*we*cm
     etha=ahpla*cm*cm
     x(1:ndim)=x(1:ndim)+etha*d(1:ndim)
   !  check if the result has converged.
     if (mod(it,10)==0 .or. rerr<5.0_dp*errf) then
        call matvec(x(1:ndim), r(1:ndim),triamesh)
        if(diag_precon) call precon(r(1:ndim),ndim)
        r(1:ndim)=bb(1:ndim) -r(1:ndim)
        rerr=sqrt(abs(dot_product(r(1:ndim),r(1:ndim))))/bmag
        !...............................................................
!        print*,'#ofiter,error:',it,rerr


        if (errf > rerr) then

           errf=rerr

           iter=it
           return
        endif
     end if
     !  make preparations for next iteration
     dum=dot_product(r0_initial(1:ndim),w(1:ndim))
     beta=dum/rho
     rho=dum
     yo(1:ndim)=w(1:ndim)+beta*ye(1:ndim)
     call matvec(yo(1:ndim), ayo(1:ndim),triamesh)
     if(diag_precon) call precon(ayo(1:ndim),ndim)
     !MAGIC
     v(1:ndim)=ayo(1:ndim)+beta*(aye(1:ndim)+beta*v(1:ndim) )
  enddo iters
  !
  call matvec(x(1:ndim), r(1:ndim),triamesh)
  if(diag_precon) call precon(r(1:ndim),ndim)
  !MAGIC
  r(1:ndim)=bb(1:ndim)-r(1:ndim)
  errf=sqrt(abs(dot_product(r(1:ndim),r(1:ndim))))/bmag
  iter=itmax
  return
end subroutine ztfqmr
