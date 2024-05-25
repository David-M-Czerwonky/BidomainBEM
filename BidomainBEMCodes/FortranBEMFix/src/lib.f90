subroutine computerhsmatlab(t2p,nt,p,np,reg,rs,js,nc,iprecision,rhs)
  use fmmutils
  use globalconstants
  implicit none

  real(kind=dp),dimension(nt),intent(out) :: rhs
 integer(kind=spint),intent(in) :: nt,np,nc
 integer(kind=spint),dimension(3,nt),intent(in) :: t2p
 real(kind=dp),dimension(2),intent(in) :: iprecision
 real(kind=dp),dimension(3,np),intent(in) :: p
 real(kind=dp),dimension(3,nc),intent(in) :: rs
 real(kind=dp),dimension(nt),intent(in) :: reg
 real(kind=dp),dimension(3,nc),intent(in) :: js

 type(coilarr) :: coil
 type(trianglemesh) :: head
 real(kind=dp),dimension(3) :: v1,v2
 real(kind=dp) :: del,relres,t2,t1
 integer(kind=spint) :: iter,i


 iprecEp=iprecision(1)
 iprecmatvec=iprecision(2)

 head%nt=nt
 head%np=np
 coil%npts=nc
 allocate(head%p(3,head%np))
 head%p=p
 allocate(head%t2p(3,head%nt))
 head%t2p=t2p
 allocate(head%reg(head%nt))
 head%reg=reg/(4.0_dp*pi)
 allocate(head%nhat(3,head%nt))
 allocate(head%vol(head%nt))
 allocate(coil%rs(3,coil%npts))
 allocate(coil%js(3,coil%npts))
 coil%rs=rs
 coil%js=js
 do i=1,head%nt
 v1=head%p(:,head%t2p(1,i))-head%p(:,head%t2p(3,i))
 v2=head%p(:,head%t2p(2,i))-head%p(:,head%t2p(3,i))
 head%nhat(1,i)=v1(2)*v2(3)-v1(3)*v2(2)
 head%nhat(2,i)=v1(3)*v2(1)-v1(1)*v2(3)
 head%nhat(3,i)=v1(1)*v2(2)-v1(2)*v2(1)
 head%vol(i)=(head%nhat(1,i)**2.0_dp+head%nhat(2,i)**2.0_dp+head%nhat(3,i)**2.0_dp)**0.5_dp
 head%nhat(:,i)=head%nhat(:,i)/head%vol(i)
 head%vol(i)=head%vol(i)*0.5_dp;
 end do
 call computerhs(coil,head,rhs)
end subroutine


subroutine computerhsmatlabks(t2p,nt,p,np,reg,rs,js,nc,iprecision,rhs)
  use fmmutils
  use globalconstants
  implicit none

  real(kind=dp),dimension(nt),intent(out) :: rhs
 integer(kind=spint),intent(in) :: nt,np,nc
 integer(kind=spint),dimension(3,nt),intent(in) :: t2p
 real(kind=dp),dimension(2),intent(in) :: iprecision
 real(kind=dp),dimension(3,np),intent(in) :: p
 real(kind=dp),dimension(3,nc),intent(in) :: rs
 real(kind=dp),dimension(nt),intent(in) :: reg
 real(kind=dp),dimension(3,nc),intent(in) :: js

 type(coilarr) :: coil
 type(trianglemesh) :: head
 real(kind=dp),dimension(3) :: v1,v2
 real(kind=dp) :: del,relres,t2,t1
 integer(kind=spint) :: iter,i


 iprecEp=iprecision(1)
 iprecmatvec=iprecision(2)

 head%nt=nt
 head%np=np
 coil%npts=nc
 allocate(head%p(3,head%np))
 head%p=p
 allocate(head%t2p(3,head%nt))
 head%t2p=t2p
 allocate(head%reg(head%nt))
 head%reg=reg/(4.0_dp*pi)
 allocate(head%nhat(3,head%nt))
 allocate(head%vol(head%nt))
 allocate(coil%rs(3,coil%npts))
 allocate(coil%js(3,coil%npts))
 coil%rs=rs
 coil%js=js
 do i=1,head%nt
 v1=head%p(:,head%t2p(1,i))-head%p(:,head%t2p(3,i))
 v2=head%p(:,head%t2p(2,i))-head%p(:,head%t2p(3,i))
 head%nhat(1,i)=v1(2)*v2(3)-v1(3)*v2(2)
 head%nhat(2,i)=v1(3)*v2(1)-v1(1)*v2(3)
 head%nhat(3,i)=v1(1)*v2(2)-v1(2)*v2(1)
 head%vol(i)=(head%nhat(1,i)**2.0_dp+head%nhat(2,i)**2.0_dp+head%nhat(3,i)**2.0_dp)**0.5_dp
 head%nhat(:,i)=head%nhat(:,i)/head%vol(i)
 head%vol(i)=head%vol(i)*0.5_dp;
 end do
 call computerhsks(coil,head,rhs)
end subroutine

subroutine fmmchargematlab(chspace,triafl,nfmm,fld,iprecision)
    use fmmutils
    use globalconstants
implicit none
real(kind=dp),dimension(3,nfmm),intent(in) :: triafl
real(kind=dp),dimension(nfmm),intent(in) :: chspace
real(kind=dp),dimension(3,nfmm),intent(out) :: fld
integer(kind=spint),intent(in) :: nfmm
real(kind=dp),intent(in),dimension(2) :: iprecision

 iprecEp=iprecision(1)
 iprecmatvec=iprecision(2)
call fmmcharge(chspace,triafl,nfmm,fld)
end subroutine


subroutine generategroupct(t2p,nt,p,np,del,ncol)
  use globalconstants
  use neargrouping
implicit none
integer(kind=spint),intent(in) :: nt,np
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp), intent(in) :: del
integer(kind=spint),intent(out) :: ncol
!internals
real(kind=dp),dimension(:,:),allocatable :: rs
integer(kind=spint),dimension(:),allocatable :: coltt,rowtt,grs
integer(kind=spint) :: ct,ct2,iloo,jloo,kloo,ns
allocate(rs(3,nt))
allocate(grs(nt))
!create groups
do iloo=1,nt
rs(:,iloo)=p(:,t2p(1,iloo))+p(:,t2p(2,iloo))+p(:,t2p(3,iloo))
end do
rs=rs/3.0_dp
ns=nt
grs=(/ (iloo, iloo = 1, ns) /)
ncol=3*nt
allocate(coltt(ncol))
allocate(rowtt(ncol))
ct=0
ct2=0
call countneighbors(t2p,nt,coltt,rowtt,ct,ct2)
!figure out how many nearfield
call groupcreatorMLct(rs,grs,ns,rs,grs,ns,del,ct)
	  ncol=ct+ct2

deallocate(coltt)
deallocate(rowtt)

end subroutine




subroutine generategroupmatlab(t2p,nt,p,np,del,colt,rowt,ncol)
use globalconstants
use neargrouping
implicit none
integer(kind=spint),intent(in) :: nt,np
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
integer(kind=spint),dimension(ncol),intent(out) :: colt,rowt
real(kind=dp),intent(in),dimension(3,np) :: p
real(kind=dp), intent(in) :: del
integer(kind=spint),intent(in) :: ncol

!internals
real(kind=dp),dimension(:,:),allocatable :: rs
integer(kind=spint),dimension(:),allocatable :: coltt,rowtt,grs
integer(kind=spint) :: ct,ct2,iloo,jloo,kloo,ns
allocate(rs(3,nt))
allocate(grs(nt))
do iloo=1,nt
rs(:,iloo)=p(:,t2p(1,iloo))+p(:,t2p(2,iloo))+p(:,t2p(3,iloo))
end do
rs=rs/3.0_dp
ns=nt
grs=(/ (iloo, iloo = 1, ns) /)
allocate(coltt(3*nt))
allocate(rowtt(3*nt))
ct2=0
ct=0
call countneighbors(t2p,nt,coltt,rowtt,ct,ct2)
ct=0
call writeneighbors(t2p,nt,coltt,rowtt,colt,rowt,ncol,ct,ct2)
!write nearfield
call groupcreatorML(rs,grs,ns,rs,grs,ns,del,ct,colt,rowt,ncol)
call QsortC(colt,rowt)
deallocate(rowtt)
deallocate(coltt)
allocate(rowtt(ncol))
allocate(coltt(ncol))
coltt=0
rowtt=0
ct=1
coltt(1)=colt(1)
rowtt(1)=rowt(1)
do iloo=2,ncol
if ((colt(iloo) .ne. colt(iloo-1)) .or. (rowt(iloo) .ne. rowt(iloo-1))) then
ct=ct+1
coltt(ct)=colt(iloo)
rowtt(ct)=rowt(iloo)
end if
end do

colt=coltt
rowt=rowtt

deallocate(rowtt)
deallocate(coltt)


deallocate(rs)
deallocate(grs)


end subroutine




subroutine createbemdatastructmatlab(t2p,nt,p,np,area,normal,row,col, &
  ncols,potout,triafl,ncolloc,colch,rowch,volch)

use globalconstants
use nearfielding

implicit none
  ! inputs:
integer(kind=spint),intent(in) :: nt,np,ncols
real(kind=dp),intent(in),dimension(3,np) :: p
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
real(kind=dp),intent(in),dimension(3,nt) :: normal
real(kind=dp),intent(in),dimension(nt) :: area
integer(kind=spint),intent(in),dimension(ncols) :: col,row
integer(kind=spint),intent(in),dimension(2) :: ncolloc

!outputs


integer(kind=spint),intent(out),dimension(ncolloc(1)*nt) :: rowch,colch
real(kind=dp),intent(out),dimension(3,ncolloc(1)*nt) :: triafl
real(kind=dp),intent(out),dimension(ncolloc(1)*nt) :: volch
real(kind=dp),intent(out),dimension(ncols) :: potout


nquad=ncolloc(1)
nquadnear=ncolloc(2)

call createbemdatastruct(t2p,nt,p,np,area,normal,row,col,ncols,potout, &
triafl,nquad,colch,rowch,volch)
end subroutine


subroutine BEMsolver(t2p,nt,p,np,reg,rs,js,nc,xval,ncolloc,iprecision,delnear,iter,relres)
 use globalconstants

  use solver
      implicit none
integer(kind=spint),intent(in) :: nt,np,nc
integer(kind=spint),dimension(3,nt),intent(in) :: t2p
integer(kind=spint),dimension(2),intent(in) :: ncolloc
real(kind=dp),dimension(2),intent(in) :: iprecision

real(kind=dp),intent(in) :: delnear
real(kind=dp),dimension(3,np),intent(in) :: p
real(kind=dp),dimension(3,nc),intent(in) :: rs
real(kind=dp),dimension(nt),intent(in) :: reg
real(kind=dp),dimension(3,nc),intent(in) :: js
real(kind=dp),dimension(nt),intent(out) :: xval

type(coilarr) :: coil
type(trianglemesh) :: head
real(kind=dp),dimension(3) :: v1,v2
real(kind=dp),dimension(nt) :: rhs
real(kind=dp) :: del,relres,t2,t1
integer(kind=spint) :: iter,i

nquad=ncolloc(1)
nquadnear=ncolloc(2)
iprecEp=iprecision(1)
iprecmatvec=iprecision(2)

call CPU_TIME(t1)
head%nt=nt
head%np=np
coil%npts=nc
allocate(head%p(3,head%np))
head%p=p
allocate(head%t2p(3,head%nt))
head%t2p=t2p
allocate(head%reg(head%nt))
head%reg=reg/(4.0_dp*pi)
allocate(head%nhat(3,head%nt))
allocate(head%vol(head%nt))
do i=1,head%nt
v1=head%p(:,head%t2p(1,i))-head%p(:,head%t2p(3,i))
v2=head%p(:,head%t2p(2,i))-head%p(:,head%t2p(3,i))
head%nhat(1,i)=v1(2)*v2(3)-v1(3)*v2(2)
head%nhat(2,i)=v1(3)*v2(1)-v1(1)*v2(3)
head%nhat(3,i)=v1(1)*v2(2)-v1(2)*v2(1)
head%vol(i)=(head%nhat(1,i)**2.0_dp+head%nhat(2,i)**2.0_dp+head%nhat(3,i)**2.0_dp)**0.5_dp
head%nhat(:,i)=head%nhat(:,i)/head%vol(i)
head%vol(i)=head%vol(i)*0.5_dp;
end do
allocate(coil%rs(3,coil%npts))
allocate(coil%js(3,coil%npts))
coil%rs=rs
coil%js=js
call CPU_TIME(t2)
!print*,'time to read data=', t2-t1
call CPU_TIME(t1)
call computerhs(coil,head,rhs)
call CPU_TIME(t2)
!print*,'time to make rhs data=', t2-t1
call CPU_TIME(t1)
del=sqrt(2.0_dp*maxval(head%vol))*delnear
call initializematvec(head,del)
call CPU_TIME(t2)
!call rwritetofile('./rownear',dble(nearmat%row),nearmat%nnz)
!call rwritetofile('./dnear',dble(nearmat%val),nearmat%nnz)


!print*,'time to make setup data=', t2-t1
xval=rhs
relres=1.0d-6
iter=100
call CPU_TIME(t1)
call ztfqmr(head%nt,rhs,xval,relres,iter,head)
call CPU_TIME(t2)
!print*,'time to solve=', t2-t1
call destroymatvec()



deallocate(head%p)
deallocate(head%t2p)
deallocate(head%reg)
deallocate(head%nhat)
deallocate(head%vol)
deallocate(coil%rs)
deallocate(coil%js)

end subroutine



subroutine BEMsolverks(t2p,nt,p,np,reg,rs,js,nc,xval,ncolloc,iprecision,delnear,iter,relres)
 use globalconstants

  use solver
      implicit none
integer(kind=spint),intent(in) :: nt,np,nc
integer(kind=spint),dimension(3,nt),intent(in) :: t2p
integer(kind=spint),dimension(2),intent(in) :: ncolloc
real(kind=dp),dimension(2),intent(in) :: iprecision

real(kind=dp),intent(in) :: delnear
real(kind=dp),dimension(3,np),intent(in) :: p
real(kind=dp),dimension(3,nc),intent(in) :: rs
real(kind=dp),dimension(nt),intent(in) :: reg
real(kind=dp),dimension(3,nc),intent(in) :: js
real(kind=dp),dimension(nt),intent(out) :: xval
integer(kind=spint),intent(out) :: iter
real(kind=dp),intent(out) :: relres
type(coilarr) :: coil
type(trianglemesh) :: head
real(kind=dp),dimension(3) :: v1,v2
real(kind=dp),dimension(nt) :: rhs
real(kind=dp) :: del,t2,t1
integer(kind=spint) :: i

nquad=ncolloc(1)
nquadnear=ncolloc(2)
iprecEp=iprecision(1)
iprecmatvec=iprecision(2)

call CPU_TIME(t1)
head%nt=nt
head%np=np
coil%npts=nc
allocate(head%p(3,head%np))
head%p=p
allocate(head%t2p(3,head%nt))
head%t2p=t2p
allocate(head%reg(head%nt))
head%reg=reg/(4.0_dp*pi)
allocate(head%nhat(3,head%nt))
allocate(head%vol(head%nt))
do i=1,head%nt
v1=head%p(:,head%t2p(1,i))-head%p(:,head%t2p(3,i))
v2=head%p(:,head%t2p(2,i))-head%p(:,head%t2p(3,i))
head%nhat(1,i)=v1(2)*v2(3)-v1(3)*v2(2)
head%nhat(2,i)=v1(3)*v2(1)-v1(1)*v2(3)
head%nhat(3,i)=v1(1)*v2(2)-v1(2)*v2(1)
head%vol(i)=(head%nhat(1,i)**2.0_dp+head%nhat(2,i)**2.0_dp+head%nhat(3,i)**2.0_dp)**0.5_dp
head%nhat(:,i)=head%nhat(:,i)/head%vol(i)
head%vol(i)=head%vol(i)*0.5_dp;
end do
allocate(coil%rs(3,coil%npts))
allocate(coil%js(3,coil%npts))
coil%rs=rs
coil%js=js
call CPU_TIME(t2)
!print*,'time to read data=', t2-t1
call CPU_TIME(t1)
!call computerhsks(coil,head,rhs)
call CPU_TIME(t2)
!print*,'time to make rhs data=', t2-t1
call CPU_TIME(t1)
del=sqrt(2.0_dp*maxval(head%vol))*delnear
!call initializematvec(head,del)
call CPU_TIME(t2)
!call rwritetofile('./rownear',dble(nearmat%row),nearmat%nnz)
!call rwritetofile('./dnear',dble(nearmat%val),nearmat%nnz)


!print*,'time to make setup data=', t2-t1
xval=rhs
relres=1.0d-6
!iter=100
call CPU_TIME(t1)
!call ztfqmr(head%nt,rhs,xval,relres,iter,head)
call CPU_TIME(t2)
!print*,'time to solve=', t2-t1
!call destroymatvec()



deallocate(head%p)
deallocate(head%t2p)
deallocate(head%reg)
deallocate(head%nhat)
deallocate(head%vol)
deallocate(coil%rs)
deallocate(coil%js)

end subroutine


subroutine evaluateEfield(t2p,nt,p,np,reg,rs,js,nc,xval,robs,nobs,Etotal,ncolloc,iprecision,delnear)
 use globalconstants
  use solver
      implicit none
integer(kind=spint),intent(in) :: nt,np,nc,nobs
integer(kind=spint),dimension(3,nt),intent(in) :: t2p
integer(kind=spint),dimension(2),intent(in) :: ncolloc
real(kind=dp),dimension(2),intent(in) :: iprecision
real(kind=dp),intent(in) :: delnear
real(kind=dp),dimension(3,np),intent(in) :: p
real(kind=dp),dimension(3,nc),intent(in) :: rs,js
real(kind=dp),dimension(nt),intent(in) :: reg
real(kind=dp),dimension(nt),intent(in) :: xval
real(kind=dp),dimension(3,nobs),intent(in) :: robs
real(kind=dp),dimension(3,nobs),intent(inout) :: Etotal

real(kind=dp),dimension(:,:),allocatable :: triaflat
integer(kind=spint),dimension(:),allocatable :: row,col


type(coilarr) :: coil
type(trianglemesh) :: head
real(kind=dp),dimension(3) :: v1,v2
real(kind=dp),dimension(nt) :: rhs
real(kind=dp) :: del,relres,t2,t1
integer(kind=spint) :: iter,i,ncol


nquad=ncolloc(1)
nquadnear=ncolloc(2)
iprecEp=iprecision(1)
iprecmatvec=iprecision(2)
!!!extract mesh and coil
call CPU_TIME(t1)
head%nt=nt
head%np=np
coil%npts=nc
allocate(head%p(3,head%np))
head%p=p
allocate(head%t2p(3,head%nt))
head%t2p=t2p
allocate(head%reg(head%nt))
head%reg=reg/(4.0_dp*pi)
allocate(head%nhat(3,head%nt))
allocate(head%vol(head%nt))
do i=1,head%nt
v1=head%p(:,head%t2p(1,i))-head%p(:,head%t2p(3,i))
v2=head%p(:,head%t2p(2,i))-head%p(:,head%t2p(3,i))
head%nhat(1,i)=v1(2)*v2(3)-v1(3)*v2(2)
head%nhat(2,i)=v1(3)*v2(1)-v1(1)*v2(3)
head%nhat(3,i)=v1(1)*v2(2)-v1(2)*v2(1)
head%vol(i)=(head%nhat(1,i)**2.0_dp+head%nhat(2,i)**2.0_dp+head%nhat(3,i)**2.0_dp)**0.5_dp
head%nhat(:,i)=head%nhat(:,i)/head%vol(i)
head%vol(i)=head%vol(i)/2.0_dp;
end do
allocate(coil%rs(3,coil%npts))
allocate(coil%js(3,coil%npts))
coil%rs=rs
coil%js=js

call CPU_TIME(t2)
!print*,'time to read data=', t2-t1
!!!!generate near group structure
del=sqrt(2.0_dp*maxval(head%vol))*delnear

call CPU_TIME(t1)
call generategroupobs(t2p,nt,p,robs,nobs,del,col,row,ncol)
call CPU_TIME(t2)
!print*,'time to generate near group=', t2-t1
Etotal=0.0_dp
!!!!add E primary
call computeEprimary(coil%rs,coil%js,robs,Etotal)
!!!!compute E secondary

call computeEsecondary(head%t2p,head%nt,head%p &
,head%np,head%vol,head%nhat,col,row,ncol &
,Etotal,robs,nobs,nquad,xval)



deallocate(head%p)
deallocate(head%t2p)
deallocate(head%reg)
deallocate(head%nhat)
deallocate(head%vol)
deallocate(coil%rs)
deallocate(coil%js)


end subroutine



subroutine evaluateEfieldks(t2p,nt,p,np,reg,rs,js,nc,xval,robs,nobs,Etotal,ncolloc,iprecision,delnear)
 use globalconstants
  use solver
      implicit none
integer(kind=spint),intent(in) :: nt,np,nc,nobs
integer(kind=spint),dimension(3,nt),intent(in) :: t2p
integer(kind=spint),dimension(2),intent(in) :: ncolloc
real(kind=dp),dimension(2),intent(in) :: iprecision
real(kind=dp),intent(in) :: delnear
real(kind=dp),dimension(3,np),intent(in) :: p
real(kind=dp),dimension(3,nc),intent(in) :: rs,js
real(kind=dp),dimension(nt),intent(in) :: reg
real(kind=dp),dimension(nt),intent(in) :: xval
real(kind=dp),dimension(3,nobs),intent(in) :: robs
real(kind=dp),dimension(3,nobs),intent(inout) :: Etotal

real(kind=dp),dimension(:,:),allocatable :: triaflat
integer(kind=spint),dimension(:),allocatable :: row,col


type(coilarr) :: coil
type(trianglemesh) :: head
real(kind=dp),dimension(3) :: v1,v2
real(kind=dp),dimension(nt) :: rhs
real(kind=dp) :: del,relres,t2,t1
integer(kind=spint) :: iter,i,ncol


nquad=ncolloc(1)
nquadnear=ncolloc(2)
iprecEp=iprecision(1)
iprecmatvec=iprecision(2)
!!!extract mesh and coil
call CPU_TIME(t1)
head%nt=nt
head%np=np
coil%npts=nc
allocate(head%p(3,head%np))
head%p=p
allocate(head%t2p(3,head%nt))
head%t2p=t2p
allocate(head%reg(head%nt))
head%reg=reg/(4.0_dp*pi)
allocate(head%nhat(3,head%nt))
allocate(head%vol(head%nt))
do i=1,head%nt
v1=head%p(:,head%t2p(1,i))-head%p(:,head%t2p(3,i))
v2=head%p(:,head%t2p(2,i))-head%p(:,head%t2p(3,i))
head%nhat(1,i)=v1(2)*v2(3)-v1(3)*v2(2)
head%nhat(2,i)=v1(3)*v2(1)-v1(1)*v2(3)
head%nhat(3,i)=v1(1)*v2(2)-v1(2)*v2(1)
head%vol(i)=(head%nhat(1,i)**2.0_dp+head%nhat(2,i)**2.0_dp+head%nhat(3,i)**2.0_dp)**0.5_dp
head%nhat(:,i)=head%nhat(:,i)/head%vol(i)
head%vol(i)=head%vol(i)/2.0_dp;
end do
allocate(coil%rs(3,coil%npts))
allocate(coil%js(3,coil%npts))
coil%rs=rs
coil%js=-js

call CPU_TIME(t2)
!print*,'time to read data=', t2-t1
!!!!generate near group structure
del=sqrt(2.0_dp*maxval(head%vol))*delnear

call CPU_TIME(t1)
call generategroupobs(t2p,nt,p,robs,nobs,del,col,row,ncol)
call CPU_TIME(t2)
!print*,'time to generate near group=', t2-t1
Etotal=0.0_dp
!!!!add E primary
call computeHprimary(coil%rs,coil%js,robs,Etotal,nobs,coil%npts,1d-3)
!!!!compute E secondary

call computeEsecondary(head%t2p,head%nt,head%p &
,head%np,head%vol,head%nhat,col,row,ncol &
,Etotal,robs,nobs,nquad,xval)



deallocate(head%p)
deallocate(head%t2p)
deallocate(head%reg)
deallocate(head%nhat)
deallocate(head%vol)
deallocate(coil%rs)
deallocate(coil%js)


end subroutine
