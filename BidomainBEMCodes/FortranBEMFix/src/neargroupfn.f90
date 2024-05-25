module neargrouping
use globalconstants
implicit none
public :: QsortC,groupcreatorMLct,groupcreatorML,countneighbors,writeneighbors
private :: Partition

contains

recursive subroutine QsortC(A,B)
  integer(kind=spint), intent(inout), dimension(:) :: A,B
  integer :: iq

  if(size(A) > 1) then
     call Partition(A,B, iq)
     call QsortC(A(:iq-1),B(:iq-1))
     call QsortC(A(iq:),B(iq:))
  endif
end subroutine QsortC

subroutine Partition(A,B, marker)
  integer(kind=spint), intent(inout), dimension(:) :: A,B
  integer, intent(out) :: marker
  integer :: i, j
  integer(kind=spint) :: temp
  integer(kind=spint) :: x,y      ! pivot point
  x = A(1)
  y = B(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if ((A(j) < x) .or. ((A(j)==x) .and. (B(j)<=y))) exit
        j = j-1
     end do
     i = i+1
     do
        if ((A(i) > x) .or. ((A(i)==x) .and. (B(i)>=y))) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
		temp = B(i)
        B(i) = B(j)
        B(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

recursive subroutine groupcreatorMLct(rs,grs,ns,ro,gro,no,del,ct)
implicit none
  integer(kind=spint),intent(in) :: ns,no
integer(kind=spint),intent(inout) :: ct
real(kind=dp),intent(in) :: del
real(kind=dp),dimension(3,ns),intent(in) :: rs
real(kind=dp),dimension(3,no),intent(in) :: ro
integer(kind=spint),dimension(ns),intent(in) :: grs
integer(kind=spint),dimension(no),intent(in) :: gro

!internals
integer(kind=spint) :: iloo,cts,cto,sloo,oloo
integer(kind=spint),dimension(:),allocatable :: keeps,keepo

real(kind=dp),dimension(3) :: mindd,maxdd,deli
real(kind=dp) :: mdeli

real(kind=dp),dimension(3,8) :: lbous,ubous


allocate(keeps(ns))
allocate(keepo(no))

do iloo=1,3
mindd(iloo)=minval(rs(iloo,:))
maxdd(iloo)=maxval(rs(iloo,:))
deli(iloo)=(maxdd(iloo)-mindd(iloo))
end do
mdeli=maxval(deli)


if (mdeli>1.5*del) then
lbous=0.0_dp
ubous=0.0_dp
!create groups
do iloo=1,8
lbous(1,iloo)=mindd(1)-1.0D-12
lbous(2,iloo)=mindd(2)-1.0D-12
lbous(3,iloo)=mindd(3)-1.0D-12
ubous(:,iloo)=mindd(:)+deli(:)/2
end do
lbous(3,5:8)=ubous(3,1)
ubous(3,5:8)=maxdd(3)+1.0D-12
lbous(2,(/3,4,7,8/))=ubous(2,1)
ubous(2,(/3,4,7,8/))=maxdd(2)+1.0D-12
lbous(1,(/2,4,6,8/))=ubous(1,1)
ubous(1,(/2,4,6,8/))=maxdd(1)+1.0D-12

!find group elements
do iloo=1,8
cto=0
cts=0
do sloo=1,ns
if   ( (rs(1,sloo) .le. ubous(1,iloo)) .AND. (rs(1,sloo) .gt. lbous(1,iloo)) &
 .AND. (rs(2,sloo) .le. ubous(2,iloo)) .AND. (rs(2,sloo) .gt. lbous(2,iloo)) &
 .AND. (rs(3,sloo) .le. ubous(3,iloo)) .AND. (rs(3,sloo) .gt. lbous(3,iloo))) then
cts=cts+1
keeps(cts)=sloo
end if
end do
do oloo=1,no
if     ((ro(1,oloo) .ge. (lbous(1,iloo)-mdeli/2)) .AND. (ro(1,oloo) .le. ubous(1,iloo)+mdeli/2 &
) .AND. (ro(2,oloo) .ge. (lbous(2,iloo)-mdeli/2)) .AND. (ro(2,oloo) .le. ubous(2,iloo)+mdeli/2 &
) .AND. (ro(3,oloo) .ge. (lbous(3,iloo)-mdeli/2)) .AND. (ro(3,oloo) .le. ubous(3,iloo)+mdeli/2)) then
cto=cto+1
keepo(cto)=oloo
end if
end do


 if ((cto .eq. 0) .OR. (cts .eq. 0)) then
 else
call groupcreatorMLct(rs(:,keeps(1:cts)),grs(keeps(1:cts)),cts, &
    ro(:,keepo(1:cto)),gro(keepo(1:cto)),cto,del,ct)
end if
end do

else
do oloo=1,no
do sloo=1,ns
if (sqrt((ro(1,oloo)-rs(1,sloo))**2.0_dp + &
         (ro(2,oloo)-rs(2,sloo))**2.0_dp + &
         (ro(3,oloo)-rs(3,sloo))**2.0_dp) < del) then
    ct=ct+1
end if
end do
end do

end if

deallocate(keepo)
deallocate(keeps)
end subroutine

recursive subroutine groupcreatorML(rs,grs,ns,ro,gro,no,del,ct,col,row,ncol)
implicit none
  integer(kind=spint),intent(in) :: ns,no,ncol
integer(kind=spint),intent(inout) :: ct
integer(kind=spint),dimension(ncol),intent(inout) :: col,row

real(kind=dp),intent(in) :: del
real(kind=dp),dimension(3,ns),intent(in) :: rs
real(kind=dp),dimension(3,no),intent(in) :: ro
integer(kind=spint),dimension(ns),intent(in) :: grs
integer(kind=spint),dimension(no),intent(in) :: gro

!internals
integer(kind=spint) :: iloo,cts,cto,sloo,oloo
integer(kind=spint),dimension(:),allocatable :: keeps,keepo

real(kind=dp),dimension(3) :: mindd,maxdd,deli
real(kind=dp) :: mdeli

real(kind=dp),dimension(3,8) :: lbous,ubous


allocate(keeps(ns))
allocate(keepo(no))

do iloo=1,3
mindd(iloo)=minval(rs(iloo,:))
maxdd(iloo)=maxval(rs(iloo,:))
deli(iloo)=(maxdd(iloo)-mindd(iloo))
end do
mdeli=maxval(deli)


if (mdeli>1.5*del) then
lbous=0.0_dp
ubous=0.0_dp
!create groups
do iloo=1,8
lbous(1,iloo)=mindd(1)-1.0D-12
lbous(2,iloo)=mindd(2)-1.0D-12
lbous(3,iloo)=mindd(3)-1.0D-12
ubous(:,iloo)=mindd(:)+deli(:)/2
end do
lbous(3,5:8)=ubous(3,1)
ubous(3,5:8)=maxdd(3)+1.0D-12
lbous(2,(/3,4,7,8/))=ubous(2,1)
ubous(2,(/3,4,7,8/))=maxdd(2)+1.0D-12
lbous(1,(/2,4,6,8/))=ubous(1,1)
ubous(1,(/2,4,6,8/))=maxdd(1)+1.0D-12

!find group elements
do iloo=1,8
cto=0
cts=0
do sloo=1,ns
if   ( (rs(1,sloo) .le. ubous(1,iloo)) .AND. (rs(1,sloo) .gt. lbous(1,iloo)) &
 .AND. (rs(2,sloo) .le. ubous(2,iloo)) .AND. (rs(2,sloo) .gt. lbous(2,iloo)) &
 .AND. (rs(3,sloo) .le. ubous(3,iloo)) .AND. (rs(3,sloo) .gt. lbous(3,iloo))) then
cts=cts+1
keeps(cts)=sloo
end if
end do
do oloo=1,no
if     ((ro(1,oloo) .ge. (lbous(1,iloo)-mdeli/2)) .AND. (ro(1,oloo) .le. ubous(1,iloo)+mdeli/2 &
) .AND. (ro(2,oloo) .ge. (lbous(2,iloo)-mdeli/2)) .AND. (ro(2,oloo) .le. ubous(2,iloo)+mdeli/2 &
) .AND. (ro(3,oloo) .ge. (lbous(3,iloo)-mdeli/2)) .AND. (ro(3,oloo) .le. ubous(3,iloo)+mdeli/2)) then
cto=cto+1
keepo(cto)=oloo
end if
end do


 if ((cto .eq. 0) .OR. (cts .eq. 0)) then
 else
call groupcreatorML(rs(:,keeps(1:cts)),grs(keeps(1:cts)),cts, &
    ro(:,keepo(1:cto)),gro(keepo(1:cto)),cto,del,ct,col,row,ncol)
end if
end do

else
do oloo=1,no
do sloo=1,ns
if (sqrt((ro(1,oloo)-rs(1,sloo))**2.0_dp + &
         (ro(2,oloo)-rs(2,sloo))**2.0_dp + &
         (ro(3,oloo)-rs(3,sloo))**2.0_dp) < del) then
    ct=ct+1
	col(ct)=grs(sloo)
	row(ct)=gro(oloo)
end if
end do
end do

end if

deallocate(keepo)
deallocate(keeps)
end subroutine

subroutine  countneighbors(t2p,nt,col,row,ct,ct2)
implicit none
integer(kind=spint),intent(in) :: nt
integer(kind=spint),intent(inout) :: ct,ct2
integer(kind=spint),dimension(3*nt),intent(out) :: col,row
integer(kind=spint),dimension(3,nt),intent(in) :: t2p
integer(kind=spint) :: iloo
do iloo=1,nt
ct=ct+1
col(ct)=t2p(1,iloo)
row(ct)=iloo
ct=ct+1
col(ct)=t2p(2,iloo)
row(ct)=iloo
ct=ct+1
col(ct)=t2p(3,iloo)
row(ct)=iloo
end do
call QsortC(col(1:ct),row(1:ct))
ct=0
ct2=1
do iloo=2,3*nt
if (col(iloo)==col(iloo-1)) then
ct2=ct2+1
else
ct=ct+ct2**2
ct2=1
end if

end do
ct2=ct+ct2**2
ct=0
end subroutine

subroutine writeneighbors(t2p,nt,col,row,colt,rowt,ncol,ct,ct2)
implicit none
integer(kind=spint),intent(in) :: nt,ncol
integer(kind=spint),intent(inout) :: ct,ct2
integer(kind=spint),dimension(3*nt),intent(in) :: col,row
integer(kind=spint),dimension(ncol),intent(out) :: colt,rowt

integer(kind=spint),dimension(3,nt),intent(in) :: t2p
integer(kind=spint) :: iloo,jloo,kloo

ct2=1
do iloo=2,3*nt
if (col(iloo)==col(iloo-1)) then
ct2=ct2+1
else
do jloo=iloo-ct2,iloo-1
do kloo=iloo-ct2,iloo-1
ct=ct+1
colt(ct)=row(jloo)
rowt(ct)=row(kloo)
end do
end do
ct2=1
end if
end do
do jloo=3*nt-ct2+1,3*nt
do kloo=3*nt-ct2+1,3*nt
ct=ct+1
colt(ct)=row(jloo)
rowt(ct)=row(kloo)
end do
end do

end subroutine




subroutine generategroup(t2p,nt,p,del,col,row,ncol)
implicit none
integer(kind=spint),intent(in) :: nt
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
integer(kind=spint),dimension(:),intent(out),allocatable :: col,row
real(kind=dp),intent(in),dimension(:,:) :: p
real(kind=dp), intent(in) :: del
integer(kind=spint),intent(out) :: ncol
!internals
real(kind=dp),dimension(:,:),allocatable :: rs
integer(kind=spint),dimension(:),allocatable :: coltt,rowtt,colt,rowt,grs
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
!write actual array
allocate(colt(ncol))
allocate(rowt(ncol))
!write neighbors
ct=0
call writeneighbors(t2p,nt,coltt,rowtt,colt,rowt,ncol,ct,ct2)
!write nearfield
call groupcreatorML(rs,grs,ns,rs,grs,ns,del,ct,colt,rowt,ncol)
!sort and remove repeats
call QsortC(colt,rowt)
ct=1
do iloo=2,ncol
if ((colt(iloo) .ne. colt(iloo-1)) .or. (rowt(iloo) .ne. rowt(iloo-1))) then
ct=ct+1
end if
end do
deallocate(coltt)
deallocate(rowtt)
allocate(col(ct))
allocate(row(ct))
ct=1
col(1)=colt(1)
row(1)=rowt(1)
do iloo=2,ncol
if ((colt(iloo) .ne. colt(iloo-1)) .or. (rowt(iloo) .ne. rowt(iloo-1))) then
ct=ct+1
col(ct)=colt(iloo)
row(ct)=rowt(iloo)
end if
end do
deallocate(rowt)
deallocate(colt)
ncol=ct
deallocate(rs)
deallocate(grs)
end subroutine


subroutine generategroupobs(t2p,nt,p,robs,nobs,del,col,row,ncol)
implicit none
integer(kind=spint),intent(in) :: nt,nobs
integer(kind=spint),intent(in),dimension(3,nt) :: t2p
real(kind=dp),intent(in),dimension(3,nobs) :: robs
integer(kind=spint),dimension(:),intent(inout),allocatable :: col,row
real(kind=dp),intent(in),dimension(:,:) :: p
real(kind=dp), intent(in) :: del
integer(kind=spint),intent(out) :: ncol
!internals
real(kind=dp),dimension(:,:),allocatable :: rs
integer(kind=spint),dimension(:),allocatable :: grs,gro
integer(kind=spint) :: ct,iloo,jloo,kloo,ns
allocate(rs(3,nt))
allocate(grs(nt))
allocate(gro(nobs))
!create groups
do iloo=1,nt
rs(:,iloo)=p(:,t2p(1,iloo))+p(:,t2p(2,iloo))+p(:,t2p(3,iloo))
end do
rs=rs/3.0_dp
ns=nt
grs=(/ (iloo, iloo = 1, ns) /)
gro=(/ (iloo, iloo = 1, nobs) /)
ncol=0
!figure out how many nearfield
call groupcreatorMLct(rs,grs,ns,robs,gro,nobs,del,ncol)
!write actual array
allocate(col(ncol))
allocate(row(ncol))
ct=0
call groupcreatorML(rs,grs,ns,robs,gro,nobs,del,ct,col,row,ncol)
!sort and remove repeats
call QsortC(col,row)

deallocate(rs)
deallocate(grs)
deallocate(gro)
end subroutine
end module neargrouping
