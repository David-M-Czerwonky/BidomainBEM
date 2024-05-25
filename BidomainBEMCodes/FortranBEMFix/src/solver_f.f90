
module solver
use globalconstants
use nearfielding
use neargrouping
use fmmutils
implicit none
real(kind=dp),dimension(:),allocatable :: r0_initial
integer(kind=spint) :: nfmm
logical :: diag_precon=.false.
type(sparse) :: chmat,nearmat
real(kind=dp),dimension(:,:),allocatable :: triafl
real(kind=dp),dimension(:,:),allocatable :: fld
real(kind=dp),dimension(:),allocatable :: chspace

contains
subroutine precon(x,ndim)
real(kind=dp),dimension(:),intent(inout) :: x
integer(kind=spint),intent(in) :: ndim


!!!!!!!!!!!!!!!!fill me later
end subroutine precon
subroutine matvec(x,y,triamesh)
real(kind=dp),dimension(:),intent(in) :: x
real(kind=dp),dimension(:),intent(out) :: y
type(trianglemesh) :: triamesh
call multiplicationroutine(x,y,triamesh)
end subroutine matvec

subroutine multiplicationroutine(x,y,triamesh)
real(kind=dp),dimension(:),intent(in) :: x
real(kind=dp),dimension(:),intent(out) :: y
type(trianglemesh),intent(in) :: triamesh
integer(kind=spint) :: i,j
real(kind=dp) :: t1,t2
!chmat has only one entry on each row
call CPU_TIME(t1)
chspace=0.0_dp
do i=1,chmat%nnz
chspace(chmat%row(i))=chspace(chmat%row(i))+x(chmat%col(i))*chmat%val(i)
end do

call CPU_TIME(t2)
!print*,'time transpose chmat mult=', t2-t1
call CPU_TIME(t1)
fld=0.0_dp
call fmmcharge(chspace,triafl,nfmm,fld)

call CPU_TIME(t2)
!print*,'time fmm mult=', t2-t1
y=0.0_dp

do j=1,3


call CPU_TIME(t1)
chspace(1:triamesh%nt)=0.0_dp
do i=1,chmat%nnz
chspace(chmat%col(i))=chspace(chmat%col(i))+ &
fld(j,chmat%row(i))*chmat%val(i)
end do

call CPU_TIME(t2)
!print*,'time chmat mult=', t2-t1

y=y+reshape(triamesh%nhat(j,:),(/triamesh%nt/))*(chspace(1:triamesh%nt))
end do
call CPU_TIME(t1)
do i=1,nearmat%nnz
y(nearmat%col(i))=y(nearmat%col(i))+ &
x(nearmat%row(i))*nearmat%val(i)
end do
call CPU_TIME(t2)
!print*,'time nearmat mult=', t2-t1

y=y/triamesh%vol
y=x*0.5_dp-triamesh%reg*y


end subroutine multiplicationroutine

subroutine initializematvec(triamesh,del)
type(trianglemesh),intent(in) :: triamesh
real(kind=dp) :: del

call generategroup(triamesh%t2p,triamesh%nt,triamesh%p, &
del,nearmat%row,nearmat%col,nearmat%nnz)
nfmm=nquad*triamesh%nt
allocate(chspace(nfmm))
allocate(fld(3,nfmm))


allocate(nearmat%val(nearmat%nnz))
allocate(triafl(3,nfmm))
allocate(chmat%col(nfmm))
allocate(chmat%row(nfmm))
allocate(chmat%val(nfmm))

call createbemdatastruct(triamesh%t2p,triamesh%nt,triamesh%p &
,triamesh%np,triamesh%vol,triamesh%nhat,nearmat%col,nearmat%row,nearmat%nnz &
,nearmat%val,triafl,nquad,chmat%col,chmat%row,chmat%val)
chmat%nnz=nfmm
call initialize_r0(triamesh%nt)


end subroutine initializematvec

subroutine destroymatvec()


deallocate(r0_initial)

deallocate(triafl)
deallocate(fld)
deallocate(chspace)

deallocate(nearmat%col)
deallocate(nearmat%row)
deallocate(nearmat%val)

deallocate(chmat%col)
deallocate(chmat%row)
deallocate(chmat%val)
end subroutine destroymatvec

subroutine initialize_r0(nt)
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  implicit none
  integer::jran,i,nt
  allocate(r0_initial(nt))
  jran=-3
  do i=1,nt
     r0_initial(i)=2.0_dp*real(ran1(jran),dp)-1.0_dp
  end do
  return
end subroutine initialize_r0


function ran1(idum)
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  implicit none
  integer,intent(inout)::idum
  integer,PARAMETER::IM=2147483647,IQ=127773,IR=2836,IA=16807,NTAB=32
  integer::NDIV
  real,parameter::EPS=1.2E-7,RNMX=1.-EPS
  real::am,ran1
  integer::j,k,iv(NTAB)=0,iy=0
  SAVE iv,iy

  NDIV=1+(IM-1)/NTAB
  AM=1./IM

  if (idum<= 0 .or. iy== 0) then
     idum=max(-idum,1)
     do j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum .lt. 0) idum=idum+IM
        if (j .le. NTAB) iv(j)=idum
     end do
     iy=iv(1)
  end if
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum .lt. 0) idum=idum+IM
  j=1+iy/NDIV
  iy=iv(j)
  iv(j)=idum
  ran1=min(AM*iy,RNMX)
  return
end function ran1
end module solver
